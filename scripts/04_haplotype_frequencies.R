# Haplotype frequencies from phased 1000 Genomes phase 3 European genotypes (via the Ensembl
# REST API), rescaled to gnomAD v4 non-Finnish-European allele-frequency margins.
# Run from the project root. Downloads are cached in data/; set REFRESH <- TRUE to refetch.

source("R/haplotypes.R")
suppressPackageStartupMessages(library(jsonlite))

REFRESH <- FALSE
geno_path  <- "data/phased_1000G_EUR_genotypes.csv"
sample_url <- "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel"

snps <- c("rs2305619", "rs3816527", "rs1840680",   # PTX3
          "rs16910526", "rs7309123",               # CLEC7A
          "rs4986790", "rs4986791",                # TLR4
          "rs4804800", "rs7248637")                # CD209

# 1. Phased genotypes for the European samples
if (REFRESH || !file.exists(geno_path)) {
  panel <- read.delim(sample_url, stringsAsFactors = FALSE)
  eur <- panel$sample[panel$super_pop == "EUR"]
  fetch_json <- function(url, tries = 5) {
    tmp <- tempfile(fileext = ".json")
    for (i in seq_len(tries)) {
      status <- system2("curl", c("-sS", "-f", "-m", "90", "-H", shQuote("Content-Type: application/json"),
                                  shQuote(url), "-o", shQuote(tmp)), stderr = FALSE)
      if (status == 0) return(fromJSON(tmp))
      Sys.sleep(3 * i)
    }
    stop("Ensembl request failed after ", tries, " tries: ", url)
  }
  fetch_one <- function(rs) {
    url <- sprintf("https://rest.ensembl.org/variation/human/%s?genotypes=1", rs)
    g <- fetch_json(url)$genotypes
    g <- g[grepl("^1000GENOMES:phase_3:", g$sample), ]
    g$sample <- sub("^1000GENOMES:phase_3:", "", g$sample)
    g <- g[g$sample %in% eur, c("sample", "genotype")]
    g$rsid <- rs
    g
  }
  geno <- do.call(rbind, lapply(snps, function(rs) { Sys.sleep(0.5); fetch_one(rs) }))
  write.csv(geno, geno_path, row.names = FALSE)
}
geno <- read.csv(geno_path, stringsAsFactors = FALSE)

# 2. Phased chromosomes: one row per chromosome, one column per SNP
stopifnot(all(grepl("^[ACGT]\\|[ACGT]$", geno$genotype)))
split_gt <- strsplit(geno$genotype, "|", fixed = TRUE)
geno$h1 <- vapply(split_gt, `[`, "", 1)
geno$h2 <- vapply(split_gt, `[`, "", 2)
samples <- sort(unique(geno$sample))
stopifnot(all(table(geno$rsid) == length(samples)))
wide <- function(col) {
  m <- sapply(snps, function(rs) {
    g <- geno[geno$rsid == rs, ]
    g[[col]][match(samples, g$sample)]
  })
  data.frame(sample = samples, m, check.names = FALSE)
}
chrom <- rbind(wide("h1"), wide("h2"))

# 3. gnomAD v4 NFE margins and the two analysed alleles per SNP
af <- read.csv("data/panel_allele_frequencies.csv", stringsAsFactors = FALSE)
nfe <- lapply(setNames(af$gnomad_nfe_freqs, af$rsid), fromJSON)
nfe <- lapply(nfe, function(x) unlist(x)[unlist(x) > 0])
stopifnot(all(sapply(nfe[snps], length) == 2))
# Order alleles so the 1000G minor allele comes first within each pair (stable, readable)
eur_minor <- sapply(snps, function(rs) {
  tb <- table(factor(chrom[[rs]], levels = names(nfe[[rs]])))
  names(tb)[which.min(tb)]
})
alleles <- lapply(setNames(snps, snps), function(rs) {
  a <- names(nfe[[rs]]); c(eur_minor[[rs]], setdiff(a, eur_minor[[rs]]))
})

blocks <- list(
  PTX3_trio  = c("rs2305619", "rs3816527", "rs1840680"),
  PTX3_a     = c("rs2305619", "rs3816527"),
  PTX3_b     = c("rs2305619", "rs1840680"),
  PTX3_c     = c("rs3816527", "rs1840680"),
  CLEC7A     = c("rs16910526", "rs7309123"),
  TLR4       = c("rs4986790", "rs4986791"),
  CD209      = c("rs4804800", "rs7248637")
)

# 4. Observed counts, rescaling, LD before and after
rows <- list(); ld_rows <- list()
for (b in names(blocks)) {
  loci <- blocks[[b]]
  tab <- haplotype_table(chrom, alleles[loci])
  margins <- lapply(loci, function(rs) nfe[[rs]])
  resc <- ipf_haplotypes(tab, margins)
  obs <- unclass(tab) / sum(tab)

  cells <- expand.grid(dimnames(tab), stringsAsFactors = FALSE)
  names(cells) <- loci
  out <- data.frame(block = b, loci = paste(loci, collapse = ";"),
                    haplotype = apply(cells, 1, paste, collapse = "-"),
                    cells, check.names = FALSE, row.names = NULL)
  out$n_chromosomes <- as.vector(unclass(tab))
  out$freq_1000G_EUR <- as.vector(obs)
  out$se_1000G_EUR <- sqrt(out$freq_1000G_EUR * (1 - out$freq_1000G_EUR) / attr(tab, "n_used"))
  out$freq_gnomAD_NFE_rescaled <- as.vector(resc)
  names(out)[names(out) %in% loci] <- paste0("allele_", seq_along(loci))
  if (length(loci) == 2) {
    dp <- rescale_pair_dprime(tab, margins[[1]][dimnames(tab)[[1]][1]], margins[[2]][dimnames(tab)[[2]][1]])
    out$freq_gnomAD_NFE_dprime <- as.vector(dp)
  }
  rows[[b]] <- out

  if (length(loci) == 2) {
    s0 <- ld_stats(tab); s1 <- ld_stats(resc)
    ld_rows[[b]] <- data.frame(block = b, loci = paste(loci, collapse = ";"),
                               n_chromosomes = attr(tab, "n_used"), n_dropped = attr(tab, "n_dropped"),
                               r2_1000G = s0[["r2"]], dprime_1000G = s0[["dprime"]],
                               r2_rescaled = s1[["r2"]], dprime_rescaled = s1[["dprime"]],
                               ipf_iterations = attr(resc, "iterations"),
                               max_margin_error = attr(resc, "max_margin_error"))
  } else {
    ld_rows[[b]] <- data.frame(block = b, loci = paste(loci, collapse = ";"),
                               n_chromosomes = attr(tab, "n_used"), n_dropped = attr(tab, "n_dropped"),
                               r2_1000G = NA, dprime_1000G = NA, r2_rescaled = NA, dprime_rescaled = NA,
                               ipf_iterations = attr(resc, "iterations"),
                               max_margin_error = attr(resc, "max_margin_error"))
  }
}
res <- dplyr::bind_rows(rows)
write.csv(res, "data/haplotype_frequencies_phased.csv", row.names = FALSE)
write.csv(do.call(rbind, ld_rows), "data/haplotype_ld_phased.csv", row.names = FALSE)
