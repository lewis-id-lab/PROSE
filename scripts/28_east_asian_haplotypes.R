# East Asian CLEC7A haplotype frequencies from phased 1000 Genomes phase 3 genotypes (via the Ensembl REST API), for the second CLEC7A
# parameter (addendum 21). Run from the project root. Downloads are cached in data/; set REFRESH <- TRUE to refetch (needs network).
#
# rs16910526 (Y238X, forward-strand C) and rs7309123 (forward-strand G/C) in the five East Asian populations: CHB (Han Beijing), CHS (Southern Han),
# CDX (Dai in Xishuangbanna, the closest 1000 Genomes population to the Thai), KHV (Kinh, Vietnam) and JPT (Japanese). Population groups used for the
# studies: Chinese Han = CHB + CHS; Thai = CDX; all East Asian = the five populations.
#
# Writes   data/phased_1000G_EAS_genotypes.csv (cache) and data/haplotype_frequencies_eas.csv

suppressPackageStartupMessages(library(jsonlite))
REFRESH <- FALSE
geno_path <- "data/phased_1000G_EAS_genotypes.csv"
sample_url <- "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel"
snps <- c("rs16910526", "rs7309123")
pops <- c("CHB", "CHS", "CDX", "KHV", "JPT")

if (REFRESH || !file.exists(geno_path)) {
  panel <- read.delim(sample_url, stringsAsFactors = FALSE)
  eas <- panel[panel$super_pop == "EAS", c("sample", "pop")]
  fetch_json <- function(url, tries = 5) {
    tmp <- tempfile(fileext = ".json")
    for (i in seq_len(tries)) {
      status <- system2("curl", c("-sS", "-f", "-m", "120", "-H", shQuote("Content-Type: application/json"), shQuote(url), "-o", shQuote(tmp)), stderr = FALSE)
      if (status == 0) return(fromJSON(tmp))
      Sys.sleep(3 * i)
    }
    stop("Ensembl request failed after ", tries, " tries: ", url)
  }
  geno <- do.call(rbind, lapply(snps, function(rs) {
    Sys.sleep(0.5)
    g <- fetch_json(sprintf("https://rest.ensembl.org/variation/human/%s?genotypes=1", rs))$genotypes
    g <- g[grepl("^1000GENOMES:phase_3:", g$sample), ]
    g$sample <- sub("^1000GENOMES:phase_3:", "", g$sample)
    g <- g[g$sample %in% eas$sample, c("sample", "genotype")]
    g$rsid <- rs
    g$pop <- eas$pop[match(g$sample, eas$sample)]
    g
  }))
  write.csv(geno, geno_path, row.names = FALSE)
}
geno <- read.csv(geno_path, stringsAsFactors = FALSE)

# Phased chromosomes, one row per chromosome
ok <- grepl("^[ACGT]\\|[ACGT]$", geno$genotype)
message("genotype strings not phased or not single bases: ", sum(!ok), " of ", nrow(geno))
geno <- geno[ok, ]
sp <- strsplit(geno$genotype, "|", fixed = TRUE)
geno$h1 <- vapply(sp, `[`, "", 1); geno$h2 <- vapply(sp, `[`, "", 2)
samples <- sort(unique(geno$sample))
stopifnot(all(table(geno$rsid) == length(samples)))
chrom <- do.call(rbind, lapply(c("h1", "h2"), function(col) {
  d <- data.frame(sample = samples, pop = geno$pop[match(samples, geno$sample)])
  for (rs in snps) { g <- geno[geno$rsid == rs, ]; d[[rs]] <- g[[col]][match(samples, g$sample)] }
  d
}))

groups <- list(`Chinese Han (CHB+CHS)` = c("CHB", "CHS"), `Thai proxy (CDX)` = "CDX", `Kinh (KHV)` = "KHV", `Japanese (JPT)` = "JPT",
               `East Asian (all five)` = pops, CHB = "CHB", CHS = "CHS")
haps <- function(x) {
  x <- x[x$rs16910526 %in% c("A", "C") & x$rs7309123 %in% c("C", "G"), ]
  tb <- table(factor(x$rs16910526, levels = c("C", "A")), factor(x$rs7309123, levels = c("G", "C")))
  d <- as.data.frame(tb, stringsAsFactors = FALSE); names(d) <- c("allele_1", "allele_2", "n_chromosomes")
  d$haplotype <- paste(d$allele_1, d$allele_2, sep = "-")
  d$freq <- d$n_chromosomes / sum(d$n_chromosomes); d$chromosomes_total <- nrow(x)
  d[, c("haplotype", "allele_1", "allele_2", "n_chromosomes", "freq", "chromosomes_total")]
}
out <- do.call(rbind, lapply(names(groups), function(g) cbind(population = g, pops = paste(groups[[g]], collapse = "+"), haps(chrom[chrom$pop %in% groups[[g]], ]))))
rownames(out) <- NULL
write.csv(out, "data/haplotype_frequencies_eas.csv", row.names = FALSE)
message("CLEC7A haplotype frequencies (rs16910526 - rs7309123):")
print(reshape(out[, c("population", "haplotype", "freq")], idvar = "population", timevar = "haplotype", direction = "wide"), digits = 3, row.names = FALSE)
message("chromosomes per group: ", paste(unique(paste0(out$population, " ", out$chromosomes_total)), collapse = "; "))
