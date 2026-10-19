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

# ---- The rs3901533 linkage (addendum 21): which CLEC7A haplotype carries Sainz's risk allele in East Asians and in Europeans? ----
# rs3901533 is the intronic SNP Sainz 2012 and Chen 2019 also typed. Ensembl labels it A/C on the forward strand; Sainz's G/T is the same SNP on
# the opposite strand (his G is the forward C, his risk allele T is the forward A). The phased genotypes are the Ensembl ones (cached); the European
# chromosomes are read from the per-haplotype rows of data/haplotype_ld_phased.csv, which come from the same source (script 04).
rs33_path <- "data/rs3901533_phased_genotypes.csv"
if (REFRESH || !file.exists(rs33_path)) {
  Sys.sleep(0.5)
  g33 <- fetch_json("https://rest.ensembl.org/variation/human/rs3901533?genotypes=1")$genotypes
  g33 <- g33[grepl("^1000GENOMES:phase_3:", g33$sample), ]
  g33$sample <- sub("^1000GENOMES:phase_3:", "", g33$sample)
  write.csv(g33[, c("sample", "genotype")], rs33_path, row.names = FALSE)
}
g33 <- read.csv(rs33_path, stringsAsFactors = FALSE)
g33 <- g33[grepl("^[AC]\\|[AC]$", g33$genotype), ]
sp33 <- strsplit(g33$genotype, "|", fixed = TRUE)
g33$a1 <- vapply(sp33, `[`, "", 1); g33$a2 <- vapply(sp33, `[`, "", 2)
alle <- rbind(data.frame(sample = g33$sample, allele = g33$a1), data.frame(sample = g33$sample, allele = g33$a2))
geno33 <- merge(alle, chrom[, c("sample", "pop", "rs16910526", "rs7309123")], by = "sample")
geno33 <- geno33[geno33$rs16910526 %in% c("A", "C") & geno33$rs7309123 %in% c("C", "G"), ]
geno33$clec7a_haplotype <- paste(geno33$rs16910526, geno33$rs7309123, sep = "-")

panel <- read.delim(sample_url, stringsAsFactors = FALSE)
eur <- panel$sample[panel$super_pop == "EUR"]
eur_g <- read.csv("data/phased_1000G_EUR_genotypes.csv", stringsAsFactors = FALSE)
eur_g <- eur_g[eur_g$rsid %in% c("rs16910526", "rs7309123") & grepl("^[ACGT]\\|[ACGT]$", eur_g$genotype), ]
sp_eur <- strsplit(eur_g$genotype, "|", fixed = TRUE)
eur_a <- do.call(rbind, lapply(1:2, function(k) { d <- data.frame(sample = eur_g$sample, allele = vapply(sp_eur, `[`, "", k), rsid = eur_g$rsid)
  reshape(d, idvar = "sample", timevar = "rsid", direction = "wide") }))
eur_x <- merge(alle[alle$sample %in% eur, ], eur_a, by = "sample")
eur_x <- eur_x[eur_x$allele.rs16910526 %in% c("A", "C") & eur_x$allele.rs7309123 %in% c("C", "G"), ]
eur_x$clec7a_haplotype <- paste(eur_x$allele.rs16910526, eur_x$allele.rs7309123, sep = "-")
ld_rows <- table(clec7a_haplotype = factor(eur_x$clec7a_haplotype, c("A-G", "A-C", "C-G")), allele = factor(eur_x$allele, c("C", "A")))
eur_h <- as.data.frame.matrix(ld_rows); eur_h$clec7a_haplotype <- rownames(eur_h)
eur_h <- cbind(data.frame(population = "European (EUR)", pops = "EUR"), transform(eur_h, chromosomes = C + A))
names(eur_h) <- sub("^(C|A)$", "rs3901533_\\1", names(eur_h))

groups33 <- list(`Chinese Han (CHB+CHS)` = c("CHB", "CHS"), `Thai proxy (CDX)` = "CDX", `East Asian (all five)` = pops)
hap33 <- do.call(rbind, Map(function(nm, pp) { x <- geno33[geno33$pop %in% pp, ]
  tb <- table(clec7a_haplotype = factor(x$clec7a_haplotype, c("A-G", "A-C", "C-G")), allele = factor(x$allele, c("C", "A")))
  d <- as.data.frame.matrix(tb); d$clec7a_haplotype <- rownames(d)
  cbind(data.frame(population = nm, pops = paste(pp, collapse = "+")), transform(d, chromosomes = C + A)) }, names(groups33), groups33))
names(hap33) <- sub("^(C|A)$", "rs3901533_\\1", names(hap33))
out33 <- rbind(hap33[, c("population", "pops", "clec7a_haplotype", "rs3901533_C", "rs3901533_A", "chromosomes")], eur_h[, c("population", "pops", "clec7a_haplotype", "rs3901533_C", "rs3901533_A", "chromosomes")])
rownames(out33) <- NULL
write.csv(out33, "data/rs3901533_linkage_eas.csv", row.names = FALSE)
message("\nrs3901533 alleles on each CLEC7A haplotype (Sainz's risk allele is the forward A):")
print(out33, row.names = FALSE)
