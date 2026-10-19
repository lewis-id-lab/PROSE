source(file.path(proj_root, "R", "haplotypes.R"))

test_that("ld_stats gives D = 0 for independence and r2 = D' = 1 for perfect association", {
  indep <- as.table(outer(c(0.3, 0.7), c(0.4, 0.6)))
  s <- ld_stats(indep)
  expect_equal(s[["D"]], 0, tolerance = 1e-12)
  expect_equal(s[["r2"]], 0, tolerance = 1e-12)

  perfect <- as.table(matrix(c(0.4, 0, 0, 0.6), 2))
  s <- ld_stats(perfect)
  expect_equal(s[["r2"]], 1)
  expect_equal(s[["dprime"]], 1)
})

test_that("ld_stats reports negative D and D' when the minor alleles are in repulsion", {
  rep_tab <- as.table(matrix(c(0.05, 0.25, 0.25, 0.45), 2))
  s <- ld_stats(rep_tab)
  expect_lt(s[["D"]], 0)
  expect_lt(s[["dprime"]], 0)
})

test_that("IPF reproduces target margins and preserves the odds ratio of a 2x2 table", {
  tab <- as.table(matrix(c(40, 10, 15, 35), 2, dimnames = list(c("a", "A"), c("b", "B"))))
  names(dimnames(tab)) <- c("l1", "l2")
  m1 <- c(a = 0.30, A = 0.70); m2 <- c(b = 0.55, B = 0.45)
  out <- ipf_haplotypes(tab, list(m1, m2))
  expect_equal(sum(out), 1)
  expect_equal(rowSums(out), m1, tolerance = 1e-9)
  expect_equal(colSums(out), m2, tolerance = 1e-9)
  expect_equal(out[1, 1] * out[2, 2] / (out[1, 2] * out[2, 1]), 40 * 35 / (10 * 15), tolerance = 1e-7)
})

test_that("IPF on three loci matches margins and preserves the three-way interaction", {
  set.seed(6408)
  arr <- array(runif(8, 1, 20), c(2, 2, 2), dimnames = list(l1 = c("a", "A"), l2 = c("b", "B"), l3 = c("c", "C")))
  m <- list(c(a = 0.2, A = 0.8), c(b = 0.5, B = 0.5), c(c = 0.35, C = 0.65))
  out <- ipf_haplotypes(as.table(arr), m)
  for (k in 1:3) expect_equal(apply(out, k, sum), m[[k]][dimnames(arr)[[k]]], tolerance = 1e-9)
  tri <- function(x) (x[1, 1, 1] * x[1, 2, 2] * x[2, 1, 2] * x[2, 2, 1]) /
                     (x[1, 1, 2] * x[1, 2, 1] * x[2, 1, 1] * x[2, 2, 2])
  expect_equal(tri(out), tri(arr), tolerance = 1e-6)
})

test_that("cells that start at zero stay at zero after IPF", {
  tab <- as.table(matrix(c(30, 0, 20, 50), 2, dimnames = list(l1 = c("a", "A"), l2 = c("b", "B"))))
  out <- ipf_haplotypes(tab, list(c(a = 0.45, A = 0.55), c(b = 0.30, B = 0.70)))
  expect_equal(out[2, 1], 0)
})

test_that("the D'-preserving rescale keeps D' and hits the new margins", {
  tab <- as.table(matrix(c(0.30, 0.05, 0.10, 0.55), 2, dimnames = list(l1 = c("a", "A"), l2 = c("b", "B"))))
  new <- rescale_pair_dprime(tab, c(a = 0.35), c(b = 0.50))
  expect_equal(sum(new), 1)
  expect_equal(sum(new[1, ]), 0.35)
  expect_equal(sum(new[, 1]), 0.50)
  expect_equal(ld_stats(new)[["dprime"]], ld_stats(tab)[["dprime"]], tolerance = 1e-9)
})

freq_path <- file.path(proj_root, "data", "haplotype_frequencies_phased.csv")
ld_path   <- file.path(proj_root, "data", "haplotype_ld_phased.csv")

test_that("phased frequencies are proper distributions on 1,006 European chromosomes", {
  skip_if_not(file.exists(freq_path), "run scripts/04_haplotype_frequencies.R first")
  h <- read.csv(freq_path)
  for (b in unique(h$block)) {
    x <- h[h$block == b, ]
    expect_equal(sum(x$n_chromosomes), 1006)
    expect_equal(sum(x$freq_1000G_EUR), 1)
    expect_equal(sum(x$freq_gnomAD_NFE_rescaled), 1, tolerance = 1e-9)
  }
})

test_that("rescaled frequencies reproduce the gnomAD NFE allele frequencies", {
  skip_if_not(file.exists(freq_path), "run scripts/04_haplotype_frequencies.R first")
  h <- read.csv(freq_path)
  af <- read.csv(file.path(proj_root, "data", "panel_allele_frequencies.csv"))
  nfe <- lapply(setNames(af$gnomad_nfe_freqs, af$rsid), jsonlite::fromJSON)
  for (b in unique(h$block)) {
    x <- h[h$block == b, ]
    loci <- strsplit(x$loci[1], ";")[[1]]
    for (k in seq_along(loci)) {
      got <- tapply(x$freq_gnomAD_NFE_rescaled, x[[paste0("allele_", k)]], sum)
      want <- unlist(nfe[[loci[k]]])[names(got)]
      expect_equal(as.numeric(got), as.numeric(want), tolerance = 1e-6, info = paste(b, loci[k]))
    }
  }
})

test_that("minor alleles are in coupling in every pair, directly from phased chromosomes", {
  skip_if_not(file.exists(ld_path), "run scripts/04_haplotype_frequencies.R first")
  h <- read.csv(freq_path)
  for (b in setdiff(unique(h$block), "PTX3_trio")) {
    x <- h[h$block == b, ]
    tab <- matrix(x$freq_1000G_EUR, 2, 2)
    expect_gt(ld_stats(as.table(tab))[["D"]], 0)
  }
})

test_that("phased LD agrees with the earlier derived LD, and the PTX3 trio has three common haplotypes", {
  skip_if_not(file.exists(ld_path), "run scripts/04_haplotype_frequencies.R first")
  ld  <- read.csv(ld_path)
  old <- read.csv(file.path(proj_root, "data", "haplotype_frequencies.csv"))
  old_r2 <- tapply(old$r2, paste(old$locus_1, old$locus_2, sep = ";"), unique)
  for (i in seq_len(nrow(ld))) {
    key <- ld$loci[i]
    if (key %in% names(old_r2) && !grepl("rs4804800", key)) {
      expect_equal(ld$r2_1000G[i], as.numeric(old_r2[[key]]), tolerance = 0.005, info = key)
    }
  }
  h <- read.csv(freq_path)
  trio <- h[h$block == "PTX3_trio", ]
  expect_equal(sum(trio$freq_1000G_EUR > 0.01), 3L)
  # The cells that are absent in 1000G stay absent after rescaling
  expect_true(all(trio$freq_gnomAD_NFE_rescaled[trio$n_chromosomes == 0] == 0))
})

test_that("the East Asian CLEC7A haplotype table is a proper frequency table with the expected structure", {
  path <- file.path(proj_root, "data", "haplotype_frequencies_eas.csv")
  skip_if_not(file.exists(path), "the East Asian haplotype table is not present")
  eas <- read.csv(path, stringsAsFactors = FALSE)
  expect_setequal(unique(eas$haplotype), c("C-G", "A-G", "C-C", "A-C"))
  expect_true(all(c("Chinese Han (CHB+CHS)", "Thai proxy (CDX)", "East Asian (all five)") %in% eas$population))
  for (g in unique(eas$population)) {
    x <- eas[eas$population == g, ]
    expect_equal(sum(x$freq), 1, tolerance = 1e-9)
    expect_equal(sum(x$n_chromosomes), x$chromosomes_total[1])
    expect_equal(x$freq, x$n_chromosomes / x$chromosomes_total, tolerance = 1e-9)
    # Y238X (C at rs16910526) is nearly absent and never found with the C allele of rs7309123
    expect_lt(sum(x$freq[x$allele_1 == "C"]), 0.02)
    expect_equal(x$n_chromosomes[x$haplotype == "C-C"], 0)
    # rs7309123 G is the common allele (0.77 to 0.88), the opposite of the European 0.43
    expect_gt(sum(x$freq[x$allele_2 == "G"]), 0.7)
  }
})

test_that("the rs3901533 linkage table shows the East Asian and European risk haplotypes on opposite rs3901533 backgrounds", {
  path <- file.path(proj_root, "data", "rs3901533_linkage_eas.csv")
  skip_if_not(file.exists(path), "the rs3901533 linkage table is not present")
  lk <- read.csv(path, stringsAsFactors = FALSE)
  expect_equal(lk$rs3901533_C + lk$rs3901533_A, lk$chromosomes)
  # Sainz's risk allele is the forward A: it sits on the A-G haplotype in Europeans and on the A-C haplotype in East Asians
  eas <- lk[lk$population == "East Asian (all five)" & lk$clec7a_haplotype == "A-G", ]
  eur <- lk[lk$population == "European (EUR)" & lk$clec7a_haplotype == "A-G", ]
  expect_gt(eas$rs3901533_A / eas$chromosomes, 0.8); expect_lt(eur$rs3901533_A / eur$chromosomes, 0.45)
  ac_eas <- lk[lk$population == "East Asian (all five)" & lk$clec7a_haplotype == "A-C", ]
  ac_eur <- lk[lk$population == "European (EUR)" & lk$clec7a_haplotype == "A-C", ]
  expect_gt(ac_eur$rs3901533_C / ac_eur$chromosomes, 0.8)
  # In East Asians rs3901533 A tags rs7309123 G almost perfectly
  cg <- lk[lk$population == "East Asian (all five)" & lk$clec7a_haplotype == "C-G", ]
  t2 <- as.table(matrix(c(ac_eas$rs3901533_C, eas$rs3901533_C + cg$rs3901533_C, ac_eas$rs3901533_A, eas$rs3901533_A + cg$rs3901533_A), 2, byrow = TRUE,
                        dimnames = list(c("C", "A"), c("C", "G"))))
  s <- ld_stats(t2); expect_gt(s[["dprime"]], 0.4); expect_gt(s[["r2"]], 0.15)  # D > 0 for C-C: rs3901533 C tags rs7309123 C, so A tags G
})
