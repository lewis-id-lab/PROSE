source(file.path(proj_root, "R", "evidence_mapping.R"))
source(file.path(proj_root, "R", "ptx3_evidence.R"))
source(file.path(proj_root, "R", "source_verification.R"))
source(file.path(proj_root, "R", "term_mapping.R"))

zip_path <- file.path(proj_root, "IMD_genetics_analysis_bundle (1).zip")
hap_path <- file.path(proj_root, "data", "haplotype_frequencies_phased.csv")
freq_path <- file.path(proj_root, "data", "panel_allele_frequencies.csv")

test_that("single-SNP haplotypes come from the allele frequencies and mark the risk allele", {
  il10 <- single_snp_haplotypes("rs1800896", "C", freq_path)
  expect_equal(sum(il10$freq), 1)
  expect_equal(sum(il10$freq[il10$risk == 1]), 0.4866, tolerance = 1e-3)
  expect_setequal(il10$a1, c("C", "T"))
  # a zero-frequency allele (IFNG rs2069705 C, 0.0 in gnomAD NFE) is dropped
  ifng <- single_snp_haplotypes("rs2069705", "G", freq_path)
  expect_setequal(ifng$a1, c("A", "G"))
  expect_equal(sum(ifng$freq[ifng$risk == 1]), 0.3379 / (0.3379 + 0.6624), tolerance = 1e-3)
  expect_equal(single_snp_haplotypes("rs2069705", "G", freq_path, source = "1000G_EUR")$freq[2], 0.327, tolerance = 1e-3)
  expect_error(single_snp_haplotypes("rs1800896", "A", freq_path), "not an observed allele")
  expect_error(single_snp_haplotypes("rs0000000", "A", freq_path), "not found")
})

test_that("the CD209 risk class is G at rs4804800 and covers 12.5% of chromosomes", {
  h <- block_haplotypes("CD209", c("G", NA), hap_path)
  expect_setequal(h$haplotype[h$risk == 1], c("G-A", "G-G"))
  expect_equal(sum(h$freq[h$risk == 1]), 0.125, tolerance = 2e-3)
})

test_that("a single-SNP term maps a study that reports the other allele with the sign reversed", {
  h <- single_snp_haplotypes("rs2069705", "G", freq_path)
  row <- function(effect, other) data.frame(variant_row_id = 1, study_id = "S", rsid = "rs2069705", usable_tier = "estimate",
    genetic_model = "dominant", effect_measure = "OR", subgroup = "whole cohort", is_primary_analysis = "yes",
    effect_allele = effect, other_allele = other, orientation_flag = "fwd_ok", comparison_as_reported = NA_character_,
    effect_estimate = 1.5, ci_lower = 1.0, ci_upper = 2.25, stringsAsFactors = FALSE)
  g <- map_extraction_row(row("G", "A"), h, loci = "rs2069705", counts_model = "carrier")
  a <- map_extraction_row(row("A", "G"), h, loci = "rs2069705", counts_model = "carrier")
  expect_gt(g$beta, 0)
  expect_lt(a$beta, 0)
  # A carriers against GG with OR 1.5 means G is protective by a comparable amount
  expect_equal(a$beta / g$beta, -1, tolerance = 0.5)
})

test_that("the whole-cohort pattern covers the Lupianez and Fisher wording added for IFNG and CD209", {
  yes <- c("phase 1 cohort (n=593), patient (recipient) genotype", "overall combined phase 1+2 cohort (n=781), patient (recipient) genotype",
           "recipient (PT) genotype")
  expect_true(all(grepl(WHOLE_COHORT_PATTERN, yes, ignore.case = TRUE)))
  expect_false(any(grepl(WHOLE_COHORT_PATTERN, c("allo-HSCT recipients, DONOR genotype (n=171)", "donor genotype, SGM1 endpoint, 2 years"), ignore.case = TRUE)))
})

test_that("an overall combined analysis is preferred to the phase analysis it contains", {
  re <- data.frame(study_id = "S", rsid = "r", cohort = "whole", whole_cohort = TRUE, genotype = "recipient", odds_ratio_scale = TRUE,
                   non_european = FALSE, priority = c(2L, 3L), primary_flag = "no", variant_row_id = 1:2, combined = c(FALSE, TRUE),
                   beta = c(0.1, 0.5), se_beta = 0.2, method = "m", flag = "", stringsAsFactors = FALSE)
  expect_equal(ptx3_study_effects(re, snp_order = "r")$beta, 0.5)
  re$combined <- FALSE
  expect_equal(ptx3_study_effects(re, snp_order = "r")$beta, 0.1)
})

test_that("the hazard-ratio fallback applies only when the odds-ratio evidence is below the threshold", {
  skip_if_not(file.exists(zip_path) && file.exists(hap_path), "the evidence bundle or the haplotype table is not present")
  vr <- read_bundle_table("extraction_variants_round1.csv", zip_path); st <- read_bundle_table("extraction_studies_round1.csv", zip_path)
  out <- tempfile("term_mapping_"); on.exit(unlink(out, recursive = TRUE))
  il10 <- run_term("il10", single_snp_haplotypes("rs1800896", "C", freq_path), "rs1800896", snp_order = "rs1800896", seed = 11L,
                   vr = vr, st = st, hr_fallback_min = 3, out_dir = out)
  expect_match(names(il10$scenarios)[1], "^primary \\(fallback\\)")
  expect_equal(il10$pooled$k[1], nrow(il10$scenarios[[1]]))
  none <- run_term("il10", single_snp_haplotypes("rs1800896", "C", freq_path), "rs1800896", snp_order = "rs1800896", seed = 11L,
                   vr = vr, st = st, hr_fallback_min = NULL, out_dir = out)
  expect_match(names(none$scenarios)[1], "^primary: ")
  expect_true(file.exists(file.path(out, "il10_pooled.csv")))
})

test_that("CD209, IL10 and IFNG on the real extraction give the documented study sets", {
  skip_if_not(file.exists(zip_path) && file.exists(hap_path), "the evidence bundle or the haplotype table is not present")
  vr <- read_bundle_table("extraction_variants_round1.csv", zip_path); st <- read_bundle_table("extraction_studies_round1.csv", zip_path)
  hc <- block_haplotypes("CD209", c("G", NA), hap_path)
  rc <- ptx3_row_effects(vr, st, hc, loci = c("rs4804800", "rs7248637"), counts_model = "carrier")
  pc <- ptx3_study_effects(rc, snp_order = c("rs4804800", "rs7248637"))
  expect_setequal(pc$study_id, c("S045", "S071", "S077"))
  expect_true(all(pc$rsid == "rs4804800"))
  expect_true(all(pc$beta > 0))
  expect_true("S080" %in% ptx3_study_effects(rc, snp_order = c("rs4804800", "rs7248637"), allow_hr = TRUE)$study_id)

  ri <- ptx3_row_effects(vr, st, single_snp_haplotypes("rs2069705", "G", freq_path), loci = "rs2069705", counts_model = "carrier")
  pi <- ptx3_study_effects(ri, snp_order = "rs2069705")
  expect_setequal(pi$study_id, c("S046", "S072"))
  expect_gt(pi$beta[pi$study_id == "S046"], 0)
  expect_lt(pi$beta[pi$study_id == "S072"], 0)                                           # the two studies disagree in direction

  rl <- ptx3_row_effects(vr, st, single_snp_haplotypes("rs1800896", "C", freq_path), loci = "rs1800896", counts_model = "carrier")
  expect_equal(ptx3_study_effects(rl, snp_order = "rs1800896")$study_id, "S077")
  expect_true(all(rl$genotype[rl$study_id == "S061"] == "unstated"))
})
