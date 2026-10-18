source(file.path(proj_root, "R", "evidence_mapping.R"))
source(file.path(proj_root, "R", "ptx3_evidence.R"))

# A small PTX3 haplotype table with the real structure: loci rs2305619, rs3816527, rs1840680; h2 = A at rs2305619
hap <- data.frame(haplotype = c("A-C-A", "A-A-A", "G-A-G"), a1 = c("A", "A", "G"), a2 = c("C", "A", "A"), a3 = c("A", "A", "G"),
                  freq = c(0.40, 0.06, 0.54), h2 = c(1L, 1L, 0L), stringsAsFactors = FALSE)

test_that("genotype tokens are read from free text and normalised", {
  expect_setequal(parse_genotype_tokens("AG+AA"), c("AG", "AA"))
  expect_setequal(parse_genotype_tokens("A/G + G/A"), "AG")
  expect_setequal(parse_genotype_tokens("G / G"), "GG")
  expect_length(parse_genotype_tokens("AML patients, reference group"), 0)
})

test_that("contrasts are read from the comparison text, including which side is the reference", {
  a <- parse_contrast("AG+AA vs GG", "dominant", "A", "G")
  expect_equal(a$method, "text"); expect_setequal(a$numerator, c("AA", "AG")); expect_equal(a$reference, "GG")
  b <- parse_contrast("GG vs AG+AA", "dominant", "A", "G")
  expect_equal(b$numerator, "GG"); expect_setequal(b$reference, c("AA", "AG"))
  c2 <- parse_contrast("AA (ref) vs AG+GG combined", "dominant", "A", "G")          # the reference is on the left
  expect_equal(c2$reference, "AA"); expect_setequal(c2$numerator, c("AG", "GG"))
  d <- parse_contrast("AC+CC vs AA (reference A)", "dominant", "C", "A")             # the reference is on the right
  expect_setequal(d$numerator, c("AC", "CC")); expect_equal(d$reference, "AA")
  e <- parse_contrast("G/G + A/G vs A/A", "dominant", "G", "A")
  expect_setequal(e$numerator, c("GG", "AG")); expect_equal(e$reference, "AA")
})

test_that("unreadable text falls back to the model label, and some cases stay unmapped", {
  dom <- parse_contrast(NA, "dominant", "A", "G")
  expect_equal(dom$method, "model_label"); expect_setequal(dom$numerator, c("AA", "AG")); expect_equal(dom$reference, "GG")
  rec <- parse_contrast(NA, "recessive", "A", "G")
  expect_equal(rec$numerator, "AA"); expect_setequal(rec$reference, c("AG", "GG"))
  expect_equal(parse_contrast(NA, "allelic", "A", "G")$type, "allele")
  expect_equal(parse_contrast("A/A + A/C genotype", "dominant", "A", "G")$method, "model_label")   # allele typo: text not usable
  expect_equal(parse_contrast("AA vs AG vs GG", "genotypic_3level", "A", "G")$flag, "more than two groups")
  expect_equal(parse_contrast(NA, "other", "A", "G")$method, "unreadable")
  expect_equal(parse_contrast("AG+AA vs GG", "dominant", NA, NA)$method, "unreadable")
})

test_that("genotype pairs form a distribution with the right h2 allele frequency", {
  p <- genotype_pairs(hap, 1)
  expect_equal(sum(p$prob), 1)
  expect_equal(sum(p$prob * p$k) / 2, 0.46)                                     # h2 frequency: 0.40 + 0.06
  expect_equal(sum(p$prob[p$genotype == "AA"]), 0.46^2)
  expect_equal(sum(p$prob[p$genotype == "AG"]), 2 * 0.46 * 0.54)
})

test_that("the implied effect of a contrast is zero at beta = 0 and matches hand calculations", {
  p1 <- genotype_pairs(hap, 1)
  expect_equal(implied_log_or_genotype(0, p1, c("AA", "AG"), "GG"), 0)
  b <- 0.7; pa <- 0.46; qa <- 0.54
  dom <- log(((pa^2 * exp(2 * b) + 2 * pa * qa * exp(b)) / qa^2) / ((pa^2 + 2 * pa * qa) / qa^2))
  expect_equal(implied_log_or_genotype(b, p1, c("AA", "AG"), "GG"), dom)
  rec <- log(((pa^2 * exp(2 * b)) / (2 * pa * qa * exp(b) + qa^2)) / ((pa^2) / (2 * pa * qa + qa^2)))
  expect_equal(implied_log_or_genotype(b, p1, "AA", c("AG", "GG")), rec)
})

test_that("the rs2305619 A allele is the h2 class, so its allelic effect is exactly beta; rs3816527 C is attenuated", {
  for (b in c(-0.8, 0.3, 1.1)) expect_equal(implied_log_or_allele(b, hap, 1, "A"), b)
  b <- 0.8
  expect_equal(implied_log_or_allele(b, hap, 3, "A"), b)                         # rs1840680 A is the same class here
  analytic <- log(0.40 * exp(b) / (0.06 * exp(b) + 0.54)) - log(0.40 / (0.06 + 0.54))
  expect_equal(implied_log_or_allele(b, hap, 2, "C"), analytic)
  expect_lt(implied_log_or_allele(b, hap, 2, "C"), b)
  expect_gt(implied_log_or_allele(b, hap, 2, "C"), 0.6 * b)                      # modest attenuation, not a different effect
})

test_that("solve_beta inverts the implied effect, refuses impossible or non-monotone cases, and gives the delta-method slope", {
  p1 <- genotype_pairs(hap, 1)
  f <- function(b) implied_log_or_genotype(b, p1, c("AA", "AG"), "GG")
  s <- solve_beta(f(0.7), f)
  expect_equal(s$beta, 0.7, tolerance = 1e-6)
  expect_equal(s$slope, (f(0.7 + 1e-4) - f(0.7 - 1e-4)) / 2e-4, tolerance = 1e-6)
  expect_true(is.na(solve_beta(50, f)$beta))
  het <- function(b) implied_log_or_genotype(b, p1, "AG", c("AA", "GG"))          # heterozygotes only: rises then falls
  expect_match(solve_beta(0.1, het)$flag, "not monotone")
  g <- function(b) implied_log_or_allele(b, hap, 1, "A")
  expect_equal(solve_beta(0.4, g)$slope, 1, tolerance = 1e-6)
})

test_that("map_to_h2 converts a standard error by the delta method", {
  ct <- list(type = "allele", allele = "A", flag = "")
  m <- map_to_h2(0.4, 0.2, ct, hap, 1)
  expect_equal(m$beta, 0.4, tolerance = 1e-6); expect_equal(m$se_beta, 0.2, tolerance = 1e-6)
  dom <- list(type = "genotype", numerator = c("AA", "AG"), reference = "GG", flag = "")
  f <- function(b) implied_log_or_genotype(b, genotype_pairs(hap, 1), dom$numerator, dom$reference)
  y <- f(0.5)
  m2 <- map_to_h2(y, 0.3, dom, hap, 1)
  expect_equal(m2$beta, 0.5, tolerance = 1e-6)
  expect_equal(m2$se_beta, 0.3 / ((f(0.5 + 1e-4) - f(0.5 - 1e-4)) / 2e-4), tolerance = 1e-6)
  expect_true(is.na(map_to_h2(0.4, 0.2, list(type = NA, flag = "no contrast"), hap, 1)$beta))
})

test_that("effects are taken from estimates and counts as documented", {
  ok <- log_or_from_estimate(1.64, 0.81, 3.48)
  expect_equal(ok$y, log(1.64)); expect_equal(ok$se, (log(3.48) - log(0.81)) / (2 * qnorm(0.975)))
  expect_true(is.na(log_or_from_estimate(1.76, 0.92, 1.76)$y))                    # estimate equal to its upper limit
  expect_true(is.na(log_or_from_estimate(1.2, 0.8, 0.5)$y))
  expect_true(is.na(log_or_from_estimate(NA, 0.8, 1.5)$y))
  expect_true(is.na(log_or_from_estimate(-1, 0.8, 1.5)$y))

  ctrl <- c(4000, 10000, 6000)
  cases <- round(ctrl * 0.05 * exp(0.5 * c(2, 1, 0)))
  tr <- log_or_from_counts_2x3(cases, ctrl)
  expect_equal(tr$y, 0.5, tolerance = 0.02)
  expect_gt(tr$se, 0)
  expect_true(is.na(log_or_from_counts_2x3(c(0, 0, 0), c(5, 5, 5))$y))

  t2 <- log_or_from_counts_2x2(11, 15, 46, 140)
  expect_equal(t2$y, log((11 / 46) / (15 / 140))); expect_equal(t2$se, sqrt(1 / 11 + 1 / 15 + 1 / 46 + 1 / 140))
  t0 <- log_or_from_counts_2x2(0, 4, 20, 30)                                      # a zero cell gets 0.5 added
  expect_equal(t0$y, log((0.5 / 20.5) / (4.5 / 30.5)))
  expect_true(is.na(log_or_from_counts_2x2(NA, 1, 2, 3)$y))
})

test_that("independent cohorts of one study combine by inverse variance", {
  y <- c(0.0, 2.0); se <- c(0.4, 0.8)
  cb <- combine_fixed(y, se)
  w <- 1 / se^2
  expect_equal(cb$y, sum(w * y) / sum(w)); expect_equal(cb$se, sqrt(1 / sum(w)))
  expect_equal(cb$q, sum(w * (y - cb$y)^2)); expect_equal(cb$k, 2L)
})

test_that("the Bayesian random-effects pooling is correct for fixed tau, recovers known parameters and is reproducible", {
  y <- c(0.2, 0.5, 0.9); se <- c(0.3, 0.25, 0.4); tau <- 0.2; sd0 <- 1
  mp <- mu_given_tau(y, se, tau, sd0)
  w <- 1 / (se^2 + tau^2)
  expect_equal(mp$mean, sum(w * y) / (sum(w) + 1 / sd0^2)); expect_equal(mp$sd, sqrt(1 / (sum(w) + 1 / sd0^2)))

  set.seed(1204)
  k <- 300; se_k <- runif(k, 0.15, 0.35); y_k <- rnorm(k, 0.4, sqrt(0.3^2 + se_k^2))
  big <- bayes_meta_grid(y_k, se_k, seed = 3)
  expect_equal(big$mu[2], 0.4, tolerance = 0.05)
  expect_equal(big$tau[2], 0.3, tolerance = 0.06)
  expect_lt(big$mu[3] - big$mu[1], 0.12)

  same <- bayes_meta_grid(rep(0.3, 6), rep(0.1, 6), seed = 4)
  expect_equal(same$mu[2], 0.3, tolerance = 0.03); expect_lt(same$tau[2], 0.1)
  a <- bayes_meta_grid(y, se, seed = 9); b <- bayes_meta_grid(y, se, seed = 9)
  expect_identical(a$draws, b$draws)
  expect_equal(colnames(a$draws), c("mu", "tau"))
  expect_gt(a$prediction[3] - a$prediction[1], a$mu[3] - a$mu[1])              # the prediction interval is wider
  # a tighter prior on tau pulls the heterogeneity down when three studies disagree
  yh <- c(0, 0.6, 1.2); sh <- rep(0.2, 3)
  expect_lt(bayes_meta_grid(yh, sh, seed = 2, prior_scale_tau = 0.05)$tau[2], bayes_meta_grid(yh, sh, seed = 2, prior_scale_tau = 0.5)$tau[2])
})

test_that("ancestry is judged by the leading descriptor, and Caucasian is not Asian", {
  expect_false(non_european_study("Caucasian"))
  expect_false(non_european_study("Caucasian (predominantly; 4% other than Caucasian including South American, Asian, African"))
  expect_true(non_european_study("Chinese (South China)"))
  expect_true(non_european_study("Thai"))
  expect_false(non_european_study("white"))
  expect_false(non_european_study("Turkish population"))
  expect_false(non_european_study(NA))
})

test_that("whole-cohort subgroups are recognised, including donor genotype in the whole cohort", {
  yes <- c("whole cohort", "All episodes (N=212)", "stratification cohort (whole cohort)", "donor genotype; whole cohort")
  no <- c("No neutropenia at presentation (N=133)", "ABPA group (n=15)", "donor genotype; SGM1 IA endpoint; 3 months", "Cohort 1 (Affymetrix 5.0)")
  expect_true(all(grepl(WHOLE_COHORT_PATTERN, yes, ignore.case = TRUE)))
  expect_false(any(grepl(WHOLE_COHORT_PATTERN, no, ignore.case = TRUE)))
})

row <- function(...) {
  base <- data.frame(variant_row_id = "V1", study_id = "S1", rsid = "rs2305619", usable_tier = "estimate", genetic_model = "dominant",
                     effect_measure = "OR", subgroup = "whole cohort", is_primary_analysis = "yes", comparison_as_reported = NA_character_,
                     effect_allele = "A", other_allele = "G", effect_estimate = 1.64, ci_lower = 0.81, ci_upper = 3.48,
                     n_cases_hom_effect = NA, n_cases_het = NA, n_cases_hom_other = NA, n_ctrl_hom_effect = NA, n_ctrl_het = NA, n_ctrl_hom_other = NA,
                     stringsAsFactors = FALSE)
  mod <- list(...); for (n in names(mod)) base[[n]] <- mod[[n]]
  base
}

test_that("reciprocal contrasts of the same data give the same per-copy effect", {
  a <- map_extraction_row(row(comparison_as_reported = "AG+AA vs GG", effect_estimate = 1.64, ci_lower = 0.81, ci_upper = 3.48), hap)
  b <- map_extraction_row(row(comparison_as_reported = "GG vs AG+AA", effect_estimate = 0.61, ci_lower = 0.29, ci_upper = 1.24), hap)
  expect_false(is.na(a$beta))
  expect_equal(a$beta, b$beta, tolerance = 0.03)
  expect_equal(a$se_beta, b$se_beta, tolerance = 0.03)
  expect_gt(a$beta, 0)
})

test_that("a protective contrast for the other allele maps to a negative h2 effect", {
  # G carriers have higher risk than A/A, so the A (h2) allele is protective
  r <- map_extraction_row(row(effect_allele = "G", other_allele = "A", comparison_as_reported = "G/G + A/G vs A/A",
                              effect_estimate = 1.253, ci_lower = 0.267, ci_upper = 5.893), hap)
  expect_lt(r$beta, 0)
})

test_that("alleles that do not belong to the locus are handled and flagged", {
  strand <- map_extraction_row(row(effect_allele = "T", other_allele = "C"), hap)
  expect_true(is.na(strand$beta)); expect_match(strand$flag, "not an allele at this locus")
  typo <- map_extraction_row(row(other_allele = "C", comparison_as_reported = "A/A + A/C genotype validation cohort",
                                 effect_estimate = 7.3, ci_lower = 1.65, ci_upper = 32.33), hap)
  expect_false(is.na(typo$beta)); expect_match(typo$flag, "replaced by G")
  expect_true(is.na(map_extraction_row(row(rsid = NA_character_), hap)$beta))
})

test_that("a 2x2 homozygote table without an effect allele uses the European minor allele and is flagged", {
  r <- map_extraction_row(row(usable_tier = "counts_2x2", genetic_model = "recessive", effect_allele = NA_character_, other_allele = NA_character_,
                              effect_estimate = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_,
                              n_cases_hom_effect = 11, n_cases_het = 15, n_ctrl_hom_effect = 46, n_ctrl_het = 140), hap)
  expect_equal(minor_allele(hap, 1), "A")
  expect_equal(minor_allele(hap, 2), "C")
  expect_match(r$flag, "inferred as the European minor allele \\(A\\)")
  expect_equal(r$y, log((11 / 46) / (15 / 140)))
  expect_gt(r$beta, 0)
  expect_true(is.na(map_extraction_row(row(usable_tier = "counts_2x2", genetic_model = "dominant", effect_allele = NA_character_, other_allele = NA_character_,
                                           n_cases_hom_effect = 1, n_cases_het = 2, n_ctrl_hom_effect = 3, n_ctrl_het = 4), hap)$beta))
})

test_that("a 2x3 table is mapped through its additive trend", {
  ctrl <- c(4000, 10000, 6000)
  cases <- round(ctrl * 0.05 * exp(0.5 * c(2, 1, 0)))
  r <- map_extraction_row(row(usable_tier = "counts_2x3", genetic_model = "genotypic_3level", comparison_as_reported = "AA vs AG vs GG",
                              n_cases_hom_effect = cases[1], n_cases_het = cases[2], n_cases_hom_other = cases[3],
                              n_ctrl_hom_effect = ctrl[1], n_ctrl_het = ctrl[2], n_ctrl_hom_other = ctrl[3]), hap)
  expect_equal(r$beta, 0.5, tolerance = 0.03)
  expect_match(r$flag, "approximation")
})

# A hand-built set of mapped rows to exercise the selection rules
mk <- function(study, rsid = "rs2305619", beta = 0.3, se_beta = 0.2, tier = "estimate", model = "dominant", primary = "yes", id = "V1",
               cohort = "whole", donor = FALSE, hr = FALSE, non_eur = FALSE, whole = TRUE) {
  data.frame(variant_row_id = id, study_id = study, rsid = rsid, tier = tier, model = model, measure = ifelse(hr, "HR", "OR"), subgroup = "whole cohort",
             primary_flag = primary, y = beta, se = se_beta, beta = beta, se_beta = se_beta, method = "estimate", flag = "", first_author = study,
             year = 2020, cohort = cohort, donor = donor, whole_cohort = whole, non_european = non_eur, odds_ratio_scale = !hr,
             priority = model_priority(tier, model), ancestry_as_reported = "x", stringsAsFactors = FALSE)
}

test_that("study selection applies the documented rules", {
  re <- rbind(
    mk("S1", beta = 0.4), mk("S2", beta = 0.1, model = "recessive", id = "V2"),
    mk("S2", beta = 0.7, model = "dominant", id = "V3"),                                    # dominant outranks recessive
    mk("S3", beta = 0.2, tier = "counts_2x3", model = "genotypic_3level", id = "V4"),
    mk("S3", beta = 0.9, id = "V5"),                                                        # counts outrank an estimate
    mk("S4", beta = 0.5, donor = TRUE), mk("S5", beta = 0.5, hr = TRUE), mk("S6", beta = 0.5, whole = FALSE),
    mk("S035", beta = -0.4), mk("S7", beta = 0.6, non_eur = TRUE),
    mk("S8", rsid = "rs1840680", beta = 0.35, id = "V8"),                                  # only the second tag SNP
    mk("S9", rsid = "rs2305619", beta = 0.3, id = "V9"), mk("S9", rsid = "rs1840680", beta = 0.8, id = "V10")  # first tag SNP wins
  )
  p <- ptx3_study_effects(re)
  expect_setequal(p$study_id, c("S1", "S2", "S3", "S7", "S8", "S9"))
  expect_equal(p$beta[p$study_id == "S2"], 0.7)
  expect_equal(p$beta[p$study_id == "S3"], 0.2)
  expect_equal(p$rsid[p$study_id == "S8"], "rs1840680")
  expect_equal(p$beta[p$study_id == "S9"], 0.3)
  expect_false(any(c("S4", "S5", "S6", "S035") %in% p$study_id))                          # donor, HR, a stratum, Herrero-Sanchez
  expect_true("S035" %in% ptx3_study_effects(re, exclude_studies = character())$study_id)
  expect_true("S5" %in% ptx3_study_effects(re, allow_hr = TRUE)$study_id)
  expect_false("S7" %in% ptx3_study_effects(re, allow_non_european = FALSE)$study_id)
  expect_equal(ptx3_study_effects(re, genotype = "donor")$study_id, "S4")
})

test_that("two cohorts in one study are combined, not counted twice", {
  re <- rbind(mk("S45", beta = 0.0, se_beta = 0.4, cohort = "stratification", id = "V1"),
              mk("S45", beta = 1.0, se_beta = 0.5, cohort = "validation", id = "V2"),
              mk("S1", beta = 0.3))
  p <- ptx3_study_effects(re)
  expect_equal(nrow(p), 2L)
  w <- 1 / c(0.4, 0.5)^2
  s45 <- p[p$study_id == "S45", ]
  expect_equal(s45$beta, sum(w * c(0, 1)) / sum(w)); expect_equal(s45$se_beta, sqrt(1 / sum(w)))
  expect_match(s45$method, "2 cohorts combined"); expect_match(s45$flag, "cohorts disagree")
  expect_equal(nrow(ptx3_study_effects(mk("S1")[0, ])), 0L)
})

zip_path <- file.path(proj_root, "IMD_genetics_analysis_bundle (1).zip")
hap_path <- file.path(proj_root, "data", "haplotype_frequencies_phased.csv")

test_that("the real PTX3 haplotype table gives the h2 class and its frequency", {
  skip_if_not(file.exists(hap_path), "the haplotype table is not present")
  g <- ptx3_haplotypes(hap_path, "gnomAD_NFE"); e <- ptx3_haplotypes(hap_path, "1000G_EUR")
  expect_equal(sum(g$freq), 1); expect_equal(sum(e$freq), 1)
  expect_equal(sum(g$freq[g$h2 == 1]), 0.482, tolerance = 0.002)                         # the gnomAD rs2305619 A frequency
  expect_equal(sum(e$freq[e$h2 == 1]), 0.456, tolerance = 0.002)
  expect_setequal(g$haplotype[g$h2 == 1], c("A-C-A", "A-A-A", "A-A-G"))
  expect_equal(g$haplotype[g$h2 == 0], "G-A-G")
})

test_that("the PTX3 pilot on the real extraction reproduces the documented primary analysis", {
  skip_if_not(file.exists(zip_path) && file.exists(hap_path), "the evidence bundle or the haplotype table is not present")
  vr <- read_bundle_table("extraction_variants_round1.csv", zip_path); st <- read_bundle_table("extraction_studies_round1.csv", zip_path)
  h <- ptx3_haplotypes(hap_path)
  re <- ptx3_row_effects(vr, st, h)
  p <- ptx3_study_effects(re)
  expect_setequal(p$study_id, c("S026", "S029", "S045", "S056", "S077"))
  expect_gt(p$beta[p$study_id == "S029"], 0)                                              # minor-allele homozygotes at higher risk
  expect_lt(p$beta[p$study_id == "S077"], 0)                                              # G carriers at higher risk, so A protective
  expect_equal(ptx3_study_effects(re, exclude_studies = character())$beta[ptx3_study_effects(re, exclude_studies = character())$study_id == "S035"] < 0, TRUE)
  expect_false(any(c("S026", "S077") %in% ptx3_study_effects(re, allow_non_european = FALSE)$study_id))
  expect_setequal(ptx3_study_effects(re, allow_non_european = FALSE)$study_id, c("S029", "S045", "S056"))
  expect_equal(ptx3_study_effects(re, genotype = "donor", allow_hr = TRUE, exclude_studies = character())$study_id, "S080")
  pooled <- bayes_meta_grid(p$beta, p$se_beta, seed = 5518)
  expect_equal(pooled$mu[2], 0.30, tolerance = 0.06)
  # the same studies through rs3816527 agree with the h2 tag SNP within the sampling error
  cross <- merge(ptx3_study_effects(re, exclude_studies = character())[, c("study_id", "beta", "se_beta")],
                 ptx3_study_effects(re, snp_order = "rs3816527", exclude_studies = character())[, c("study_id", "beta")], by = "study_id")
  expect_lt(max(abs(cross$beta.x - cross$beta.y)), 0.25)
})
