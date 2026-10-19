source(file.path(proj_root, "R", "decision_curve.R"))
source(file.path(proj_root, "R", "evidence_mapping.R"))
source(file.path(proj_root, "R", "genetic_layer.R"))

lp0 <- profile_logit(pars)
q0 <- independent_profile_probs(pars$marginal)
p4 <- c(PTX3 = 0.48, CLEC7A = 0.08, TLR4 = 0.06, CD209 = 0.125)
cells0 <- genetic_cells(q0, lp0, dosage_grid(p4))

test_that("the panel is the four core terms and the optional ones are not in it", {
  expect_equal(names(GENETIC_TERMS), c("PTX3", "CLEC7A", "TLR4", "CD209"))
})

test_that("risk-class frequencies match the documented values", {
  hap <- file.path(proj_root, "data", "haplotype_frequencies_phased.csv")
  skip_if_not(file.exists(hap))
  expect_equal(genetic_term_frequencies(path = hap), c(PTX3 = 0.482, CLEC7A = 0.077, TLR4 = 0.061, CD209 = 0.125), tolerance = 0.01)
})

test_that("the dosage grid is a Hardy-Weinberg distribution with the right margins", {
  g <- dosage_grid(p4)
  expect_equal(dim(g$G), c(81L, 4L))
  expect_equal(sum(g$prob), 1)
  expect_equal(as.numeric(crossprod(g$G, g$prob)), 2 * unname(p4))
  expect_equal(sum(g$prob[g$G[, "PTX3"] == 2]), p4[["PTX3"]]^2)
})

test_that("cells are profile by dosage with weights that sum to one and the baseline prevalence", {
  expect_equal(length(cells0$w), 128L * 81L)
  expect_equal(sum(cells0$w), 1)
  expect_equal(cells0$prevalence, sum(q0 * plogis(lp0)))
  # genotype is independent of the profile
  expect_equal(as.numeric(tapply(cells0$w, cells0$profile, sum)), as.numeric(q0))
})

test_that("recentring holds the mean risk at the baseline prevalence (D5)", {
  beta <- c(PTX3 = 0.4, CLEC7A = -0.3, TLR4 = 0.9, CD209 = 0.6)
  lp <- cells0$lp_base + drop(cells0$G %*% beta)
  a <- recentre_shift(lp, cells0$w, cells0$prevalence)
  expect_equal(sum(cells0$w * plogis(lp + a)), cells0$prevalence, tolerance = 1e-9)
  expect_equal(recentre_shift(cells0$lp_base, cells0$w, cells0$prevalence), 0, tolerance = 1e-9)
})

test_that("with no genetic effect every score equals the baseline", {
  zero <- setNames(rep(0, 4), names(p4))
  r <- genetic_nb(cells0, zero, zero, c(0.05, 0.2))
  b <- r[r$scenario == "baseline", ]
  for (s in c("ceiling", "realised")) {
    x <- r[r$scenario == s, ]
    expect_equal(x$nb, b$nb, tolerance = 1e-9); expect_equal(x$auc, b$auc, tolerance = 1e-9)
  }
})

test_that("a score with the true effects never has lower net benefit than the baseline, and a mis-specified one can", {
  set.seed(4172)
  for (i in 1:10) {
    beta <- setNames(rnorm(4, 0.3, 0.6), names(p4))
    r <- genetic_nb(cells0, beta, beta * 0, c(0.03, 0.05, 0.1, 0.2))
    nb <- function(s) r$nb[r$scenario == s]
    expect_gte(min(nb("ceiling") - nb("baseline")), -1e-10)          # the true-risk rule is optimal at every threshold
    expect_equal(nb("realised"), nb("baseline"), tolerance = 1e-9)    # score with zero weights = baseline
  }
  flipped <- genetic_nb(cells0, c(PTX3 = 0.5, CLEC7A = 0.5, TLR4 = 0.5, CD209 = 0.5), c(PTX3 = -1, CLEC7A = -1, TLR4 = -1, CD209 = -1), 0.05)
  expect_lt(flipped$nb[flipped$scenario == "realised"], flipped$nb[flipped$scenario == "baseline"])
})

test_that("a larger true effect gives a larger gain for the ceiling score", {
  gain <- function(b) { r <- genetic_nb(cells0, setNames(rep(b, 4), names(p4)), setNames(rep(b, 4), names(p4)), 0.05)
    r$nb[r$scenario == "ceiling"] - r$nb[r$scenario == "baseline"] }
  expect_lt(gain(0.2), gain(0.5)); expect_lt(gain(0.5), gain(1.0))
})

test_that("effects are drawn from the posteriors, reproducibly, and the predictive reading is more spread out", {
  draws <- lapply(setNames(1:4, names(p4)), function(j) data.frame(mu = rnorm(500, j / 10, 0.1), tau = rep(0.3, 500)))
  a <- draw_term_effects(draws, 200, "mean", seed = 11L); b <- draw_term_effects(draws, 200, "mean", seed = 11L)
  expect_equal(a, b)
  expect_equal(dim(a), c(200L, 4L)); expect_equal(colnames(a), names(p4))
  expect_equal(a[, "CD209"], draws$CD209$mu[1:200])
  pr <- draw_term_effects(draws, 200, "predictive", seed = 11L)
  expect_gt(sd(pr[, 1]), sd(a[, 1]))
  expect_equal(pooled_estimates(draws), vapply(draws, function(d) median(d$mu), numeric(1)))
  expect_error(draw_term_effects(draws, 600, "mean", seed = 1L))
})

test_that("run_genetic_layer tests the gain against the MCID and gives one row per draw, score and threshold", {
  eff <- rbind(c(PTX3 = 0.3, CLEC7A = 0, TLR4 = 0.4, CD209 = 0.6), c(PTX3 = 1, CLEC7A = 1, TLR4 = 1, CD209 = 1))
  res <- run_genetic_layer(cells0, eff, colMeans(eff), c(0.05, 0.2), mcid = c(`0.05` = 0.002, `0.2` = 0.001))
  expect_equal(nrow(res), 2 * 3 * 2)
  expect_true(all(res$delta_nb[res$scenario == "baseline"] == 0))
  expect_equal(res$meets_mcid, res$delta_nb >= res$mcid)
  expect_equal(res$mcid_truth, 0.2 * res$tp[res$scenario == "baseline"][match(paste(res$draw, res$threshold), paste(res$draw[res$scenario == "baseline"], res$threshold[res$scenario == "baseline"]))], tolerance = 1e-12)
  expect_true(any(res$meets_mcid) && any(!res$meets_mcid))
})

test_that("genotype is fixed within a patient and has the Hardy-Weinberg dosage distribution", {
  patient <- rep(paste0("P", 1:20000), times = sample(1:3, 20000, replace = TRUE))
  g <- assign_patient_genotypes(patient, p4, seed = 5311L)
  expect_equal(nrow(g), length(patient))
  expect_true(all(vapply(split(seq_along(patient), patient)[1:500], function(i) nrow(unique(g[i, , drop = FALSE])) == 1L, logical(1))))
  one_per_patient <- g[!duplicated(patient), ]
  expect_equal(unname(colMeans(one_per_patient)), 2 * unname(p4), tolerance = 0.03)
  expect_equal(assign_patient_genotypes(patient, p4, seed = 5311L), g)
})

test_that("a replicated cohort has distinct patients in every repetition", {
  co <- list(rows = data.frame(patient = c("A", "A", "B"), case = c(0, 1, 0), stringsAsFactors = FALSE), X = matrix(c(0, 1, 1), 3, 1, dimnames = list(NULL, "x")))
  rc <- replicate_cohort(co, 4)
  expect_equal(nrow(rc$X), 12L)
  expect_equal(length(unique(rc$rows$patient)), 8L)
  expect_equal(sum(rc$rows$case), 4)
})

test_that("the exact engine agrees with a patient-level simulation", {
  set.seed(7826)
  n <- 1500
  X <- matrix(rbinom(n * 7, 1, rep(pars$marginal, each = n)), n, 7, dimnames = list(NULL, PREDICTORS))
  co <- list(rows = data.frame(patient = rep(paste0("P", 1:600), length.out = n), case = 0L, stringsAsFactors = FALSE), X = X)
  rc <- replicate_cohort(co, 60)
  beta <- c(PTX3 = 1, CLEC7A = 0.5, TLR4 = 1, CD209 = 1)
  q <- tabulate(profile_index(X), 128) / n
  cells <- genetic_cells(q, lp0, dosage_grid(p4))
  ex <- genetic_nb(cells, beta, beta, c(0.05, 0.2))
  sim <- simulate_nb(rc, lp0, p4, beta, cells$prevalence, c(0.05, 0.2), seed = 9031L)
  for (k in 1:2) {
    exact_gain <- ex$nb[ex$scenario == "ceiling"][k] - ex$nb[ex$scenario == "baseline"][k]
    expect_equal(sim$ceiling[k] - sim$baseline[k], exact_gain, tolerance = 0.002, scale = 1)
  }
})

# ---- attribution, ensemble and sampling variability -------------------------------------------------------------------------

test_that("the lean rule evaluation agrees with dca_exact", {
  r <- plogis(cells0$lp_base + drop(cells0$G %*% c(0.3, 0, 0.5, 0.6)))
  rt <- plogis(cells0$lp_base)
  th <- c(0.03, 0.05, 0.2)
  d <- dca_exact(cells0$w, rt, r, thresholds = th)
  x <- tp_fp(cells0$w, rt, r, th)
  expect_equal(unname(x[, "tp"]), d$tp); expect_equal(unname(x[, "fp"]), d$fp)
  expect_equal(nb_from(x, th), d$nb)
})

test_that("Shapley contributions sum to the gain of the full score, for the ceiling and the realised score", {
  beta <- c(PTX3 = 0.3, CLEC7A = -0.1, TLR4 = 0.5, CD209 = 0.7); bhat <- c(PTX3 = 0.2, CLEC7A = 0, TLR4 = 0.3, CD209 = 0.4)
  th <- c(0.05, 0.2)
  full <- genetic_nb(cells0, beta, bhat, th)
  gain <- function(s) full$nb[full$scenario == s] - full$nb[full$scenario == "baseline"]
  ceil <- shapley_gain(cells0, beta, beta, th)
  real <- shapley_gain(cells0, beta, bhat, th)
  expect_equal(unname(rowSums(ceil$phi)), ceil$total); expect_equal(ceil$total, gain("ceiling"), ignore_attr = TRUE)
  expect_equal(unname(rowSums(real$phi)), real$total); expect_equal(real$total, gain("realised"), ignore_attr = TRUE)
  expect_equal(dim(ceil$phi), c(2L, 4L)); expect_equal(colnames(ceil$alone), names(p4))
})

test_that("a term with no effect gets no credit, and a lone term gets all of it", {
  th <- 0.05
  lone <- shapley_gain(cells0, c(PTX3 = 0, CLEC7A = 0, TLR4 = 0, CD209 = 0.8), c(PTX3 = 0, CLEC7A = 0, TLR4 = 0, CD209 = 0.8), th)
  expect_equal(unname(lone$phi[1, c("PTX3", "CLEC7A", "TLR4")]), c(0, 0, 0), tolerance = 1e-9)
  expect_equal(unname(lone$phi[1, "CD209"]), lone$total, tolerance = 1e-9)
  expect_equal(unname(lone$alone[1, "CD209"]), lone$total, tolerance = 1e-9)
  expect_equal(unname(lone$leave_out[1, "CD209"]), lone$total, tolerance = 1e-9)
})

test_that("two terms with the same frequency and effect get the same contribution", {
  pp <- c(PTX3 = 0.3, CLEC7A = 0.3, TLR4 = 0.1, CD209 = 0.1)
  cc <- genetic_cells(q0, lp0, dosage_grid(pp))
  s <- shapley_gain(cc, c(PTX3 = 0.6, CLEC7A = 0.6, TLR4 = 0.3, CD209 = 0.3), c(PTX3 = 0.6, CLEC7A = 0.6, TLR4 = 0.3, CD209 = 0.3), c(0.05, 0.2))
  expect_equal(s$phi[, "PTX3"], s$phi[, "CLEC7A"], tolerance = 1e-9)
  expect_equal(s$phi[, "TLR4"], s$phi[, "CD209"], tolerance = 1e-9)
})

test_that("an ensemble of identical members reproduces the single-distribution layer, with the member's own MCID", {
  eff <- rbind(c(PTX3 = 0.3, CLEC7A = 0, TLR4 = 0.4, CD209 = 0.6), c(PTX3 = 0.8, CLEC7A = 0.2, TLR4 = 0.1, CD209 = 0.2), c(PTX3 = 0, CLEC7A = 0, TLR4 = 0, CD209 = 0))
  bh <- colMeans(eff); th <- c(0.05, 0.2)
  qm <- rbind(q0, q0)
  ens <- run_genetic_layer_ensemble(qm, c(1, 2, 1), lp0, dosage_grid(p4), eff, bh, th)
  single <- run_genetic_layer(cells0, eff, bh, th, mcid = c(`0.05` = 0.002, `0.2` = 0.001))
  expect_equal(ens$delta_nb, single$delta_nb, tolerance = 1e-12)
  expect_equal(unique(ens$mcid[ens$threshold == 0.05]), baseline_mcid(q0, lp0, 0.05)[1], tolerance = 1e-12)
  expect_equal(ens$meets_mcid, ens$delta_nb >= ens$mcid)
  expect_setequal(ens$member, c(1, 2)); expect_equal(nrow(ens), 3 * 3 * 2)
  # a different member changes the answer
  q_alt <- independent_profile_probs(pars$marginal * 1.5)
  q_alt <- q_alt / sum(q_alt)
  alt <- run_genetic_layer_ensemble(rbind(q0, q_alt), c(1, 2, 1), lp0, dosage_grid(p4), eff, bh, th)
  expect_gt(max(abs(alt$delta_nb[alt$draw == 2] - ens$delta_nb[ens$draw == 2])), 1e-5)
})

test_that("a simulated study has the cohort's size, is reproducible, and averages to the exact gain", {
  set.seed(6105)
  n <- 1500
  X <- matrix(rbinom(n * 7, 1, rep(pars$marginal, each = n)), n, 7, dimnames = list(NULL, PREDICTORS))
  co <- list(rows = data.frame(patient = rep(paste0("P", 1:700), length.out = n), case = 0L, stringsAsFactors = FALSE), X = X)
  cg <- cohort_groups(co)
  expect_equal(length(cg$groups), 700L)
  beta <- c(PTX3 = 1, CLEC7A = 0.5, TLR4 = 1, CD209 = 1); th <- c(0.05, 0.2)
  cells <- genetic_cells(tabulate(profile_index(X), 128) / n, lp0, dosage_grid(p4))
  sh <- genetic_scores(cells, beta, beta)
  mc <- c(`0.05` = 0.004, `0.2` = 0.001)
  set.seed(1); a <- simulate_study(cg, lp0, p4, sh, beta, beta, th, mc, resample = FALSE)
  set.seed(1); b <- simulate_study(cg, lp0, p4, sh, beta, beta, th, mc, resample = FALSE)
  expect_equal(a, b)
  expect_equal(nrow(a), 3 * 2); expect_true(all(a$delta_nb[a$scenario == "baseline"] == 0))
  expect_equal(a$mcid_sample[1], 0.2 * a$tp[a$scenario == "baseline" & a$threshold == 0.05])
  expect_equal(a$mcid[a$threshold == 0.2][1], 0.001)
  ex <- genetic_nb(cells, beta, beta, th)
  exact_gain <- ex$nb[ex$scenario == "ceiling"] - ex$nb[ex$scenario == "baseline"]
  sims <- do.call(rbind, lapply(1:300, function(i) { s <- simulate_study(cg, lp0, p4, sh, beta, beta, th, mc, resample = FALSE); s[s$scenario == "ceiling", ] }))
  expect_equal(as.numeric(tapply(sims$delta_nb, sims$threshold, mean)), exact_gain, tolerance = 0.0015, scale = 1)
  # resampling patients changes the clinical mix, so it adds variance
  rs <- do.call(rbind, lapply(1:300, function(i) { s <- simulate_study(cg, lp0, p4, sh, beta, beta, th, mc, resample = TRUE); s[s$scenario == "ceiling" & s$threshold == 0.05, ] }))
  expect_gt(sd(rs$delta_nb), 0)
})

# ---- allo-HSCT stratum, donor genotype, and studies drawn from a profile distribution -----------------------------------------------

test_that("cells carry a true-risk offset that the scores do not", {
  expect_equal(cells0$pred_mean, cells0$prevalence)                     # no offset: truth and baseline score agree
  ca <- genetic_cells(q0, lp0, dosage_grid(p4), offset = 1.5)
  expect_gt(ca$prevalence, ca$pred_mean)
  expect_equal(ca$pred_mean, cells0$pred_mean)
  zero <- setNames(rep(0, 4), names(p4))
  r <- genetic_nb(ca, zero, zero, c(0.05, 0.2))
  b <- r[r$scenario == "baseline", ]
  for (s in c("ceiling", "realised")) expect_equal(r$nb[r$scenario == s], b$nb, tolerance = 1e-9)   # no genetic information, no gain
  sc <- genetic_scores(ca, c(PTX3 = 0.4, CLEC7A = 0, TLR4 = 0.5, CD209 = 0.6), zero)
  expect_equal(sum(ca$w * plogis(sc$lp_true)), ca$prevalence, tolerance = 1e-9)        # the truth keeps its event rate
  expect_equal(sum(ca$w * plogis(sc$lp_ceiling)), ca$pred_mean, tolerance = 1e-9)      # the score keeps the baseline's mean risk
})

test_that("the baseline MCID accounts for a true-risk offset", {
  th <- c(0.05, 0.2)
  expect_equal(baseline_mcid(q0, lp0, th), 0.2 * dca_exact(q0, plogis(lp0), thresholds = th)$tp)
  d <- dca_exact(q0, plogis(lp0 + 1.2), plogis(lp0), thresholds = th)
  expect_equal(baseline_mcid(q0, lp0, th, offset = 1.2), d$mcid_b)
  expect_gt(baseline_mcid(q0, lp0, th, offset = 1.2)[1], baseline_mcid(q0, lp0, th)[1])
})

test_that("the allo stratum reproduces the Table 1 margins it is built from", {
  st <- allo_stratum(q0, lp0)
  r <- plogis(lp0)
  expect_equal(sum(st$q), 1)
  expect_equal(st$offset, log((44 / 133) / (248 / 3994)))
  expect_equal(sum(st$share * st$q * plogis(lp0 + st$offset)) / sum(q0 * r), 44 / 133)                  # P(allo | case)
  expect_equal(sum(st$share * st$q * (1 - plogis(lp0 + st$offset))) / sum(q0 * (1 - r)), 248 / 3994)    # P(allo | no case)
  expect_equal(st$event_rate, sum(st$q * plogis(lp0 + st$offset)))
  # the stratum leans to higher-risk profiles than the cohort
  expect_gt(sum(st$q * r), sum(q0 * r))
})

test_that("donor draws load where they exist and a term without donor evidence is held at zero", {
  d <- tempfile("donor_"); dir.create(d); on.exit(unlink(d, recursive = TRUE))
  write.csv(data.frame(mu = rnorm(300, 0.2, 0.1), tau = 0.2), file.path(d, "ptx3_h2_donor_posterior_draws.csv"), row.names = FALSE)
  write.csv(data.frame(mu = rnorm(300, -0.3, 0.1), tau = 0.2), file.path(d, "clec7a_donor_posterior_draws.csv"), row.names = FALSE)
  dd <- load_donor_draws(dir = d)
  expect_equal(names(dd), names(GENETIC_TERMS))
  expect_null(dd$TLR4); expect_null(dd$CD209)
  e <- draw_donor_effects(dd, 100, "mean", seed = 3L)
  expect_equal(dim(e), c(100L, 4L)); expect_true(all(e[, c("TLR4", "CD209")] == 0))
  expect_equal(e[, "PTX3"], dd$PTX3$mu[1:100])
  expect_equal(draw_donor_effects(dd, 100, "predictive", seed = 3L), draw_donor_effects(dd, 100, "predictive", seed = 3L))
  est <- donor_pooled_estimates(dd)
  expect_equal(est[["TLR4"]], 0); expect_equal(est[["PTX3"]], median(dd$PTX3$mu))
  expect_error(draw_donor_effects(dd, 600, "mean", seed = 1L))
})

test_that("a study drawn from a profile distribution is reproducible and averages to the exact gain, with and without an offset", {
  beta <- c(PTX3 = 1, CLEC7A = 0.5, TLR4 = 1, CD209 = 1); th <- c(0.05, 0.2)
  counts <- c(rep(1L, 600), rep(2L, 300), rep(4L, 40))
  for (off in c(0, 1.5)) {
    cells <- genetic_cells(q0, lp0, dosage_grid(p4), offset = off)
    sh <- genetic_scores(cells, beta, beta)
    mc <- setNames(baseline_mcid(q0, lp0, th, offset = off), as.character(th))
    set.seed(11); a <- simulate_study_q(counts, q0, lp0, p4, sh, beta, beta, th, mc, offset = off)
    set.seed(11); b <- simulate_study_q(counts, q0, lp0, p4, sh, beta, beta, th, mc, offset = off)
    expect_equal(a, b)
    expect_equal(nrow(a), 6L); expect_true(all(a$delta_nb[a$scenario == "baseline"] == 0))
    expect_equal(a$mcid[a$threshold == 0.05][1], unname(mc[1]))
    set.seed(12)
    sims <- do.call(rbind, lapply(1:400, function(i) { s <- simulate_study_q(counts, q0, lp0, p4, sh, beta, beta, th, mc, offset = off); s[s$scenario == "ceiling", ] }))
    ex <- genetic_nb(cells, beta, beta, th)
    exact_gain <- ex$nb[ex$scenario == "ceiling"] - ex$nb[ex$scenario == "baseline"]
    se <- tapply(sims$delta_nb, sims$threshold, sd) / sqrt(400)
    expect_lt(max(abs(as.numeric(tapply(sims$delta_nb, sims$threshold, mean)) - exact_gain) / as.numeric(se)), 4)
  }
  expect_lt(abs(unique(a$admissions) - sum(counts)), 0.2 * sum(counts))   # resampled patients change the total a little
  fixed <- simulate_study_q(counts, q0, lp0, p4, sh, beta, beta, th, mc, offset = 1.5, resample_counts = FALSE)
  expect_equal(unique(fixed$admissions), sum(counts))
})

test_that("the data-based allo stratum is the flagged admissions' own distribution, recalibrated to their observed event rate", {
  set.seed(8821)
  n <- 3000
  X <- matrix(rbinom(n * 7, 1, rep(pars$marginal, each = n)), n, 7, dimnames = list(NULL, PREDICTORS))
  allo <- rbinom(n, 1, 0.1) == 1
  case <- rbinom(n, 1, ifelse(allo, 0.2, 0.03))
  co <- list(rows = data.frame(case = case, allo_hsct = allo), X = X)
  st <- allo_stratum_data(co, lp0)
  expect_equal(sum(st$q), 1)
  expect_equal(st$n, sum(allo)); expect_equal(st$events, sum(case[allo])); expect_equal(st$share, mean(allo))
  expect_equal(st$q, tabulate(profile_index(X[allo, ]), 128) / sum(allo))
  expect_equal(sum(st$q * plogis(lp0 + st$offset)), sum(case[allo]) / sum(allo), tolerance = 1e-9)
  expect_equal(st$predicted_rate_baseline, sum(st$q * plogis(lp0)))
  expect_gt(st$offset, 0)                                                # the baseline under-predicts in a stratum that is really riskier
  co$rows$allo_hsct <- NA
  expect_error(allo_stratum_data(co, lp0), "no allo-HSCT flag")
  # a stratum built this way works with the cells and the MCID like the others
  cells <- genetic_cells(st$q, lp0, dosage_grid(p4), offset = st$offset)
  expect_equal(cells$prevalence, sum(allo) / sum(allo) * sum(case[allo]) / sum(allo), tolerance = 1e-9)
  expect_gt(baseline_mcid(st$q, lp0, 0.05, offset = st$offset), 0)
})

test_that("a gene with several haplotype classes has multinomial genotypes, and the joint grid is a product", {
  hap <- data.frame(haplotype = c("C-G", "A-G", "A-C"), freq = c(0.08, 0.35, 0.57), stringsAsFactors = FALSE)
  part <- haplotype_gene_part(hap, c(CLEC7A = "C-G", CLEC7A_AG = "A-G"))
  expect_equal(nrow(part$G), 6L); expect_equal(sum(part$prob), 1)
  expect_equal(sum(part$G[, "CLEC7A"] * part$prob), 2 * 0.08); expect_equal(sum(part$G[, "CLEC7A_AG"] * part$prob), 2 * 0.35)
  expect_true(all(rowSums(part$G) <= 2))
  expect_equal(sum(part$prob[part$G[, "CLEC7A"] == 2]), 0.08^2)
  expect_equal(sum(part$prob[part$G[, "CLEC7A"] == 1 & part$G[, "CLEC7A_AG"] == 1]), 2 * 0.08 * 0.35)
  jg <- joint_grid(list(dosage_grid(c(PTX3 = 0.48)), part, dosage_grid(c(TLR4 = 0.06))))
  expect_equal(dim(jg$G), c(3L * 6L * 3L, 4L)); expect_equal(colnames(jg$G), c("PTX3", "CLEC7A", "CLEC7A_AG", "TLR4"))
  expect_equal(sum(jg$prob), 1)
  expect_equal(sum(jg$prob[jg$G[, "TLR4"] == 2]), 0.06^2)
  # with the second parameter at zero the joint grid gives the same layer as independent binomials for the other terms
  q <- q0; p3 <- c(PTX3 = 0.48, CLEC7A = 0.08, TLR4 = 0.06)
  indep <- genetic_nb(genetic_cells(q, lp0, dosage_grid(p3)), c(PTX3 = 0.3, CLEC7A = 0.1, TLR4 = 0.4), c(PTX3 = 0.3, CLEC7A = 0.1, TLR4 = 0.4), c(0.05, 0.2))
  jg2 <- joint_grid(list(dosage_grid(c(PTX3 = 0.48)), haplotype_gene_part(hap[c(1, 2, 3), ] |> transform(freq = c(0.08, 0.35, 0.57)), c(CLEC7A = "C-G", CLEC7A_AG = "A-G")), dosage_grid(c(TLR4 = 0.06))))
  b <- c(PTX3 = 0.3, CLEC7A = 0.1, CLEC7A_AG = 0, TLR4 = 0.4)
  joint <- genetic_nb(genetic_cells(q, lp0, jg2), b, b, c(0.05, 0.2))
  expect_equal(joint$nb[joint$scenario == "ceiling"], indep$nb[indep$scenario == "ceiling"], tolerance = 1e-10)
})

test_that("the smoothed allo stratum keeps the observed share and event rate and spreads over more profiles than the data", {
  set.seed(3394)
  n <- 3000
  X <- matrix(rbinom(n * 7, 1, rep(pars$marginal, each = n)), n, 7, dimnames = list(NULL, PREDICTORS))
  allo <- rbinom(n, 1, plogis(-3.5 + 3 * X[, "high_risk_chemo"] + 1.2 * X[, "severe_lymphopenia"] + 0.6 * X[, "prolonged_neutropenia"])) == 1
  case <- rbinom(n, 1, ifelse(allo, 0.2, 0.03))
  co <- list(rows = data.frame(case = case, allo_hsct = allo), X = X)
  sm <- allo_stratum_smooth(co, lp0); dt <- allo_stratum_data(co, lp0)
  expect_equal(sum(sm$q), 1); expect_equal(sm$share, mean(allo), tolerance = 1e-9)
  expect_equal(sum(sm$q * plogis(lp0 + sm$offset)), sum(case[allo]) / sum(allo), tolerance = 1e-9)
  expect_gte(sm$profiles_present, dt$profiles_present)
  expect_lt(0.5 * sum(abs(sm$q - dt$q)), 0.25)              # smoothing, not a different distribution
  expect_length(sm$coefficients, 8L)
  co$rows$allo_hsct <- NA
  expect_error(allo_stratum_smooth(co, lp0), "no allo-HSCT flag")
})
