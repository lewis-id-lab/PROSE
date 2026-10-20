source(file.path(proj_root, "R", "decision_curve.R"))
source(file.path(proj_root, "R", "evidence_mapping.R"))
source(file.path(proj_root, "R", "genetic_layer.R"))
source(file.path(proj_root, "R", "aim4_sweep.R"))
source(file.path(proj_root, "R", "nri.R"))
source(file.path(proj_root, "R", "calibration_check.R"))

lp_c <- profile_logit(pars)
q_c <- independent_profile_probs(pars$marginal)
p_c <- c(PTX3 = 0.48, CLEC7A = 0.08, TLR4 = 0.06, CD209 = 0.125)
grid_c <- dosage_grid(p_c)
b_true <- c(PTX3 = 0.4, CLEC7A = 0.2, TLR4 = 0.5, CD209 = 0.7)
b_score <- c(PTX3 = 0.3, CLEC7A = 0, TLR4 = 0.35, CD209 = 0.6)

test_that("profile shifts hold each profile's clinical risk", {
  gl <- drop(grid_c$G %*% b_true)
  R <- plogis(outer(lp_c + profile_shifts(lp_c, gl, grid_c$prob), gl, "+"))
  expect_equal(drop(R %*% grid_c$prob), plogis(lp_c), tolerance = 1e-10)
})

test_that("the current comparison reproduces the genetic-layer engine", {
  cells <- genetic_cells(q_c, lp_c, grid_c)
  r <- genetic_nb(cells, b_true, b_score, c(0.05, 0.2))
  ref_ceil <- r$nb[r$scenario == "ceiling"] - r$nb[r$scenario == "baseline"]
  ref_real <- r$nb[r$scenario == "realised"] - r$nb[r$scenario == "baseline"]
  g <- calibration_gains(q_c, lp_c, grid_c, b_true, b_score, c(0.05, 0.2), comparisons = "current")
  expect_equal(g$gain_oracle, ref_ceil, tolerance = 1e-10)
  expect_equal(g$gain_score, ref_real, tolerance = 1e-10)
})

test_that("with no genetic effect every comparison gives zero gain", {
  z <- setNames(rep(0, 4), names(p_c))
  g <- calibration_gains(q_c, lp_c, grid_c, z, z, c(0.05, 0.2))
  expect_equal(g$gain_oracle, rep(0, 6), tolerance = 1e-12)
  expect_equal(g$gain_score, rep(0, 6), tolerance = 1e-12)
})

test_that("the marginal clinical risk is the best clinical-only rule under the global truth", {
  g <- calibration_gains(q_c, lp_c, grid_c, b_true * 3, b_score, c(0.05, 0.1, 0.2), comparisons = c("current", "marginal"))
  cur <- g[g$comparison == "current", ]; mar <- g[g$comparison == "marginal", ]
  expect_true(all(mar$nb_clinical >= cur$nb_clinical - 1e-12))
})

test_that("calibration diagnostics tolerate risks that round to exactly 0 or 1", {
  x <- c(0, 0.2, 0.5, 1); y <- c(0.01, 0.2, 0.6, 0.99)
  f <- calibration_fit(y, x, rep(0.25, 4))
  expect_true(all(is.finite(f)))
})

test_that("calibration diagnostics: the clinical score is exactly calibrated under the profile-preserving truth", {
  g <- calibration_gains(q_c, lp_c, grid_c, b_true * 3, b_score, 0.05, calibration = TRUE)
  pr <- g[g$comparison == "profile", ]
  expect_equal(pr$clin_slope, 1, tolerance = 1e-6)
  expect_equal(pr$clin_citl, 0, tolerance = 1e-6)
  expect_lt(pr$clin_max_abs_err, 1e-9)
  mar <- g[g$comparison == "marginal", ]
  expect_equal(mar$clin_slope, 1, tolerance = 1e-6)
  expect_gt(g$clin_max_abs_err[g$comparison == "current"], 0)
})

test_that("the realised sweep reproduces sweep_cell_realised for the current comparison", {
  th <- c(0.05, 0.2)
  ref <- sweep_cell_realised(q_c, lp_c, 3, 1.6, 0.3, 0.15, th, R = 3L, N = 300L, seed = 11L)
  g <- calibration_gains_realised(q_c, lp_c, 3, 1.6, 0.3, 0.15, th, R = 3L, N = 300L, seed = 11L)
  cur <- g[g$comparison == "current", ]
  expect_equal(cur$gain, ref$gain, tolerance = 1e-12)
  expect_equal(cur$replicate, ref$replicate)
})

test_that("the realised sweep: comparisons share draws, and no genetic effect gives zero gain", {
  th <- c(0.05, 0.2)
  g <- calibration_gains_realised(q_c, lp_c, 4, 1, 0.3, 0, th, R = 2L, N = 200L, seed = 5L)
  expect_equal(g$gain, rep(0, nrow(g)), tolerance = 1e-9)
  # a pure-noise score (no effect, imprecise weights) cannot help under any comparison
  gn <- calibration_gains_realised(q_c, lp_c, 4, 1, 0.3, 0.2, th, R = 2L, N = 200L, seed = 5L)
  expect_true(all(gn$gain <= 1e-9))
  g2 <- calibration_gains_realised(q_c, lp_c, 4, 1.5, 0.3, 0.2, th, R = 2L, N = 200L, seed = 5L)
  expect_equal(sort(unique(g2$comparison)), c("current", "marginal", "profile"))
  g3 <- calibration_gains_realised(q_c, lp_c, 4, 1.5, 0.3, 0.2, th, comparisons = "profile", R = 2L, N = 200L, seed = 5L)
  expect_equal(g3$gain, g2$gain[g2$comparison == "profile"], tolerance = 1e-12)
})

test_that("per-profile shifts on a simulated sample hold each profile's clinical risk", {
  set.seed(8)
  S <- rbinom(150, 8, 0.3) * 0.4 + rnorm(150, 0, 0.1)
  a <- profile_shifts(lp_c, S, rep(1 / 150, 150))
  r <- plogis(outer(lp_c + a, S, "+"))
  expect_equal(rowMeans(r), plogis(lp_c), tolerance = 1e-9)
})

test_that("saturated profiles do not break the profile shifts", {
  a <- profile_shifts(c(-3, 2), c(0, 45, 90), rep(1 / 3, 3))
  expect_true(all(is.finite(a)))
})

test_that("NRI by scheme: current reproduces the genetic-layer engine and a null effect gives no reclassification", {
  cells <- genetic_cells(q_c, lp_c, grid_c)
  s <- genetic_scores(cells, b_true, b_score)
  ref <- reclassification_metrics(cells$w, s$lp_true, cells$lp_base, s$lp_real)
  got <- nri_by_scheme(q_c, lp_c, grid_c, b_true, b_score)
  expect_equal(got["current", ], ref, tolerance = 1e-8)
  z <- setNames(rep(0, 4), names(p_c))
  nul <- nri_by_scheme(q_c, lp_c, grid_c, z, z)
  expect_equal(unname(nul[, c("nri_cat5", "nri_cont", "dnb_05", "dnb_20")]), matrix(0, 3, 4), tolerance = 1e-9)
})

test_that("NRI by scheme: under the global truth the gain over the marginal risk cannot exceed the gain over the clinical score", {
  got <- nri_by_scheme(q_c, lp_c, grid_c, b_true * 3, b_true * 3)
  expect_true(all(is.finite(got)))
  expect_gte(got["current", "dnb_05"], got["marginal", "dnb_05"] - 1e-12)
})

test_that("detail columns: the current comparison reproduces the engine's true positives, flagged share and aROC", {
  th <- c(0.05, 0.2)
  cells <- genetic_cells(q_c, lp_c, grid_c)
  r <- genetic_nb(cells, b_true, b_score, th)
  g <- calibration_gains(q_c, lp_c, grid_c, b_true, b_score, th, comparisons = "current", detail = TRUE)
  base <- r[r$scenario == "baseline", ]; ceil <- r[r$scenario == "ceiling", ]; real <- r[r$scenario == "realised", ]
  expect_equal(g$extra_tp_oracle, ceil$tp - base$tp, tolerance = 1e-10)
  expect_equal(g$extra_tp_score, real$tp - base$tp, tolerance = 1e-10)
  expect_equal(g$flagged_oracle, ceil$flagged, tolerance = 1e-10)
  expect_equal(g$auc_clinical, base$auc, tolerance = 1e-8)
  expect_equal(g$auc_oracle, ceil$auc, tolerance = 1e-8)
  expect_equal(g$auc_score, real$auc, tolerance = 1e-8)
})
