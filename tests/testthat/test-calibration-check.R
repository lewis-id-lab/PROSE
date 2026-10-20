source(file.path(proj_root, "R", "decision_curve.R"))
source(file.path(proj_root, "R", "evidence_mapping.R"))
source(file.path(proj_root, "R", "genetic_layer.R"))
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
