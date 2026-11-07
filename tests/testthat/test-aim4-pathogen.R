source(file.path(proj_root, "R", "decision_curve.R"))
source(file.path(proj_root, "R", "evidence_mapping.R"))
source(file.path(proj_root, "R", "genetic_layer.R"))
source(file.path(proj_root, "R", "aim4_sweep.R"))
source(file.path(proj_root, "R", "calibration_check.R"))
source(file.path(proj_root, "R", "aim4_pathogen.R"))

lp_t <- profile_logit(pars)
q_t <- independent_profile_probs(pars$marginal)
at <- c(0.05, 0.20)

test_that("with every case genotype-sensitive (f = 1) the scenario is the exact sweep, and the fitted weight is the true one", {
  a <- sweep_cell(q_t, lp_t, 8, 1.6, 0.3, at)
  b <- pathogen_cell(q_t, lp_t, 8, 1.6, 0.3, 1, at)
  expect_equal(b$gain_specific, a$gain, tolerance = 1e-10)
  expect_equal(b$gain_diluted, a$gain, tolerance = 1e-8)
  expect_equal(b$or_diluted, rep(1.6, 2), tolerance = 1e-6)
})

test_that("dilution lowers the fitted odds ratio and the gain, and the event rate is held", {
  b <- lapply(c(1, 0.8, 0.6, 0.4), function(f) pathogen_cell(q_t, lp_t, 10, 2, 0.3, f, at))
  or <- vapply(b, function(x) x$or_diluted[1], numeric(1)); g <- vapply(b, function(x) x$gain_diluted[1], numeric(1))
  expect_true(all(diff(or) < 0)); expect_true(all(diff(g) < 0)); expect_gt(or[4], 1)
  pmf <- dosage_pmf(10, 0.3)
  cells <- genetic_cells(q_t, lp_t, list(G = matrix(seq_along(pmf) - 1, ncol = 1, dimnames = list(NULL, "D")), prob = pmf))
  r <- pathogen_risk(cells$lp_base, cells$G[, 1] * log(2), cells$w, cells$prevalence, 0.6)
  expect_equal(sum(cells$w * r), cells$prevalence, tolerance = 1e-10)
})

test_that("the realised version without estimation error reproduces the exact gain, and noise lowers it", {
  ex <- pathogen_cell(q_t, lp_t, 10, 2, 0.3, 0.8, at)
  r0 <- pathogen_cell_realised(q_t, lp_t, 10, 2, 0.3, 0.8, 0, at, R = 8L, seed = 7L)
  m0 <- aggregate(gain ~ score + threshold, r0, mean)
  expect_equal(m0$gain[m0$score == "specific" & m0$threshold == 0.05], ex$gain_specific[1], tolerance = 0.15)
  expect_equal(m0$gain[m0$score == "all_mould" & m0$threshold == 0.05], ex$gain_diluted[1], tolerance = 0.15)
  r2 <- pathogen_cell_realised(q_t, lp_t, 10, 2, 0.3, 0.8, 0.5, at, R = 8L, seed = 7L)
  m2 <- aggregate(gain ~ score + threshold, r2, mean)
  expect_lt(m2$gain[m2$score == "specific" & m2$threshold == 0.05], m0$gain[m0$score == "specific" & m0$threshold == 0.05])
  expect_equal(nrow(r2), 8L * 2L * 2L)
})

test_that("under the profile-preserving truth f = 1 reproduces the profile comparison of script 31, and the gain is lower than under the global truth", {
  pmf <- dosage_pmf(8, 0.3)
  grid <- list(G = matrix(seq_along(pmf) - 1, ncol = 1, dimnames = list(NULL, "D")), prob = pmf)
  ref <- calibration_gains(q_t, lp_t, grid, c(D = log(1.6)), thresholds = at, comparisons = "profile")
  b <- pathogen_cell(q_t, lp_t, 8, 1.6, 0.3, 1, at, scheme = "profile")
  expect_equal(b$gain_specific, ref$gain_oracle, tolerance = 1e-10)
  g <- pathogen_cell(q_t, lp_t, 8, 1.6, 0.3, 1, at, scheme = "global")
  expect_lt(b$gain_specific[1], g$gain_specific[1])
})

test_that("the realised version under the profile truth reproduces the exact gain without error, and a dilution still lowers the gain", {
  ex <- pathogen_cell(q_t, lp_t, 10, 2, 0.3, 0.8, at, scheme = "profile")
  r0 <- pathogen_cell_realised(q_t, lp_t, 10, 2, 0.3, 0.8, 0, at, R = 8L, seed = 7L, scheme = "profile")
  m0 <- aggregate(gain ~ score + threshold, r0, mean)
  expect_equal(m0$gain[m0$score == "specific" & m0$threshold == 0.05], ex$gain_specific[1], tolerance = 0.15)
  expect_equal(m0$gain[m0$score == "all_mould" & m0$threshold == 0.05], ex$gain_diluted[1], tolerance = 0.15)
  g <- vapply(c(1, 0.8, 0.6), function(f) pathogen_cell(q_t, lp_t, 10, 2, 0.3, f, at, scheme = "profile")$gain_diluted[1], numeric(1))
  expect_true(all(diff(g) < 0))
})

test_that("the scheme follows PROSE_SCHEME by default", {
  old <- Sys.getenv("PROSE_SCHEME", NA); on.exit(if (is.na(old)) Sys.unsetenv("PROSE_SCHEME") else Sys.setenv(PROSE_SCHEME = old))
  Sys.setenv(PROSE_SCHEME = "profile"); a <- pathogen_cell(q_t, lp_t, 6, 1.8, 0.3, 0.8, at)
  Sys.setenv(PROSE_SCHEME = "global"); b <- pathogen_cell(q_t, lp_t, 6, 1.8, 0.3, 0.8, at)
  expect_equal(a, pathogen_cell(q_t, lp_t, 6, 1.8, 0.3, 0.8, at, scheme = "profile")); expect_false(isTRUE(all.equal(a$gain_specific, b$gain_specific)))
})

test_that("the stress version reproduces the exact gain without error, and each departure lowers the realised gain", {
  old <- Sys.getenv("PROSE_SCHEME", NA); on.exit(if (is.na(old)) Sys.unsetenv("PROSE_SCHEME") else Sys.setenv(PROSE_SCHEME = old))
  Sys.setenv(PROSE_SCHEME = "profile")
  ex <- pathogen_cell(q_t, lp_t, 8, 2.5, 0.3, 0.8, at)
  r0 <- pathogen_cell_stress(q_t, lp_t, 8, 2.5, 0.3, 0.8, 0, at, R = 8L, seed = 11L)
  m0 <- aggregate(gain ~ score + threshold, r0, mean)
  expect_equal(m0$gain[m0$score == "specific" & m0$threshold == 0.05], ex$gain_specific[1], tolerance = 0.15)
  expect_true(all(r0$true_selected == 8))
  g <- function(...) { x <- pathogen_cell_stress(q_t, lp_t, 8, 2.5, 0.3, 0.8, 0.3, at, R = 12L, seed = 3L, ...); median(x$gain[x$score == "specific" & x$threshold == 0.05]) }
  base <- g()
  expect_lt(g(tag_r = 0.5), base)
  cu <- pathogen_cell_stress(q_t, lp_t, 8, 2.5, 0.3, 0.8, 0.3, at, R = 12L, seed = 3L, winners_curse = TRUE)
  expect_lt(mean(cu$true_selected), 8)
  Sys.setenv(PROSE_SCHEME = "global")
  expect_error(pathogen_cell_stress(q_t, lp_t, 8, 2.5, 0.3, 0.8, 0.3, at), "profile")
})
