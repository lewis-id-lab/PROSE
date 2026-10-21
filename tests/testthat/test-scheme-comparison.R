# The stratum and whole-cohort gains used to compare the global and the profile-preserving truth in Aim 3b (R/scheme_comparison.R)
source(file.path(proj_root, "R", "decision_curve.R"))
source(file.path(proj_root, "R", "evidence_mapping.R"))
source(file.path(proj_root, "R", "genetic_layer.R"))
source(file.path(proj_root, "R", "aim3b.R"))
source(file.path(proj_root, "R", "calibration_check.R"))
source(file.path(proj_root, "R", "scheme_comparison.R"))

lp_c2 <- profile_logit(pars)
q_c2 <- independent_profile_probs(pars$marginal)
p_c2 <- c(PTX3 = 0.48, CLEC7A = 0.08, TLR4 = 0.06, CD209 = 0.125)
grid_c2 <- dosage_grid(p_c2)
st_c2 <- neutropenia_strata(q_c2, lp_c2)
cells_c2 <- list(N = trim_cells(genetic_cells(st_c2$N$q, lp_c2, grid_c2)), NN = trim_cells(genetic_cells(st_c2$NN$q, lp_c2, grid_c2)))
share_c2 <- c(N = st_c2$N$share, NN = st_c2$NN$share)
b_c2 <- c(PTX3 = 0.4, CLEC7A = 0.2, TLR4 = 0.5, CD209 = 0.7); bs_c2 <- c(PTX3 = 0.3, CLEC7A = 0, TLR4 = 0.35, CD209 = 0.6)
th_c2 <- c(0.05, 0.2)

test_that("with no effect every gain is zero under both schemes", {
  z <- b_c2 * 0
  for (sch in c("global", "profile")) {
    g <- stratum_gains(cells_c2$NN, z, z, sch, th_c2)
    expect_equal(unname(unlist(g[c("ceil_clin", "ceil_marg", "real_clin", "real_marg")])), rep(0, 8), tolerance = 1e-9)
  }
})

test_that("under the profile scheme the clinical score is calibrated, so the comparator makes no difference", {
  g <- stratum_gains(cells_c2$NN, b_c2 * 3, bs_c2, "profile", th_c2)
  expect_equal(g$ceil_clin, g$ceil_marg, tolerance = 1e-9)
  expect_equal(g$real_clin, g$real_marg, tolerance = 1e-9)
  expect_equal(g$slope, 1, tolerance = 1e-6)
  gg <- stratum_gains(cells_c2$NN, b_c2 * 3, bs_c2, "global", th_c2)
  expect_lt(gg$slope, 1)
  expect_true(all(gg$ceil_marg <= gg$ceil_clin + 1e-12))
})

test_that("a realised score that carries the true effects is the ceiling score, under both schemes", {
  for (sch in c("global", "profile")) {
    g <- stratum_gains(cells_c2$N, b_c2 * 2, b_c2 * 2, sch, th_c2)
    expect_equal(g$real_clin, g$ceil_clin, tolerance = 1e-9)
    expect_equal(g$real_marg, g$ceil_marg, tolerance = 1e-9)
  }
})

test_that("the ceiling gain agrees with multi_score_nb on the same stratum", {
  for (sch in c("global", "profile")) {
    g <- stratum_gains(cells_c2$N, b_c2, bs_c2, sch, th_c2)
    m <- multi_score_nb(cells_c2$N, b_c2, list(baseline = b_c2 * 0, ceiling = b_c2, realised = bs_c2), th_c2, scheme = sch)
    gain <- function(sc) m$nb[m$score == sc] - m$nb[m$score == "baseline"]
    expect_equal(g$ceil_clin, gain("ceiling"), tolerance = 1e-9)
    expect_equal(g$real_clin, gain("realised"), tolerance = 1e-9)
  }
})

test_that("the whole-cohort gain is the share-weighted sum of the strata's", {
  tr <- list(N = b_c2 * 0, NN = b_c2)
  w <- whole_gains(cells_c2, share_c2, tr, bs_c2, "global", th_c2)
  nn <- stratum_gains(cells_c2$NN, b_c2, bs_c2, "global", th_c2)
  n0 <- stratum_gains(cells_c2$N, b_c2 * 0, bs_c2, "global", th_c2)
  expect_equal(w$ceil_clin, share_c2[["NN"]] * nn$ceil_clin + share_c2[["N"]] * n0$ceil_clin, tolerance = 1e-12)
  expect_equal(w$slope_NN, nn$slope)
})
