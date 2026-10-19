source(file.path(proj_root, "R", "decision_curve.R"))
source(file.path(proj_root, "R", "evidence_mapping.R"))
source(file.path(proj_root, "R", "genetic_layer.R"))
source(file.path(proj_root, "R", "aim4_sweep.R"))

lp0 <- profile_logit(pars)
q0 <- independent_profile_probs(pars$marginal)

test_that("the exact convolution of identical loci is the binomial, and of different loci matches enumeration", {
  expect_equal(sum(locus_pmf(0.3)), 1)
  expect_equal(convolve_pmf(locus_pmf(0.2), locus_pmf(0.2)), dbinom(0:4, 4, 0.2))
  for (k in c(1, 5, 20)) expect_equal(dosage_pmf(k, 0.35), dbinom(0:(2 * k), 2 * k, 0.35))
  ps <- c(0.1, 0.4, 0.25)
  g <- as.matrix(expand.grid(0:2, 0:2, 0:2))
  pr <- apply(g, 1, function(x) prod(dbinom(x, 2, ps)))
  brute <- as.numeric(tapply(pr, rowSums(g), sum))
  expect_equal(dosage_pmf(3, ps), brute)
  expect_equal(sum(dosage_pmf(20, 0.07)), 1)
})

test_that("a sweep cell equals the explicit multi-locus calculation, and has no gain without an effect", {
  th <- c(0.05, 0.2)
  x <- sweep_cell(q0, lp0, 3, 1.6, 0.3, th)
  p3 <- c(T1 = 0.3, T2 = 0.3, T3 = 0.3); b3 <- setNames(rep(log(1.6), 3), names(p3))
  ex <- genetic_nb(genetic_cells(q0, lp0, dosage_grid(p3)), b3, b3, th)
  expect_equal(x$gain, ex$nb[ex$scenario == "ceiling"] - ex$nb[ex$scenario == "baseline"], tolerance = 1e-10)
  expect_equal(x$baseline_nb, ex$nb[ex$scenario == "baseline"], tolerance = 1e-10)
  expect_equal(x$auc, ex$auc[ex$scenario == "ceiling"], tolerance = 1e-10)
  expect_equal(sweep_cell(q0, lp0, 5, 1, 0.3, th)$gain, c(0, 0), tolerance = 1e-9)
  expect_true(all(x$gain >= -1e-12))
})

test_that("the gain grows with the number of loci and with the odds ratio", {
  th <- 0.05
  g <- function(k, or, p) sweep_cell(q0, lp0, k, or, p, th)$gain
  ks <- vapply(c(1, 2, 4, 8, 16), function(k) g(k, 1.5, 0.3), numeric(1)); expect_true(all(diff(ks) > 0))
  os <- vapply(c(1.1, 1.3, 1.6, 2, 3), function(o) g(6, o, 0.3), numeric(1)); expect_true(all(diff(os) > 0))
  expect_gt(g(10, 1.5, 0.3), g(10, 1.5, 0.05))
})

test_that("first_reaching finds the smallest value on the axis, or NA", {
  expect_equal(first_reaching(1:5, c(0, 0.1, 0.5, 0.9, 1.2), 0.5), 3)
  expect_true(is.na(first_reaching(1:3, c(0, 0.1, 0.2), 0.5)))
})

test_that("the discovery-sample formulas are inverses and give the expected scale", {
  expect_equal(sample_for_sigma(sigma_for_sample(1000, 0.5, 0.3), 0.5, 0.3), 1000)
  expect_equal(sigma_for_sample(4078, 0.0326, 0.3), 0.142, tolerance = 0.01)
  expect_gt(sigma_for_sample(1000, 0.5, 0.05), sigma_for_sample(1000, 0.5, 0.3))        # rarer alleles are estimated less precisely
})

test_that("the realised-score simulation reproduces the exact ceiling at zero error, and loses gain as the error grows", {
  th <- c(0.05, 0.2)
  ex <- sweep_cell(q0, lp0, 5, 2.0, 0.3, th)$gain
  r0 <- sweep_cell_realised(q0, lp0, 5, 2.0, 0.3, 0, th, R = 10, N = 2500, seed = 21L)
  m0 <- tapply(r0$gain, r0$threshold, mean)
  expect_equal(as.numeric(m0), ex, tolerance = 0.06)
  expect_equal(nrow(r0), 10L * 2L); expect_setequal(names(r0), c("n_loci", "odds_ratio", "raf", "sigma", "replicate", "threshold", "gain"))
  r1 <- sweep_cell_realised(q0, lp0, 5, 2.0, 0.3, 0.2, th, R = 10, N = 2500, seed = 21L)
  r2 <- sweep_cell_realised(q0, lp0, 5, 2.0, 0.3, 1.0, th, R = 10, N = 2500, seed = 21L)
  mg <- function(r) tapply(r$gain, r$threshold, mean)
  expect_true(all(mg(r2) < mg(r1))); expect_true(all(mg(r1) < m0 * 1.05))
  expect_true(all(mg(r2) < 0.7 * m0))                                                     # large errors cost a lot
  expect_equal(sweep_cell_realised(q0, lp0, 5, 2.0, 0.3, 0.2, th, R = 3, N = 1500, seed = 5L), sweep_cell_realised(q0, lp0, 5, 2.0, 0.3, 0.2, th, R = 3, N = 1500, seed = 5L))
})
