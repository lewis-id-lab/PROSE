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
