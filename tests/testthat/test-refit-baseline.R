source(file.path(proj_root, "R", "bs_data.R"))
source(file.path(proj_root, "R", "decision_curve.R"))
source(file.path(proj_root, "R", "boscore_data.R"))
source(file.path(proj_root, "R", "outcome_risk.R"))
source(file.path(proj_root, "R", "refit_baseline.R"))

true_coef <- c("(Intercept)" = -5.45, pars$coef[PREDICTORS])

make_cohort <- function(n = 8000, seed = 4402, clustered = FALSE) {
  set.seed(seed)
  X <- sapply(PREDICTORS, function(p) rbinom(n, 1, pars$marginal[[p]]))
  y <- rbinom(n, 1, plogis(true_coef[[1]] + X %*% true_coef[PREDICTORS]))
  patient <- if (clustered) sample(sprintf("P%04d", 1:(n / 2)), n, replace = TRUE) else sprintf("P%05d", seq_len(n))
  list(rows = data.frame(patient = patient, case = y), X = X)
}

test_that("fit_baseline recovers known coefficients and returns usable parameters", {
  cohort <- make_cohort(n = 80000)
  f <- fit_baseline(cohort)
  expect_equal(names(f$coef), c("(Intercept)", PREDICTORS))
  expect_lt(max(abs(f$coef - true_coef)), 0.25)
  expect_equal(sum(f$fitted), sum(cohort$rows$case), tolerance = 1e-6)   # calibrated in the large
  rp <- refit_params(f$coef)
  expect_equal(profile_logit(rp, all_profiles()), as.numeric(f$coef[[1]] + all_profiles() %*% f$coef[PREDICTORS]))
  expect_equal(profile_risk_from_coef(f$coef), plogis(profile_logit(rp, all_profiles())))
})

test_that("with one row per patient and no small-sample factor the cluster-robust variance is the HC0 sandwich", {
  cohort <- make_cohort(n = 3000)
  D <- design_matrix(cohort$X); y <- cohort$rows$case
  fit <- glm(y ~ cohort$X, family = binomial); p <- fitted(fit)
  bread <- solve(t(D) %*% diag(p * (1 - p)) %*% D)
  meat <- t(D) %*% diag((y - p)^2) %*% D
  hc0 <- bread %*% meat %*% bread
  expect_equal(unname(cluster_robust_vcov(D, y, p, cohort$rows$patient, small_sample = FALSE)), unname(hc0), tolerance = 1e-8)
})

test_that("clustering inflates the robust variance only when outcomes are correlated within patients", {
  set.seed(7712)
  n <- 6000; g <- 600
  ids <- sprintf("P%03d", 1:g)
  patient <- rep(ids, each = n / g)
  frailty <- rnorm(g, 0, 1.2)[match(patient, ids)]                               # shared patient effect
  # Covariates that are constant within a patient are the case where clustering matters most
  Xp <- sapply(PREDICTORS, function(p) rbinom(g, 1, pars$marginal[[p]]))
  X <- Xp[match(patient, ids), ]
  y <- rbinom(n, 1, plogis(-2.5 + X %*% true_coef[PREDICTORS] + frailty))
  f <- fit_baseline(list(rows = data.frame(patient = patient, case = y), X = X))
  expect_gt(mean(f$se_robust / f$se_model), 1.3)
  indep <- make_cohort(n = 6000, seed = 12)
  fi <- fit_baseline(indep)
  expect_equal(mean(fi$se_robust / fi$se_model), 1, tolerance = 0.12)
})

test_that("the optimism bootstrap shows positive optimism for a small sample and little for a large one", {
  th <- c(0.05, 0.20)
  small <- make_cohort(n = 700, seed = 31)
  large <- make_cohort(n = 30000, seed = 32)
  o_s <- optimism_bootstrap(small, th, B = 150, seed = 5)
  o_l <- optimism_bootstrap(large, th, B = 60, seed = 5)
  expect_equal(colnames(o_s), c("auc_apparent", "auc_test", "slope_test", "nb_apparent_0.05", "nb_apparent_0.2", "nb_test_0.05", "nb_test_0.2"))
  opt <- function(o) mean(o[, "auc_apparent"] - o[, "auc_test"])
  expect_gt(opt(o_s), 0.005)
  expect_lt(opt(o_l), 0.004)
  expect_lt(mean(o_s[, "slope_test"]), 1)                       # an overfitted model's slope on new data is below 1
  expect_true(all(is.finite(o_s)))
})

test_that("the optimism bootstrap is reproducible for a seed", {
  cohort <- make_cohort(n = 1500, seed = 9)
  a <- optimism_bootstrap(cohort, 0.05, B = 20, seed = 77)
  b <- optimism_bootstrap(cohort, 0.05, B = 20, seed = 77)
  expect_identical(a, b)
})

test_that("observed_nb_rows matches a hand calculation", {
  p <- c(0.01, 0.10, 0.30, 0.60); y <- c(0, 1, 0, 1)
  expect_equal(observed_nb_rows(p, y, 0.20), (1 / 4) - (1 / 4) * 0.20 / 0.80)   # rows 3 and 4 flagged: one case, one non-case
  expect_equal(observed_nb_rows(p, y, c(0.05, 0.99)), c(2 / 4 - 1 / 4 * 0.05 / 0.95, 0))
})

restricted <- file.path(proj_root, "data", "restricted", "base_cohort_pseudonymised.csv")

test_that("on the real base cohort the refit matches the earlier sensitivity analysis and has the expected form", {
  skip_if_not(file.exists(restricted), "the pseudonymised base cohort is gitignored and not present")
  base <- select_cohort(cohort_from_pseudonymised(restricted))
  f <- fit_baseline(base)
  expect_equal(unname(f$coef), c(-5.591, 0.764, 0.862, 0.803, 0.933, 1.421, 1.575, 1.985), tolerance = 0.02)
  expect_true(all(f$se_robust > 0))
  expect_equal(sum(f$fitted), 133, tolerance = 1e-6)
  q <- tabulate(profile_index(base$X), 128L) / nrow(base$X)
  r <- profile_risk_from_coef(f$coef)
  expect_equal(dca_exact(q, r, r, 0.05)$mcid_b, 0.00458, tolerance = 1e-4)   # scenario C of addendum 4
})
