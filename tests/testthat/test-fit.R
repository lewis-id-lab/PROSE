fit_path <- file.path(proj_root, "output", "fitted_theta.csv")

load_fit <- function() {
  d <- read.csv(fit_path)
  setNames(d$value, d$parameter)[THETA_NAMES]
}

test_that("the saved fit is feasible and reproduces the published targets on fresh draws", {
  skip_if_not(file.exists(fit_path), "run scripts/01_fit_correlation.R first")
  th <- load_fit()
  expect_false(anyNA(th))
  expect_true(all(th >= THETA_LOWER - 1e-9 & th <= THETA_UPPER + 1e-9))

  pd <- pd_adjust(build_corr(th))
  expect_equal(pd$shift, 0)

  tt <- target_table(pars, targets)
  m <- profile_metrics(profile_probs_mc(pd$R, pars$marginal, 5e6, seed = 8153), pars)
  z <- (metrics_vector(m) - tt$published) / tt$se
  names(z) <- tt$target

  # Global targets: prevalence within 1 SE, aROC inside the published 95% CI
  expect_lt(abs(z[["prevalence"]]), 1)
  expect_gt(m$auc, 0.84); expect_lt(m$auc, 0.86)
  # Per-predictor group prevalences: none beyond 3 SE
  expect_lt(max(abs(z[grepl("prev:", names(z))])), 3)
  expect_lt(sum(z^2), 15)
})

test_that("the saved matrix and profile table match the saved parameters", {
  skip_if_not(file.exists(fit_path), "run scripts/01_fit_correlation.R first")
  th <- load_fit()
  Rsaved <- as.matrix(read.csv(file.path(proj_root, "output", "fitted_correlation_matrix.csv"), row.names = 1))
  expect_equal(Rsaved, pd_adjust(build_corr(th))$R, ignore_attr = TRUE, tolerance = 1e-8)

  pt <- read.csv(file.path(proj_root, "output", "profile_probabilities.csv"))
  expect_equal(nrow(pt), 128L)
  expect_equal(sum(pt$q), 1, tolerance = 1e-6)
  expect_equal(pt$risk, plogis(profile_logit(pars, profiles)))
  tt <- target_table(pars, targets)
  m <- profile_metrics(pt$q, pars)
  z <- (metrics_vector(m) - tt$published) / tt$se
  expect_lt(max(abs(z)), 2)
})

test_that("the fit beats independence and the exchangeable benchmark, and version 2 beats version 1", {
  skip_if_not(file.exists(fit_path), "run scripts/01_fit_correlation.R first")
  loss <- read.csv(file.path(proj_root, "output", "reproduction_loss.csv"))
  chi2 <- setNames(loss$chi2, loss$model)
  expect_lt(chi2[["fitted"]], chi2[["fitted_v1"]])
  expect_lt(chi2[["fitted_v1"]], chi2[["exchangeable_0.30"]])
  expect_lt(chi2[["exchangeable_0.30"]], chi2[["independent"]])
})
