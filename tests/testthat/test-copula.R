theta_ref <- c(rho_A = 0.24, d_malig_chemo = 0.05, d_chemo_neut = 0.30, d_neut_lymph = 0.20,
               d_chemo_lymph = 0.10, rho_A_steroid = 0.05, rho_A_cmv = 0.25, rho_steroid_cmv = 0.15,
               rho_prior = 0.15)

test_that("structured correlation matrix is valid and placed correctly", {
  R <- build_corr(theta_ref)
  expect_true(isSymmetric(R))
  expect_equal(unname(diag(R)), rep(1, 7))
  expect_false(is.null(tryCatch(chol(R), error = function(e) NULL)))
  expect_equal(R["high_risk_chemo", "prolonged_neutropenia"], 0.24 + 0.30)
  expect_equal(R["uncontrolled_malignancy", "high_risk_chemo"], 0.24 + 0.05)
  expect_equal(R["prolonged_neutropenia", "severe_lymphopenia"], 0.24 + 0.20)
  expect_equal(R["high_risk_chemo", "severe_lymphopenia"], 0.24 + 0.10)
  expect_equal(R["uncontrolled_malignancy", "severe_lymphopenia"], 0.24)
  expect_equal(R["uncontrolled_malignancy", "prolonged_neutropenia"], 0.24)
  expect_equal(R["high_dose_steroids", "cmv"], 0.15)
  expect_equal(R["prior_mould_disease", "cmv"], 0.15)
  expect_equal(R["cmv", "prolonged_neutropenia"], 0.25)
})

test_that("version 1 parameter vectors still build, with chemotherapy-lymphopenia at the baseline", {
  R1 <- build_corr(theta_ref[THETA_NAMES_V1])
  expect_equal(R1["high_risk_chemo", "severe_lymphopenia"], 0.24)
  R2 <- build_corr(replace(theta_ref, "d_chemo_lymph", 0))
  expect_equal(R1, R2)
})

test_that("parameter bounds are named consistently", {
  expect_setequal(names(THETA_LOWER), THETA_NAMES)
  expect_setequal(names(THETA_UPPER), THETA_NAMES)
  expect_setequal(names(THETA_UPPER_CAPPED), THETA_NAMES)
  expect_setequal(names(THETA_LOWER_V1), THETA_NAMES_V1)
  expect_true(all(THETA_LOWER <= THETA_UPPER))
  expect_true(all(THETA_UPPER_CAPPED <= THETA_UPPER))
})

test_that("a non-positive-definite matrix is rejected by the sampler", {
  bad <- exchangeable_corr(-0.5)
  set.seed(4390)
  Z0 <- matrix(rnorm(1000 * 7), 1000)
  expect_null(profile_probs_crn(bad, Z0, pars$marginal))
})

test_that("pd_adjust leaves valid matrices alone and repairs invalid ones continuously", {
  good <- build_corr(theta_ref)
  expect_identical(pd_adjust(good)$R, good)
  expect_equal(pd_adjust(good)$penalty, 0)

  bad <- exchangeable_corr(-0.5)
  fixed <- pd_adjust(bad)
  expect_gt(fixed$penalty, 0)
  expect_equal(unname(diag(fixed$R)), rep(1, 7))
  expect_false(is.null(tryCatch(chol(fixed$R), error = function(e) NULL)))

  # Just inside and just outside the boundary give nearly identical matrices
  eps <- 1e-3
  near_in  <- exchangeable_corr(-1 / 6 + 0.9 * eps / 1)
  near_out <- exchangeable_corr(-1 / 6 - 1e-6)
  expect_lt(max(abs(pd_adjust(near_in)$R - pd_adjust(near_out)$R)), 0.01)
})

test_that("the fitting objective is finite and continuous across the feasibility boundary", {
  set.seed(2915)
  Z0 <- matrix(rnorm(5e4 * 7), 5e4)
  obj <- make_objective(Z0, pars, targets)
  th <- replace(theta_ref, c("d_chemo_neut", "d_neut_lymph", "d_chemo_lymph"), c(0.7, 0.7, 0))
  at <- function(d) replace(th, "d_chemo_lymph", d)
  min_eig <- function(d) min(eigen(build_corr(at(d)), symmetric = TRUE, only.values = TRUE)$values)

  # d_chemo_lymph = 0 is infeasible here (chemotherapy and lymphopenia are too weakly linked
  # given how strongly neutropenia is tied to each); a larger value restores feasibility
  expect_lt(min_eig(0), 1e-3)
  expect_gt(min_eig(0.7), 1e-3)
  expect_true(is.finite(obj(at(0))))

  boundary <- uniroot(function(d) min_eig(d) - 1e-3, c(0, 0.7), tol = 1e-10)$root
  expect_lt(abs(obj(at(boundary - 1e-4)) - obj(at(boundary + 1e-4))), 5)
})

test_that("the copula preserves the published predictor margins", {
  q <- profile_probs_mc(build_corr(theta_ref), pars$marginal, 1e6, seed = 8412)
  expect_equal(sum(q), 1)
  expect_lt(max(abs(as.numeric(crossprod(profiles, q)) - pars$marginal)), 0.003)
})

test_that("an identity latent matrix reproduces the independent distribution", {
  q <- profile_probs_mc(diag(7), pars$marginal, 2e6, seed = 1936)
  expect_lt(max(abs(q - independent_profile_probs(pars$marginal))), 0.002)
})

test_that("exchangeable rho = 0.30 reproduces the documented 3.3% prevalence and aROC of 0.85", {
  q <- profile_probs_mc(exchangeable_corr(0.30), pars$marginal, 2e6, seed = 6127)
  m <- profile_metrics(q, pars)
  expect_lt(abs(m$prevalence - 0.033), 0.003)
  expect_lt(abs(m$auc - 0.852), 0.01)
})

test_that("numerical integration reproduces the independent distribution and agrees with sampling", {
  qe <- profile_probs_exact(diag(7), pars$marginal, abseps = 1e-6, seed = 3301)
  expect_equal(sum(qe), 1)
  expect_lt(max(abs(qe - independent_profile_probs(pars$marginal))), 1e-5)

  R <- build_corr(theta_ref)
  q1 <- profile_probs_exact(R, pars$marginal, abseps = 1e-6, seed = 3302)
  q2 <- profile_probs_exact(R, pars$marginal, abseps = 1e-6, seed = 7714)
  expect_lt(max(abs(q1 - q2)), 2e-5)
  qm <- profile_probs_mc(R, pars$marginal, 4e6, seed = 5128)
  expect_lt(max(abs(q1 - qm)), 0.0015)
})

test_that("more latent correlation raises both prevalence and discrimination", {
  m <- lapply(c(0, 0.2, 0.4), function(r) {
    profile_metrics(profile_probs_mc(exchangeable_corr(r), pars$marginal, 1e6, seed = 3358), pars)
  })
  expect_true(all(diff(sapply(m, `[[`, "prevalence")) > 0))
  expect_true(all(diff(sapply(m, `[[`, "auc")) > 0))
})

test_that("target table has the 16 published targets with positive standard errors", {
  tt <- target_table(pars, targets)
  expect_equal(nrow(tt), 16L)
  expect_true(all(tt$se > 0))
})
