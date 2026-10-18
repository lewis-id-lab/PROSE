source(file.path(proj_root, "R", "decision_curve.R"))
source(file.path(proj_root, "R", "boscore_data.R"))
source(file.path(proj_root, "R", "outcome_risk.R"))

risk_pub <- plogis(profile_logit(pars, all_profiles()))

# Synthetic cohort whose outcomes follow known coefficients, drawn from the independent profile distribution
make_cohort <- function(n = 20000, coefs, seed = 6150) {
  set.seed(seed)
  X <- sapply(PREDICTORS, function(p) rbinom(n, 1, pars$marginal[[p]]))
  lp <- coefs[["(Intercept)"]] + X %*% coefs[PREDICTORS]
  list(rows = data.frame(patient = sample(sprintf("P%04d", 1:6000), n, replace = TRUE), case = rbinom(n, 1, plogis(lp))), X = X)
}
true_coef <- c("(Intercept)" = -5.45, pars$coef[PREDICTORS])

test_that("the refit recovers known coefficients and returns a valid risk for every profile", {
  cohort <- make_cohort(coefs = true_coef, n = 60000)
  cf <- refit_coef(cohort)
  expect_equal(names(cf), c("(Intercept)", PREDICTORS))
  expect_lt(max(abs(cf - true_coef)), 0.3)
  r <- profile_risk_from_coef(cf)
  expect_length(r, 128L)
  expect_true(all(r > 0 & r < 1))
  expect_equal(r, plogis(profile_logit(list(intercept = cf[["(Intercept)"]], coef = cf[PREDICTORS]), all_profiles())))
})

test_that("an in-sample logistic refit is calibrated in the large", {
  cohort <- make_cohort(coefs = true_coef)
  q <- tabulate(profile_index(cohort$X), 128L) / nrow(cohort$X)
  expect_equal(sum(q * profile_risk_from_coef(refit_coef(cohort))), mean(cohort$rows$case), tolerance = 1e-6)
})

test_that("when the refit equals the published equation, all three parametric scenarios coincide", {
  q <- independent_profile_probs(pars$marginal)
  s <- mcid_scenarios(q, risk_pub, risk_pub, c(0.05, 0.20))
  for (col in c("nb", "tp", "fp", "flagged", "mcid_b")) {
    a <- s[[col]][startsWith(s$scenario, "A")]
    expect_equal(s[[col]][startsWith(s$scenario, "B")], a)
    expect_equal(s[[col]][startsWith(s$scenario, "C")], a)
  }
})

test_that("a higher true risk raises the true positives and the MCID when the rule is unchanged", {
  q <- independent_profile_probs(pars$marginal)
  higher <- pmin(0.99, risk_pub * 1.3)
  s <- mcid_scenarios(q, risk_pub, higher, c(0.05, 0.20))
  a <- s[startsWith(s$scenario, "A"), ]; b <- s[startsWith(s$scenario, "B"), ]
  expect_true(all(b$tp > a$tp))
  expect_equal(b$mcid_b, 0.2 * b$tp)
  expect_equal(b$flagged, a$flagged)                          # same rule, so the same people are flagged
})

test_that("the observed-outcome scenario matches a hand calculation", {
  X <- matrix(0, 6, 7, dimnames = list(NULL, PREDICTORS)); X[1:3, "cmv"] <- 1       # three high-risk rows, three low-risk
  cohort <- list(rows = data.frame(case = c(1, 0, 0, 1, 0, 0), patient = 1:6), X = X)
  rp <- risk_pub
  d <- observed_outcome_scenario(cohort, rp, 0.05)
  p <- rp[profile_index(X)]
  flagged <- p >= 0.05
  expect_equal(d$flagged, mean(flagged))
  expect_equal(d$tp, sum(flagged & cohort$rows$case == 1) / 6)
  expect_equal(d$fp, sum(flagged & cohort$rows$case == 0) / 6)
  expect_equal(d$nb, d$tp - d$fp * 0.05 / 0.95)
  expect_equal(d$mcid_b, 0.2 * d$tp)
  # a row subset is respected
  d2 <- observed_outcome_scenario(cohort, rp, 0.05, rows = 1:3)
  expect_equal(d2$flagged, mean(p[1:3] >= 0.05))
})

test_that("on the real base cohort the refit truth gives a higher MCID than the published truth", {
  restricted <- file.path(proj_root, "data", "restricted", "base_cohort_pseudonymised.csv")
  skip_if_not(file.exists(restricted), "the pseudonymised base cohort is gitignored and not present")
  base <- select_cohort(cohort_from_pseudonymised(restricted))
  q <- tabulate(profile_index(base$X), 128L) / nrow(base$X)
  s <- mcid_scenarios(q, risk_pub, profile_risk_from_coef(refit_coef(base)), 0.05)
  a <- s$mcid_b[startsWith(s$scenario, "A")]; b <- s$mcid_b[startsWith(s$scenario, "B")]
  expect_equal(a, 0.004043, tolerance = 5e-4)                  # the adopted value in addendum 3
  expect_gt(b, a)
})
