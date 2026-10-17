test_that("published coefficients reproduce the Table 2 adjusted odds ratios", {
  expect_lt(max(abs(exp(pars$coef) - pars$odds_ratio)), 0.005)
  expect_equal(pars$intercept, -5.45)
})

test_that("published prevalences agree with the published group counts", {
  expect_equal(c(pars$n_case, pars$n_non), c(133, 3994))
  expect_lt(max(abs(pars$count_case / pars$n_case - pars$prev_case)), 5e-4)
  expect_lt(max(abs(pars$count_non / pars$n_non - pars$prev_non)), 5e-4)
  overall <- (pars$count_case + pars$count_non) / (pars$n_case + pars$n_non)
  expect_lt(max(abs(overall - pars$marginal)), 5e-4)
})

test_that("validation targets load as published", {
  expect_equal(targets$prevalence, 0.033)
  expect_equal(targets$aroc, 0.85)
  expect_equal(targets$aroc_boot, 0.83)
  # 133/4127 is 3.22%; the 3.3% target is kept as pre-specified (D6) but is within 0.1 pp
  expect_lt(abs(133 / 4127 - targets$prevalence), 0.001)
})

test_that("profile enumeration is complete and indexed consistently", {
  expect_equal(dim(profiles), c(128L, 7L))
  expect_equal(nrow(unique(profiles)), 128L)
  expect_equal(profile_index(profiles), 1:128)
})

test_that("independent profile distribution is a proper distribution with the right margins", {
  q <- independent_profile_probs(pars$marginal)
  expect_equal(sum(q), 1)
  expect_equal(as.numeric(crossprod(profiles, q)), as.numeric(pars$marginal[colnames(profiles)]))
})

test_that("grouped expected AUC equals brute-force pairwise AUC, including ties", {
  set.seed(5207)
  lp <- round(rnorm(40), 1)
  w1 <- runif(40); w0 <- runif(40)
  brute <- 0
  for (i in 1:40) for (j in 1:40) {
    brute <- brute + w1[i] * w0[j] * ((lp[i] > lp[j]) + 0.5 * (lp[i] == lp[j]))
  }
  brute <- brute / (sum(w1) * sum(w0))
  expect_equal(expected_auc(lp, w1, w0), brute)
})

test_that("profile metrics satisfy the case / non-case mixture identity", {
  q <- profile_probs_mc(exchangeable_corr(0.2), pars$marginal, 5e5, seed = 2764)
  m <- profile_metrics(q, pars)
  expect_equal(m$prevalence * m$case_prev + (1 - m$prevalence) * m$noncase_prev, m$marginal)
})

test_that("independence reproduces the documented 2.27% prevalence and aROC of 0.775", {
  m <- profile_metrics(independent_profile_probs(pars$marginal), pars)
  expect_lt(abs(m$prevalence - 0.0227), 0.001)
  expect_lt(abs(m$auc - 0.775), 0.005)
})
