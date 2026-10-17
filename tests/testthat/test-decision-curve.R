source(file.path(proj_root, "R", "decision_curve.R"))

risk_all <- plogis(profile_logit(pars, profiles))
q_ind <- independent_profile_probs(pars$marginal, profiles)

test_that("dca_exact agrees with net_benefit() at several thresholds", {
  curve <- dca_exact(q_ind, risk_all, thresholds = c(0.01, 0.05, 0.10, 0.20, 0.30))
  for (i in seq_len(nrow(curve))) {
    ref <- net_benefit(q_ind, pars, curve$threshold[i])
    expect_equal(curve$nb[i], ref$nb)
    expect_equal(curve$nb_all[i], ref$nb_all)
    expect_equal(curve$tp[i], ref$tp)
    expect_equal(curve$fp[i], ref$fp)
    expect_equal(curve$flagged[i], ref$flagged)
  }
})

test_that("hand calculation, with decisions and outcomes using different risks", {
  w <- c(0.5, 0.3, 0.2)
  rt <- c(0.02, 0.10, 0.40)    # what actually happens
  rp <- c(0.01, 0.30, 0.06)    # what the rule sees
  cv <- dca_exact(w, rt, rp, thresholds = c(0.05, 0.20))
  # t = 0.05 flags cells 2 and 3; t = 0.20 flags cell 2 only
  expect_equal(cv$tp, c(0.3 * 0.10 + 0.2 * 0.40, 0.3 * 0.10))
  expect_equal(cv$fp, c(0.3 * 0.90 + 0.2 * 0.60, 0.3 * 0.90))
  expect_equal(cv$nb[1], cv$tp[1] - cv$fp[1] * 0.05 / 0.95)
  expect_equal(cv$flagged, c(0.5, 0.3))
  expect_equal(cv$mcid_b, 0.2 * cv$tp)
  prev <- sum(w * rt)
  expect_equal(cv$prevalence, rep(prev, 2))
  expect_equal(cv$snb, cv$nb / prev)
})

test_that("true positives, false positives and flagged share fall as the threshold rises", {
  cv <- dca_exact(q_ind, risk_all)
  expect_true(all(diff(cv$tp) <= 1e-12))
  expect_true(all(diff(cv$fp) <= 1e-12))
  expect_true(all(diff(cv$flagged) <= 1e-12))
})

test_that("a threshold below every risk reproduces treat-all; above every risk, treat-none", {
  cv <- dca_exact(q_ind, risk_all, thresholds = c(min(risk_all) / 2, 0.999))
  expect_equal(cv$nb[1], cv$nb_all[1])
  expect_equal(cv$nb[2], 0)
  expect_equal(cv$flagged, c(1, 0))
})

test_that("comparing a curve with itself gives zero gain, which does not meet a positive MCID", {
  cv <- dca_exact(q_ind, risk_all)
  cmp <- dca_compare(cv, cv)
  expect_true(all(cmp$delta_nb == 0))
  expect_false(any(cmp$meets_mcid[cmp$mcid_b > 0]))
})

test_that("a better-informed rule shows a positive gain over the baseline", {
  # Outcomes follow a richer risk than the baseline sees: decisions on the richer risk do better
  set.seed(8836)
  rt <- pmin(0.95, risk_all * exp(rnorm(128, 0, 0.5)))
  base <- dca_exact(q_ind, rt, risk_all, thresholds = c(0.05, 0.10))
  better <- dca_exact(q_ind, rt, rt, thresholds = c(0.05, 0.10))
  expect_true(all(dca_compare(base, better)$delta_nb > 0))
})

test_that("dca_summary picks out the requested thresholds", {
  s <- dca_summary(dca_exact(q_ind, risk_all), at = c(0.05, 0.20))
  expect_equal(s$threshold, c(0.05, 0.20))
  expect_equal(s$tp_per_1000, 1000 * dca_exact(q_ind, risk_all, thresholds = c(0.05, 0.20))$tp)
})

test_that("plot helpers return ggplot objects", {
  cv <- dca_exact(q_ind, risk_all, thresholds = seq(0.01, 0.30, by = 0.01))
  expect_s3_class(plot_dca(cv), "ggplot")
  expect_s3_class(plot_dca_rates(cv), "ggplot")
  expect_s3_class(plot_mcid(cv), "ggplot")
})
