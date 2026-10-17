test_that("net benefit matches a hand calculation on a two-profile distribution", {
  risk <- c(0.02, 0.30)
  q <- c(0.9, 0.1)
  nb <- net_benefit(q, pars, 0.10, risk = risk)
  # Only the second profile is flagged at t = 0.10
  tp <- 0.1 * 0.30; fp <- 0.1 * 0.70
  expect_equal(nb$tp, tp)
  expect_equal(nb$fp, fp)
  expect_equal(nb$nb, tp - fp * 0.10 / 0.90)
  prev <- 0.9 * 0.02 + 0.1 * 0.30
  expect_equal(nb$nb_all, prev - (1 - prev) * 0.10 / 0.90)
})

test_that("the net-benefit code reproduces the exchangeable-0.30 reference table", {
  ref <- read.csv(file.path(proj_root, "data", "baseline_net_benefit_reference.csv"))
  q <- profile_probs_mc(exchangeable_corr(0.30), pars$marginal, 4e6, seed = 6127)
  for (thr in c(0.05, 0.10, 0.20)) {
    r <- ref[abs(ref$threshold - thr) < 1e-9, ]
    nb <- net_benefit(q, pars, thr)
    # Reference cohort has 3.25% prevalence against 3.27% here, so a small tolerance
    expect_lt(abs(nb$nb - r$nb_model), 0.001)
    expect_lt(abs(nb$tp * 1000 - r$tp_per_1000), 1.5)
    expect_lt(abs(nb$flagged - r$pct_flagged), 0.01)
  }
})

test_that("treat-all net benefit is negative at thresholds above prevalence", {
  q <- independent_profile_probs(pars$marginal)
  expect_lt(net_benefit(q, pars, 0.05)$nb_all, 0)
})

test_that("a screening rule cannot do worse than treat-none when threshold equals its cut-off", {
  q <- profile_probs_mc(exchangeable_corr(0.30), pars$marginal, 1e6, seed = 4821)
  expect_gt(net_benefit(q, pars, 0.05)$nb, 0)
  expect_gt(net_benefit(q, pars, 0.05)$nb, net_benefit(q, pars, 0.05)$nb_all)
})
