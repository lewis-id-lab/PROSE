source(file.path(proj_root, "R", "decision_curve.R"))
source(file.path(proj_root, "R", "evidence_mapping.R"))
source(file.path(proj_root, "R", "genetic_layer.R"))
source(file.path(proj_root, "R", "prs_benchmark.R"))

lp_t <- profile_logit(pars)
q_t <- independent_profile_probs(pars$marginal)

test_that("the truth is recentred to the baseline prevalence", {
  b <- log(2)
  a <- prs_shift(q_t, lp_t, b)
  z <- prs_nodes(400)
  expect_equal(sum(outer(q_t, rep(1 / 400, 400)) * plogis(outer(lp_t, b * z, "+") + a)), sum(q_t * plogis(lp_t)), tolerance = 1e-10)
})

test_that("a score with no accuracy gives no gain, and gain rises with accuracy", {
  b <- log(2.5)
  g <- vapply(c(0, 0.5, 0.8, 1), function(r) prs_gain(q_t, lp_t, b, r, 0.2)$gain, numeric(1))
  expect_equal(g[1], 0, tolerance = 1e-12)
  expect_true(all(diff(g) > 0))
})

test_that("the ceiling equals the genetic-layer engine on the same nodes", {
  b <- log(3)
  z <- prs_nodes(400)
  cells <- genetic_cells(q_t, lp_t, list(G = matrix(z, dimnames = list(NULL, "PRS")), prob = rep(1 / 400, 400)))
  r <- genetic_nb(cells, c(PRS = b), c(PRS = b), c(0.05, 0.2))
  ref <- r$nb[r$scenario == "ceiling"] - r$nb[r$scenario == "baseline"]
  expect_equal(prs_gain(q_t, lp_t, b, 1, c(0.05, 0.2))$gain, ref, tolerance = 1e-8)
})

test_that("discovery size and accuracy are inverse", {
  b <- log(2); m <- 1e4; n <- 20000
  rho <- prs_accuracy(n, 0.5, b, m)
  expect_equal(prs_n_effective(rho, b, m), n * 0.25)
  expect_lt(prs_accuracy(n, 0.5, b, 6e4), rho)
})
