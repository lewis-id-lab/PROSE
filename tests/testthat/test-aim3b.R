source(file.path(proj_root, "R", "decision_curve.R"))
source(file.path(proj_root, "R", "evidence_mapping.R"))
source(file.path(proj_root, "R", "genetic_layer.R"))
source(file.path(proj_root, "R", "aim3b.R"))

lp0 <- profile_logit(pars)
q0 <- independent_profile_probs(pars$marginal)
p4 <- c(PTX3 = 0.48, CLEC7A = 0.08, TLR4 = 0.06, CD209 = 0.125)
st0 <- neutropenia_strata(q0, lp0)
b4 <- c(PTX3 = 0.3, CLEC7A = 0.1, TLR4 = 0.4, CD209 = 0.6)

test_that("the neutropenia strata partition the cohort and the baseline's events", {
  expect_equal(st0$N$share + st0$NN$share, 1)
  expect_equal(sum(st0$N$q), 1); expect_equal(sum(st0$NN$q), 1)
  neut <- all_profiles()[, "prolonged_neutropenia"] == 1
  expect_true(all(st0$N$q[!neut] == 0)); expect_true(all(st0$NN$q[neut] == 0))
  expect_equal(st0$N$share, sum(q0[neut]), tolerance = 1e-12)
  expect_equal(st0$N$events + st0$NN$events, sum(q0 * plogis(lp0)))
  expect_equal(st0$w, st0$N$events / (st0$N$events + st0$NN$events))
  expect_gt(st0$N$rate, st0$NN$rate)                       # neutropenia is a risk factor
})

test_that("modified stratum effects keep the pooled effect as the event-weighted average", {
  for (r in c(4, 1, 0.25, 0, 1.7)) {
    m <- modify_effects(b4, r, st0$w)
    expect_equal(st0$w * m$N + (1 - st0$w) * m$NN, b4, tolerance = 1e-12)
    if (r > 0) expect_equal(m$N, r * m$NN, tolerance = 1e-12)
  }
  expect_equal(modify_effects(b4, 1, st0$w)$N, b4)
  fo <- modify_effects(b4, Inf, st0$w); expect_equal(fo$NN, b4 * 0); expect_equal(fo$N, b4 / st0$w)
  bo <- modify_effects(b4, 0, st0$w); expect_equal(bo$N, b4 * 0); expect_equal(bo$NN, b4 / (1 - st0$w))
  expect_error(modify_effects(b4, -5, st0$w), "not compatible")
})

test_that("scenarios given as fractions are applied directly, and the difference matches the ratio form", {
  f <- scenario_effects(b4, c(N = -0.2, NN = 1.9), st0$w)
  expect_equal(f$N, b4 * -0.2); expect_equal(f$NN, b4 * 1.9)
  expect_equal(scenario_difference(c(N = -0.2, NN = 1.9), st0$w), -2.1)
  for (r in c(Inf, 4, 1, 0.25, 0)) {
    m <- modify_effects(rep(1, 4), r, st0$w)
    expect_equal(scenario_difference(r, st0$w), m$N[[1]] - m$NN[[1]], tolerance = 1e-12)
  }
})

cells_N <- trim_cells(genetic_cells(st0$N$q, lp0, dosage_grid(p4)))
cells_NN <- trim_cells(genetic_cells(st0$NN$q, lp0, dosage_grid(p4)))

test_that("trimmed cells keep only the weighted cells and their total weight", {
  full <- genetic_cells(st0$N$q, lp0, dosage_grid(p4))
  expect_lt(length(cells_N$w), length(full$w)); expect_true(all(cells_N$w > 0))
  expect_equal(sum(cells_N$w), 1); expect_equal(nrow(cells_N$G), length(cells_N$w))
  expect_equal(cells_N$prevalence, st0$N$rate)
})

test_that("multi_score_nb gives the baseline for a zero score, never favours a score over the true model, and rewards the right effects", {
  zero <- b4 * 0
  r <- multi_score_nb(cells_N, b4, list(baseline = zero, ceiling = b4, wrong = -b4, none = zero), c(0.05, 0.2))
  nb <- function(s) r$nb[r$score == s]
  expect_equal(nb("none"), nb("baseline"))
  expect_equal(nb("baseline"), dca_exact(cells_N$w, plogis(cells_N$lp_base + drop(cells_N$G %*% b4) + recentre_shift(cells_N$lp_base + drop(cells_N$G %*% b4), cells_N$w, cells_N$prevalence)),
                                         plogis(cells_N$lp_base), thresholds = c(0.05, 0.2))$nb)
  expect_gte(min(nb("ceiling") - nb("baseline")), -1e-10)
  expect_gt(min(nb("ceiling") - nb("wrong")), 0)
})

test_that("run_aim3b combines the strata by their shares and puts all of a one-sided truth's gain in that stratum", {
  eff <- rbind(b4, b4 * 1.5)
  cells <- list(N = cells_N, NN = cells_NN)
  res <- run_aim3b(st0, cells, eff, b4, list(homogeneous = 1, `H_F only` = Inf, `Brunel only` = 0), c(0.05, 0.2))
  expect_equal(nrow(res), 2 * 3 * 3 * 5 * 2)
  w <- res[res$level == "whole" & res$score == "ceiling", ]; n <- res[res$level == "N" & res$score == "ceiling", ]; nn <- res[res$level == "NN" & res$score == "ceiling", ]
  key <- function(d) paste(d$draw, d$scenario, d$threshold)
  expect_equal(w$delta_nb, st0$N$share * n$delta_nb[match(key(w), key(n))] + st0$NN$share * nn$delta_nb[match(key(w), key(nn))], tolerance = 1e-12)
  expect_true(all(abs(res$delta_nb[res$level == "N" & res$scenario == "Brunel only" & res$score == "ceiling"]) < 1e-12))   # no effect in N, so the oracle is the baseline
  expect_true(all(abs(res$delta_nb[res$level == "NN" & res$scenario == "H_F only" & res$score == "ceiling"]) < 1e-12))
  expect_true(all(res$delta_nb[res$score == "baseline"] == 0))
  expect_true(all(res$delta_nb[res$score == "ceiling"] >= -1e-10))
})

test_that("the interaction fit recovers the true difference on a large cohort and is calibrated when there is none", {
  set.seed(5530)
  n <- 2000
  X <- matrix(rbinom(n * 7, 1, rep(pars$marginal, each = n)), n, 7, dimnames = list(NULL, PREDICTORS))
  co <- list(rows = data.frame(patient = rep(paste0("P", 1:900), length.out = n), case = 0L, stringsAsFactors = FALSE), X = X)
  rc <- replicate_cohort(co, 25); cg <- cohort_groups(rc); neut <- rc$X[, "prolonged_neutropenia"] == 1
  q <- tabulate(profile_index(X), 128) / n; st <- neutropenia_strata(q, lp0)
  cn <- trim_cells(genetic_cells(st$N$q, lp0, dosage_grid(p4))); cnn <- trim_cells(genetic_cells(st$NN$q, lp0, dosage_grid(p4)))
  for (r in c(Inf, 0.25, 1)) {
    tr <- modify_effects(b4, r, st$w)
    sh <- list(N = genetic_scores(cn, tr$N, tr$N)$a_true, NN = genetic_scores(cnn, tr$NN, tr$NN)$a_true)
    set.seed(41)
    est <- t(vapply(1:3, function(i) simulate_interaction_study(cg, neut, lp0, p4, sh, tr$N, tr$NN, b4, resample = FALSE), numeric(5)))
    expect_equal(mean(est[, "delta"]), scenario_difference(r, st$w), tolerance = 0.35, scale = 1)
    expect_true(all(est[, "events_neutropenic"] <= est[, "events"]))
  }
  tr <- modify_effects(b4, 4, st$w); sh <- list(N = 0, NN = 0)
  set.seed(9); a <- simulate_interaction_study(cg, neut, lp0, p4, sh, tr$N, tr$NN, b4)
  set.seed(9); b <- simulate_interaction_study(cg, neut, lp0, p4, sh, tr$N, tr$NN, b4)
  expect_equal(a, b)
})
