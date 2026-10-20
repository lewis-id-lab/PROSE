# The profile-preserving construction of the genetic truth on cells (addendum 24, D27): R/genetic_layer.R and R/aim3b.R
source(file.path(proj_root, "R", "decision_curve.R"))
source(file.path(proj_root, "R", "evidence_mapping.R"))
source(file.path(proj_root, "R", "genetic_layer.R"))
source(file.path(proj_root, "R", "aim3b.R"))
source(file.path(proj_root, "R", "calibration_check.R"))

lp_s <- profile_logit(pars)
q_s <- independent_profile_probs(pars$marginal)
p_s <- c(PTX3 = 0.48, CLEC7A = 0.08, TLR4 = 0.06, CD209 = 0.125)
grid_s <- dosage_grid(p_s)
cells_s <- genetic_cells(q_s, lp_s, grid_s)
b_s <- c(PTX3 = 0.4, CLEC7A = 0.2, TLR4 = 0.5, CD209 = 0.7)
bs_s <- c(PTX3 = 0.3, CLEC7A = 0, TLR4 = 0.35, CD209 = 0.6)
th_s <- c(0.05, 0.2)

# Mean of a cell-level quantity within each clinical profile, over genotype
by_profile_mean <- function(x, cells) tapply(cells$w * x, cells$profile, sum) / tapply(cells$w, cells$profile, sum)

test_that("the scheme comes from PROSE_SCHEME, defaults to global, and rejects anything else", {
  withr::local_envvar(PROSE_SCHEME = NA)
  expect_equal(genetic_scheme(), "global")
  expect_equal(scheme_suffix(), "")
  withr::local_envvar(PROSE_SCHEME = "profile")
  expect_equal(genetic_scheme(), "profile")
  expect_equal(scheme_suffix(), "_profile")
  withr::local_envvar(PROSE_SCHEME = "nonsense")
  expect_error(genetic_scheme())
})

test_that("the global scheme is unchanged: one shift that holds the prevalence", {
  s <- genetic_scores(cells_s, b_s, bs_s, scheme = "global")
  a <- recentre_shift(cells_s$lp_base + drop(cells_s$G %*% b_s), cells_s$w, cells_s$prevalence)
  expect_equal(s$a_true, a)
  expect_length(s$a_true, 1)
  expect_equal(sum(cells_s$w * plogis(s$lp_true)), cells_s$prevalence, tolerance = 1e-9)
})

test_that("the profile scheme holds each profile's risk, for the truth and for both scores", {
  s <- genetic_scores(cells_s, b_s, bs_s, scheme = "profile")
  for (lp in list(s$lp_true, s$lp_ceiling, s$lp_real))
    expect_equal(as.numeric(by_profile_mean(plogis(lp), cells_s)), plogis(lp_s), tolerance = 1e-9)
  expect_length(s$a_true, length(q_s))
})

test_that("the cell-level profile shifts agree with the matrix implementation in calibration_check.R", {
  s <- genetic_scores(cells_s, b_s, bs_s, scheme = "profile")
  gl_t <- drop(grid_s$G %*% b_s); gl_s <- drop(grid_s$G %*% bs_s)
  expect_equal(plogis(s$lp_true), as.numeric(scheme_risk(q_s, lp_s, gl_t, grid_s$prob, "profile")), tolerance = 1e-8)
  expect_equal(plogis(s$lp_real), as.numeric(scheme_risk(q_s, lp_s, gl_s, grid_s$prob, "profile")), tolerance = 1e-8)
  g <- genetic_scores(cells_s, b_s, bs_s, scheme = "global")
  expect_equal(plogis(g$lp_true), as.numeric(scheme_risk(q_s, lp_s, gl_t, grid_s$prob, "global")), tolerance = 1e-8)
})

test_that("trimmed stratum cells keep the property, and the shift is a per-profile vector that reproduces the cell logits", {
  st <- neutropenia_strata(q_s, lp_s)
  cN <- trim_cells(genetic_cells(st$N$q, lp_s, grid_s))
  s <- genetic_scores(cN, b_s, b_s, scheme = "profile")
  expect_equal(as.numeric(by_profile_mean(plogis(s$lp_true), cN)), plogis(lp_s[sort(unique(cN$profile))]), tolerance = 1e-9)
  gl <- drop(cN$G %*% b_s)
  expect_equal(cN$lp_base + gl + sh_at(s$a_true, cN$profile), s$lp_true, tolerance = 1e-9)
})

test_that("a profile with no weight gets no shift and does not disturb the others", {
  q0 <- q_s; q0[c(3, 17)] <- 0; q0 <- q0 / sum(q0)
  cz <- genetic_cells(q0, lp_s, grid_s)
  s <- genetic_scores(cz, b_s, bs_s, scheme = "profile")
  expect_true(all(is.finite(s$lp_true)))
  keep <- setdiff(seq_along(q0), c(3, 17))
  expect_equal(as.numeric(by_profile_mean(plogis(s$lp_true), cz)[as.character(keep)]), plogis(lp_s[keep]), tolerance = 1e-9)
  expect_equal(s$a_true[c(3, 17)], c(0, 0))
})

test_that("with a true-risk offset the truth holds the offset risk and the scores hold the clinical risk", {
  co <- genetic_cells(q_s, lp_s, grid_s, offset = 0.7)
  s <- genetic_scores(co, b_s, bs_s, scheme = "profile")
  expect_equal(as.numeric(by_profile_mean(plogis(s$lp_true), co)), plogis(lp_s + 0.7), tolerance = 1e-9)
  expect_equal(as.numeric(by_profile_mean(plogis(s$lp_ceiling), co)), plogis(lp_s), tolerance = 1e-9)
})

test_that("with no genetic effect every scheme gives the baseline", {
  z <- setNames(rep(0, 4), names(p_s))
  for (sch in c("global", "profile")) {
    s <- genetic_scores(cells_s, z, z, scheme = sch)
    expect_equal(s$lp_true, cells_s$lp_base, tolerance = 1e-9)
    expect_equal(s$lp_real, cells_s$lp_base, tolerance = 1e-9)
  }
})

test_that("Shapley contributions sum to the gain of the full score under both schemes", {
  for (sch in c("global", "profile")) {
    r <- genetic_nb(cells_s, b_s, bs_s, th_s, scheme = sch)
    gain <- function(sc) r$nb[r$scenario == sc] - r$nb[r$scenario == "baseline"]
    sh <- shapley_gain(cells_s, b_s, bs_s, th_s, scheme = sch)
    expect_equal(unname(sh$total), gain("realised"), tolerance = 1e-9)
    expect_equal(unname(rowSums(sh$phi)), gain("realised"), tolerance = 1e-9)
    sc <- shapley_gain(cells_s, b_s, b_s, th_s, scheme = sch)
    expect_equal(unname(sc$total), gain("ceiling"), tolerance = 1e-9)
  }
})

test_that("multi_score_nb gives the engine's gain for the same scores under both schemes", {
  zero <- b_s * 0
  for (sch in c("global", "profile")) {
    m <- multi_score_nb(cells_s, b_s, list(baseline = zero, ceiling = b_s, realised = bs_s), th_s, scheme = sch)
    r <- genetic_nb(cells_s, b_s, bs_s, th_s, scheme = sch)
    for (sc in c("ceiling", "realised"))
      expect_equal(m$nb[m$score == sc] - m$nb[m$score == "baseline"], r$nb[r$scenario == sc] - r$nb[r$scenario == "baseline"], tolerance = 1e-9)
  }
})

test_that("under the profile scheme the baseline is calibrated: its true positives are those of the clinical score", {
  r <- genetic_nb(cells_s, b_s * 3, bs_s, th_s, scheme = "profile")
  b <- r[r$scenario == "baseline", ]
  expect_equal(b$tp, vapply(th_s, function(t) sum(q_s * plogis(lp_s) * (plogis(lp_s) >= t)), numeric(1)), tolerance = 1e-9)
})

test_that("a simulated study under the profile scheme averages to the exact gain", {
  skip_on_cran()
  withr::local_envvar(PROSE_SCHEME = "profile")
  cohort <- list(X = all_profiles()[sample.int(128, 1500, replace = TRUE, prob = q_s), ], rows = data.frame(patient = seq_len(1500)))
  cg <- cohort_groups(cohort)
  q_emp <- tabulate(profile_index(cohort$X), 128L) / 1500
  cells <- genetic_cells(q_emp, lp_s, grid_s)
  b <- b_s * 2
  sh <- genetic_scores(cells, b, b)
  set.seed(4471)
  sims <- do.call(rbind, lapply(1:300, function(i) { s <- simulate_study(cg, lp_s, p_s, sh, b, b, c(0.05), c(`0.05` = 0.00458), resample = FALSE); s[s$scenario == "ceiling", ] }))
  ex <- genetic_nb(cells, b, b, 0.05)
  exact <- ex$nb[ex$scenario == "ceiling"] - ex$nb[ex$scenario == "baseline"]
  expect_lt(abs(mean(sims$delta_nb) - exact), 4 * sd(sims$delta_nb) / sqrt(nrow(sims)) + 1e-4)
})
