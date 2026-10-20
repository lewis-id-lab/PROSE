# Aim 4: the genetic-architecture requirements sweep (addendum 16). Requires R/baseline.R, R/decision_curve.R and R/genetic_layer.R.
#
# k independent loci, each with the same risk-allele frequency p and the same per-allele odds ratio, log-additive, in Hardy-Weinberg equilibrium and unlinked
# (an assumption; the real panel is four terms). The genetic score is the number of risk alleles D, whose distribution is the exact convolution of the
# k per-locus genotype distributions, which for identical loci is Binomial(2k, p). The true logit is the baseline logit + log(OR) D, recentred so that
# the event rate is unchanged (D5). The score is the ceiling one (true weights), so the sweep gives the best a genetic score of that architecture could do.

locus_pmf <- function(p) c((1 - p)^2, 2 * p * (1 - p), p^2)

# Exact convolution of two probability mass functions on 0, 1, 2, ...
convolve_pmf <- function(a, b) {
  out <- numeric(length(a) + length(b) - 1L)
  for (i in seq_along(a)) out[i:(i + length(b) - 1L)] <- out[i:(i + length(b) - 1L)] + a[i] * b
  out
}

# Distribution of the number of risk alleles over n_loci loci with frequencies p (one value for all loci, or one per locus), by exact convolution
dosage_pmf <- function(n_loci, p) {
  p <- rep_len(p, n_loci)
  Reduce(convolve_pmf, lapply(p, locus_pmf))
}

# Net-benefit gain of the ceiling score for one architecture on a profile distribution q with baseline logits lp
sweep_cell <- function(q, lp_profile, n_loci, odds_ratio, raf, thresholds) {
  pmf <- dosage_pmf(n_loci, raf)
  grid <- list(G = matrix(seq_along(pmf) - 1, ncol = 1, dimnames = list(NULL, "D")), prob = pmf)
  cells <- genetic_cells(q, lp_profile, grid)
  b <- c(D = log(odds_ratio))
  r <- genetic_nb(cells, b, b, thresholds)
  base <- r[r$scenario == "baseline", ]; ceil <- r[r$scenario == "ceiling", ]
  data.frame(n_loci = n_loci, odds_ratio = odds_ratio, raf = raf, threshold = thresholds,
             baseline_nb = base$nb, gain = ceil$nb - base$nb, extra_cases_per_1000 = 1000 * (ceil$tp - base$tp),
             extra_false_positives_per_1000 = 1000 * (ceil$fp - base$fp), flagged = ceil$flagged, auc = ceil$auc, baseline_auc = base$auc)
}

# The smallest k (or OR) on a sorted axis at which the gain reaches the MCID; NA if never
first_reaching <- function(axis, gain, mcid) { i <- which(gain >= mcid); if (length(i)) axis[min(i)] else NA_real_ }

# ---- The sweep with a realised score (addendum 17) ---------------------------------------------------------------------------------------------
# The weights of the score are estimates: each locus's log odds ratio carries an independent error with standard deviation sigma (its stated precision), so the score
# is sum_j (log OR + e_j) g_j while the truth is log OR times the number of risk alleles. The error is unbiased: the loci are fixed in advance (no winner's curse) and the
# estimates come from a discovery study of fixed size. With estimated weights the joint distribution of (number of alleles, score) has no small exact form, so the
# genotypes are simulated: N individuals per replicate crossed with the exact profile distribution (the truth's recentring and the baseline are exact). Each replicate is one
# discovery study's weights; the spread over replicates is the study-to-study variation in how good the score is.
#
# Per-locus standard error for a discovery sample of n subjects with case fraction pi: sigma = 1 / sqrt(n pi (1 - pi) 2 p (1 - p)).
sigma_for_sample <- function(n, case_fraction, raf) 1 / sqrt(n * case_fraction * (1 - case_fraction) * 2 * raf * (1 - raf))
sample_for_sigma <- function(sigma, case_fraction, raf) 1 / (sigma^2 * case_fraction * (1 - case_fraction) * 2 * raf * (1 - raf))

# Gain of the realised score in each of R replicates (rows: replicate x threshold)
sweep_cell_realised <- function(q, lp_profile, n_loci, odds_ratio, raf, sigma, thresholds, R = 12L, N = 2000L, seed = 1L) {
  # The global-intercept construction only; calibration_gains_realised() (R/calibration_check.R) does all three comparisons
  if (genetic_scheme() != "global") stop("sweep_cell_realised() builds the global-intercept truth; use calibration_gains_realised()", call. = FALSE)
  set.seed(seed)
  beta <- log(odds_ratio)
  pmf <- dosage_pmf(n_loci, raf)
  grid <- list(G = matrix(seq_along(pmf) - 1, ncol = 1, dimnames = list(NULL, "D")), prob = pmf)
  cells <- genetic_cells(q, lp_profile, grid)
  a_true <- genetic_scores(cells, c(D = beta), c(D = beta))$a_true
  base_nb <- nb_from(tp_fp(cells$w, stats::plogis(cells$lp_base + cells$G[, 1] * beta + a_true), stats::plogis(cells$lp_base), thresholds), thresholds)
  npro <- length(q)
  wt <- rep(q, times = N) / N
  lpb <- rep(lp_profile, times = N)
  gains <- vapply(seq_len(R), function(r) {
    G <- matrix(stats::rbinom(N * n_loci, 2, raf), N, n_loci)
    D <- rowSums(G); S <- drop(G %*% (beta + stats::rnorm(n_loci, 0, sigma)))
    rt <- stats::plogis(lpb + beta * rep(D, each = npro) + a_true)
    lps <- lpb + rep(S, each = npro)
    lps <- lps + recentre_shift(lps, wt, cells$pred_mean)
    nb_from(tp_fp(wt, rt, stats::plogis(lps), thresholds), thresholds) - base_nb
  }, numeric(length(thresholds)))
  gains <- matrix(gains, nrow = length(thresholds))
  data.frame(n_loci = n_loci, odds_ratio = odds_ratio, raf = raf, sigma = sigma, replicate = rep(seq_len(R), each = length(thresholds)),
             threshold = rep(thresholds, times = R), gain = as.numeric(gains))
}
