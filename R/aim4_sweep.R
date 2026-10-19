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
