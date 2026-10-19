# Polygenic score benchmark (addendum 23): how strong a polygenic score would have to be, and how large its discovery study, for
# the score to reach the MCID on the refit baseline.
#
# Model. A polygenic score G ~ N(0, 1), independent of the clinical profile (an assumption, as for the panel), with true effect
# beta on the log-odds scale per standard deviation:
#   true logit = baseline logit + beta * G + a,   a chosen so that prevalence is unchanged (D5).
# Ceiling: the score knows G and beta. Realised: the score is an estimate S from a discovery study, with S ~ N(0, 1) and
# correlation rho with G, so G = rho * S + sqrt(1 - rho^2) * U. It carries the weight beta * rho (the regression of G on S),
# recentred to the baseline mean risk: the best linear use of an imperfect score, recalibrated in the target cohort, which is
# optimistic. Everything is exact over the 128 clinical profiles by equal-probability normal quantile nodes, with the truth
# integrated over the part of G the score does not capture.
#
# Discovery size. For M_e independent causal loci of equal small effect, each estimated by logistic regression in n subjects with
# case fraction pi, the accuracy of the summed score is
#   rho^2 = n pi (1 - pi) beta^2 / (n pi (1 - pi) beta^2 + M_e)
# (Daetwyler et al. 2008, written on the log-odds scale; the per-locus precision is the one in sigma_for_sample(), R/aim4_sweep.R).
# Loci are taken as known in advance, so there is no selection and no winner's curse (optimistic).
#
# Requires R/baseline.R (expected_auc), R/decision_curve.R (dca_exact) and R/genetic_layer.R (recentre_shift).

# Equal-probability nodes of the standard normal
prs_nodes <- function(n) stats::qnorm(stats::ppoints(n))

# The shift a that holds the prevalence for a true polygenic effect beta (depends on beta only, not on the score's accuracy)
prs_shift <- function(q, lp, beta, prevalence = sum(q * stats::plogis(lp)), z = prs_nodes(400)) {
  w <- as.numeric(outer(q, rep(1 / length(z), length(z))))
  recentre_shift(as.numeric(outer(lp, beta * z, "+")), w, prevalence)
}

# Net-benefit gain over the baseline of a polygenic score with accuracy rho (rho = 1 is the ceiling), with the expected aROC of the
# score alone and of the combined score, under the true model.
prs_gain <- function(q, lp, beta, rho, thresholds, a = NULL, s = prs_nodes(400), u = prs_nodes(64)) {
  stopifnot(rho >= 0, rho <= 1)
  pred_mean <- sum(q * stats::plogis(lp))
  if (is.null(a)) a <- prs_shift(q, lp, beta, pred_mean)
  ns <- length(s)
  lp_c <- rep(lp, times = ns)
  w <- rep(q, times = ns) / ns
  g_known <- beta * rho * rep(s, each = length(lp))
  sd_u <- beta * sqrt(1 - rho^2)
  r_true <- if (sd_u == 0) stats::plogis(lp_c + g_known + a) else
    rowMeans(stats::plogis(outer(lp_c + g_known + a, sd_u * u, "+")))
  lp_score <- lp_c + g_known
  lp_score <- lp_score + recentre_shift(lp_score, w, pred_mean)
  d1 <- dca_exact(w, r_true, stats::plogis(lp_score), thresholds = thresholds)
  d0 <- dca_exact(w, r_true, stats::plogis(lp_c), thresholds = thresholds)
  data.frame(threshold = thresholds, gain = d1$nb - d0$nb,
             auc_prs = expected_auc(g_known, w * r_true, w * (1 - r_true)),
             auc_baseline = expected_auc(lp_c, w * r_true, w * (1 - r_true)),
             auc_combined = expected_auc(lp_score, w * r_true, w * (1 - r_true)))
}

# Smallest accuracy rho at which the score reaches the MCID at one threshold; NA if even the ceiling (rho = 1) falls short
prs_rho_needed <- function(q, lp, beta, mcid, threshold, a = prs_shift(q, lp, beta)) {
  f <- function(rho) prs_gain(q, lp, beta, rho, threshold, a = a)$gain - mcid
  if (f(1) < 0) return(NA_real_)
  stats::uniroot(f, c(0, 1), tol = 1e-5)$root
}

# Effective discovery size n pi (1 - pi) needed for accuracy rho, and the accuracy a given discovery study achieves
prs_n_effective <- function(rho, beta, m_eff) m_eff * rho^2 / (beta^2 * (1 - rho^2))
prs_accuracy <- function(n, case_fraction, beta, m_eff) {
  ne <- n * case_fraction * (1 - case_fraction)
  sqrt(ne * beta^2 / (ne * beta^2 + m_eff))
}
