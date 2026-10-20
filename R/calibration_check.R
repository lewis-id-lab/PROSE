# Does adding genetics preserve the clinical baseline fairly? (addendum 24)
#
# The genetic layer (D5) adds beta * g to the clinical logit and moves one global intercept to hold the overall prevalence:
#   true logit = eta_X + beta * g + a.
# Within a clinical profile X the true risk averaged over genotype, m(X) = sum_g p_g expit(eta_X + beta * g + a), is not expit(eta_X)
# (Jensen's inequality), so the clinical score is miscalibrated profile by profile under the simulated truth even though the overall
# prevalence is right. The best clinical-only rule at threshold t flags m(X) >= t, so NB(m) >= NB(expit(eta)) under that truth, and
# part of the reported gain could be the correction of an induced miscalibration rather than genetic information.
#
# Three comparisons, each on profile x genotype cells (profile varying fastest, as in genetic_cells()):
#   current     global-intercept truth; clinical comparator expit(eta_X); genetic score recentred once (as genetic_scores()).
#   marginal    the same truth; clinical comparator the actual clinical-only risk m(X) under it.
#   profile     profile-preserving truth: an intercept a_X per profile so that the mean risk over genotype is expit(eta_X) in every
#               profile, with the conditional genetic effects unchanged; clinical comparator expit(eta_X), which is then exactly
#               calibrated; the genetic score is built the same way (an intercept per profile that preserves its clinical risk).
# In each, "oracle" is the true risk as the score (the ceiling score when the score's weights are the true ones), and "score" is a
# score carrying beta_score (the realised score).
#
# Requires R/decision_curve.R (dca_exact) and R/genetic_layer.R (recentre_shift).

# Per-profile shifts a_X with sum_g p_g expit(lp_X + gl_g + a_X) = target_X, by vectorised Newton steps
profile_shifts <- function(lp, gl, p, target = stats::plogis(lp), tol = 1e-12, max_iter = 200L) {
  a <- rep(-sum(p * gl), length(lp))
  for (it in seq_len(max_iter)) {
    r <- stats::plogis(outer(lp + a, gl, "+"))
    f <- drop(r %*% p) - target
    d <- drop((r * (1 - r)) %*% p)
    step <- pmax(pmin(f / d, 2), -2)
    a <- a - step
    if (max(abs(step)) < tol) break
  }
  a
}

# True-risk matrix (profiles x genotypes) under one scheme for genetic linear predictor gl
scheme_risk <- function(q, lp, gl, p, scheme = c("global", "profile")) {
  scheme <- match.arg(scheme)
  if (scheme == "global") {
    W <- outer(q, p)
    a <- recentre_shift(as.numeric(outer(lp, gl, "+")), as.numeric(W), sum(q * stats::plogis(lp)))
    stats::plogis(outer(lp, gl, "+") + a)
  } else {
    stats::plogis(outer(lp + profile_shifts(lp, gl, p), gl, "+"))
  }
}

# Weighted logistic recalibration of true probabilities y on predicted probabilities x: slope, and calibration-in-the-large
# (intercept with the predicted logit as an offset). Exact for an infinite population.
calibration_fit <- function(y, x, w) {
  keep <- w > 0
  y <- y[keep]; lx <- stats::qlogis(x[keep]); w <- w[keep] / sum(w[keep])
  fam <- stats::quasibinomial()
  slope <- if (stats::sd(lx) < 1e-12) NA_real_ else stats::glm.fit(cbind(1, lx), y, weights = w, family = fam)$coefficients[2]
  citl <- stats::glm.fit(matrix(1, length(y)), y, weights = w, offset = lx, family = fam)$coefficients[1]
  c(slope = unname(slope), citl = unname(citl))
}

# Gains under the three comparisons for one set of true effects and score weights.
#   grid  list(G = genotype matrix with named columns, prob = genotype probabilities), independent of the profile
calibration_gains <- function(q, lp, grid, beta_true, beta_score = beta_true, thresholds,
                              comparisons = c("current", "marginal", "profile"), calibration = FALSE) {
  p <- grid$prob
  gl_t <- drop(grid$G %*% beta_true[colnames(grid$G)])
  gl_s <- drop(grid$G %*% beta_score[colnames(grid$G)])
  W <- outer(q, p)
  np <- length(q); ng <- length(p)
  clin <- stats::plogis(lp)
  nb <- function(R, S) dca_exact(as.numeric(W), as.numeric(R), as.numeric(S), thresholds = thresholds)
  one <- function(comparison, R, C, S) {
    m_true <- drop(R %*% p)
    d_c <- nb(R, matrix(C, np, ng)); d_o <- nb(R, R); d_s <- nb(R, S)
    out <- data.frame(comparison = comparison, threshold = thresholds, nb_clinical = d_c$nb,
                      gain_oracle = d_o$nb - d_c$nb, gain_score = d_s$nb - d_c$nb,
                      tp_clinical = d_c$tp, extra_tp_oracle = d_o$tp - d_c$tp,
                      flag_disagree = vapply(thresholds, function(t) sum(q[(C >= t) != (m_true >= t)]), numeric(1)),
                      stringsAsFactors = FALSE)
    if (calibration) {
      cc <- calibration_fit(m_true, C, q)
      sc <- calibration_fit(as.numeric(R), as.numeric(S), as.numeric(W))
      out$clin_slope <- cc[["slope"]]; out$clin_citl <- cc[["citl"]]
      out$clin_mean_abs_err <- sum(q * abs(m_true - C)); out$clin_max_abs_err <- max(abs(m_true - C)[q > 0])
      out$score_slope <- sc[["slope"]]; out$score_citl <- sc[["citl"]]
    }
    out
  }
  res <- list()
  if (any(c("current", "marginal") %in% comparisons)) {
    R0 <- scheme_risk(q, lp, gl_t, p, "global")
    S0 <- scheme_risk(q, lp, gl_s, p, "global")
    if ("current" %in% comparisons) res$current <- one("current", R0, clin, S0)
    if ("marginal" %in% comparisons) res$marginal <- one("marginal", R0, drop(R0 %*% p), S0)
  }
  if ("profile" %in% comparisons) {
    RA <- scheme_risk(q, lp, gl_t, p, "profile")
    SA <- scheme_risk(q, lp, gl_s, p, "profile")
    res$profile <- one("profile", RA, clin, SA)
  }
  do.call(rbind, unname(res))
}
