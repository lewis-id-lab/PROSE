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
# Requires R/decision_curve.R (dca_exact) and R/genetic_layer.R (recentre_shift, tp_fp, nb_from). calibration_gains_realised also needs
# R/aim4_sweep.R (dosage_pmf); nri_by_scheme needs R/nri.R (reclassification_metrics).

# Per-profile shifts a_X with sum_g p_g expit(lp_X + gl_g + a_X) = target_X, by vectorised Newton steps
profile_shifts <- function(lp, gl, p, target = stats::plogis(lp), tol = 1e-12, max_iter = 200L) {
  a <- rep(-sum(p * gl), length(lp))
  for (it in seq_len(max_iter)) {
    r <- stats::plogis(outer(lp + a, gl, "+"))
    f <- drop(r %*% p) - target
    # Floor keeps a saturated profile (every risk rounds to 0 or 1) from giving NaN; the step is clamped below
    d <- pmax(drop((r * (1 - r)) %*% p), 1e-300)
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
  # Clamp so that risks that round to exactly 0 or 1 (very large genetic logits) do not give an infinite logit
  eps <- 1e-12
  y <- y[keep]; lx <- stats::qlogis(pmin(pmax(x[keep], eps), 1 - eps)); w <- w[keep] / sum(w[keep])
  fam <- stats::quasibinomial()
  slope <- if (stats::sd(lx) < 1e-12) NA_real_ else stats::glm.fit(cbind(1, lx), y, weights = w, family = fam)$coefficients[2]
  citl <- stats::glm.fit(matrix(1, length(y)), y, weights = w, offset = lx, family = fam)$coefficients[1]
  c(slope = unname(slope), citl = unname(citl))
}

# Gains under the three comparisons for one set of true effects and score weights.
#   grid  list(G = genotype matrix with named columns, prob = genotype probabilities), independent of the profile
#   calibration  add calibration slope and calibration-in-the-large of the clinical and genetic scores
#   detail       add the score's extra true cases, the flagged shares and the expected aROC of the clinical, ceiling and realised scores
#                (needs expected_auc from R/baseline.R); fractions of the cohort, as in dca_exact
calibration_gains <- function(q, lp, grid, beta_true, beta_score = beta_true, thresholds,
                              comparisons = c("current", "marginal", "profile"), calibration = FALSE, detail = FALSE) {
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
    if (detail) {
      Cm <- matrix(C, np, ng)
      auc <- function(P) expected_auc(stats::qlogis(pmin(pmax(as.numeric(P), 1e-12), 1 - 1e-12)), as.numeric(W * R), as.numeric(W * (1 - R)))
      out$extra_tp_score <- d_s$tp - d_c$tp
      out$flagged_clinical <- d_c$flagged; out$flagged_oracle <- d_o$flagged; out$flagged_score <- d_s$flagged
      out$auc_clinical <- auc(Cm); out$auc_oracle <- auc(R); out$auc_score <- auc(S)
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

# Gains of a realised (estimated-weight) score in each of R replicate discovery studies, under the three comparisons (the analogue of
# sweep_cell_realised in R/aim4_sweep.R, which it reproduces for comparison "current" at the same seed). The genotypes and weight errors of
# a replicate are shared by the comparisons, so they are paired. The truth and the baseline are exact; the score is simulated on N
# genotypes crossed with the exact profile distribution.
#   current   global-intercept truth; comparator expit(eta_X); score recentred once
#   marginal  the same truth; comparator m(X) = E[p | X] under it
#   profile   profile-preserving truth; comparator expit(eta_X); the score gets an intercept per profile that holds its mean risk over
#             the N simulated genotypes at expit(eta_X)
calibration_gains_realised <- function(q, lp, n_loci, odds_ratio, raf, sigma, thresholds,
                                       comparisons = c("current", "marginal", "profile"), R = 12L, N = 2000L, seed = 1L) {
  set.seed(seed)
  beta <- log(odds_ratio)
  pmf <- dosage_pmf(n_loci, raf)
  gl <- beta * (seq_along(pmf) - 1)
  np <- length(q); ng <- length(pmf)
  clin <- stats::plogis(lp)
  Wc <- as.numeric(outer(q, pmf))
  need_global <- any(c("current", "marginal") %in% comparisons)
  truth <- list(global = if (need_global) scheme_risk(q, lp, gl, pmf, "global"),
                profile = if ("profile" %in% comparisons) scheme_risk(q, lp, gl, pmf, "profile"))
  scheme_of <- c(current = "global", marginal = "global", profile = "profile")
  comparator <- list(current = clin, marginal = if (need_global) drop(truth$global %*% pmf), profile = clin)
  base_nb <- lapply(setNames(comparisons, comparisons), function(cmp) {
    nb_from(tp_fp(Wc, as.numeric(truth[[scheme_of[[cmp]]]]), as.numeric(matrix(comparator[[cmp]], np, ng)), thresholds), thresholds)
  })
  wt <- rep(q, times = N) / N
  lpb <- rep(lp, times = N)
  gains <- lapply(setNames(comparisons, comparisons), function(cmp) matrix(NA_real_, length(thresholds), R))
  for (r in seq_len(R)) {
    G <- matrix(stats::rbinom(N * n_loci, 2, raf), N, n_loci)
    D <- rowSums(G); S <- drop(G %*% (beta + stats::rnorm(n_loci, 0, sigma)))
    if (need_global) {
      lps <- lpb + rep(S, each = np)
      pred_global <- stats::plogis(lps + recentre_shift(lps, wt, sum(q * clin)))
    }
    if ("profile" %in% comparisons) {
      pred_profile <- as.numeric(stats::plogis(outer(lp + profile_shifts(lp, S, rep(1 / N, N)), S, "+")))
    }
    for (cmp in comparisons) {
      rt <- as.numeric(truth[[scheme_of[[cmp]]]][, D + 1L])
      pred <- if (cmp == "profile") pred_profile else pred_global
      gains[[cmp]][, r] <- nb_from(tp_fp(wt, rt, pred, thresholds), thresholds) - base_nb[[cmp]]
    }
  }
  do.call(rbind, lapply(comparisons, function(cmp) {
    data.frame(comparison = cmp, n_loci = n_loci, odds_ratio = odds_ratio, raf = raf, sigma = sigma,
               replicate = rep(seq_len(R), each = length(thresholds)), threshold = rep(thresholds, times = R),
               gain = as.numeric(gains[[cmp]]), stringsAsFactors = FALSE)
  }))
}

# Reclassification metrics (R/nri.R) of a genetic score against the clinical comparator under each comparison, for one genotype grid.
# Returns a matrix with one row per comparison. Truth, score and comparator are built as in calibration_gains.
nri_by_scheme <- function(q, lp, grid, beta_true, beta_score = beta_true, comparisons = c("current", "marginal", "profile"),
                          thresholds = c(0.05, 0.20)) {
  p <- grid$prob
  gl_t <- drop(grid$G %*% beta_true[colnames(grid$G)])
  gl_s <- drop(grid$G %*% beta_score[colnames(grid$G)])
  W <- as.numeric(outer(q, p))
  np <- length(q); ng <- length(p)
  clin <- matrix(stats::plogis(lp), np, ng)
  if (any(c("current", "marginal") %in% comparisons)) { Rg <- scheme_risk(q, lp, gl_t, p, "global"); Sg <- scheme_risk(q, lp, gl_s, p, "global") }
  if ("profile" %in% comparisons) { Rp <- scheme_risk(q, lp, gl_t, p, "profile"); Sp <- scheme_risk(q, lp, gl_s, p, "profile") }
  lg <- function(x) stats::qlogis(as.numeric(x))
  out <- lapply(comparisons, function(cmp) {
    Rm <- if (cmp == "profile") Rp else Rg
    Sm <- if (cmp == "profile") Sp else Sg
    Cm <- if (cmp == "marginal") matrix(drop(Rg %*% p), np, ng) else clin
    reclassification_metrics(W, lg(Rm), lg(Cm), lg(Sm), thresholds)
  })
  out <- do.call(rbind, out)
  rownames(out) <- comparisons
  out
}
