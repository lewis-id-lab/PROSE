# Aim 4: the pathogen-specificity scenario (scope decision 1; not in addenda 16 and 17). Requires R/baseline.R, R/decision_curve.R, R/genetic_layer.R and R/aim4_sweep.R.
#
# The sweep of script 23 treats every case of invasive mould disease as genotype-sensitive. Here only a fraction f of cases (the target pathogen, for example
# Aspergillus) arises through the pathway the loci act on; the other 1 - f arise independently of the genotype. The outcome is still "any proven or probable IMD".
#
#   true risk = (1 - f) * plogis(lp_base) + f * plogis(lp_base + log(OR) D + a)
#
# a is the shift that holds the target-pathogen component, and so the total event rate, at its baseline value (D5), so f = 1 reproduces script 23 exactly.
# Two scores are compared, both recentred to the baseline score's mean risk:
#   specific  carries the target-pathogen effect log(OR) (what a pathogen-specific discovery study would estimate), applied to every admission
#   diluted   carries the effect fitted to the all-mould outcome (what a discovery study on any IMD would estimate): the best single-weight score for this truth,
#             found by quasi-binomial regression of the true risk on the baseline logit (offset) and D
# The MCID is the baseline's, as in script 23 (20% more true cases at no extra false positives of the baseline under its own predicted risk).

pathogen_risk <- function(lp_base, gl, w, prevalence, f) {
  a <- recentre_shift(lp_base + gl, w, prevalence)
  (1 - f) * stats::plogis(lp_base) + f * stats::plogis(lp_base + gl + a)
}

# Effect of D on the all-mould outcome: the weighted fractional-logit fit with the baseline logit as an offset
diluted_beta <- function(lp_base, D, w, r_true) {
  if (all(D == D[1])) return(0)
  fit <- suppressWarnings(stats::glm(r_true ~ D + offset(lp_base), family = stats::quasibinomial(), weights = w, start = c(0, 0)))
  unname(stats::coef(fit)["D"])
}

pathogen_cell <- function(q, lp_profile, n_loci, odds_ratio, raf, f, thresholds) {
  pmf <- dosage_pmf(n_loci, raf)
  cells <- genetic_cells(q, lp_profile, list(G = matrix(seq_along(pmf) - 1, ncol = 1, dimnames = list(NULL, "D")), prob = pmf))
  D <- cells$G[, 1]; beta <- log(odds_ratio)
  r_true <- pathogen_risk(cells$lp_base, D * beta, cells$w, cells$prevalence, f)
  b_dil <- diluted_beta(cells$lp_base, D, cells$w, r_true)
  score <- function(b) { lp <- cells$lp_base + b * D; lp + recentre_shift(lp, cells$w, cells$pred_mean) }
  one <- function(lp) { r <- tp_fp(cells$w, r_true, stats::plogis(lp), thresholds); list(r = r, nb = nb_from(r, thresholds)) }
  base <- one(cells$lp_base); spec <- one(score(beta)); dil <- one(score(b_dil))
  data.frame(n_loci = n_loci, odds_ratio = odds_ratio, raf = raf, case_fraction = f, threshold = thresholds,
             beta_diluted = b_dil, or_diluted = exp(b_dil), baseline_nb = base$nb,
             gain_specific = spec$nb - base$nb, gain_diluted = dil$nb - base$nb,
             extra_cases_per_1000_specific = 1000 * (spec$r[, "tp"] - base$r[, "tp"]),
             extra_cases_per_1000_diluted = 1000 * (dil$r[, "tp"] - base$r[, "tp"]))
}

# ---- With estimation error in the weights (addendum 39; the counterpart of sweep_cell_realised in R/aim4_sweep.R) -------------------------------
# The score's weights are estimates with an independent error of standard deviation sigma on each locus (no winner's curse), around one of two centres:
#   specific  the target-pathogen effect log(OR)         (a discovery study that knows each case's pathogen)
#   all-mould the effect fitted to the all-mould outcome (a discovery study on any proven or probable mould disease)
# The truth is the pathogen mixture of pathogen_risk(); genotypes are simulated (N individuals crossed with the exact profile distribution), R replicate
# discovery studies, the same simulated genotypes and errors for both centres. Returns one row per replicate, score and threshold.
pathogen_cell_realised <- function(q, lp_profile, n_loci, odds_ratio, raf, f, sigma, thresholds, R = 12L, N = 2000L, seed = 1L) {
  set.seed(seed)
  beta <- log(odds_ratio)
  pmf <- dosage_pmf(n_loci, raf)
  cells <- genetic_cells(q, lp_profile, list(G = matrix(seq_along(pmf) - 1, ncol = 1, dimnames = list(NULL, "D")), prob = pmf))
  D0 <- cells$G[, 1]
  a <- recentre_shift(cells$lp_base + beta * D0, cells$w, cells$prevalence)
  r_cells <- (1 - f) * stats::plogis(cells$lp_base) + f * stats::plogis(cells$lp_base + beta * D0 + a)
  centres <- c(specific = beta, all_mould = diluted_beta(cells$lp_base, D0, cells$w, r_cells))
  base_nb <- nb_from(tp_fp(cells$w, r_cells, stats::plogis(cells$lp_base), thresholds), thresholds)
  npro <- length(q)
  wt <- rep(q, times = N) / N
  lpb <- rep(lp_profile, times = N)
  out <- lapply(seq_len(R), function(r) {
    G <- matrix(stats::rbinom(N * n_loci, 2, raf), N, n_loci)
    D <- rowSums(G)
    eps <- stats::rnorm(n_loci, 0, sigma)
    rt <- (1 - f) * stats::plogis(lpb) + f * stats::plogis(lpb + beta * rep(D, each = npro) + a)
    do.call(rbind, lapply(names(centres), function(sc) {
      lps <- lpb + rep(drop(G %*% (centres[[sc]] + eps)), each = npro)
      lps <- lps + recentre_shift(lps, wt, cells$pred_mean)
      data.frame(replicate = r, score = sc, threshold = thresholds, gain = nb_from(tp_fp(wt, rt, stats::plogis(lps), thresholds), thresholds) - base_nb)
    }))
  })
  cbind(n_loci = n_loci, odds_ratio = odds_ratio, raf = raf, case_fraction = f, sigma = sigma, do.call(rbind, out))
}
