# Aim 3b: the genetic layer split by neutropenia status (addendum 15). Requires R/baseline.R, R/decision_curve.R and R/genetic_layer.R.
#
# The strata are the neutropenic and non-neutropenic admissions of the base cohort (the prolonged-neutropenia predictor). The true per-copy
# effects can differ between them; the event-weighted average of the two stratum effects is held at the pooled effect, so the whole-cohort
# evidence stays consistent: w beta_N + (1 - w) beta_NN = beta, with w the share of the baseline's events in the neutropenic stratum.
# The modification ratio is r = beta_N / beta_NN (infinity: all of the effect in neutropenic admissions; 0: all in non-neutropenic ones).

# Profile distribution, share and baseline event rate of each stratum; w is the neutropenic stratum's share of the baseline's events
neutropenia_strata <- function(q, lp_profile, profiles = all_profiles()) {
  neut <- profiles[, "prolonged_neutropenia"] == 1
  risk <- stats::plogis(lp_profile)
  one <- function(sel) {
    qs <- ifelse(sel, q, 0)
    list(q = qs / sum(qs), share = sum(qs), rate = sum(qs * risk) / sum(qs), events = sum(qs * risk))
  }
  N <- one(neut); NN <- one(!neut)
  list(N = N, NN = NN, w = N$events / (N$events + NN$events))
}

# Stratum effects with the pooled effect held as the event-weighted average
modify_effects <- function(beta, r, w) {
  if (is.infinite(r)) return(list(N = beta / w, NN = beta * 0))
  denom <- r * w + 1 - w
  if (denom <= 0) stop("ratio ", r, " is not compatible with the pooled effect and an event share of ", round(w, 3), call. = FALSE)
  nn <- beta / denom
  list(N = r * nn, NN = nn)
}

# A scenario is either a modification ratio (the pooled effect held as the event-weighted average) or stratum fractions of the pooled effect,
# c(N = f_N, NN = f_NN), applied directly (used for the data-based Brunel scenario, whose stratum effects cannot be reconciled with the
# Bologna event mix: a ratio of -0.11 would need an effect nine times the pooled one in the non-neutropenic stratum)
scenario_effects <- function(beta, spec, w) {
  if (length(spec) == 2L && !is.null(names(spec))) return(list(N = beta * spec[["N"]], NN = beta * spec[["NN"]]))
  modify_effects(beta, spec, w)
}

# Difference of the two stratum fractions (the true interaction in units of the pooled-estimate score)
scenario_difference <- function(spec, w) {
  if (length(spec) == 2L && !is.null(names(spec))) return(spec[["N"]] - spec[["NN"]])
  if (is.infinite(spec)) return(1 / w)
  f_nn <- 1 / (spec * w + 1 - w)
  (spec - 1) * f_nn
}

# Drop the cells that carry no weight (profiles outside the stratum)
trim_cells <- function(cells) {
  k <- cells$w > 0
  cells$w <- cells$w[k]; cells$lp_base <- cells$lp_base[k]; cells$G <- cells$G[k, , drop = FALSE]; cells$profile <- cells$profile[k]
  cells
}

# True positives, false positives and net benefit of several fixed scores under one true model, in one stratum's cells. The truth and
# each score carry their own per-copy effects and are recentred to the stratum's baseline event rate (D5); an all-zero score is the baseline.
#   score_betas: named list of effect vectors (named by term)
multi_score_nb <- function(cells, beta_true, score_betas, thresholds, scheme = genetic_scheme()) {
  w <- cells$w; terms <- colnames(cells$G)
  clin <- stats::plogis(cells$lp_base)
  lp_t <- cells$lp_base + drop(cells$G %*% beta_true[terms])
  r_true <- stats::plogis(lp_t + cell_shift(lp_t, w, cells$profile, clin, cells$prevalence, scheme))
  do.call(rbind, lapply(names(score_betas), function(s) {
    b <- score_betas[[s]][terms]
    lp <- cells$lp_base
    if (any(b != 0)) { lp <- lp + drop(cells$G %*% b); lp <- lp + cell_shift(lp, w, cells$profile, clin, cells$pred_mean, scheme) }
    x <- tp_fp(w, r_true, stats::plogis(lp), thresholds)
    data.frame(score = s, threshold = thresholds, tp = x[, "tp"], fp = x[, "fp"], nb = nb_from(x, thresholds), stringsAsFactors = FALSE)
  }))
}

# The whole run: for each draw of the pooled effects and each modification scenario, the net-benefit gain of each score over the baseline in
# each stratum (per stratum admission) and in the whole cohort. `cells` is list(N =, NN =) of trimmed cells, `strata` from neutropenia_strata().
#   r_values: named list of scenarios (ratios, or c(N =, NN =) fractions; see scenario_effects()); assumed: the two ratios used by the assumed-H_F and assumed-Brunel scores
run_aim3b <- function(strata, cells, effects, beta_hat, r_values, thresholds, assumed = c(feys = 4, brunel = 0.25), draw_id = seq_len(nrow(effects))) {
  w <- strata$w
  af <- modify_effects(beta_hat, assumed[["feys"]], w); ab <- modify_effects(beta_hat, assumed[["brunel"]], w)
  zero <- beta_hat * 0
  out <- vector("list", nrow(effects) * length(r_values)); k <- 0L
  for (i in seq_len(nrow(effects))) {
    for (sc in names(r_values)) {
      tr <- scenario_effects(effects[i, ], r_values[[sc]], w)
      per <- lapply(c("N", "NN"), function(s) {
        scores <- list(baseline = zero, ceiling = tr[[s]], common = beta_hat,
                       assumed_feys = af[[s]], assumed_brunel = ab[[s]])
        multi_score_nb(cells[[s]], tr[[s]], scores, thresholds)
      })
      names(per) <- c("N", "NN")
      key <- paste(per$N$score, per$N$threshold)
      whole <- per$N; whole$nb <- strata$N$share * per$N$nb + strata$NN$share * per$NN$nb
      whole$tp <- strata$N$share * per$N$tp + strata$NN$share * per$NN$tp
      whole$fp <- strata$N$share * per$N$fp + strata$NN$share * per$NN$fp
      lev <- list(N = per$N, NN = per$NN, whole = whole)
      d <- do.call(rbind, lapply(names(lev), function(l) {
        x <- lev[[l]]
        b <- x[x$score == "baseline", ]
        j <- match(x$threshold, b$threshold)
        data.frame(level = l, score = x$score, threshold = x$threshold, nb = x$nb, delta_nb = x$nb - b$nb[j], delta_tp = x$tp - b$tp[j], stringsAsFactors = FALSE)
      }))
      k <- k + 1L
      out[[k]] <- cbind(draw = draw_id[i], scenario = sc, d)
    }
  }
  do.call(rbind, out)
}

# MCID of the baseline in a stratum: 20% more true cases at no extra false positives, per threshold
stratum_mcid <- function(strata_one, lp_profile, thresholds) baseline_mcid(strata_one$q, lp_profile, thresholds)

# One simulated study of the cohort's size with stratum-specific true effects, and the interaction fit that tests whether the effect is
# larger in neutropenic admissions: outcome ~ N + S + S:N + offset(baseline logit), with S the pooled-estimate genetic score.
#   shifts: list(N =, NN =) recentring shifts of the truth (genetic_scores()$a_true of each stratum's cells)
#   sampler: optional function(n) returning an n x terms matrix of dosages drawn jointly (for terms that are not independent binomials, such as the two CLEC7A parameters)
simulate_interaction_study <- function(cg, neut, lp_profile, p, shifts, beta_true_N, beta_true_NN, beta_hat, resample = TRUE, sampler = NULL) {
  by <- cg$groups
  idx <- if (resample) sample.int(length(by), replace = TRUE) else seq_along(by)
  rows <- unlist(by[idx], use.names = FALSE)
  pat <- rep(seq_along(idx), lengths(by[idx]))
  g <- if (is.null(sampler)) matrix(vapply(p, function(pj) stats::rbinom(length(idx), 2, pj), numeric(length(idx))), nrow = length(idx), dimnames = list(NULL, names(p))) else sampler(length(idx))
  G <- g[pat, , drop = FALSE]
  n <- neut[rows]
  prof <- cg$prof[rows]
  lp_b <- lp_profile[prof]
  gl <- ifelse(n, drop(G %*% beta_true_N[colnames(g)]), drop(G %*% beta_true_NN[colnames(g)]))
  y <- stats::rbinom(length(rows), 1, stats::plogis(lp_b + gl + ifelse(n, sh_at(shifts$N, prof), sh_at(shifts$NN, prof))))
  S <- drop(G %*% beta_hat[colnames(g)])
  # The neutropenia main effect is in the model: the truth is recentred within each stratum, so the strata differ in intercept
  nn <- as.numeric(n); SN <- S * nn
  fit <- tryCatch(suppressWarnings(stats::glm(y ~ nn + S + SN + offset(lp_b), family = stats::binomial)), error = function(e) NULL)
  if (is.null(fit) || !fit$converged) return(c(delta = NA_real_, se = NA_real_, p = NA_real_, events = sum(y), events_neutropenic = sum(y[n])))
  cf <- summary(fit)$coefficients
  c(delta = cf["SN", 1], se = cf["SN", 2], p = cf["SN", 4], events = sum(y), events_neutropenic = sum(y[n]))
}
