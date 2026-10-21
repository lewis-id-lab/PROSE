# Gains of the ceiling and the realised score in the Aim 3b strata under the global or the profile-preserving truth, against the clinical score and
# against the best clinical-only rule m(X) (addenda 27 and 28). Requires R/genetic_layer.R (cell_shift, tp_fp, nb_from), R/decision_curve.R and
# R/calibration_check.R (calibration_fit). The scheme is an argument here, not the environment variable, so both can be computed side by side.

# Mean of a cell-level risk within each clinical profile, mapped back to the cells: the marginal risk m(X) of the cells' profile
profile_mean <- function(x, cl) as.numeric((tapply(cl$w * x, cl$profile, sum) / tapply(cl$w, cl$profile, sum))[as.character(cl$profile)])

# One stratum (trimmed cells `cl`): the truth has per-copy effects b; the ceiling score carries b and the realised score carries bs. Returns the net-benefit
# gain of each over the clinical score and over m(X), one value per threshold, and the clinical score's calibration slope within the stratum.
stratum_gains <- function(cl, b, bs, scheme, th) {
  nbf <- function(r, pred) nb_from(tp_fp(cl$w, r, pred, th), th)
  clin <- stats::plogis(cl$lp_base)
  lp_t <- cl$lp_base + drop(cl$G %*% b[colnames(cl$G)])
  r <- stats::plogis(lp_t + cell_shift(lp_t, cl$w, cl$profile, clin, cl$prevalence, scheme))
  lp_s <- cl$lp_base + drop(cl$G %*% bs[colnames(cl$G)])
  s <- stats::plogis(lp_s + cell_shift(lp_s, cl$w, cl$profile, clin, cl$pred_mean, scheme))
  m <- profile_mean(r, cl)
  nb_clin <- nbf(r, clin); nb_marg <- nbf(r, m); nb_ceil <- nbf(r, r); nb_real <- nbf(r, s)
  first <- !duplicated(cl$profile); o <- order(cl$profile[first])
  cal <- calibration_fit(m[first][o], clin[first][o], as.numeric(tapply(cl$w, cl$profile, sum)))
  list(ceil_clin = nb_ceil - nb_clin, ceil_marg = nb_ceil - nb_marg, real_clin = nb_real - nb_clin, real_marg = nb_real - nb_marg, slope = unname(cal["slope"]))
}

# The whole cohort from the two strata: `cells` is list(N =, NN =) of trimmed cells, `share` the strata's shares of admissions, `tr` list(N =, NN =) of
# per-copy effects in the truth, `bs` the weights the realised score carries (the same in both strata)
whole_gains <- function(cells, share, tr, bs, scheme, th) {
  g <- list(N = stratum_gains(cells$N, tr$N, bs, scheme, th), NN = stratum_gains(cells$NN, tr$NN, bs, scheme, th))
  comb <- function(k) share[["N"]] * g$N[[k]] + share[["NN"]] * g$NN[[k]]
  list(ceil_clin = comb("ceil_clin"), ceil_marg = comb("ceil_marg"), real_clin = comb("real_clin"), real_marg = comb("real_marg"),
       slope_N = g$N$slope, slope_NN = g$NN$slope)
}
