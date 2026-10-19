# Five-term Aim 3b and Shapley attribution for one profile distribution (a sensitivity cohort or a correlation-ensemble member), addendum 18.
# Requires R/baseline.R, R/decision_curve.R, R/genetic_layer.R and R/aim3b.R.

# Aim 3b rows (run_aim3b) for the draws `idx` of `eff` in the profile distribution q, with the MCID of the baseline in that distribution
# (whole cohort and each stratum) and the neutropenic stratum's share of admissions attached to every row.
aim3b_member_rows <- function(q, lp, grid, eff, idx, beta_hat, r_values, thresholds, assumed = c(feys = 4, brunel = 0.25)) {
  st <- neutropenia_strata(q, lp)
  cells <- list(N = trim_cells(genetic_cells(st$N$q, lp, grid)), NN = trim_cells(genetic_cells(st$NN$q, lp, grid)))
  r <- run_aim3b(st, cells, eff[idx, , drop = FALSE], beta_hat, r_values, thresholds, assumed, draw_id = idx)
  m <- list(whole = baseline_mcid(q, lp, thresholds), N = baseline_mcid(st$N$q, lp, thresholds), NN = baseline_mcid(st$NN$q, lp, thresholds))
  r$mcid <- mapply(function(l, th) m[[l]][match(th, thresholds)], r$level, r$threshold)
  r$share_N <- st$N$share
  r
}

# Shapley rows (shapley_gain) for the draws `idx`, ceiling and realised scores, in the profile distribution q
attribution_member_rows <- function(q, lp, grid, eff, idx, beta_hat, thresholds) {
  cells <- genetic_cells(q, lp, grid)
  terms <- colnames(cells$G)
  do.call(rbind, lapply(idx, function(i) do.call(rbind, lapply(c("ceiling", "realised"), function(sc) {
    s <- shapley_gain(cells, eff[i, ], if (sc == "ceiling") eff[i, ] else beta_hat, thresholds)
    cbind(draw = i, score = sc, threshold = thresholds, total = s$total,
          stats::setNames(as.data.frame(s$phi), paste0("phi_", terms)),
          stats::setNames(as.data.frame(s$alone), paste0("alone_", terms)),
          stats::setNames(as.data.frame(s$leave_out), paste0("loo_", terms)))
  }))))
}

# Summary of Aim 3b rows over draws (and members): gain, range, extra cases, share reaching the MCID of the row, and the neutropenic stratum's share of the whole-cohort gain
summarise_aim3b_rows <- function(lay) {
  s <- by_key(lay[lay$score != "baseline", ], c("truth", "scenario", "level", "score", "threshold"), function(d) data.frame(
    gain = stats::median(d$delta_nb), lo = stats::quantile(d$delta_nb, 0.05, names = FALSE), hi = stats::quantile(d$delta_nb, 0.95, names = FALSE),
    extra_cases_per_1000 = 1000 * stats::median(d$delta_tp), p_meets_mcid = mean(d$delta_nb >= d$mcid), mcid = stats::median(d$mcid)))
  key <- c("truth", "scenario", "score", "threshold", "draw")
  whole <- lay[lay$level == "whole" & lay$score != "baseline", c(key, "delta_nb")]
  n <- lay[lay$level == "N", c(key, "delta_nb", "share_N")]
  whole$gain_N <- n$share_N[match(do.call(paste, whole[key]), do.call(paste, n[key]))] * n$delta_nb[match(do.call(paste, whole[key]), do.call(paste, n[key]))]
  whole <- whole[whole$delta_nb > 0, ]
  sh <- by_key(whole, c("truth", "scenario", "score", "threshold"), function(d) data.frame(share_of_gain_from_neutropenic = stats::median(d$gain_N / d$delta_nb)))
  s <- merge(s, sh, by = c("truth", "scenario", "score", "threshold"), all.x = TRUE)
  s$share_of_gain_from_neutropenic[s$level != "whole"] <- NA
  s
}

# Summary of Shapley rows: medians, shares of the total, and the share of draws in which the two CLEC7A terms together contribute more than any other single term
summarise_attribution_rows <- function(a, terms) {
  s <- summarise_attribution(a, terms)
  others <- setdiff(terms, c("CLEC7A", "CLEC7A_AG"))
  clec <- a[["phi_CLEC7A"]] + a[["phi_CLEC7A_AG"]]
  a$clec_top <- clec > do.call(pmax, lapply(others, function(t) a[[paste0("phi_", t)]]))
  top <- by_key(a, c("truth", "score", "threshold"), function(d) data.frame(p_clec7a_combined_largest = mean(d$clec_top)))
  s$share_CLEC7A_combined <- s$share_CLEC7A + s$share_CLEC7A_AG
  merge(s, top, by = c("truth", "score", "threshold"))
}
