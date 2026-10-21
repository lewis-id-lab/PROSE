# Why the profile truth gives the realised score a larger gain than the global truth when all of the effect is in non-neutropenic admissions (addendum 30).
# Run from the project root, after scripts 08, 11 and 37.
#
# Addenda 27 and 28 found that in this scenario the global intercept makes the CEILING score look better, but the REALISED score (which carries the pooled estimates in both
# strata) worse. The scenario: the truth has no effect in the neutropenic stratum (N) and an effect 4.93 times the pooled one in the non-neutropenic stratum (NN); the realised
# score carries the pooled estimates (1 times) in both. Pooled point estimates, 5% threshold, per 1,000 admissions of the stratum.
#
#   1. A 2 x 2: how the truth is built (global or profile) crossed with how the score is recentred (global or profile), by stratum: is it the truth or the score?
#   2. The flows in NN: the flags the score adds and removes against the clinical score, how risky those patients are under each truth, and the net benefit of each flow.
#   3. The same by the clinical-risk band of the patient's profile, for the realised and the ceiling score, and what the global intercept does to true risk in each band.
#
# Writes   output/aim3b_sign_reversal_factorial.csv, aim3b_sign_reversal_flows.csv, aim3b_sign_reversal_bands.csv

for (f in c("baseline", "bs_data", "decision_curve", "boscore_data", "refit_baseline", "evidence_mapping", "genetic_layer", "aim3b", "calibration_check", "scheme_comparison")) source(file.path("R", paste0(f, ".R")))

base <- select_cohort(load_base_cohort(), include_twins = FALSE); lp <- profile_logit(load_refit_params())
q <- tabulate(profile_index(base$X), 2L^length(PREDICTORS)) / nrow(base$X); strata <- neutropenia_strata(q, lp)
cells <- lapply(list(N = strata$N, NN = strata$NN), function(s) trim_cells(genetic_cells(s$q, lp, dosage_grid(genetic_term_frequencies()))))
bhat <- pooled_estimates(load_term_draws()); w_ev <- strata$w; share <- c(N = strata$N$share, NN = strata$NN$share); th <- 0.05
tr <- scenario_effects(bhat, 0, w_ev)

# Truth built under scheme st, score recentred under scheme ss; the comparator is the clinical score
fact <- function(cl, b, bs, st, ss) {
  clin <- plogis(cl$lp_base); lp_t <- cl$lp_base + drop(cl$G %*% b[colnames(cl$G)])
  r <- plogis(lp_t + cell_shift(lp_t, cl$w, cl$profile, clin, cl$prevalence, st))
  lp_s <- cl$lp_base + drop(cl$G %*% bs[colnames(cl$G)]); s <- plogis(lp_s + cell_shift(lp_s, cl$w, cl$profile, clin, cl$pred_mean, ss))
  f <- function(pred) { fl <- pred >= th; c(flagged = sum(cl$w[fl]), tp = sum((cl$w * r)[fl]), fp = sum((cl$w * (1 - r))[fl])) * 1000 }
  a <- f(clin); b2 <- f(s); nb <- function(v) v[["tp"]] - v[["fp"]] * th / (1 - th)
  c(gain = nb(b2) - nb(a), d_flagged = b2[["flagged"]] - a[["flagged"]], d_tp = b2[["tp"]] - a[["tp"]], d_fp = b2[["fp"]] - a[["fp"]])
}
combos <- expand.grid(truth_built = c("global", "profile"), score_recentred = c("global", "profile"), stringsAsFactors = FALSE)
factorial <- do.call(rbind, lapply(seq_len(nrow(combos)), function(i) do.call(rbind, lapply(c("N", "NN"), function(s) {
  data.frame(combos[i, ], stratum = s, share = share[[s]], t(fact(cells[[s]], tr[[s]], bhat, combos$truth_built[i], combos$score_recentred[i])), row.names = NULL) }))))
factorial$contribution_to_whole <- factorial$gain * factorial$share

# The NN stratum
cl <- cells$NN; clin <- plogis(cl$lp_base)
band <- cut(clin, c(0, 0.01, 0.03, 0.05, 0.10, 1), labels = c("<1%", "1-3%", "3-5%", "5-10%", ">10%"))
truth_r <- function(st) { lp_t <- cl$lp_base + drop(cl$G %*% tr$NN[colnames(cl$G)]); plogis(lp_t + cell_shift(lp_t, cl$w, cl$profile, clin, cl$prevalence, st)) }
score_pred <- list(realised = function() { lp_s <- cl$lp_base + drop(cl$G %*% bhat[colnames(cl$G)]); plogis(lp_s + cell_shift(lp_s, cl$w, cl$profile, clin, cl$pred_mean, "global")) },
                   ceiling = function(st) truth_r(st))
flows <- do.call(rbind, lapply(c("realised", "ceiling"), function(sc) do.call(rbind, lapply(c("global", "profile"), function(st) {
  r <- truth_r(st); p <- if (sc == "realised") score_pred$realised() else score_pred$ceiling(st)
  g <- ifelse(p >= th & clin < th, "added", ifelse(clin >= th & p < th, "removed", "other"))
  do.call(rbind, lapply(c("added", "removed"), function(k) do.call(rbind, lapply(c("all", levels(band)), function(b) {
    i <- g == k & (b == "all" | band == b); if (!any(i)) return(NULL)
    data.frame(score = sc, truth = st, group = k, clinical_risk_band = b, per_1000 = 1000 * sum(cl$w[i]), mean_true_risk = sum((cl$w * r)[i]) / sum(cl$w[i]),
               mean_clinical_risk = sum((cl$w * clin)[i]) / sum(cl$w[i]), nb = (if (k == "added") 1 else -1) * 1000 * sum((cl$w * (r - (1 - r) * th / (1 - th)))[i])) }))))
}))))
bands <- do.call(rbind, lapply(levels(band), function(b) { i <- band == b; rg <- truth_r("global")
  data.frame(clinical_risk_band = b, share_of_stratum = sum(cl$w[i]), mean_clinical_risk = sum((cl$w * clin)[i]) / sum(cl$w[i]), true_risk_global = sum((cl$w * rg)[i]) / sum(cl$w[i]),
             global_over_clinical = sum((cl$w * rg)[i]) / sum((cl$w * clin)[i])) }))

dir.create("output", showWarnings = FALSE)
write.csv(factorial, "output/aim3b_sign_reversal_factorial.csv", row.names = FALSE)
write.csv(flows, "output/aim3b_sign_reversal_flows.csv", row.names = FALSE)
write.csv(bands, "output/aim3b_sign_reversal_bands.csv", row.names = FALSE)
message("1. Whole-cohort gain of the realised score, per 1,000 admissions (truth built x score recentred):")
print(aggregate(contribution_to_whole ~ truth_built + score_recentred, factorial, sum), row.names = FALSE, digits = 4)
message("\n2. Stratum NN, realised score, flows added and removed (all bands):"); print(flows[flows$score == "realised" & flows$clinical_risk_band == "all", ], row.names = FALSE, digits = 3)
message("\n3. What the global intercept does to true risk, by the clinical-risk band of the profile:"); print(bands, row.names = FALSE, digits = 3)
