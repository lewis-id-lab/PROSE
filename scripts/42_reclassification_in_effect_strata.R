# Is the realised score's gain in a stratum where the truth DOES have an effect governed by the same threshold-crossing mechanism as its loss in the null stratum? (addendum 34)
# Run from the project root, after scripts 22, 40 and 41.
#
# Addendum 33 showed that in the neutropenic stratum, where the truth has no effect, the realised score's loss is the cost of reclassifying a few large profiles that lie close to
# the 5% line. Which patients a score reclassifies depends only on the score (its weights and the clinical profile), not on the truth; the truth only values them. For ANY truth the
# net-benefit change is the sum over cells of weight x (true risk - 5%) / 0.95 x (flag by the score - flag by the clinical score), so a score gains when the patients it adds have a
# true risk above the threshold and the patients it removes have one below it.
#
# This takes the six scenarios as the paper uses them (script 21; the Brunel-estimate stratum fractions are N -0.21 and NN 1.87), for the four- and five-term panels, and for each
# stratum (N neutropenic, NN non-neutropenic) reports, under the profile-preserving truth with the score recentred globally (as scripts 39 to 41):
#   the realised score's gain, and what it would be if the truth had no effect (the null cost of the same reclassification);
#   the flags it adds and removes (the same patients in every scenario), the mean true risk of those patients, and the share that are "right" (added with true risk >= 5%,
#     removed with true risk < 5%);
#   how concentrated the reclassification and the gain are in a few profiles.
# And, for the homogeneous scenario, the gain and the null cost against the scale k of the score's weights (the truth fixed).
# Gates: the null gains of N reproduce script 41; the NN gain in "all in non-neutropenic" reproduces script 40's profile-truth gain.
#
# Writes   output/reclassification_effect_strata_flows.csv, reclassification_effect_strata_scaling.csv

for (f in c("baseline", "bs_data", "decision_curve", "boscore_data", "refit_baseline", "evidence_mapping", "genetic_layer", "aim3b", "calibration_check", "scheme_comparison")) source(file.path("R", paste0(f, ".R")))

base <- select_cohort(load_base_cohort(), include_twins = FALSE); lp <- profile_logit(load_refit_params())
q <- tabulate(profile_index(base$X), 2L^length(PREDICTORS)) / nrow(base$X); strata <- neutropenia_strata(q, lp); w_ev <- strata$w; th <- 0.05
s5 <- clec7a_second_setup(n_draw = 10L)
models <- list(`four terms` = list(grid = dosage_grid(genetic_term_frequencies()), bhat = pooled_estimates(load_term_draws())), `five terms` = list(grid = s5$grid, bhat = s5$beta_hat))
r_values <- list(`all in neutropenic (H_F only)` = Inf, `ratio 4` = 4, homogeneous = 1, `Brunel direction, ratio 1/4` = 0.25, `all in non-neutropenic` = 0, `Brunel estimate` = c(N = -0.21, NN = 1.87))
u <- function(r) (r - th) / (1 - th)

stratum_setup <- function(md, S) {
  cl <- trim_cells(genetic_cells(strata[[S]]$q, lp, md$grid)); clin <- plogis(cl$lp_base)
  list(cl = cl, clin = clin, flag_c = clin >= th)
}
score_flags <- function(st, bs) { lp_s <- st$cl$lp_base + drop(st$cl$G %*% bs[colnames(st$cl$G)]); plogis(lp_s + cell_shift(lp_s, st$cl$w, st$cl$profile, st$clin, st$cl$pred_mean, "global")) >= th }
truth_r <- function(st, b) { lp_t <- st$cl$lp_base + drop(st$cl$G %*% b[colnames(st$cl$G)]); plogis(lp_t + cell_shift(lp_t, st$cl$w, st$cl$profile, st$clin, st$cl$prevalence, "profile")) }
analyse <- function(st, flag_s, r) {
  w <- st$cl$w; fc <- st$flag_c; add <- flag_s & !fc; rem <- fc & !flag_s; d <- as.numeric(flag_s) - as.numeric(fc)
  prof <- tapply(w * u(r) * d, st$cl$profile, sum) * 1000; rw <- tapply(w * abs(d), st$cl$profile, sum)
  top3 <- function(x) sum(sort(abs(x), decreasing = TRUE)[1:3]) / sum(abs(x))
  list(gain = 1000 * sum(w * u(r) * d), null_gain = 1000 * sum(w * u(st$clin) * d),
       added_per_1000 = 1000 * sum(w[add]), removed_per_1000 = 1000 * sum(w[rem]),
       added_mean_true_risk = if (any(add)) sum((w * r)[add]) / sum(w[add]) else NA_real_, removed_mean_true_risk = if (any(rem)) sum((w * r)[rem]) / sum(w[rem]) else NA_real_,
       added_share_right = if (any(add)) sum(w[add & r >= th]) / sum(w[add]) else NA_real_, removed_share_right = if (any(rem)) sum(w[rem & r < th]) / sum(w[rem]) else NA_real_,
       nb_added = 1000 * sum((w * u(r))[add]), nb_removed = -1000 * sum((w * u(r))[rem]),
       profiles_reclassified = sum(rw > 0), top3_profiles_share_of_reclassified_weight = sum(sort(rw, decreasing = TRUE)[1:3]) / sum(rw), top3_profiles_share_of_abs_gain = top3(prof))
}

flows <- do.call(rbind, lapply(names(models), function(m) do.call(rbind, lapply(c("N", "NN"), function(S) {
  md <- models[[m]]; st <- stratum_setup(md, S); fs <- score_flags(st, md$bhat)
  do.call(rbind, lapply(names(r_values), function(sc) { tr <- scenario_effects(md$bhat, r_values[[sc]], w_ev)[[S]]
    cbind(data.frame(model = m, stratum = S, scenario = sc, truth_effect_multiple = unname(tr[1] / md$bhat[1])), as.data.frame(analyse(st, fs, truth_r(st, tr)))) })) }))))
rownames(flows) <- NULL
flows$value_of_the_truth_effect <- flows$gain - flows$null_gain

# Scaling of the score's weights, the truth fixed at the homogeneous scenario (the pooled effect in both strata)
scaling <- do.call(rbind, lapply(names(models), function(m) do.call(rbind, lapply(c("N", "NN"), function(S) { md <- models[[m]]; st <- stratum_setup(md, S); r <- truth_r(st, md$bhat)
  do.call(rbind, lapply(c(0.5, 0.75, 1, 1.25, 1.5, 2), function(k) { fs <- score_flags(st, k * md$bhat); a <- analyse(st, fs, r)
    data.frame(model = m, stratum = S, k = k, gain = a$gain, null_gain = a$null_gain, value_of_the_truth_effect = a$gain - a$null_gain, reclassified_per_1000 = a$added_per_1000 + a$removed_per_1000,
               added_share_right = a$added_share_right) })) }))))

# Gates
ng <- read.csv("output/null_stratum_loss_terms.csv")
for (m in names(models)) { g <- flows$null_gain[flows$model == m & flows$stratum == "N" & flows$scenario == "homogeneous"]
  stopifnot(abs(g - ng$gain[ng$panel == m & grepl("^all", ng$score)]) < 1e-8) }
iv <- read.csv("output/aim3b_sign_reversal_intervention.csv")
for (m in names(models)) stopifnot(abs(flows$gain[flows$model == m & flows$stratum == "NN" & flows$scenario == "all in non-neutropenic"] - iv$realised_gain[iv$model == m & grepl("^every band \\(the profile", iv$intervention)]) < 1e-8)

dir.create("output", showWarnings = FALSE)
write.csv(flows, "output/reclassification_effect_strata_flows.csv", row.names = FALSE)
write.csv(scaling, "output/reclassification_effect_strata_scaling.csv", row.names = FALSE)
message("Realised score, per 1,000 stratum admissions at 5%, profile truth:"); print(flows[, c("model", "stratum", "scenario", "truth_effect_multiple", "gain", "null_gain", "value_of_the_truth_effect", "added_share_right")], row.names = FALSE, digits = 3)
