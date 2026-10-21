# Intervention test of the explanation in addendum 30, for the four- and the five-term panel (addenda 31 and 32). Run from the project root, after scripts 22 and 39.
#
# Addendum 30 explained why the profile truth gives the realised score a larger gain than the global truth when all of the effect is in non-neutropenic admissions: the
# global intercept lowers true risk in the profiles with a clinical risk of 3% to 5%, where the realised score's added flags come from. This changes the truth in chosen
# clinical-risk bands of the non-neutropenic stratum (NN) and asks whether that alone reproduces the gap. Pooled point estimates, 5%, per 1,000 NN admissions.
#
# The truth in a band is either the global truth's risk or the profile truth's, or the global truth's risk rescaled to the clinical score's mean in each profile (keeps the
# global truth's pattern across genotypes and removes only its shift of the profile's mean). For each intervention:
#   realised score  the pooled estimates, recentred globally, as in script 39; the comparator is the clinical score
#   ceiling score   the truth itself (the oracle for that truth)
#   gap recovered   (gain with the intervention - gain under the global truth) / (gain under the profile truth - gain under the global truth)
# The five-term panel (PTX3, CLEC7A, CLEC7A_AG, TLR4, CD209, with the joint CLEC7A genotypes of script 22) is run the same way; its effect scaling is as in script 37.
# Gates: the four-term rows must reproduce script 39 (the no-band and every-band cases); for both panels the no-band case must reproduce stratum_gains() of
# R/scheme_comparison.R for the realised and the ceiling score, and the every-band case its ceiling gain under the profile truth.
#
# Also writes, for each panel, the flags the realised score adds and removes against the clinical score by band, under each truth (what addendum 30 used to explain the gap).
#
# Writes   output/aim3b_sign_reversal_intervention.csv, output/aim3b_sign_reversal_intervention_flows.csv

for (f in c("baseline", "bs_data", "decision_curve", "boscore_data", "refit_baseline", "evidence_mapping", "genetic_layer", "aim3b", "calibration_check", "scheme_comparison")) source(file.path("R", paste0(f, ".R")))

base <- select_cohort(load_base_cohort(), include_twins = FALSE); lp <- profile_logit(load_refit_params())
q <- tabulate(profile_index(base$X), 2L^length(PREDICTORS)) / nrow(base$X); strata <- neutropenia_strata(q, lp); w_ev <- strata$w; th <- 0.05
s5 <- clec7a_second_setup(n_draw = 10L)
models <- list(`four terms` = list(grid = dosage_grid(genetic_term_frequencies()), bhat = pooled_estimates(load_term_draws())), `five terms` = list(grid = s5$grid, bhat = s5$beta_hat))

run_model <- function(name, grid, bhat) {
  tr <- scenario_effects(bhat, 0, w_ev)$NN                               # all of the effect in the non-neutropenic stratum
  cl <- trim_cells(genetic_cells(strata$NN$q, lp, grid)); clin <- plogis(cl$lp_base)
  band <- cut(clin, c(0, 0.01, 0.03, 0.05, 0.10, 1), labels = c("<1%", "1-3%", "3-5%", "5-10%", ">10%")); bands <- levels(band)
  lp_t <- cl$lp_base + drop(cl$G %*% tr[colnames(cl$G)])
  r_global <- plogis(lp_t + cell_shift(lp_t, cl$w, cl$profile, clin, cl$prevalence, "global"))
  r_profile <- plogis(lp_t + cell_shift(lp_t, cl$w, cl$profile, clin, cl$prevalence, "profile"))
  r_rescaled <- pmin(r_global * clin / profile_mean(r_global, cl), 1)     # the global truth's pattern, with each profile's mean put back to the clinical risk
  lp_s <- cl$lp_base + drop(cl$G %*% bhat[colnames(cl$G)]); s_real <- plogis(lp_s + cell_shift(lp_s, cl$w, cl$profile, clin, cl$pred_mean, "global"))
  nb_gain <- function(r, pred) { f <- function(p) { fl <- p >= th; 1000 * (sum((cl$w * r)[fl]) - sum((cl$w * (1 - r))[fl]) * th / (1 - th)) }; f(pred) - f(clin) }
  build <- function(changed, to) { r <- r_global; i <- band %in% changed; r[i] <- (if (to == "profile") r_profile else r_rescaled)[i]; r }
  res_one <- function(label, r) data.frame(intervention = label, realised_gain = nb_gain(r, s_real), ceiling_gain = nb_gain(r, r))
  ivs <- list(
    list("no band (the global truth)", character(0), "profile"), list("every band (the profile truth)", bands, "profile"),
    list("3-5% only, to the profile truth", "3-5%", "profile"), list("every band except 3-5%, to the profile truth", setdiff(bands, "3-5%"), "profile"),
    list("<1% only, to the profile truth", "<1%", "profile"), list("1-3% only, to the profile truth", "1-3%", "profile"),
    list("5-10% only, to the profile truth", "5-10%", "profile"), list(">10% only, to the profile truth", ">10%", "profile"),
    list("1-3% and 3-5%, to the profile truth", c("1-3%", "3-5%"), "profile"),
    list("3-5% only, mean rescaled (global pattern kept)", "3-5%", "rescaled"), list("1-3% and 3-5%, mean rescaled (global pattern kept)", c("1-3%", "3-5%"), "rescaled"),
    list("<1% only, mean rescaled (global pattern kept)", "<1%", "rescaled"), list("every band, mean rescaled (global pattern kept)", bands, "rescaled"))
  out <- do.call(rbind, lapply(ivs, function(v) cbind(res_one(v[[1]], build(v[[2]], v[[3]])), bands_changed = paste(v[[2]], collapse = " + "))))
  g0 <- out[out$intervention == "no band (the global truth)", ]; g1 <- out[out$intervention == "every band (the profile truth)", ]
  out$realised_gap_recovered <- (out$realised_gain - g0$realised_gain) / (g1$realised_gain - g0$realised_gain)
  out$ceiling_gap_recovered <- (out$ceiling_gain - g0$ceiling_gain) / (g1$ceiling_gain - g0$ceiling_gain)
  r_rev <- r_profile; i <- band == "3-5%"; r_rev[i] <- r_global[i]                      # the reverse: the profile truth with the global values back in 3-5%
  rev <- cbind(res_one("profile truth with the global values in 3-5% only", r_rev), bands_changed = "3-5% back to global")
  rev$realised_gap_recovered <- (rev$realised_gain - g1$realised_gain) / (g0$realised_gain - g1$realised_gain); rev$ceiling_gap_recovered <- (rev$ceiling_gain - g1$ceiling_gain) / (g0$ceiling_gain - g1$ceiling_gain)
  out <- cbind(model = name, rbind(out, rev))
  # Gates
  sg0 <- stratum_gains(cl, tr, bhat, "global", th); sg1 <- stratum_gains(cl, tr, bhat, "profile", th)
  stopifnot(abs(g0$realised_gain - 1000 * sg0$real_clin) < 1e-9, abs(g0$ceiling_gain - 1000 * sg0$ceil_clin) < 1e-9, abs(g1$ceiling_gain - 1000 * sg1$ceil_clin) < 1e-9)
  # The flags the realised score adds and removes against the clinical score, by band, under each truth
  flows <- do.call(rbind, lapply(c("global", "profile"), function(st) { r <- if (st == "global") r_global else r_profile
    g <- ifelse(s_real >= th & clin < th, "added", ifelse(clin >= th & s_real < th, "removed", "other"))
    do.call(rbind, lapply(c("added", "removed"), function(k) do.call(rbind, lapply(bands, function(b) { i <- g == k & band == b; if (!any(i)) return(NULL)
      data.frame(model = name, truth = st, group = k, clinical_risk_band = b, per_1000 = 1000 * sum(cl$w[i]), mean_true_risk = sum((cl$w * r)[i]) / sum(cl$w[i]),
                 nb = (if (k == "added") 1 else -1) * 1000 * sum((cl$w * (r - (1 - r) * th / (1 - th)))[i])) })))) }))
  list(out = out, flows = flows)
}
res <- Map(function(nm, md) run_model(nm, md$grid, md$bhat), names(models), models)
out <- do.call(rbind, lapply(res, `[[`, "out")); flows <- do.call(rbind, lapply(res, `[[`, "flows")); rownames(out) <- NULL; rownames(flows) <- NULL

# The four-term rows reproduce script 39
ref <- read.csv("output/aim3b_sign_reversal_factorial.csv"); ref <- ref[ref$stratum == "NN" & ref$score_recentred == "global", ]
f4 <- out[out$model == "four terms", ]
stopifnot(abs(f4$realised_gain[f4$intervention == "no band (the global truth)"] - ref$gain[ref$truth_built == "global"]) < 1e-9,
          abs(f4$realised_gain[f4$intervention == "every band (the profile truth)"] - ref$gain[ref$truth_built == "profile"]) < 1e-9)

dir.create("output", showWarnings = FALSE)
write.csv(out, "output/aim3b_sign_reversal_intervention.csv", row.names = FALSE)
write.csv(flows, "output/aim3b_sign_reversal_intervention_flows.csv", row.names = FALSE)
for (m in names(models)) { message("\n", m, ": stratum NN, 5%, per 1,000 NN admissions:"); x <- out[out$model == m, ]
  print(x[, c("intervention", "realised_gain", "realised_gap_recovered", "ceiling_gain", "ceiling_gap_recovered")], row.names = FALSE, digits = 3) }
