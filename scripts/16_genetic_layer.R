# First pass of the genetic layer (addendum 11). Run from the project root, after scripts 08, 11, 12, 14 and 15.
#
# Adds the PTX3, CLEC7A, TLR4 and CD209 terms (D14) to the refit baseline on the base cohort (4,078 admissions in 1,945 patients) and
# computes, exactly, the net-benefit gain of a genotype-informed score over the baseline at 5% (primary) and 20% (secondary) for each
# draw of the true effects from the pooled posteriors. The gain is compared with the pre-specified MCID (20% more true cases at no
# extra false positives, from the refit baseline: output/dca_baseline_refit_summary.csv).
#
# Two readings of the true effect: "mean" (the pooled mean effect, a posterior draw of mu) and "predictive" (adds between-study
# heterogeneity, theta ~ N(mu, tau), as for a new population). Scores: ceiling (true effects in the score) and realised (posterior
# medians in the score). The exact engine is checked against a patient-level simulation of a 100-fold replicated cohort.
#
# Not yet carried through: the predictor-correlation ensemble (the empirical profile distribution is used), the Aim 3b strata,
# donor genotype, and the CLEC7A second parameter (held at zero).
#
# Writes   output/genetic_layer_effects.csv       the true effects drawn for each world
#          output/genetic_layer_draws.csv         per draw, score and threshold: net benefit, gain, MCID test
#          output/genetic_layer_summary.csv       summaries over draws
#          output/genetic_layer_engine_check.csv  exact engine against simulation

source("R/baseline.R")
source("R/bs_data.R")
source("R/decision_curve.R")
source("R/boscore_data.R")
source("R/refit_baseline.R")
source("R/evidence_mapping.R")
source("R/genetic_layer.R")

at <- c(0.05, 0.20)
n_draw <- 2000L

base <- select_cohort(load_base_cohort(), include_twins = FALSE)
pars <- load_refit_params()
lp_profile <- profile_logit(pars)
q <- tabulate(profile_index(base$X), 2L^length(PREDICTORS)) / nrow(base$X)

mcid_tab <- read.csv("output/dca_baseline_refit_summary.csv")
mcid <- setNames(mcid_tab$mcid_b[match(at, mcid_tab$threshold)], as.character(at))

p <- genetic_term_frequencies()
cells <- genetic_cells(q, lp_profile, dosage_grid(p))
draws <- load_term_draws()
beta_hat <- pooled_estimates(draws)

message("base cohort: ", nrow(base$X), " admissions, ", length(unique(base$rows$patient)), " patients; baseline prevalence ", sprintf("%.4f", cells$prevalence))
message("risk-class frequencies: ", paste(names(p), sprintf("%.3f", p), collapse = ", "))
message("pooled estimates (log OR per copy): ", paste(names(beta_hat), sprintf("%.3f", beta_hat), collapse = ", "))
message("MCID: ", paste(at, sprintf("%.5f", mcid), sep = " = ", collapse = ", "))

effects <- list(mean = draw_term_effects(draws, n_draw, "mean", seed = 6417L),
                predictive = draw_term_effects(draws, n_draw, "predictive", seed = 3052L))
res <- do.call(rbind, lapply(names(effects), function(h) cbind(truth = h, run_genetic_layer(cells, effects[[h]], beta_hat, at, mcid))))

summarise <- function(d) {
  data.frame(
    nb = median(d$nb), delta_nb = median(d$delta_nb), delta_lo = quantile(d$delta_nb, 0.05, names = FALSE), delta_hi = quantile(d$delta_nb, 0.95, names = FALSE),
    extra_cases_per_1000 = 1000 * median(d$delta_tp), p_meets_mcid = mean(d$meets_mcid), mcid = d$mcid[1],
    auc = median(d$auc), flagged = median(d$flagged))
}
base_tp <- res[res$scenario == "baseline", c("truth", "draw", "threshold", "tp")]
res$delta_tp <- res$tp - base_tp$tp[match(paste(res$truth, res$draw, res$threshold), paste(base_tp$truth, base_tp$draw, base_tp$threshold))]
keys <- unique(res[, c("truth", "scenario", "threshold")])
summary_tab <- do.call(rbind, lapply(seq_len(nrow(keys)), function(i) {
  d <- res[res$truth == keys$truth[i] & res$scenario == keys$scenario[i] & res$threshold == keys$threshold[i], ]
  cbind(keys[i, ], summarise(d))
}))

# Check of the exact engine: simulate patient genotypes and outcomes on a 100-fold replicated cohort at the pooled estimates
rep_cohort <- replicate_cohort(base, 100L)
sim <- simulate_nb(rep_cohort, lp_profile, p, beta_hat, cells$prevalence, at, seed = 8264L)
ex <- genetic_nb(cells, beta_hat, beta_hat, at)
check <- data.frame(threshold = at,
                    exact_baseline = ex$nb[ex$scenario == "baseline"], simulated_baseline = sim$baseline,
                    exact_ceiling = ex$nb[ex$scenario == "ceiling"], simulated_ceiling = sim$ceiling)
check$exact_gain <- check$exact_ceiling - check$exact_baseline
check$simulated_gain <- check$simulated_ceiling - check$simulated_baseline

dir.create("output", showWarnings = FALSE)
write.csv(do.call(rbind, lapply(names(effects), function(h) cbind(truth = h, draw = seq_len(n_draw), as.data.frame(effects[[h]])))), "output/genetic_layer_effects.csv", row.names = FALSE)
write.csv(res[res$scenario != "baseline", c("truth", "draw", "scenario", "threshold", "nb", "delta_nb", "delta_tp", "mcid", "meets_mcid", "auc", "flagged")], "output/genetic_layer_draws.csv", row.names = FALSE)
write.csv(summary_tab, "output/genetic_layer_summary.csv", row.names = FALSE)
write.csv(check, "output/genetic_layer_engine_check.csv", row.names = FALSE)

message("\nSummary over ", n_draw, " draws (gain = net-benefit gain over the baseline under the same truth):")
print(cbind(summary_tab[1:3], round(summary_tab[-(1:3)], 5)), row.names = FALSE)
message("\nExact engine against simulation:")
print(round(check, 5), row.names = FALSE)
