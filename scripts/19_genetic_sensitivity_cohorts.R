# The genetic layer on the sensitivity cohorts (addendum 12, extension). Run from the project root, after scripts 16 and 17.
#
# Settings (the primary one, the refit baseline on the 4,078-admission base cohort, is in scripts 16 and 17):
#   published_4078  the published Stanzani equation as the baseline model and as the clinical part of the truth, on the base cohort
#   published_4127  the published equation on the 4,127 rows of the paper (twin rows kept): the D8 baseline
#   refit_4127      a refit of the seven predictors on the 4,127 rows, on those rows
# For each: the exact layer (2,000 draws per reading), Shapley attribution (1,000 draws), and the sampling variability of a study of
# the cohort's size (four runs of 2,000 studies). The effects, the seeds and the pooled estimates are those of the primary analysis, so
# the settings differ only in the cohort and the baseline model. The MCID is the setting's own (20% more true cases of its baseline).
#
# Set SENS_ONLY to a subset of the setting names to run them one at a time; per-setting files are cached in output/ and merged.
# Writes   output/genetic_sensitivity_{layer,attribution,sampling}_<setting>.csv and the merged genetic_sensitivity_*_summary.csv

source("R/baseline.R"); source("R/bs_data.R"); source("R/decision_curve.R"); source("R/boscore_data.R")
source("R/refit_baseline.R"); source("R/evidence_mapping.R"); source("R/genetic_layer.R"); source("R/genetic_helpers.R")

# Set Sys.setenv(PROSE_SCHEME = "profile") before sourcing to run under the profile-preserving truth (addendum 25); files then end in _profile
sfx <- scheme_suffix()
at <- c(0.05, 0.20)
all_settings <- c("published_4078", "published_4127", "refit_4127")
run_these <- if (exists("SENS_ONLY")) SENS_ONLY else all_settings

full <- load_base_cohort()
cohorts <- list(`4078` = select_cohort(full, include_twins = FALSE), `4127` = select_cohort(full, include_twins = TRUE))
published <- load_baseline_params()
spec <- list(
  published_4078 = list(cohort = cohorts[["4078"]], lp = profile_logit(published)),
  published_4127 = list(cohort = cohorts[["4127"]], lp = profile_logit(published)),
  refit_4127 = list(cohort = cohorts[["4127"]], lp = profile_logit(refit_params(fit_baseline(cohorts[["4127"]])$coef))))

p <- genetic_term_frequencies(); grid <- dosage_grid(p)
draws <- load_term_draws(); beta_hat <- pooled_estimates(draws)
beta_hi <- vapply(draws, function(d) unname(stats::quantile(d$mu, 0.95)), numeric(1))
effects_file <- read.csv("output/genetic_layer_effects.csv"); terms <- names(GENETIC_TERMS)
n_sim <- 2000L
drawn <- function(h) as.matrix(effects_file[effects_file$truth == h & effects_file$draw <= n_sim, terms])
fixed <- function(b) matrix(b, n_sim, length(b), byrow = TRUE, dimnames = list(NULL, terms))
effects <- list(mean = drawn("mean"), predictive = drawn("predictive"))
runs <- list(fixed_medians = list(eff = fixed(beta_hat), seed = 5290L), fixed_upper = list(eff = fixed(beta_hi), seed = 6841L),
             draws_mean = list(eff = drawn("mean"), seed = 2736L), draws_predictive = list(eff = drawn("predictive"), seed = 9417L))

cl <- start_workers()
for (s in run_these) {
  co <- spec[[s]]$cohort; lp <- spec[[s]]$lp
  q <- tabulate(profile_index(co$X), 2L^length(PREDICTORS)) / nrow(co$X)
  cells <- genetic_cells(q, lp, grid)
  mcid <- setNames(baseline_mcid(q, lp, at), as.character(at))
  message(s, ": ", nrow(co$X), " admissions, ", length(unique(co$rows$patient)), " patients; baseline prevalence ",
          sprintf("%.4f", cells$prevalence), "; MCID ", paste(sprintf("%.5f", mcid), collapse = " and "))
  lay <- exact_layer(cl, cells, effects, beta_hat, at, mcid)
  att <- attribution_layer(cl, cells, effects, beta_hat, at)
  smp <- sampling_layer(cl, cohort_groups(co), lp, p, cells, runs, beta_hat, at, mcid)
  out <- list(layer = cbind(setting = s, summarise_layer(lay)),
              attribution = cbind(setting = s, summarise_attribution(att, terms)),
              sampling = cbind(setting = s, by_key(smp, c("analysis", "scenario", "threshold"), summarise_sampling)))
  exact_fixed <- do.call(rbind, lapply(c("fixed_medians", "fixed_upper"), function(r) {
    x <- genetic_nb(cells, runs[[r]]$eff[1, ], beta_hat, at); b <- x[x$scenario == "baseline", ]
    g <- x[x$scenario != "baseline", ]
    data.frame(analysis = r, scenario = g$scenario, threshold = g$threshold, exact_gain = g$nb - b$nb[match(g$threshold, b$threshold)])
  }))
  out$sampling$exact_gain <- exact_fixed$exact_gain[match(paste(out$sampling$analysis, out$sampling$scenario, out$sampling$threshold),
                                                           paste(exact_fixed$analysis, exact_fixed$scenario, exact_fixed$threshold))]
  for (n in names(out)) write.csv(out[[n]], sprintf("output/genetic_sensitivity_%s_%s%s.csv", n, s, sfx), row.names = FALSE)
  message("\n", s, ": exact layer"); show_df(out$layer)
}
parallel::stopCluster(cl)

have <- all_settings[file.exists(sprintf("output/genetic_sensitivity_layer_%s%s.csv", all_settings, sfx))]
for (n in c("layer", "attribution", "sampling"))
  write.csv(do.call(rbind, lapply(have, function(s) read.csv(sprintf("output/genetic_sensitivity_%s_%s%s.csv", n, s, sfx)))),
            sprintf("output/genetic_sensitivity_%s_summary%s.csv", n, sfx), row.names = FALSE)
