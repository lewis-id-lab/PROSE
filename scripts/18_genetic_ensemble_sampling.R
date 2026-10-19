# The correlation ensemble combined with the sampling variability of a study (addendum 12, extension). Run from the project root, after
# scripts 16 and 17 (script 17 caches the 150 member profile distributions).
#
# Each simulated study draws one copula member and builds its clinical mix from that member: patients are resampled from the real
# admissions-per-patient counts (so genotype clustering is kept), and each admission's profile is drawn independently from the member's
# distribution (which ignores covariate correlation within a patient). The truth, the recentring and the MCID are the member's own.
# Controls: the same simulation with the empirical profile distribution in place of a member (member 0), to separate the effect of the
# independent-profile approximation from the effect of the copula; and the real-row resampling of script 17.
#
# Writes   output/genetic_ensemble_sampling_draws.csv, output/genetic_ensemble_sampling_summary.csv

source("R/baseline.R"); source("R/bs_data.R"); source("R/decision_curve.R"); source("R/boscore_data.R")
source("R/refit_baseline.R"); source("R/evidence_mapping.R"); source("R/genetic_layer.R"); source("R/genetic_helpers.R")

at <- c(0.05, 0.20)
base <- select_cohort(load_base_cohort(), include_twins = FALSE)
lp_profile <- profile_logit(load_refit_params())
q_emp <- tabulate(profile_index(base$X), 2L^length(PREDICTORS)) / nrow(base$X)
q_members <- as.matrix(read.csv("output/correlation_ensemble_profile_probs.csv")[, -1])
counts <- unname(lengths(cohort_groups(base)$groups))
p <- genetic_term_frequencies(); grid <- dosage_grid(p)
draws <- load_term_draws(); beta_hat <- pooled_estimates(draws)
beta_hi <- vapply(draws, function(d) unname(stats::quantile(d$mu, 0.95)), numeric(1))
effects_file <- read.csv("output/genetic_layer_effects.csv"); terms <- names(GENETIC_TERMS)
n_sim <- 2000L
fixed <- function(b) matrix(b, n_sim, length(b), byrow = TRUE, dimnames = list(NULL, terms))
drawn <- function(h) as.matrix(effects_file[effects_file$truth == h & effects_file$draw <= n_sim, terms])

runs <- list(
  fixed_medians = list(eff = fixed(beta_hat), seed = 3318L, member = TRUE),
  fixed_upper = list(eff = fixed(beta_hi), seed = 7064L, member = TRUE),
  draws_mean = list(eff = drawn("mean"), seed = 1529L, member = TRUE),
  draws_predictive = list(eff = drawn("predictive"), seed = 8843L, member = TRUE),
  empirical_iid_fixed_medians = list(eff = fixed(beta_hat), seed = 4207L, member = FALSE),
  empirical_iid_fixed_upper = list(eff = fixed(beta_hi), seed = 6396L, member = FALSE))

cl <- start_workers()
samp <- do.call(rbind, lapply(names(runs), function(r) {
  ctx <- list(counts = counts, qm = q_members, q_emp = q_emp, lp = lp_profile, grid = grid, p = p, eff = runs[[r]]$eff, bh = beta_hat, at = at,
              seed = runs[[r]]$seed, use_member = runs[[r]]$member, members = nrow(q_members))
  cbind(analysis = r, parallel_rows(cl, chunk_ids(n_sim), function(ch, ctx) {
    set.seed(ctx$seed + ch[1])
    cache <- list()
    do.call(rbind, lapply(ch, function(i) {
      m <- if (ctx$use_member) ((i - 1L) %% ctx$members) + 1L else 0L
      key <- as.character(m)
      if (is.null(cache[[key]])) {
        q <- if (m == 0L) ctx$q_emp else ctx$qm[m, ]
        cc <- genetic_cells(q, ctx$lp, ctx$grid)
        cache[[key]] <<- list(q = q, cells = cc, mcid = setNames(baseline_mcid(q, ctx$lp, ctx$at), as.character(ctx$at)))
      }
      cm <- cache[[key]]
      sh <- genetic_scores(cm$cells, ctx$eff[i, ], ctx$bh)
      s <- simulate_study_q(ctx$counts, cm$q, ctx$lp, ctx$p, sh, ctx$eff[i, ], ctx$bh, ctx$at, cm$mcid)
      cbind(draw = i, member = m, s[s$scenario != "baseline", ])
    }))
  }, ctx))
}))
parallel::stopCluster(cl)

summ <- by_key(samp, c("analysis", "scenario", "threshold"), summarise_sampling)
ref <- read.csv("output/genetic_sampling_summary.csv")
rk <- function(d) paste(d$analysis, d$scenario, d$threshold)
ref_map <- c(fixed_medians = "fixed_medians", fixed_upper = "fixed_upper", draws_mean = "draws_mean", draws_predictive = "draws_predictive",
             empirical_iid_fixed_medians = "fixed_medians", empirical_iid_fixed_upper = "fixed_upper")
j <- match(paste(ref_map[summ$analysis], summ$scenario, summ$threshold), rk(ref))
summ$real_rows_mean_gain <- ref$mean_gain[j]; summ$real_rows_sd_gain <- ref$sd_gain[j]; summ$real_rows_p_reaches_mcid <- ref$p_reaches_mcid[j]

dir.create("output", showWarnings = FALSE)
write.csv(samp[, c("analysis", "draw", "member", "scenario", "threshold", "tp", "fp", "nb", "delta_nb", "delta_tp", "mcid_sample", "mcid", "events", "admissions")],
          "output/genetic_ensemble_sampling_draws.csv", row.names = FALSE)
write.csv(summ, "output/genetic_ensemble_sampling_summary.csv", row.names = FALSE)
message("Sampling variability with the clinical mix drawn from a copula member (observed gain of fixed scores):")
show_df(summ)
