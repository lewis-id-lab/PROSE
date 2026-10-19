# Donor genotype in the allogeneic HSCT stratum: the D1 secondary analysis (addendum 12, extension). Run from the project root, after
# scripts 12, 14, 15 and 16 (they write the donor posterior draws).
#
# Two strata. "data" (primary, addendum 14): the base cohort's admission-level allo-HSCT flag (script 08) gives the stratum's own profile
# distribution, and a recalibration offset on the baseline logit sets its predicted event rate to the observed one (allo_stratum_data()).
# "margins" (addendum 13, kept as a sensitivity analysis): built from the published Table 1 margins alone, under the assumption that allo
# status depends on the outcome only, given the clinical profile (allo_stratum()). In both the baseline score has no allo item, so its
# decisions are unchanged and it under-predicts in the stratum.
#
# In the stratum the genetic score uses donor genotype and the donor effects (D1: after allo-HSCT the innate immune system is donor-derived).
# Donor dosages are Hardy-Weinberg with the same class frequencies, independent of the recipient's profile. CD209 has no donor-genotype
# study, so its donor effect is not estimable and is held at zero. Donor evidence is one study each for PTX3 (Fisher 2017, hazard ratio)
# and CLEC7A (Ceesay 2016), two for TLR4, and so is weaker than the recipient evidence.
#
# Gains are per allo admission. The whole-cohort equivalent is the gain times the stratum's share (decisions elsewhere stay with the
# baseline score). The stratum MCID is 20% more true cases among its admissions at no extra false positives, from its own baseline; the
# whole-cohort MCID of D10 is also shown.
#
# Writes   output/genetic_donor_evidence.csv, genetic_donor_strata.csv, genetic_donor_stratum_summary.csv, genetic_donor_baseline_stratum.csv,
#          genetic_donor_sampling_summary.csv (the last three have a `stratum` column)

source("R/baseline.R"); source("R/bs_data.R"); source("R/decision_curve.R"); source("R/boscore_data.R")
source("R/refit_baseline.R"); source("R/evidence_mapping.R"); source("R/genetic_layer.R"); source("R/genetic_helpers.R")

at <- c(0.05, 0.20)
base <- select_cohort(load_base_cohort(), include_twins = FALSE)
lp <- profile_logit(load_refit_params())
q <- tabulate(profile_index(base$X), 2L^length(PREDICTORS)) / nrow(base$X)
mcid_tab <- read.csv("output/dca_baseline_refit_summary.csv")
mcid_whole <- setNames(mcid_tab$mcid_b[match(at, mcid_tab$threshold)], as.character(at))

p <- genetic_term_frequencies(); grid <- dosage_grid(p); terms <- names(GENETIC_TERMS)
strata <- list(margins = allo_stratum(q, lp))
strata$margins$n <- round(strata$margins$share * nrow(base$X))
if (!all(is.na(base$rows$allo_hsct))) strata <- c(list(data = allo_stratum_data(base, lp)), strata) else message("no allo-HSCT flag in the cohort: only the Table 1 construction is run")

strata_tab <- do.call(rbind, lapply(names(strata), function(n) { s <- strata[[n]]
  data.frame(stratum = n, share = s$share, admissions = s$n, event_rate = s$event_rate, odds_multiplier = exp(s$offset),
             profiles_present = sum(s$q > 0), baseline_predicted_rate = sum(s$q * plogis(lp)),
             tv_distance_to_cohort = 0.5 * sum(abs(s$q - q))) }))
message("strata:"); show_df(strata_tab)

donor <- load_donor_draws(); beta_hat <- donor_pooled_estimates(donor)
beta_hi <- vapply(donor, function(d) if (is.null(d)) 0 else unname(stats::quantile(d$mu, 0.95)), numeric(1))
n_draw <- 2000L
effects <- list(mean = draw_donor_effects(donor, n_draw, "mean", seed = 2905L), predictive = draw_donor_effects(donor, n_draw, "predictive", seed = 5173L))

evidence <- data.frame(term = terms, donor_studies = vapply(donor, function(d) if (is.null(d)) 0L else 1L, integer(1)),
                       log_or_median = unname(beta_hat), or_median = exp(unname(beta_hat)),
                       or_lo = vapply(donor, function(d) if (is.null(d)) NA_real_ else exp(unname(stats::quantile(d$mu, 0.05))), numeric(1)),
                       or_hi = vapply(donor, function(d) if (is.null(d)) NA_real_ else exp(unname(stats::quantile(d$mu, 0.95))), numeric(1)),
                       recipient_or_median = exp(unname(pooled_estimates(load_term_draws()))), row.names = NULL)
evidence$donor_studies <- c(PTX3 = 1L, CLEC7A = 1L, TLR4 = 2L, CD209 = 0L)[terms]

fixed <- function(b) matrix(b, n_draw, length(b), byrow = TRUE, dimnames = list(NULL, terms))
runs <- list(fixed_medians = list(eff = fixed(beta_hat), seed = 4630L), fixed_upper = list(eff = fixed(beta_hi), seed = 9082L),
             draws_mean = list(eff = effects$mean, seed = 1847L), draws_predictive = list(eff = effects$predictive, seed = 7315L))

cl <- start_workers()
out <- lapply(names(strata), function(sn) {
  st <- strata[[sn]]
  cells <- genetic_cells(st$q, lp, grid, offset = st$offset)
  mcid <- setNames(baseline_mcid(st$q, lp, at, offset = st$offset), as.character(at))
  message(sprintf("\n%s stratum: share %.4f (%d admissions), event rate %.4f, odds multiplier %.2f; stratum MCID %s",
                  sn, st$share, st$n, st$event_rate, exp(st$offset), paste(sprintf("%.5f", mcid), collapse = " and ")))
  lay <- exact_layer(cl, cells, effects, beta_hat, at, mcid)
  key <- paste(lay$truth, lay$draw, lay$threshold)
  bse <- lay[lay$scenario == "baseline", ]
  lay$delta_tp <- lay$tp - bse$tp[match(key, paste(bse$truth, bse$draw, bse$threshold))]
  lay$whole_gain <- st$share * lay$delta_nb
  lay$mcid_whole <- unname(mcid_whole[as.character(lay$threshold)])
  summ <- by_key(lay[lay$scenario != "baseline", ], c("truth", "scenario", "threshold"), function(d) data.frame(
    stratum_gain = median(d$delta_nb), lo = quantile(d$delta_nb, 0.05, names = FALSE), hi = quantile(d$delta_nb, 0.95, names = FALSE),
    extra_cases_per_1000_allo = 1000 * median(d$delta_tp), p_meets_stratum_mcid = mean(d$meets_mcid), stratum_mcid = d$mcid[1],
    whole_cohort_gain = median(d$whole_gain), p_meets_whole_mcid = mean(d$whole_gain >= d$mcid_whole), whole_mcid = d$mcid_whole[1],
    auc = median(d$auc)))
  ref <- genetic_nb(cells, beta_hat, beta_hat, at); b0 <- ref[ref$scenario == "baseline", ]
  baseline_stratum <- data.frame(threshold = at, nb = b0$nb, tp_per_1000_allo = 1000 * b0$tp, fp_per_1000_allo = 1000 * b0$fp, flagged = b0$flagged, auc = b0$auc, mcid = unname(mcid))
  # Sampling variability of a study with a stratum's worth of admissions (each admission its own donor)
  samp <- do.call(rbind, lapply(names(runs), function(r) {
    ctx <- list(n = st$n, q = st$q, lp = lp, p = p, cells = cells, eff = runs[[r]]$eff, bh = beta_hat, at = at, mcid = mcid, offset = st$offset, seed = runs[[r]]$seed)
    cbind(analysis = r, parallel_rows(cl, chunk_ids(n_draw), function(ch, ctx) {
      set.seed(ctx$seed + ch[1])
      do.call(rbind, lapply(ch, function(i) {
        sh <- genetic_scores(ctx$cells, ctx$eff[i, ], ctx$bh)
        s <- simulate_study_q(rep(1L, ctx$n), ctx$q, ctx$lp, ctx$p, sh, ctx$eff[i, ], ctx$bh, ctx$at, ctx$mcid, offset = ctx$offset, resample_counts = FALSE)
        cbind(draw = i, s[s$scenario != "baseline", ])
      }))
    }, ctx))
  }))
  samp_sum <- by_key(samp, c("analysis", "scenario", "threshold"), summarise_sampling)
  exact_fixed <- do.call(rbind, lapply(c("fixed_medians", "fixed_upper"), function(r) {
    x <- genetic_nb(cells, runs[[r]]$eff[1, ], beta_hat, at); b <- x[x$scenario == "baseline", ]; g <- x[x$scenario != "baseline", ]
    data.frame(analysis = r, scenario = g$scenario, threshold = g$threshold, exact_gain = g$nb - b$nb[match(g$threshold, b$threshold)])
  }))
  samp_sum$exact_gain <- exact_fixed$exact_gain[match(paste(samp_sum$analysis, samp_sum$scenario, samp_sum$threshold), paste(exact_fixed$analysis, exact_fixed$scenario, exact_fixed$threshold))]
  samp_sum$study_admissions <- st$n
  list(summ = cbind(stratum = sn, summ), baseline = cbind(stratum = sn, baseline_stratum), sampling = cbind(stratum = sn, samp_sum))
})
parallel::stopCluster(cl)
collect <- function(n) do.call(rbind, lapply(out, `[[`, n))

dir.create("output", showWarnings = FALSE)
write.csv(evidence, "output/genetic_donor_evidence.csv", row.names = FALSE)
write.csv(strata_tab, "output/genetic_donor_strata.csv", row.names = FALSE)
write.csv(collect("summ"), "output/genetic_donor_stratum_summary.csv", row.names = FALSE)
write.csv(collect("baseline"), "output/genetic_donor_baseline_stratum.csv", row.names = FALSE)
write.csv(collect("sampling"), "output/genetic_donor_sampling_summary.csv", row.names = FALSE)
message("\nDonor evidence:"); show_df(evidence)
message("\nBaseline in the strata:"); show_df(collect("baseline"))
message("\nDonor-genotype gain per allo admission and whole-cohort equivalent:"); show_df(collect("summ"))
message("\nSampling variability of a stratum-sized study:"); show_df(collect("sampling"))
