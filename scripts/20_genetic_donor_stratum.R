# Donor genotype in the allogeneic HSCT stratum: the D1 secondary analysis (addendum 12, extension). Run from the project root, after
# scripts 12, 14, 15 and 16 (they write the donor posterior draws).
#
# The base cohort has no admission-level allo-HSCT flag, so the stratum is built from the published Table 1 margins (allo is the reason
# for admission in 44 of 133 cases and 248 of 3,994 non-cases) under the assumption that allo status depends on the outcome only, given
# the clinical profile (allo_stratum() in R/genetic_layer.R). That gives a share of 7.1% of admissions, a stratum event rate of 15.2% and
# a true risk with odds 5.3 times the baseline score's at every profile. The baseline score does not know this (it has no allo item),
# so its decisions are unchanged and it under-predicts in the stratum.
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
# Writes   output/genetic_donor_evidence.csv, genetic_donor_stratum_summary.csv, genetic_donor_sampling_summary.csv

source("R/baseline.R"); source("R/bs_data.R"); source("R/decision_curve.R"); source("R/boscore_data.R")
source("R/refit_baseline.R"); source("R/evidence_mapping.R"); source("R/genetic_layer.R"); source("R/genetic_helpers.R")

at <- c(0.05, 0.20)
base <- select_cohort(load_base_cohort(), include_twins = FALSE)
lp <- profile_logit(load_refit_params())
q <- tabulate(profile_index(base$X), 2L^length(PREDICTORS)) / nrow(base$X)
mcid_tab <- read.csv("output/dca_baseline_refit_summary.csv")
mcid_whole <- setNames(mcid_tab$mcid_b[match(at, mcid_tab$threshold)], as.character(at))

p <- genetic_term_frequencies(); grid <- dosage_grid(p); terms <- names(GENETIC_TERMS)
st <- allo_stratum(q, lp)
cells <- genetic_cells(st$q, lp, grid, offset = st$offset)
mcid <- setNames(baseline_mcid(st$q, lp, at, offset = st$offset), as.character(at))
n_stratum <- round(st$share * nrow(base$X))
message(sprintf("allo-HSCT stratum: share %.4f (%d of %d admissions), event rate %.4f, odds multiplier %.2f; stratum MCID %s",
                st$share, n_stratum, nrow(base$X), st$event_rate, exp(st$offset), paste(sprintf("%.5f", mcid), collapse = " and ")))

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

cl <- start_workers()
lay <- exact_layer(cl, cells, effects, beta_hat, at, mcid)
parallel::stopCluster(cl)

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
shifts_of <- function(eff) genetic_scores(cells, eff, beta_hat)
fixed <- function(b) matrix(b, n_draw, length(b), byrow = TRUE, dimnames = list(NULL, terms))
runs <- list(fixed_medians = list(eff = fixed(beta_hat), seed = 4630L), fixed_upper = list(eff = fixed(beta_hi), seed = 9082L),
             draws_mean = list(eff = effects$mean, seed = 1847L), draws_predictive = list(eff = effects$predictive, seed = 7315L))
cl <- start_workers()
samp <- do.call(rbind, lapply(names(runs), function(r) {
  ctx <- list(n = n_stratum, q = st$q, lp = lp, p = p, cells = cells, eff = runs[[r]]$eff, bh = beta_hat, at = at, mcid = mcid, offset = st$offset, seed = runs[[r]]$seed)
  cbind(analysis = r, parallel_rows(cl, chunk_ids(n_draw), function(ch, ctx) {
    set.seed(ctx$seed + ch[1])
    do.call(rbind, lapply(ch, function(i) {
      sh <- genetic_scores(ctx$cells, ctx$eff[i, ], ctx$bh)
      s <- simulate_study_q(rep(1L, ctx$n), ctx$q, ctx$lp, ctx$p, sh, ctx$eff[i, ], ctx$bh, ctx$at, ctx$mcid, offset = ctx$offset, resample_counts = FALSE)
      cbind(draw = i, s[s$scenario != "baseline", ])
    }))
  }, ctx))
}))
parallel::stopCluster(cl)
samp_sum <- by_key(samp, c("analysis", "scenario", "threshold"), summarise_sampling)
exact_fixed <- do.call(rbind, lapply(c("fixed_medians", "fixed_upper"), function(r) {
  x <- genetic_nb(cells, runs[[r]]$eff[1, ], beta_hat, at); b <- x[x$scenario == "baseline", ]; g <- x[x$scenario != "baseline", ]
  data.frame(analysis = r, scenario = g$scenario, threshold = g$threshold, exact_gain = g$nb - b$nb[match(g$threshold, b$threshold)])
}))
samp_sum$exact_gain <- exact_fixed$exact_gain[match(paste(samp_sum$analysis, samp_sum$scenario, samp_sum$threshold), paste(exact_fixed$analysis, exact_fixed$scenario, exact_fixed$threshold))]
samp_sum$study_admissions <- n_stratum

dir.create("output", showWarnings = FALSE)
write.csv(evidence, "output/genetic_donor_evidence.csv", row.names = FALSE)
write.csv(summ, "output/genetic_donor_stratum_summary.csv", row.names = FALSE)
write.csv(baseline_stratum, "output/genetic_donor_baseline_stratum.csv", row.names = FALSE)
write.csv(samp_sum, "output/genetic_donor_sampling_summary.csv", row.names = FALSE)
message("\nDonor evidence:"); show_df(evidence)
message("\nBaseline in the stratum:"); show_df(baseline_stratum)
message("\nDonor-genotype gain in the stratum (per allo admission) and whole-cohort equivalent:"); show_df(summ)
message("\nSampling variability of a stratum-sized study:"); show_df(samp_sum)
