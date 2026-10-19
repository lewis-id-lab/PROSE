# Aim 3b: the genetic layer split by neutropenia status (addendum 15). Run from the project root, after scripts 08, 11, 14 and 16.
# The prediction and the design were committed before this script was run.
#
# Writes   output/aim3b_strata.csv                 shares, event rates and the neutropenic stratum's event share w
#          output/aim3b_brunel_modification.csv    Brunel 2018's stratified per-copy effects and the difference between strata
#          output/aim3b_layer_summary.csv          net-benefit gain by truth, scenario, score, level (stratum or whole cohort) and threshold
#          output/aim3b_interaction_summary.csv    the falsification test: the interaction fit in simulated studies of the cohort's size

source("R/baseline.R"); source("R/bs_data.R"); source("R/decision_curve.R"); source("R/boscore_data.R")
source("R/refit_baseline.R"); source("R/evidence_mapping.R"); source("R/ptx3_evidence.R"); source("R/source_verification.R")
source("R/genetic_layer.R"); source("R/genetic_helpers.R"); source("R/aim3b.R")

at <- c(0.05, 0.20)
base <- select_cohort(load_base_cohort(), include_twins = FALSE)
lp <- profile_logit(load_refit_params())
q <- tabulate(profile_index(base$X), 2L^length(PREDICTORS)) / nrow(base$X)
strata <- neutropenia_strata(q, lp)
# Set AIM3B_SECOND_CLEC7A <- TRUE before sourcing to run the analysis with the second CLEC7A parameter (addendum 17): five terms, joint CLEC7A genotypes,
# outputs with the suffix _second_clec7a.
second <- exists("AIM3B_SECOND_CLEC7A") && isTRUE(AIM3B_SECOND_CLEC7A)
suffix <- if (second) "_second_clec7a" else ""
p <- genetic_term_frequencies(); grid <- dosage_grid(p); terms <- names(GENETIC_TERMS)
if (second) { setup <- clec7a_second_setup(n_draw = 1000L); grid <- setup$grid; terms <- setup$terms }
cells <- lapply(list(N = strata$N, NN = strata$NN), function(s) trim_cells(genetic_cells(s$q, lp, grid)))
mcid_tab <- read.csv("output/dca_baseline_refit_summary.csv")
mcid_whole <- setNames(mcid_tab$mcid_b[match(at, mcid_tab$threshold)], as.character(at))
mcid_stratum <- list(N = setNames(stratum_mcid(strata$N, lp, at), as.character(at)), NN = setNames(stratum_mcid(strata$NN, lp, at), as.character(at)))

strata_tab <- data.frame(stratum = c("neutropenic", "non-neutropenic"), share = c(strata$N$share, strata$NN$share),
                         baseline_event_rate = c(strata$N$rate, strata$NN$rate), event_share = c(strata$w, 1 - strata$w),
                         observed_events = c(sum(base$rows$case[base$X[, "prolonged_neutropenia"] == 1]), sum(base$rows$case[base$X[, "prolonged_neutropenia"] == 0])),
                         observed_admissions = c(sum(base$X[, "prolonged_neutropenia"] == 1), sum(base$X[, "prolonged_neutropenia"] == 0)),
                         mcid_5 = c(mcid_stratum$N[["0.05"]], mcid_stratum$NN[["0.05"]]), mcid_20 = c(mcid_stratum$N[["0.2"]], mcid_stratum$NN[["0.2"]]))
message("strata:"); show_df(strata_tab)

# Brunel 2018: stratified per-copy effects (rows for the whole cohort, no neutropenia at presentation, neutropenia at presentation)
vr <- read_bundle_table("extraction_variants_round1.csv"); st <- read_bundle_table("extraction_studies_round1.csv")
vr <- apply_corrections(vr, read.csv("data/extraction_corrections.csv", stringsAsFactors = FALSE))$data
hap <- block_haplotypes("PTX3_trio", risk_alleles = c("A", NA, NA), source = "gnomAD_NFE")
brow <- vr[vr$study_id == "S029" & vr$rsid %in% c("rs2305619", "rs3816527") & grepl("^All episodes \\(N=212\\)|^No neutropenia at presentation \\(N=133\\)|^Neutropenia at presentation \\(N=79\\)", vr$subgroup), ]
bm <- do.call(rbind, lapply(seq_len(nrow(brow)), function(i) {
  m <- map_extraction_row(brow[i, ], hap, PTX3_LOCI, "trend")
  data.frame(rsid = brow$rsid[i], stratum = ifelse(grepl("^All", brow$subgroup[i]), "all", ifelse(grepl("^No neutropenia", brow$subgroup[i]), "non-neutropenic", "neutropenic")),
             log_or = m$y, se = m$se, beta_per_copy = m$beta, se_beta = m$se_beta)
}))
brunel <- do.call(rbind, lapply(unique(bm$rsid), function(s) {
  x <- bm[bm$rsid == s, ]; n <- x[x$stratum == "neutropenic", ]; nn <- x[x$stratum == "non-neutropenic", ]
  d <- n$beta_per_copy - nn$beta_per_copy; se <- sqrt(n$se_beta^2 + nn$se_beta^2)
  data.frame(rsid = s, beta_neutropenic = n$beta_per_copy, se_neutropenic = n$se_beta, beta_non_neutropenic = nn$beta_per_copy, se_non_neutropenic = nn$se_beta,
             difference = d, se_difference = se, z = d / se, p_two_sided = 2 * pnorm(-abs(d / se)), lo95 = d - 1.96 * se, hi95 = d + 1.96 * se,
             lo90 = d - 1.645 * se, hi90 = d + 1.645 * se, ratio_point = n$beta_per_copy / nn$beta_per_copy)
}))
message("\nBrunel 2018, stratified per-copy effects (PTX3 h2):"); show_df(brunel)
# The data-based scenario applies Brunel's stratum effects directly, as fractions of its whole-cohort effect (rs2305619). The planned
# ratio-constrained version is infeasible: the ratio is -0.11 and the Bologna event share w = 0.80 would need an effect about nine times
# the pooled one in the non-neutropenic stratum (deviation from addendum 15, recorded in addendum 15's results).
b_all <- bm$beta_per_copy[bm$rsid == "rs2305619" & bm$stratum == "all"]
brunel_fractions <- c(N = brunel$beta_neutropenic[brunel$rsid == "rs2305619"] / b_all, NN = brunel$beta_non_neutropenic[brunel$rsid == "rs2305619"] / b_all)
r_values <- list(`H_F only` = Inf, `H_F, ratio 4` = 4, homogeneous = 1, `Brunel direction, ratio 1/4` = 0.25, `Brunel only` = 0, `Brunel estimate` = brunel_fractions)
message("\nscenarios: ", paste(names(r_values), vapply(r_values, function(x) if (length(x) == 2) sprintf("fractions N %.2f, NN %.2f", x[["N"]], x[["NN"]]) else sprintf("ratio %.2f", x), ""), sep = " = ", collapse = "; "), "; event share in neutropenic admissions w = ", sprintf("%.3f", strata$w))

pooled <- load_term_draws(); beta_hat <- pooled_estimates(pooled)
beta_hi <- vapply(pooled, function(d) unname(stats::quantile(d$mu, 0.95)), numeric(1))
effects_file <- read.csv("output/genetic_layer_effects.csv")
n_draw <- 1000L
effects <- list(mean = as.matrix(effects_file[effects_file$truth == "mean" & effects_file$draw <= n_draw, c("PTX3", "CLEC7A", "TLR4", "CD209")]),
                predictive = as.matrix(effects_file[effects_file$truth == "predictive" & effects_file$draw <= n_draw, c("PTX3", "CLEC7A", "TLR4", "CD209")]))
if (second) { effects <- setup$effects; beta_hat <- setup$beta_hat; beta_hi <- setup$beta_hi }

cl <- start_workers(files = c("R/baseline.R", "R/decision_curve.R", "R/genetic_layer.R", "R/aim3b.R"))
lay <- do.call(rbind, lapply(names(effects), function(h) {
  ctx <- list(strata = strata, cells = cells, eff = effects[[h]], bh = beta_hat, r = r_values, at = at)
  cbind(truth = h, parallel_rows(cl, chunk_ids(n_draw), function(ch, ctx) run_aim3b(ctx$strata, ctx$cells, ctx$eff[ch, , drop = FALSE], ctx$bh, ctx$r, ctx$at, draw_id = ch), ctx))
}))

mc <- function(level, th) switch(level, whole = mcid_whole[[as.character(th)]], N = mcid_stratum$N[[as.character(th)]], NN = mcid_stratum$NN[[as.character(th)]])
lay$mcid <- mapply(mc, lay$level, lay$threshold)
layer_sum <- by_key(lay[lay$score != "baseline", ], c("truth", "scenario", "level", "score", "threshold"), function(d) data.frame(
  gain = median(d$delta_nb), lo = quantile(d$delta_nb, 0.05, names = FALSE), hi = quantile(d$delta_nb, 0.95, names = FALSE),
  extra_cases_per_1000 = 1000 * median(d$delta_tp), p_meets_mcid = mean(d$delta_nb >= d$mcid), mcid = d$mcid[1]))
# The neutropenic stratum's share of the whole-cohort gain (median over draws with a positive whole-cohort gain)
w_id <- paste(lay$truth, lay$scenario, lay$score, lay$threshold, lay$draw)
wide <- lay[lay$level == "whole", c("truth", "scenario", "score", "threshold", "draw", "delta_nb")]
ndf <- lay[lay$level == "N", ]; ndf$key <- paste(ndf$truth, ndf$scenario, ndf$score, ndf$threshold, ndf$draw)
wide$gain_N <- strata$N$share * ndf$delta_nb[match(paste(wide$truth, wide$scenario, wide$score, wide$threshold, wide$draw), ndf$key)]
wide <- wide[wide$score != "baseline" & wide$delta_nb > 0, ]
share_N <- by_key(wide, c("truth", "scenario", "score", "threshold"), function(d) data.frame(share_of_gain_from_neutropenic = median(d$gain_N / d$delta_nb), n_positive = nrow(d)))
layer_sum <- merge(layer_sum, share_N, by = c("truth", "scenario", "score", "threshold"), all.x = TRUE)
layer_sum$share_of_gain_from_neutropenic[layer_sum$level != "whole"] <- NA

# ---- The falsification test in simulated studies -------------------------------------------------------------------------------------
cg <- cohort_groups(base); neut <- base$X[, "prolonged_neutropenia"] == 1
n_sim <- 1000L
levels_tab <- list(medians = beta_hat, upper = beta_hi)
inter <- do.call(rbind, lapply(names(levels_tab), function(lv) do.call(rbind, lapply(names(r_values), function(sc) {
  tr <- scenario_effects(levels_tab[[lv]], r_values[[sc]], strata$w)
  sh <- list(N = genetic_scores(cells$N, tr$N, tr$N)$a_true, NN = genetic_scores(cells$NN, tr$NN, tr$NN)$a_true)
  ctx <- list(cg = cg, neut = neut, lp = lp, p = p, grid = if (second) grid else NULL, shifts = sh, tN = tr$N, tNN = tr$NN, bh = beta_hat, seed = 3706L + nchar(lv) * 100L + match(sc, names(r_values)))
  res <- parallel::parLapply(cl, chunk_ids(n_sim, 15L), function(ch, ctx) {
    set.seed(ctx$seed + ch[1])
    smp <- if (!is.null(ctx$grid)) function(n) ctx$grid$G[sample.int(nrow(ctx$grid$G), n, replace = TRUE, prob = ctx$grid$prob), , drop = FALSE] else NULL
    t(vapply(ch, function(i) simulate_interaction_study(ctx$cg, ctx$neut, ctx$lp, ctx$p, ctx$shifts, ctx$tN, ctx$tNN, ctx$bh, sampler = smp), numeric(5)))
  }, ctx)
  r <- as.data.frame(do.call(rbind, res))
  data.frame(effect_level = lv, scenario = sc, true_difference_in_score_units = scenario_difference(r_values[[sc]], strata$w),
             studies = nrow(r), failed_fits = sum(is.na(r$delta)),
             median_delta = median(r$delta, na.rm = TRUE), sd_delta = sd(r$delta, na.rm = TRUE), median_se = median(r$se, na.rm = TRUE),
             p_delta_positive = mean(r$delta > 0, na.rm = TRUE),
             p_sig_positive = mean(r$p < 0.05 & r$delta > 0, na.rm = TRUE), p_sig_negative = mean(r$p < 0.05 & r$delta < 0, na.rm = TRUE),
             median_events = median(r$events), median_events_neutropenic = median(r$events_neutropenic))
}))))
parallel::stopCluster(cl)

dir.create("output", showWarnings = FALSE)
if (!second) {
  write.csv(strata_tab, "output/aim3b_strata.csv", row.names = FALSE)
  write.csv(brunel, "output/aim3b_brunel_modification.csv", row.names = FALSE)
}
write.csv(layer_sum, paste0("output/aim3b_layer_summary", suffix, ".csv"), row.names = FALSE)
write.csv(inter, paste0("output/aim3b_interaction_summary", suffix, ".csv"), row.names = FALSE)
message("\nGain by scenario (whole cohort, mean reading, 5%):")
show_df(layer_sum[layer_sum$truth == "mean" & layer_sum$level == "whole" & layer_sum$threshold == 0.05, c("scenario", "score", "gain", "lo", "hi", "extra_cases_per_1000", "p_meets_mcid", "share_of_gain_from_neutropenic")])
message("\nThe interaction test:"); show_df(inter)
