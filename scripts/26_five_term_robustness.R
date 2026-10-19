# Does CLEC7A's lead hold? The five-term Aim 3b layer and the Shapley attribution on the sensitivity cohorts and on the correlation ensemble (addendum 18).
# Run from the project root, after scripts 16, 19, 22 and 24. Set ROBUST_ONLY to a subset of the setting names to run them one at a time; per-setting files are cached in output/ and merged.
#
# Settings: published_4078, published_4127, refit_4127 (as script 19: each with its own baseline and its own MCIDs) and ensemble (the refit baseline on the base cohort, with a
# clinical correlation structure drawn from the 150 near-optimal members, 20 or so draws of the effects per member, as script 17).
# 1,000 draws of the true effects per reading (pooled mean, predictive), the second CLEC7A parameter included, as in script 24 and the five-term run of script 21.
# The interaction test of Aim 3b is not repeated here (it depends on the cohort rows and was reported for the primary cohort).
#
# Writes   output/robust5_{aim3b,attribution}_<setting>.csv and the merged output/robust5_aim3b_summary.csv, robust5_attribution_summary.csv

source("R/baseline.R"); source("R/bs_data.R"); source("R/copula.R"); source("R/decision_curve.R"); source("R/boscore_data.R"); source("R/refit_baseline.R")
source("R/evidence_mapping.R"); source("R/ptx3_evidence.R"); source("R/source_verification.R")
source("R/genetic_layer.R"); source("R/genetic_helpers.R"); source("R/aim3b.R"); source("R/robustness.R")

at <- c(0.05, 0.20)
all_settings <- c("published_4078", "published_4127", "refit_4127", "ensemble")
run_these <- if (exists("ROBUST_ONLY")) ROBUST_ONLY else all_settings
n_draw <- 1000L

full <- load_base_cohort()
cohorts <- list(`4078` = select_cohort(full, include_twins = FALSE), `4127` = select_cohort(full, include_twins = TRUE))
published <- load_baseline_params()
lp_refit <- profile_logit(load_refit_params())
spec <- list(
  published_4078 = list(cohort = cohorts[["4078"]], lp = profile_logit(published)),
  published_4127 = list(cohort = cohorts[["4127"]], lp = profile_logit(published)),
  refit_4127 = list(cohort = cohorts[["4127"]], lp = profile_logit(refit_params(fit_baseline(cohorts[["4127"]])$coef))),
  ensemble = list(cohort = cohorts[["4078"]], lp = lp_refit))
profile_q <- function(co) tabulate(profile_index(co$X), 2L^length(PREDICTORS)) / nrow(co$X)

setup <- clec7a_second_setup(n_draw = n_draw)
terms <- setup$terms

# Brunel's stratum effects as fractions of its whole-cohort effect (as script 21)
vr <- read_bundle_table("extraction_variants_round1.csv"); st <- read_bundle_table("extraction_studies_round1.csv")
vr <- apply_corrections(vr, read.csv("data/extraction_corrections.csv", stringsAsFactors = FALSE))$data
hap <- block_haplotypes("PTX3_trio", risk_alleles = c("A", NA, NA), source = "gnomAD_NFE")
brow <- vr[which(vr$study_id == "S029" & vr$rsid %in% "rs2305619" & grepl("^All episodes \\(N=212\\)|^No neutropenia at presentation \\(N=133\\)|^Neutropenia at presentation \\(N=79\\)", vr$subgroup)), ]
bb <- vapply(seq_len(nrow(brow)), function(i) map_extraction_row(brow[i, ], hap, PTX3_LOCI, "trend")$beta, numeric(1))
names(bb) <- ifelse(grepl("^All", brow$subgroup), "all", ifelse(grepl("^No neutropenia", brow$subgroup), "non", "neut"))
r_values <- list(`H_F only` = Inf, `H_F, ratio 4` = 4, homogeneous = 1, `Brunel direction, ratio 1/4` = 0.25, `Brunel only` = 0,
                 `Brunel estimate` = c(N = unname(bb["neut"] / bb["all"]), NN = unname(bb["non"] / bb["all"])))

q_members <- as.matrix(read.csv("output/correlation_ensemble_profile_probs.csv")[, -1])
cl <- start_workers(files = c("R/baseline.R", "R/decision_curve.R", "R/genetic_layer.R", "R/aim3b.R", "R/robustness.R"))

for (s in run_these) {
  sp <- spec[[s]]
  qs <- if (s == "ensemble") q_members else matrix(profile_q(sp$cohort), nrow = 1)
  member_of_draw <- rep_len(seq_len(nrow(qs)), n_draw)
  tasks <- if (s == "ensemble") lapply(seq_len(nrow(qs)), function(m) list(member = m, idx = which(member_of_draw == m))) else
    lapply(split(seq_len(n_draw), rep_len(seq_len(30), n_draw)), function(ix) list(member = 1L, idx = ix))
  message(s, ": ", nrow(qs), " profile distribution(s), ", length(tasks), " tasks")
  lay <- do.call(rbind, lapply(names(setup$effects), function(h) {
    ctx <- list(q = qs, lp = sp$lp, grid = setup$grid, eff = setup$effects[[h]], bh = setup$beta_hat, r = r_values, at = at)
    cbind(truth = h, parallel_rows(cl, tasks, function(tk, ctx) { r <- aim3b_member_rows(ctx$q[tk$member, ], ctx$lp, ctx$grid, ctx$eff, tk$idx, ctx$bh, ctx$r, ctx$at); r$member <- tk$member; r }, ctx))
  }))
  att <- do.call(rbind, lapply(names(setup$effects), function(h) {
    ctx <- list(q = qs, lp = sp$lp, grid = setup$grid, eff = setup$effects[[h]], bh = setup$beta_hat, at = at)
    cbind(truth = h, parallel_rows(cl, tasks, function(tk, ctx) { r <- attribution_member_rows(ctx$q[tk$member, ], ctx$lp, ctx$grid, ctx$eff, tk$idx, ctx$bh, ctx$at); r$member <- tk$member; r }, ctx))
  }))
  write.csv(cbind(setting = s, summarise_aim3b_rows(lay)), sprintf("output/robust5_aim3b_%s.csv", s), row.names = FALSE)
  write.csv(cbind(setting = s, summarise_attribution_rows(att, terms)), sprintf("output/robust5_attribution_%s.csv", s), row.names = FALSE)
  a <- read.csv(sprintf("output/robust5_attribution_%s.csv", s))
  message("\n", s, ": Shapley shares (CLEC7A combined, CD209, PTX3, TLR4) and the share of draws in which CLEC7A is the largest contributor")
  show_df(a[, c("truth", "score", "threshold", "total", "share_CLEC7A_combined", "share_CD209", "share_PTX3", "share_TLR4", "p_clec7a_combined_largest")])
}
parallel::stopCluster(cl)

have <- all_settings[file.exists(sprintf("output/robust5_attribution_%s.csv", all_settings))]
for (n in c("aim3b", "attribution")) write.csv(do.call(rbind, lapply(have, function(s) read.csv(sprintf("output/robust5_%s_%s.csv", n, s)))), sprintf("output/robust5_%s_summary.csv", n), row.names = FALSE)
