# Aim 4: the genetic-architecture requirements sweep (addendum 16). Run from the project root, after scripts 11 and 16.
#
# 1 to 20 loci x odds ratio 1.1 to 3.0 (by 0.1) x risk-allele frequency 0.05 to 0.50 (by 0.05), identical independent loci, exact convolution of the genotype
# distributions, on the base cohort's empirical profile distribution with the refit baseline. Net-benefit gain of the ceiling score at 5% (primary) and 20%
# (secondary) against the MCID (20% more true cases at no extra false positives of the baseline). Robustness: the published equation as the baseline, and the
# 150-member correlation ensemble on a coarse subset of the grid.
#
# Writes   output/aim4_sweep.csv, aim4_frontier.csv, aim4_sweep_published.csv, aim4_sweep_ensemble.csv, output/figures/aim4_sweep_5pct.png, aim4_sweep_20pct.png

source("R/baseline.R"); source("R/bs_data.R"); source("R/decision_curve.R"); source("R/boscore_data.R"); source("R/refit_baseline.R")
source("R/genetic_layer.R"); source("R/genetic_helpers.R"); source("R/aim4_sweep.R")
suppressPackageStartupMessages(library(ggplot2))

at <- c(0.05, 0.20)
base <- select_cohort(load_base_cohort(), include_twins = FALSE)
q <- tabulate(profile_index(base$X), 2L^length(PREDICTORS)) / nrow(base$X)
lp_refit <- profile_logit(load_refit_params()); lp_pub <- profile_logit(load_baseline_params())
grid <- expand.grid(n_loci = 1:20, odds_ratio = round(seq(1.1, 3.0, by = 0.1), 2), raf = round(seq(0.05, 0.50, by = 0.05), 2))
message("grid: ", nrow(grid), " architectures")

cl <- start_workers(files = c("R/baseline.R", "R/decision_curve.R", "R/genetic_layer.R", "R/aim4_sweep.R"))
run_grid <- function(g, q, lp) {
  ctx <- list(g = g, q = q, lp = lp, at = at)
  parallel_rows(cl, split(seq_len(nrow(g)), rep_len(seq_len(30), nrow(g))), function(ch, ctx)
    do.call(rbind, lapply(ch, function(i) sweep_cell(ctx$q, ctx$lp, ctx$g$n_loci[i], ctx$g$odds_ratio[i], ctx$g$raf[i], ctx$at))), ctx)
}
mcid_of <- function(q, lp) setNames(baseline_mcid(q, lp, at), as.character(at))
add_mcid <- function(d, m) { d$mcid <- unname(m[as.character(d$threshold)]); d$ratio <- d$gain / d$mcid; d$meets_mcid <- d$gain >= d$mcid; d }

sw <- add_mcid(run_grid(grid, q, lp_refit), mcid_of(q, lp_refit))
swp <- add_mcid(run_grid(grid, q, lp_pub), mcid_of(q, lp_pub))
message("MCID (refit): ", paste(sprintf("%.5f", mcid_of(q, lp_refit)), collapse = " and "), "; (published): ", paste(sprintf("%.5f", mcid_of(q, lp_pub)), collapse = " and "))

# Frontiers: the smallest number of loci at each OR and frequency, and the smallest OR at each number of loci and frequency
frontier <- function(d, label) {
  k <- do.call(rbind, lapply(split(d, list(d$odds_ratio, d$raf, d$threshold), drop = TRUE), function(x) { x <- x[order(x$n_loci), ]
    data.frame(baseline = label, kind = "min_loci", odds_ratio = x$odds_ratio[1], raf = x$raf[1], n_loci = NA_real_, threshold = x$threshold[1], value = first_reaching(x$n_loci, x$gain, x$mcid[1])) }))
  o <- do.call(rbind, lapply(split(d, list(d$n_loci, d$raf, d$threshold), drop = TRUE), function(x) { x <- x[order(x$odds_ratio), ]
    data.frame(baseline = label, kind = "min_odds_ratio", odds_ratio = NA_real_, raf = x$raf[1], n_loci = x$n_loci[1], threshold = x$threshold[1], value = first_reaching(x$odds_ratio, x$gain, x$mcid[1])) }))
  rbind(k, o)
}
fr <- rbind(frontier(sw, "refit"), frontier(swp, "published"))
rownames(fr) <- NULL

# The correlation ensemble on a coarse subset
q_members <- as.matrix(read.csv("output/correlation_ensemble_profile_probs.csv")[, -1])
coarse <- expand.grid(n_loci = c(1, 2, 5, 10, 20), odds_ratio = c(1.2, 1.5, 2.0, 3.0), raf = c(0.1, 0.3, 0.5))
ens <- do.call(rbind, parallel::parLapply(cl, split(seq_len(nrow(q_members)), rep_len(seq_len(30), nrow(q_members))), function(ms, ctx) {
  do.call(rbind, lapply(ms, function(m) {
    mc <- setNames(baseline_mcid(ctx$q[m, ], ctx$lp, ctx$at), as.character(ctx$at))
    d <- do.call(rbind, lapply(seq_len(nrow(ctx$g)), function(i) sweep_cell(ctx$q[m, ], ctx$lp, ctx$g$n_loci[i], ctx$g$odds_ratio[i], ctx$g$raf[i], ctx$at)))
    d$member <- m; d$mcid <- unname(mc[as.character(d$threshold)]); d
  }))
}, list(q = q_members, lp = lp_refit, g = coarse, at = at)))
ens$ratio <- ens$gain / ens$mcid
ens_sum <- by_key(ens, c("n_loci", "odds_ratio", "raf", "threshold"), function(d) data.frame(
  ratio_min = min(d$ratio), ratio_median = median(d$ratio), ratio_max = max(d$ratio), share_of_members_reaching_mcid = mean(d$gain >= d$mcid)))
parallel::stopCluster(cl)

dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)
write.csv(sw, "output/aim4_sweep.csv", row.names = FALSE)
write.csv(fr, "output/aim4_frontier.csv", row.names = FALSE)
write.csv(swp, "output/aim4_sweep_published.csv", row.names = FALSE)
write.csv(ens_sum, "output/aim4_sweep_ensemble.csv", row.names = FALSE)

plot_sweep <- function(d, th, file) {
  x <- d[abs(d$threshold - th) < 1e-9, ]
  p <- ggplot(x, aes(n_loci, odds_ratio)) +
    geom_raster(aes(fill = ratio)) +
    geom_contour(aes(z = ratio), breaks = 1, colour = "white", linewidth = 0.5) +
    facet_wrap(~ paste0("risk-allele frequency ", raf), nrow = 2) +
    scale_fill_viridis_c(trans = "log10", name = "gain / MCID") +
    labs(x = "Number of loci", y = "Odds ratio per risk allele")
  ggsave(file, p, width = 11, height = 5.5, dpi = 150)
}
plot_sweep(sw, 0.05, "output/figures/aim4_sweep_5pct.png")
plot_sweep(sw, 0.20, "output/figures/aim4_sweep_20pct.png")

message("\nGain at the corners (5%): "); show_df(sw[sw$threshold == 0.05 & sw$n_loci %in% c(1, 20) & sw$odds_ratio %in% c(1.1, 1.5, 3.0) & sw$raf %in% c(0.05, 0.3, 0.5), c("n_loci", "odds_ratio", "raf", "gain", "extra_cases_per_1000", "auc", "ratio")])
message("\nSmallest number of loci reaching the MCID (refit baseline, 5%):")
w <- fr[fr$baseline == "refit" & fr$kind == "min_loci" & fr$threshold == 0.05, ]
print(xtabs(value ~ odds_ratio + raf, data = transform(w, value = ifelse(is.na(value), 0, value))))
