# The NRI comparison under the three comparisons (addendum 24, continued from script 31). Run from the project root, after script 29.
#
# As script 29 (the four-term panel, ceiling and realised, and the PTX3 structures, over the same 2,000 draws of the true effects, pooled-mean
# reading, seed 6417), but with the truth, the comparator and the scores built under each comparison of R/calibration_check.R:
#   current    global-intercept truth; comparator expit(eta_X)      (reproduces script 29)
#   marginal   the same truth; comparator m(X) = E[p | X]
#   profile    profile-preserving truth; comparator expit(eta_X); scores with an intercept per profile
# Under "marginal" and "profile" the comparator is calibrated against the simulated truth, so reclassification counts only what the
# genetic score adds to a calibrated clinical score.
#
# Writes   output/calibration_check_nri.csv        per comparison, scenario and metric: value at the pooled estimates, median and 90% range
#          output/calibration_check_nri_grid.csv   per comparison and structure, per-copy log OR: NRI and net-benefit gain

source("R/baseline.R"); source("R/bs_data.R"); source("R/decision_curve.R"); source("R/boscore_data.R"); source("R/refit_baseline.R")
source("R/evidence_mapping.R"); source("R/genetic_layer.R"); source("R/genetic_helpers.R"); source("R/nri.R"); source("R/calibration_check.R")

n_draw <- 2000L
comparisons <- c("current", "marginal", "profile")
base <- select_cohort(load_base_cohort(), include_twins = FALSE)
lp <- profile_logit(load_refit_params())
q <- tabulate(profile_index(base$X), 2L^length(PREDICTORS)) / nrow(base$X)

draws <- load_term_draws()
beta_hat <- pooled_estimates(draws)
effects <- draw_term_effects(draws, n_draw, "mean", seed = 6417L)

grid4 <- dosage_grid(genetic_term_frequencies())
gP <- ptx3_diplotype_grid()
snp_freq <- setNames(as.numeric(crossprod(gP$G, gP$prob))[2:4] / 2, colnames(gP$G)[2:4])
gI <- dosage_grid(snp_freq)
zero <- c(h2 = 0, rs2305619 = 0, rs3816527 = 0, rs1840680 = 0)
snps <- c("rs2305619", "rs3816527", "rs1840680")

# All scenarios, all comparisons, for one vector of true effects b: rows "comparison | scenario"
scenario_metrics <- function(b, ctx) {
  bp <- b[["PTX3"]]
  parts <- list(
    "four-term panel, ceiling"     = nri_by_scheme(ctx$q, ctx$lp, ctx$grid4, b, b),
    "four-term panel, realised"    = nri_by_scheme(ctx$q, ctx$lp, ctx$grid4, b, ctx$beta_hat),
    "PTX3 h2 haplotype"            = nri_by_scheme(ctx$q, ctx$lp, ctx$gP, replace(ctx$zero, "h2", bp), replace(ctx$zero, "h2", bp)),
    "PTX3 three SNPs, real LD"     = nri_by_scheme(ctx$q, ctx$lp, ctx$gP, replace(ctx$zero, "h2", bp), replace(ctx$zero, ctx$snps, bp)),
    "PTX3 three SNPs, independent" = nri_by_scheme(ctx$q, ctx$lp, ctx$gI, setNames(rep(bp, 3), ctx$snps), setNames(rep(bp, 3), ctx$snps)))
  out <- do.call(rbind, lapply(names(parts), function(s) { m <- parts[[s]]; rownames(m) <- paste(rownames(m), s, sep = " | "); m }))
  out
}
ctx <- list(q = q, lp = lp, grid4 = grid4, gP = gP, gI = gI, beta_hat = beta_hat, zero = zero, snps = snps)

point <- scenario_metrics(beta_hat, ctx)

cl <- start_workers(cores = 12, files = c("R/baseline.R", "R/decision_curve.R", "R/genetic_layer.R", "R/nri.R", "R/calibration_check.R"))
parallel::clusterExport(cl, "scenario_metrics")
t0 <- Sys.time()
per_draw <- unlist(parallel::parLapply(cl, split(seq_len(n_draw), rep_len(seq_len(60), n_draw)), function(ix, ctx, effects) {
  lapply(ix, function(i) scenario_metrics(effects[i, ], ctx))
}, ctx, effects), recursive = FALSE)
message(sprintf("draws: %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
# Restore draw order (chunks interleave the draws)
per_draw <- per_draw[order(unlist(split(seq_len(n_draw), rep_len(seq_len(60), n_draw))))]
arr <- simplify2array(per_draw)

summary_tab <- do.call(rbind, lapply(rownames(point), function(rn) {
  parts <- strsplit(rn, " | ", fixed = TRUE)[[1]]
  data.frame(comparison = parts[1], scenario = parts[2], metric = colnames(point), at_pooled = point[rn, ],
             median = apply(arr[rn, , ], 1, median),
             lo = apply(arr[rn, , ], 1, quantile, probs = 0.05, names = FALSE),
             hi = apply(arr[rn, , ], 1, quantile, probs = 0.95, names = FALSE), row.names = NULL)
}))

# The current comparison must reproduce script 29
ref <- read.csv("output/nri_comparison.csv")
cur <- summary_tab[summary_tab$comparison == "current", ]
j <- match(paste(cur$scenario, cur$metric), paste(ref$scenario, ref$metric))
stopifnot(!anyNA(j), nrow(cur) == nrow(ref),
          max(abs(cur$at_pooled - ref$at_pooled[j])) < 1e-8, max(abs(cur$median - ref$median[j])) < 1e-8)

# NRI and net benefit against the per-copy PTX3 effect
grid_beta <- seq(0, 1.6, by = 0.02)
effect_grid <- do.call(rbind, parallel::parLapply(cl, grid_beta, function(bp, ctx) {
  mk <- function(structure, m) data.frame(comparison = rownames(m), structure = structure, log_or = bp, m, row.names = NULL)
  rbind(mk("PTX3 h2 haplotype", nri_by_scheme(ctx$q, ctx$lp, ctx$gP, replace(ctx$zero, "h2", bp), replace(ctx$zero, "h2", bp))),
        mk("PTX3 three SNPs, independent", nri_by_scheme(ctx$q, ctx$lp, ctx$gI, setNames(rep(bp, 3), ctx$snps), setNames(rep(bp, 3), ctx$snps))))
}, ctx))
parallel::stopCluster(cl)

ref_grid <- read.csv("output/nri_effect_grid.csv")
cg <- effect_grid[effect_grid$comparison == "current", ]
stopifnot(nrow(cg) == nrow(ref_grid), max(abs(cg$nri_cat5 - ref_grid$nri_cat5)) < 1e-8, max(abs(cg$dnb_05 - ref_grid$dnb_05)) < 1e-8)

dir.create("output", showWarnings = FALSE)
write.csv(summary_tab, "output/calibration_check_nri.csv", row.names = FALSE)
write.csv(effect_grid, "output/calibration_check_nri_grid.csv", row.names = FALSE)

message("\nAt the pooled estimates, by comparison:")
print(summary_tab[summary_tab$metric %in% c("nri_cat5", "nri_cat3", "nri_cont", "dnb_05", "dnb_20"),
                  c("scenario", "comparison", "metric", "at_pooled", "median", "lo", "hi")], row.names = FALSE, digits = 3)
