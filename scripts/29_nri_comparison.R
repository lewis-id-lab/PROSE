# Tertiary NRI comparison (addendum 22): the haplotype panel against the predecessor's three-SNP independent model. Run from the
# project root, after scripts 08, 11, 12, 14, 15 and 16.
#
# The predecessor simulated the three PTX3 SNPs (rs2305619, rs3816527, rs1840680) as independent multiplicative risks and reported a
# net reclassification improvement of 22.5%. Its per-SNP effects are not on record (only a homozygote OR of 15.2 for rs1840680), so
# this script isolates the linkage error by giving every PTX3 structure the same per-copy effect, the pooled PTX3 h2 effect, and
# asks what each structure does to NRI and to net benefit. All quantities are exact over profile x genotype cells; truths are
# recentred to the baseline prevalence and scores to the baseline mean risk (D5), as in the genetic layer.
#
# Scenarios (refit baseline, base cohort):
#   four-term panel, ceiling      the primary panel, true effects in the score (as script 16)
#   four-term panel, realised     pooled medians in the score
#   PTX3 h2 haplotype             truth and score on h2, the correct structure
#   PTX3 three SNPs, real LD      truth on h2; the score adds the effect once per SNP, on real (linked) chromosomes: what the
#                                 three-SNP model does when applied to real patients
#   PTX3 three SNPs, independent  the predecessor's simulation: the three SNPs drawn independently, each with the effect, truth and
#                                 score alike
#
# NRI: categorical at 5% (two categories, the decision threshold), three categories (cuts at 5% and 20%), and continuous.
# Uncertainty: the same 2,000 draws of the true effects as script 16 (pooled-mean reading, seed 6417).
# Effect grid: NRI and net-benefit gain against the per-copy PTX3 effect for the h2 and the independent three-SNP structures.
#
# Writes   output/nri_comparison.csv    per scenario and metric: value at the pooled estimates, median and 90% range over draws
#          output/nri_effect_grid.csv   per structure and per-copy log OR: NRI and net-benefit gain

source("R/baseline.R")
source("R/bs_data.R")
source("R/decision_curve.R")
source("R/boscore_data.R")
source("R/refit_baseline.R")
source("R/evidence_mapping.R")
source("R/genetic_layer.R")
source("R/nri.R")

n_draw <- 2000L

base <- select_cohort(load_base_cohort(), include_twins = FALSE)
lp_profile <- profile_logit(load_refit_params())
q <- tabulate(profile_index(base$X), 2L^length(PREDICTORS)) / nrow(base$X)

draws <- load_term_draws()
beta_hat <- pooled_estimates(draws)
effects <- draw_term_effects(draws, n_draw, "mean", seed = 6417L)

cells4 <- genetic_cells(q, lp_profile, dosage_grid(genetic_term_frequencies()))
gP <- ptx3_diplotype_grid()
cellsP <- genetic_cells(q, lp_profile, gP)
snp_freq <- setNames(as.numeric(crossprod(gP$G, gP$prob))[2:4] / 2, colnames(gP$G)[2:4])
cellsI <- genetic_cells(q, lp_profile, dosage_grid(snp_freq))

message("PTX3 risk-allele frequencies (gnomAD NFE): ", paste(names(snp_freq), sprintf("%.3f", snp_freq), collapse = ", "))

# All scenarios for one vector of true effects b (named by term)
scenario_metrics <- function(b) {
  s4 <- genetic_scores(cells4, b, beta_hat)
  bp <- b[["PTX3"]]
  zero <- c(h2 = 0, rs2305619 = 0, rs3816527 = 0, rs1840680 = 0)
  sP <- genetic_scores(cellsP, replace(zero, "h2", bp), replace(zero, c("rs2305619", "rs3816527", "rs1840680"), bp))
  sI <- genetic_scores(cellsI, setNames(rep(bp, 3), names(snp_freq)), setNames(rep(bp, 3), names(snp_freq)))
  rbind(
    "four-term panel, ceiling"     = reclassification_metrics(cells4$w, s4$lp_true, cells4$lp_base, s4$lp_ceiling),
    "four-term panel, realised"    = reclassification_metrics(cells4$w, s4$lp_true, cells4$lp_base, s4$lp_real),
    "PTX3 h2 haplotype"            = reclassification_metrics(cellsP$w, sP$lp_true, cellsP$lp_base, sP$lp_ceiling),
    "PTX3 three SNPs, real LD"     = reclassification_metrics(cellsP$w, sP$lp_true, cellsP$lp_base, sP$lp_real),
    "PTX3 three SNPs, independent" = reclassification_metrics(cellsI$w, sI$lp_true, cellsI$lp_base, sI$lp_ceiling)
  )
}

point <- scenario_metrics(beta_hat)
per_draw <- lapply(seq_len(n_draw), function(i) scenario_metrics(effects[i, ]))
arr <- simplify2array(per_draw)

summary_tab <- do.call(rbind, lapply(rownames(point), function(sc) {
  data.frame(scenario = sc, metric = colnames(point), at_pooled = point[sc, ],
             median = apply(arr[sc, , ], 1, median),
             lo = apply(arr[sc, , ], 1, quantile, probs = 0.05, names = FALSE),
             hi = apply(arr[sc, , ], 1, quantile, probs = 0.95, names = FALSE),
             row.names = NULL)
}))

# NRI against the per-copy PTX3 effect, for the correct (h2) and the predecessor's (independent three-SNP) structure
grid_beta <- seq(0, 1.6, by = 0.02)
effect_grid <- do.call(rbind, lapply(grid_beta, function(bp) {
  zero <- c(h2 = 0, rs2305619 = 0, rs3816527 = 0, rs1840680 = 0)
  sP <- genetic_scores(cellsP, replace(zero, "h2", bp), zero)
  sI <- genetic_scores(cellsI, setNames(rep(bp, 3), names(snp_freq)), setNames(rep(bp, 3), names(snp_freq)))
  rbind(
    data.frame(structure = "PTX3 h2 haplotype", log_or = bp,
               t(reclassification_metrics(cellsP$w, sP$lp_true, cellsP$lp_base, sP$lp_ceiling))),
    data.frame(structure = "PTX3 three SNPs, independent", log_or = bp,
               t(reclassification_metrics(cellsI$w, sI$lp_true, cellsI$lp_base, sI$lp_ceiling)))
  )
}))

dir.create("output", showWarnings = FALSE)
write.csv(summary_tab, "output/nri_comparison.csv", row.names = FALSE)
write.csv(effect_grid, "output/nri_effect_grid.csv", row.names = FALSE)

message("\nAt the pooled estimates (PTX3 per-copy log OR ", sprintf("%.3f", beta_hat[["PTX3"]]), "):")
print(round(point[, c("nri_cat5", "nri_cat3", "nri_cont", "dnb_05", "dauc")], 4))
