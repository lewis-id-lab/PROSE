# Shapley attribution of the genetic-layer gain with the second CLEC7A parameter (addendum 17). Run from the project root, after scripts 16 and 22.
# Five terms: PTX3, CLEC7A (Y238X copies), CLEC7A_AG (A-G copies), TLR4 and CD209, with the joint genotype structure of CLEC7A. Same design as the attribution in
# script 17 (1,000 draws per reading, ceiling and realised scores, empirical cohort). The combined CLEC7A share is the sum of the two CLEC7A terms (Shapley values add).
#
# Writes   output/genetic_attribution_second_clec7a_summary.csv

source("R/baseline.R"); source("R/bs_data.R"); source("R/decision_curve.R"); source("R/boscore_data.R"); source("R/refit_baseline.R")
source("R/evidence_mapping.R"); source("R/genetic_layer.R"); source("R/genetic_helpers.R")

at <- c(0.05, 0.20)
base <- select_cohort(load_base_cohort(), include_twins = FALSE)
lp <- profile_logit(load_refit_params())
q <- tabulate(profile_index(base$X), 2L^length(PREDICTORS)) / nrow(base$X)
setup <- clec7a_second_setup(n_draw = 1000L)
cells <- genetic_cells(q, lp, setup$grid)

cl <- start_workers()
a <- attribution_layer(cl, cells, setup$effects, setup$beta_hat, at, n_attr = 1000L)
parallel::stopCluster(cl)
summ <- summarise_attribution(a, setup$terms)
summ$share_CLEC7A_combined <- summ$share_CLEC7A + summ$share_CLEC7A_AG
summ$phi_CLEC7A_combined <- summ$phi_CLEC7A + summ$phi_CLEC7A_AG

old <- read.csv("output/genetic_attribution_summary.csv")
key <- c("truth", "score", "threshold")
summ <- merge(summ, old[c(key, "total", paste0("share_", c("PTX3", "CLEC7A", "TLR4", "CD209")))], by = key, suffixes = c("", "_four_terms"))
write.csv(summ, "output/genetic_attribution_second_clec7a_summary.csv", row.names = FALSE)
message("Shapley shares with the second CLEC7A parameter (four-term shares in brackets via the _four_terms columns of the file):")
show_df(summ[, c(key, "total", paste0("share_", setup$terms), "share_CLEC7A_combined")])
