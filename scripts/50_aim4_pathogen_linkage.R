# Linkage among the causal loci in the pathogen-specificity scenario (addendum 43). Run from the project root, after scripts 11, 31 and 43, with
# Sys.setenv(PROSE_SCHEME = "profile") (the only truth this script builds).
#
# As script 47 (f = 1 and 0.8; 12 replicate discovery studies of 2,000 simulated genotypes; the mixture truth of addendum 41), with the causal loci in blocks of 4
# adjacent loci whose alleles are correlated rho = 0, 0.3 and 0.6 on each chromosome (R/aim4_pathogen.R, pathogen_cell_linked). The true dosage distribution is exact.
# The weights are estimated two ways: "marginal" (each locus on its own, as in a single-SNP discovery study: the estimate carries the effect of the linked neighbours
# and the errors are correlated within a block) and "joint" (each block estimated together). sigma = 0.10 and 0.20 is the standard error of one locus's marginal estimate.
# Grid: odds ratio 2.0, 2.5, 3.0 x frequency 0.1, 0.3, 0.5 x 1 to 20 loci. rho = 0 is the independent loci of script 48 ("base").
#
# Writes   output/aim4_pathogen_linkage_sweep_profile.csv, aim4_pathogen_linkage_frontier_profile.csv, output/figures/aim4_pathogen_linkage_frontier_profile.png

Sys.setenv(PROSE_SCHEME = "profile")
source("R/baseline.R"); source("R/bs_data.R"); source("R/decision_curve.R"); source("R/boscore_data.R"); source("R/refit_baseline.R")
source("R/genetic_layer.R"); source("R/genetic_helpers.R"); source("R/aim4_sweep.R"); source("R/calibration_check.R"); source("R/aim4_pathogen.R")
suppressPackageStartupMessages(library(ggplot2))

at <- c(0.05, 0.20)
base <- select_cohort(load_base_cohort(), include_twins = FALSE)
q <- tabulate(profile_index(base$X), 2L^length(PREDICTORS)) / nrow(base$X)
lp <- profile_logit(load_refit_params())
mcid <- setNames(baseline_mcid(q, lp, at), as.character(at))

grid <- expand.grid(n_loci = 1:20, odds_ratio = c(2.0, 2.5, 3.0), raf = c(0.1, 0.3, 0.5), case_fraction = c(1, 0.8), rho = c(0, 0.3, 0.6), sigma = c(0.10, 0.20))
grid$seed <- 16000L + seq_len(nrow(grid))
message("grid: ", nrow(grid), " cells, 12 replicates each")

cl <- start_workers(files = c("R/baseline.R", "R/decision_curve.R", "R/genetic_layer.R", "R/aim4_sweep.R", "R/calibration_check.R", "R/aim4_pathogen.R"))
rep_gains <- parallel_rows(cl, split(seq_len(nrow(grid)), rep_len(seq_len(90), nrow(grid))), function(ch, ctx)
  do.call(rbind, lapply(ch, function(i) pathogen_cell_linked(ctx$q, ctx$lp, ctx$g$n_loci[i], ctx$g$odds_ratio[i], ctx$g$raf[i], ctx$g$case_fraction[i], ctx$g$sigma[i], ctx$at,
    block_size = 4L, rho = ctx$g$rho[i], seed = ctx$g$seed[i]))), list(g = grid, q = q, lp = lp, at = at))
parallel::stopCluster(cl)

keys <- c("rho", "n_loci", "odds_ratio", "raf", "case_fraction", "sigma", "estimation", "score", "threshold")
rep_gains$mcid <- unname(mcid[as.character(rep_gains$threshold)])
summ <- by_key(rep_gains, keys, function(d) data.frame(gain_median = median(d$gain), gain_lo = quantile(d$gain, 0.05, names = FALSE), gain_hi = quantile(d$gain, 0.95, names = FALSE),
                                                        share_reaching_mcid = mean(d$gain >= d$mcid), gain_exact = d$gain_exact[1], mcid = d$mcid[1]))
summ$fraction_of_exact <- summ$gain_median / summ$gain_exact
frontier <- do.call(rbind, lapply(split(summ, list(summ$rho, summ$odds_ratio, summ$raf, summ$case_fraction, summ$sigma, summ$estimation, summ$score, summ$threshold), drop = TRUE), function(x) { x <- x[order(x$n_loci), ]
  data.frame(rho = x$rho[1], odds_ratio = x$odds_ratio[1], raf = x$raf[1], case_fraction = x$case_fraction[1], sigma = x$sigma[1], estimation = x$estimation[1], score = x$score[1], threshold = x$threshold[1],
             exact_reaches = first_reaching(x$n_loci, x$gain_exact, x$mcid[1]), median_reaches = first_reaching(x$n_loci, x$gain_median, x$mcid[1])) }))
rownames(frontier) <- NULL
write.csv(summ, "output/aim4_pathogen_linkage_sweep_profile.csv", row.names = FALSE)
write.csv(frontier, "output/aim4_pathogen_linkage_frontier_profile.csv", row.names = FALSE)

pl <- frontier[frontier$threshold == 0.05 & frontier$raf == 0.3 & frontier$sigma == 0.10 & frontier$score == "specific", ]
pl$loci <- ifelse(is.na(pl$median_reaches), 21, pl$median_reaches)
pl$panel <- paste0("f = ", pl$case_fraction, ", frequency 0.3, SE 0.10, 5% threshold")
dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)
ggsave("output/figures/aim4_pathogen_linkage_frontier_profile.png", ggplot(pl, aes(odds_ratio, loci, colour = factor(rho), linetype = estimation)) + geom_line() + geom_point() + facet_wrap(~ panel) +
         scale_y_continuous(breaks = c(1, 5, 10, 15, 20, 21), labels = c(1, 5, 10, 15, 20, ">20")) +
         labs(x = "Odds ratio per risk allele (target-pathogen pathway)", y = "Fewest loci that reach the MCID (median replicate)", colour = "Allele correlation in a block of 4", linetype = "Weights estimated"), width = 9, height = 4.5, dpi = 150)
