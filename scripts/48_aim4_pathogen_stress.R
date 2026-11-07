# Stress tests of the pathogen-specificity scenario with estimated weights (addendum 42). Run from the project root, after scripts 11, 31 and 43, with
# Sys.setenv(PROSE_SCHEME = "profile") (the only truth this script builds).
#
# As script 47 (f = 1 and 0.8; sigma = 0.10 and 0.20; 12 replicate discovery studies of 2,000 simulated genotypes), with the three departures of the planned
# analysis that the architecture sweep leaves out, each alone and all together (R/aim4_pathogen.R, pathogen_cell_stress):
#   base      the identical independent loci of script 47
#   locus     loci differ in frequency and precision (the error of a weight is inversely proportional to sqrt(2 p (1 - p)))
#   curse     the winner's curse: the k loci are picked from 5 k candidates (the rest null) by their z-scores
#   tag       incomplete tagging: the scored genotype is a tag correlated r = 0.7 with the causal one
#   all       the three together
# Grid: odds ratio 2.0, 2.5, 3.0 x frequency 0.1, 0.3, 0.5 x 1 to 20 loci. A check at sigma = 0 compares the base scenario with script 43's exact gains.
# Linkage among the causal loci themselves is not modelled.
#
# Writes   output/aim4_pathogen_stress_sweep_profile.csv, aim4_pathogen_stress_frontier_profile.csv, aim4_pathogen_stress_check_profile.csv,
#          output/figures/aim4_pathogen_stress_frontier_profile.png

Sys.setenv(PROSE_SCHEME = "profile")
source("R/baseline.R"); source("R/bs_data.R"); source("R/decision_curve.R"); source("R/boscore_data.R"); source("R/refit_baseline.R")
source("R/genetic_layer.R"); source("R/genetic_helpers.R"); source("R/aim4_sweep.R"); source("R/calibration_check.R"); source("R/aim4_pathogen.R")
suppressPackageStartupMessages(library(ggplot2))

at <- c(0.05, 0.20)
base <- select_cohort(load_base_cohort(), include_twins = FALSE)
q <- tabulate(profile_index(base$X), 2L^length(PREDICTORS)) / nrow(base$X)
lp <- profile_logit(load_refit_params())
mcid <- setNames(baseline_mcid(q, lp, at), as.character(at))
exact <- read.csv("output/aim4_pathogen_sweep_profile.csv")

scen <- data.frame(scenario = c("base", "locus", "curse", "tag", "all"), locus_specific = c(FALSE, TRUE, FALSE, FALSE, TRUE),
                   winners_curse = c(FALSE, FALSE, TRUE, FALSE, TRUE), tag_r = c(1, 1, 1, 0.7, 0.7), stringsAsFactors = FALSE)
arch <- expand.grid(n_loci = 1:20, odds_ratio = c(2.0, 2.5, 3.0), raf = c(0.1, 0.3, 0.5), case_fraction = c(1, 0.8), sigma = c(0.10, 0.20))
grid <- merge(arch, scen, by = NULL)
grid$seed <- 11000L + seq_len(nrow(grid))
check_grid <- merge(expand.grid(n_loci = c(2, 5, 10, 20), odds_ratio = c(2.0, 2.5, 3.0), raf = c(0.1, 0.3, 0.5), case_fraction = c(1, 0.8), sigma = 0), scen[1, ], by = NULL)
check_grid$seed <- 12500L + seq_len(nrow(check_grid))
message("grid: ", nrow(grid), " architecture x scenario x precision cells, 12 replicates each; check: ", nrow(check_grid), " cells")

cl <- start_workers(files = c("R/baseline.R", "R/decision_curve.R", "R/genetic_layer.R", "R/aim4_sweep.R", "R/calibration_check.R", "R/aim4_pathogen.R"))
run <- function(g) parallel_rows(cl, split(seq_len(nrow(g)), rep_len(seq_len(90), nrow(g))), function(ch, ctx)
  do.call(rbind, lapply(ch, function(i) cbind(scenario = ctx$g$scenario[i], pathogen_cell_stress(ctx$q, ctx$lp, ctx$g$n_loci[i], ctx$g$odds_ratio[i], ctx$g$raf[i], ctx$g$case_fraction[i], ctx$g$sigma[i], ctx$at,
    locus_specific = ctx$g$locus_specific[i], winners_curse = ctx$g$winners_curse[i], tag_r = ctx$g$tag_r[i], seed = ctx$g$seed[i])))),
  list(g = g, q = q, lp = lp, at = at))
rep_gains <- run(grid)
chk <- run(check_grid)
parallel::stopCluster(cl)

keys <- c("scenario", "n_loci", "odds_ratio", "raf", "case_fraction", "sigma", "score", "threshold")
rep_gains$mcid <- unname(mcid[as.character(rep_gains$threshold)])
summ <- by_key(rep_gains, keys, function(d) data.frame(
  gain_median = median(d$gain), gain_lo = quantile(d$gain, 0.05, names = FALSE), gain_hi = quantile(d$gain, 0.95, names = FALSE),
  share_reaching_mcid = mean(d$gain >= d$mcid), true_selected = mean(d$true_selected), mcid = d$mcid[1]))
ex_col <- c(specific = "gain_specific", all_mould = "gain_diluted")
key6 <- function(x) paste(x$n_loci, x$odds_ratio, x$raf, x$case_fraction, x$threshold)
j <- match(key6(summ), key6(exact))
summ$gain_exact <- vapply(seq_len(nrow(summ)), function(i) exact[[ex_col[[summ$score[i]]]]][j[i]], numeric(1))
bs <- summ[summ$scenario == "base", ]
jb <- match(paste(summ$n_loci, summ$odds_ratio, summ$raf, summ$case_fraction, summ$sigma, summ$score, summ$threshold), paste(bs$n_loci, bs$odds_ratio, bs$raf, bs$case_fraction, bs$sigma, bs$score, bs$threshold))
summ$gain_base <- bs$gain_median[jb]
summ$fraction_of_base <- summ$gain_median / summ$gain_base

check <- by_key(chk, keys, function(d) data.frame(mc_gain_mean = mean(d$gain), mc_gain_sd = sd(d$gain)))
jc <- match(key6(check), key6(exact))
check$exact_gain <- vapply(seq_len(nrow(check)), function(i) exact[[ex_col[[check$score[i]]]]][jc[i]], numeric(1))
check$difference <- check$mc_gain_mean - check$exact_gain

frontier <- do.call(rbind, lapply(split(summ, list(summ$scenario, summ$odds_ratio, summ$raf, summ$case_fraction, summ$sigma, summ$score, summ$threshold), drop = TRUE), function(x) { x <- x[order(x$n_loci), ]
  data.frame(scenario = x$scenario[1], odds_ratio = x$odds_ratio[1], raf = x$raf[1], case_fraction = x$case_fraction[1], sigma = x$sigma[1], score = x$score[1], threshold = x$threshold[1],
             median_reaches = first_reaching(x$n_loci, x$gain_median, x$mcid[1]), reliable_reaches = { i <- which(x$share_reaching_mcid >= 0.8); if (length(i)) x$n_loci[min(i)] else NA_real_ }) }))
rownames(frontier) <- NULL

write.csv(summ, "output/aim4_pathogen_stress_sweep_profile.csv", row.names = FALSE)
write.csv(frontier, "output/aim4_pathogen_stress_frontier_profile.csv", row.names = FALSE)
write.csv(check, "output/aim4_pathogen_stress_check_profile.csv", row.names = FALSE)

pl <- frontier[frontier$threshold == 0.05 & frontier$raf == 0.3 & frontier$sigma == 0.10 & frontier$score == "specific", ]
pl$loci <- ifelse(is.na(pl$median_reaches), 21, pl$median_reaches)
pl$panel <- paste0("f = ", pl$case_fraction, ", frequency 0.3, SE 0.10, 5% threshold")
pl$scenario <- factor(pl$scenario, levels = scen$scenario)
dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)
ggsave("output/figures/aim4_pathogen_stress_frontier_profile.png", ggplot(pl, aes(odds_ratio, loci, colour = scenario)) + geom_line() + geom_point() + facet_wrap(~ panel) +
         scale_y_continuous(breaks = c(1, 5, 10, 15, 20, 21), labels = c(1, 5, 10, 15, 20, ">20")) +
         labs(x = "Odds ratio per risk allele (target-pathogen pathway)", y = "Fewest loci that reach the MCID (median replicate)", colour = NULL), width = 9, height = 4.5, dpi = 150)

message("\nSimulation against the exact gain at sigma = 0 (base scenario):"); show_df(check[, c("n_loci", "odds_ratio", "raf", "case_fraction", "score", "threshold", "mc_gain_mean", "exact_gain", "difference")][1:12, ])
message("\nShare of architectures whose median replicate reaches the 5% MCID:")
s5 <- summ[summ$threshold == 0.05, ]; s5$meets <- s5$gain_median >= s5$mcid
print(reshape(aggregate(meets ~ scenario + case_fraction + sigma + score, s5, mean), idvar = c("scenario", "case_fraction", "score"), timevar = "sigma", direction = "wide"), row.names = FALSE)
