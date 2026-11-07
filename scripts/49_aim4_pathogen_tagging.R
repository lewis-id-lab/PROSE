# Where does incomplete tagging start to cost? The tagging stress test of script 48 over r = 0.5, 0.7, 0.9 and 1.0 (addendum 42). Run from the project root,
# after scripts 11, 31, 43 and 48.
#
# As script 48's "tag" scenario alone (identical loci, no winner's curse; the scored variant is a tag correlated r with the causal one, and its estimated weight is
# centred on r times the causal effect), at f = 1 and 0.8 and sigma = 0.10, on odds ratios 2.0, 2.5, 3.0 x frequencies 0.1, 0.3, 0.5 x 1 to 20 loci (12 replicate
# discovery studies of 2,000 simulated genotypes), under the profile-preserving truth. r = 1 is the causal variant itself.
#
# Writes   output/aim4_pathogen_tagging_sweep_profile.csv, aim4_pathogen_tagging_frontier_profile.csv, output/figures/aim4_pathogen_tagging_frontier_profile.png

Sys.setenv(PROSE_SCHEME = "profile")
source("R/baseline.R"); source("R/bs_data.R"); source("R/decision_curve.R"); source("R/boscore_data.R"); source("R/refit_baseline.R")
source("R/genetic_layer.R"); source("R/genetic_helpers.R"); source("R/aim4_sweep.R"); source("R/calibration_check.R"); source("R/aim4_pathogen.R")
suppressPackageStartupMessages(library(ggplot2))

at <- c(0.05, 0.20)
base <- select_cohort(load_base_cohort(), include_twins = FALSE)
q <- tabulate(profile_index(base$X), 2L^length(PREDICTORS)) / nrow(base$X)
lp <- profile_logit(load_refit_params())
mcid <- setNames(baseline_mcid(q, lp, at), as.character(at))

grid <- expand.grid(n_loci = 1:20, odds_ratio = c(2.0, 2.5, 3.0), raf = c(0.1, 0.3, 0.5), case_fraction = c(1, 0.8), tag_r = c(0.5, 0.7, 0.9, 1.0), sigma = 0.10)
grid$seed <- 14000L + seq_len(nrow(grid))
message("grid: ", nrow(grid), " cells, 12 replicates each")

cl <- start_workers(files = c("R/baseline.R", "R/decision_curve.R", "R/genetic_layer.R", "R/aim4_sweep.R", "R/calibration_check.R", "R/aim4_pathogen.R"))
rep_gains <- parallel_rows(cl, split(seq_len(nrow(grid)), rep_len(seq_len(60), nrow(grid))), function(ch, ctx)
  do.call(rbind, lapply(ch, function(i) cbind(tag_r = ctx$g$tag_r[i], pathogen_cell_stress(ctx$q, ctx$lp, ctx$g$n_loci[i], ctx$g$odds_ratio[i], ctx$g$raf[i], ctx$g$case_fraction[i], ctx$g$sigma[i], ctx$at,
    tag_r = ctx$g$tag_r[i], seed = ctx$g$seed[i])))), list(g = grid, q = q, lp = lp, at = at))
parallel::stopCluster(cl)

keys <- c("tag_r", "n_loci", "odds_ratio", "raf", "case_fraction", "score", "threshold")
rep_gains$mcid <- unname(mcid[as.character(rep_gains$threshold)])
summ <- by_key(rep_gains, keys, function(d) data.frame(gain_median = median(d$gain), gain_lo = quantile(d$gain, 0.05, names = FALSE), gain_hi = quantile(d$gain, 0.95, names = FALSE),
                                                        share_reaching_mcid = mean(d$gain >= d$mcid), mcid = d$mcid[1]))
r1 <- summ[summ$tag_r == 1, ]
jb <- match(paste(summ$n_loci, summ$odds_ratio, summ$raf, summ$case_fraction, summ$score, summ$threshold), paste(r1$n_loci, r1$odds_ratio, r1$raf, r1$case_fraction, r1$score, r1$threshold))
summ$gain_causal <- r1$gain_median[jb]
summ$fraction_of_causal <- summ$gain_median / summ$gain_causal
frontier <- do.call(rbind, lapply(split(summ, list(summ$tag_r, summ$odds_ratio, summ$raf, summ$case_fraction, summ$score, summ$threshold), drop = TRUE), function(x) { x <- x[order(x$n_loci), ]
  data.frame(tag_r = x$tag_r[1], odds_ratio = x$odds_ratio[1], raf = x$raf[1], case_fraction = x$case_fraction[1], score = x$score[1], threshold = x$threshold[1],
             median_reaches = first_reaching(x$n_loci, x$gain_median, x$mcid[1])) }))
rownames(frontier) <- NULL
write.csv(summ, "output/aim4_pathogen_tagging_sweep_profile.csv", row.names = FALSE)
write.csv(frontier, "output/aim4_pathogen_tagging_frontier_profile.csv", row.names = FALSE)

pl <- frontier[frontier$threshold == 0.05 & frontier$raf == 0.3 & frontier$score == "specific", ]
pl$loci <- ifelse(is.na(pl$median_reaches), 21, pl$median_reaches)
pl$panel <- paste0("f = ", pl$case_fraction, ", frequency 0.3, SE 0.10, 5% threshold")
ggsave("output/figures/aim4_pathogen_tagging_frontier_profile.png", ggplot(pl, aes(odds_ratio, loci, colour = factor(tag_r))) + geom_line() + geom_point() + facet_wrap(~ panel) +
         scale_y_continuous(breaks = c(1, 5, 10, 15, 20, 21), labels = c(1, 5, 10, 15, 20, ">20")) +
         labs(x = "Odds ratio per risk allele (target-pathogen pathway)", y = "Fewest loci that reach the MCID (median replicate)", colour = "Tag-causal correlation r"), width = 9, height = 4.5, dpi = 150)
