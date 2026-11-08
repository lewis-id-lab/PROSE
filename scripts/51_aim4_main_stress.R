# The stress tests in the main architecture sweep (addendum 44). Run from the project root, after scripts 11, 31 and 43, with PROSE_SCHEME = "profile" (set below).
#
# The main sweep (scripts 23, 25 and 31 to 34) treats every case as genotype-sensitive (f = 1) and, in its realised version, gives the score weights with an independent
# error and identical independent loci scored at the causal variant. Here the realised score of the same grid (odds ratio 1.2 to 3.0, 10 values, x frequency 0.1, 0.3, 0.5 x 1 to
# 20 loci; 600 architectures) meets the departures of addenda 42 and 43, at f = 1 and sigma = 0.10 and 0.40 (12 replicate discovery studies of 2,000 simulated genotypes):
#   base      identical independent loci, as in the main sweep (the check)
#   locus     unequal frequencies and precision
#   curse     the winner's curse (k loci picked from 5 k candidates)
#   tag0.9    the scored variant a tag correlated 0.9 with the causal one
#   tag0.7    correlated 0.7
#   all       locus, curse and tag at 0.7 together
#   linked    blocks of 4 causal loci with allele correlation 0.3 and 0.6, weights estimated per locus (marginal) or per block (joint)
# Profile-preserving truth. The weights are the true ones (the target-pathogen and the all-mould centres coincide at f = 1 up to 0.1 in the log odds ratio; the target-pathogen score is used).
#
# Writes   output/aim4_main_stress_sweep_profile.csv, aim4_main_stress_frontier_profile.csv, output/figures/aim4_main_stress_frontier_profile.png

Sys.setenv(PROSE_SCHEME = "profile")
source("R/baseline.R"); source("R/bs_data.R"); source("R/decision_curve.R"); source("R/boscore_data.R"); source("R/refit_baseline.R")
source("R/genetic_layer.R"); source("R/genetic_helpers.R"); source("R/aim4_sweep.R"); source("R/calibration_check.R"); source("R/aim4_pathogen.R")
suppressPackageStartupMessages(library(ggplot2))

at <- c(0.05, 0.20)
base <- select_cohort(load_base_cohort(), include_twins = FALSE)
q <- tabulate(profile_index(base$X), 2L^length(PREDICTORS)) / nrow(base$X)
lp <- profile_logit(load_refit_params())
mcid <- setNames(baseline_mcid(q, lp, at), as.character(at))
exact <- read.csv("output/aim4_pathogen_sweep_profile.csv"); exact <- exact[exact$case_fraction == 1, ]

scen <- data.frame(scenario = c("base", "locus", "curse", "tag0.9", "tag0.7", "all", "linked0.3", "linked0.6"), kind = c(rep("stress", 6), "linked", "linked"),
                   locus_specific = c(FALSE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE), winners_curse = c(FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE),
                   tag_r = c(1, 1, 1, 0.9, 0.7, 0.7, 1, 1), rho = c(0, 0, 0, 0, 0, 0, 0.3, 0.6), stringsAsFactors = FALSE)
arch <- expand.grid(n_loci = 1:20, odds_ratio = c(1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 2.0, 2.5, 3.0), raf = c(0.1, 0.3, 0.5), sigma = c(0.10, 0.40))
grid <- merge(arch, scen, by = NULL)
grid$seed <- 18000L + seq_len(nrow(grid))
message("grid: ", nrow(grid), " cells, 12 replicates each")

cl <- start_workers(files = c("R/baseline.R", "R/decision_curve.R", "R/genetic_layer.R", "R/aim4_sweep.R", "R/calibration_check.R", "R/aim4_pathogen.R"))
rep_gains <- parallel_rows(cl, split(seq_len(nrow(grid)), rep_len(seq_len(120), nrow(grid))), function(ch, ctx)
  do.call(rbind, lapply(ch, function(i) { g <- ctx$g[i, ]
    r <- if (g$kind == "linked") pathogen_cell_linked(ctx$q, ctx$lp, g$n_loci, g$odds_ratio, g$raf, 1, g$sigma, ctx$at, block_size = 4L, rho = g$rho, seed = g$seed)
         else { x <- pathogen_cell_stress(ctx$q, ctx$lp, g$n_loci, g$odds_ratio, g$raf, 1, g$sigma, ctx$at, locus_specific = g$locus_specific, winners_curse = g$winners_curse, tag_r = g$tag_r, seed = g$seed); x$estimation <- "per_locus"; x }
    r <- r[r$score == "specific", c("n_loci", "odds_ratio", "raf", "sigma", "estimation", "replicate", "threshold", "gain")]
    cbind(scenario = g$scenario, r) })), list(g = grid, q = q, lp = lp, at = at))
parallel::stopCluster(cl)

rep_gains$mcid <- unname(mcid[as.character(rep_gains$threshold)])
rep_gains$scenario <- ifelse(rep_gains$scenario %in% c("linked0.3", "linked0.6"), paste0(rep_gains$scenario, "_", rep_gains$estimation), rep_gains$scenario)
keys <- c("scenario", "n_loci", "odds_ratio", "raf", "sigma", "threshold")
summ <- by_key(rep_gains, keys, function(d) data.frame(gain_median = median(d$gain), gain_lo = quantile(d$gain, 0.05, names = FALSE), gain_hi = quantile(d$gain, 0.95, names = FALSE),
                                                        share_reaching_mcid = mean(d$gain >= d$mcid), mcid = d$mcid[1]))
key5 <- function(x) paste(x$n_loci, x$odds_ratio, x$raf, x$threshold)
summ$gain_exact <- exact$gain_specific[match(key5(summ), key5(exact))]
summ$fraction_of_exact <- summ$gain_median / summ$gain_exact
frontier <- do.call(rbind, lapply(split(summ, list(summ$scenario, summ$odds_ratio, summ$raf, summ$sigma, summ$threshold), drop = TRUE), function(x) { x <- x[order(x$n_loci), ]
  data.frame(scenario = x$scenario[1], odds_ratio = x$odds_ratio[1], raf = x$raf[1], sigma = x$sigma[1], threshold = x$threshold[1],
             exact_reaches = first_reaching(x$n_loci, x$gain_exact, x$mcid[1]), median_reaches = first_reaching(x$n_loci, x$gain_median, x$mcid[1])) }))
rownames(frontier) <- NULL
write.csv(summ, "output/aim4_main_stress_sweep_profile.csv", row.names = FALSE)
write.csv(frontier, "output/aim4_main_stress_frontier_profile.csv", row.names = FALSE)

pl <- frontier[frontier$threshold == 0.05 & frontier$raf == 0.3 & frontier$odds_ratio >= 1.4 & frontier$scenario %in% c("base", "tag0.9", "tag0.7", "all", "linked0.3_marginal", "linked0.3_joint"), ]
pl$loci <- ifelse(is.na(pl$median_reaches), 21, pl$median_reaches)
pl$panel <- paste0("SE ", pl$sigma, ", frequency 0.3, 5% threshold")
ex <- unique(data.frame(odds_ratio = pl$odds_ratio, loci = ifelse(is.na(pl$exact_reaches), 21, pl$exact_reaches)))
dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)
ggsave("output/figures/aim4_main_stress_frontier_profile.png", ggplot(pl, aes(odds_ratio, loci, colour = scenario)) + geom_line() + geom_point() + facet_wrap(~ panel) +
         scale_y_continuous(breaks = c(1, 5, 10, 15, 20, 21), labels = c(1, 5, 10, 15, 20, ">20")) +
         labs(x = "Odds ratio per risk allele", y = "Fewest loci that reach the MCID (median replicate)", colour = NULL), width = 10, height = 4.5, dpi = 150)
