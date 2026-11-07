# The pathogen-specificity scenario with estimation error in the weights (addendum 39). Run from the project root, after scripts 11, 23 and 43.
#
# As script 43 (only a fraction f of cases arises through the pathway the loci act on; MCID +0.00458 at 5% and +0.00160 at 20%) and as script 25 (the
# weights of the score are estimates with an independent error of sd sigma on each locus's log odds ratio, 12 replicate discovery studies of 2,000
# simulated genotypes each), with the two centres of the weights of script 43: the target-pathogen effect and the all-mould fit. Grid: odds ratio 1.5,
# 1.8, 2.0, 2.5, 3.0 x frequency 0.1, 0.3, 0.5 x 1 to 20 loci x f = 0.8 and 0.6 x sigma = 0.10, 0.20, 0.40 (f = 1 is script 25; at f = 0.4 nothing
# reaches the MCID even with exact weights). A check at sigma = 0 compares the simulation with the exact gains of script 43.
#
# Truth: PROSE_SCHEME = "profile" gives the primary (profile-preserving) truth, D27; the default is the global-intercept truth (addendum 39); files then end in
# _profile, and the sigma = 0 check is against the same-scheme output of script 43.
#
# Writes   output/aim4_pathogen_realised_sweep.csv (per architecture, score and threshold: median gain, 5% and 95% of replicates, share reaching the MCID),
#          output/aim4_pathogen_realised_frontier.csv, output/aim4_pathogen_realised_check.csv, output/figures/aim4_pathogen_realised_frontier.png

source("R/baseline.R"); source("R/bs_data.R"); source("R/decision_curve.R"); source("R/boscore_data.R"); source("R/refit_baseline.R")
source("R/genetic_layer.R"); source("R/genetic_helpers.R"); source("R/aim4_sweep.R"); source("R/calibration_check.R"); source("R/aim4_pathogen.R")
suppressPackageStartupMessages(library(ggplot2))

sfx <- scheme_suffix()
at <- c(0.05, 0.20)
base <- select_cohort(load_base_cohort(), include_twins = FALSE)
q <- tabulate(profile_index(base$X), 2L^length(PREDICTORS)) / nrow(base$X)
lp <- profile_logit(load_refit_params())
mcid <- setNames(baseline_mcid(q, lp, at), as.character(at))
exact <- read.csv(paste0("output/aim4_pathogen_sweep", sfx, ".csv"))

grid <- expand.grid(n_loci = 1:20, odds_ratio = c(1.5, 1.8, 2.0, 2.5, 3.0), raf = c(0.1, 0.3, 0.5), case_fraction = c(0.8, 0.6), sigma = c(0.10, 0.20, 0.40))
grid$seed <- 8000L + seq_len(nrow(grid))
check_grid <- expand.grid(n_loci = c(2, 5, 10, 20), odds_ratio = c(1.5, 2.0, 3.0), raf = c(0.1, 0.3, 0.5), case_fraction = c(0.8, 0.6), sigma = 0)
check_grid$seed <- 9500L + seq_len(nrow(check_grid))
message("grid: ", nrow(grid), " architectures x precisions, 12 replicates each; check: ", nrow(check_grid), " cells at sigma = 0")

cl <- start_workers(files = c("R/baseline.R", "R/decision_curve.R", "R/genetic_layer.R", "R/aim4_sweep.R", "R/calibration_check.R", "R/aim4_pathogen.R"))
run <- function(g) parallel_rows(cl, split(seq_len(nrow(g)), rep_len(seq_len(60), nrow(g))), function(ch, ctx)
  do.call(rbind, lapply(ch, function(i) pathogen_cell_realised(ctx$q, ctx$lp, ctx$g$n_loci[i], ctx$g$odds_ratio[i], ctx$g$raf[i], ctx$g$case_fraction[i], ctx$g$sigma[i], ctx$at, seed = ctx$g$seed[i]))),
  list(g = g, q = q, lp = lp, at = at))
rep_gains <- run(grid)
chk <- run(check_grid)
parallel::stopCluster(cl)

keys <- c("n_loci", "odds_ratio", "raf", "case_fraction", "sigma", "score", "threshold")
rep_gains$mcid <- unname(mcid[as.character(rep_gains$threshold)])
summ <- by_key(rep_gains, keys, function(d) data.frame(
  gain_median = median(d$gain), gain_lo = quantile(d$gain, 0.05, names = FALSE), gain_hi = quantile(d$gain, 0.95, names = FALSE),
  share_reaching_mcid = mean(d$gain >= d$mcid), mcid = d$mcid[1]))
ex_col <- c(specific = "gain_specific", all_mould = "gain_diluted")
j <- match(paste(summ$n_loci, summ$odds_ratio, summ$raf, summ$case_fraction, summ$threshold), paste(exact$n_loci, exact$odds_ratio, exact$raf, exact$case_fraction, exact$threshold))
summ$gain_exact <- vapply(seq_len(nrow(summ)), function(i) exact[[ex_col[[summ$score[i]]]]][j[i]], numeric(1))
summ$fraction_of_exact <- summ$gain_median / summ$gain_exact

check <- by_key(chk, keys, function(d) data.frame(mc_gain_mean = mean(d$gain), mc_gain_sd = sd(d$gain)))
jc <- match(paste(check$n_loci, check$odds_ratio, check$raf, check$case_fraction, check$threshold), paste(exact$n_loci, exact$odds_ratio, exact$raf, exact$case_fraction, exact$threshold))
check$exact_gain <- vapply(seq_len(nrow(check)), function(i) exact[[ex_col[[check$score[i]]]]][jc[i]], numeric(1))
check$difference <- check$mc_gain_mean - check$exact_gain

# Frontier: the fewest loci whose median replicate reaches the MCID, and the fewest at which at least 80% of replicates do
frontier <- do.call(rbind, lapply(split(summ, list(summ$odds_ratio, summ$raf, summ$case_fraction, summ$sigma, summ$score, summ$threshold), drop = TRUE), function(x) { x <- x[order(x$n_loci), ]
  data.frame(odds_ratio = x$odds_ratio[1], raf = x$raf[1], case_fraction = x$case_fraction[1], sigma = x$sigma[1], score = x$score[1], threshold = x$threshold[1],
             median_reaches = first_reaching(x$n_loci, x$gain_median, x$mcid[1]), reliable_reaches = { i <- which(x$share_reaching_mcid >= 0.8); if (length(i)) x$n_loci[min(i)] else NA_real_ }) }))
rownames(frontier) <- NULL

write.csv(summ, paste0("output/aim4_pathogen_realised_sweep", sfx, ".csv"), row.names = FALSE)
write.csv(frontier, paste0("output/aim4_pathogen_realised_frontier", sfx, ".csv"), row.names = FALSE)
write.csv(check, paste0("output/aim4_pathogen_realised_check", sfx, ".csv"), row.names = FALSE)

pl <- frontier[frontier$threshold == 0.20 | frontier$threshold == 0.05, ]
pl$loci <- ifelse(is.na(pl$median_reaches), 21, pl$median_reaches)
pl$panel <- paste0("threshold ", round(pl$threshold * 100), "%, f = ", pl$case_fraction, ", frequency 0.3")
pl <- pl[pl$raf == 0.3, ]
pl$weights <- paste0(ifelse(pl$score == "specific", "specific", "all-mould"), ", SE ", pl$sigma)
dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)
ggsave(paste0("output/figures/aim4_pathogen_realised_frontier", sfx, ".png"), ggplot(pl, aes(odds_ratio, loci, colour = factor(sigma), linetype = score)) + geom_line() + geom_point() +
         facet_wrap(~ panel, nrow = 2) + scale_y_continuous(breaks = c(1, 5, 10, 15, 20, 21), labels = c(1, 5, 10, 15, 20, ">20")) +
         labs(x = "Odds ratio per risk allele (target-pathogen pathway)", y = "Fewest loci that reach the MCID (median replicate)", colour = "SE of the weights", linetype = "Weights centred on"),
       width = 12, height = 6, dpi = 150)

message("\nSimulation against the exact gain at sigma = 0 (difference in gain):"); show_df(check[, c("n_loci", "odds_ratio", "raf", "case_fraction", "score", "threshold", "mc_gain_mean", "exact_gain", "difference")])
message("\nFewest loci reaching the 5% MCID (median replicate; NA = none within 20), frequency 0.3, f = 0.8:")
w <- frontier[frontier$threshold == 0.05 & frontier$raf == 0.3 & frontier$case_fraction == 0.8, ]
print(reshape(w[c("odds_ratio", "sigma", "score", "median_reaches")], idvar = c("odds_ratio", "score"), timevar = "sigma", direction = "wide"), row.names = FALSE)
