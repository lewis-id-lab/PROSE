# The Aim 4 sweep with a realised score (addendum 17). Run from the project root, after scripts 11, 16 and 23.
#
# As script 23 (1 to 20 loci, identical independent loci, base cohort, refit baseline, MCID +0.00458 at 5% and +0.00160 at 20%), but the weights of the score are estimates with
# an independent error of standard deviation sigma on each locus's log odds ratio: sigma = 0.05, 0.10, 0.20 and 0.40 (a discovery sample of about 4,000, 1,000, 240 or 60 subjects
# for a balanced design at frequency 0.3). The grid is odds ratio 1.2 to 3.0 (10 values) x frequency 0.1, 0.3, 0.5 x 1 to 20 loci. 12 replicate discovery studies per architecture
# (2,000 simulated genotypes each). The ceiling score of script 23 is the sigma = 0 case, and a check of the simulation against it is part of the output.
#
# Writes   output/aim4_realised_sweep.csv (per architecture, precision and threshold: median gain, 5% and 95% of replicates, share of replicates reaching the MCID),
#          output/aim4_realised_frontier.csv, output/aim4_realised_check.csv, output/figures/aim4_realised_frontier.png

source("R/baseline.R"); source("R/bs_data.R"); source("R/decision_curve.R"); source("R/boscore_data.R"); source("R/refit_baseline.R")
source("R/genetic_layer.R"); source("R/genetic_helpers.R"); source("R/aim4_sweep.R")
suppressPackageStartupMessages(library(ggplot2))

at <- c(0.05, 0.20)
base <- select_cohort(load_base_cohort(), include_twins = FALSE)
q <- tabulate(profile_index(base$X), 2L^length(PREDICTORS)) / nrow(base$X)
lp <- profile_logit(load_refit_params())
mcid <- setNames(baseline_mcid(q, lp, at), as.character(at))
exact <- read.csv("output/aim4_sweep.csv")

grid <- expand.grid(n_loci = 1:20, odds_ratio = c(1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 2.0, 2.5, 3.0), raf = c(0.1, 0.3, 0.5), sigma = c(0.05, 0.10, 0.20, 0.40))
grid$seed <- 7000L + seq_len(nrow(grid))
check_grid <- expand.grid(n_loci = c(2, 5, 10, 20), odds_ratio = c(1.5, 2.0, 3.0), raf = c(0.1, 0.3, 0.5), sigma = 0)
check_grid$seed <- 9000L + seq_len(nrow(check_grid))
message("grid: ", nrow(grid), " architectures x precisions, 12 replicates each; check: ", nrow(check_grid), " cells at sigma = 0")

cl <- start_workers(files = c("R/baseline.R", "R/decision_curve.R", "R/genetic_layer.R", "R/aim4_sweep.R"))
run <- function(g) parallel_rows(cl, split(seq_len(nrow(g)), rep_len(seq_len(60), nrow(g))), function(ch, ctx)
  do.call(rbind, lapply(ch, function(i) sweep_cell_realised(ctx$q, ctx$lp, ctx$g$n_loci[i], ctx$g$odds_ratio[i], ctx$g$raf[i], ctx$g$sigma[i], ctx$at, seed = ctx$g$seed[i]))),
  list(g = g, q = q, lp = lp, at = at))
rep_gains <- run(grid)
chk <- run(check_grid)
parallel::stopCluster(cl)

rep_gains$mcid <- unname(mcid[as.character(rep_gains$threshold)])
summ <- by_key(rep_gains, c("n_loci", "odds_ratio", "raf", "sigma", "threshold"), function(d) data.frame(
  gain_median = median(d$gain), gain_lo = quantile(d$gain, 0.05, names = FALSE), gain_hi = quantile(d$gain, 0.95, names = FALSE),
  share_reaching_mcid = mean(d$gain >= d$mcid), mcid = d$mcid[1]))
j <- match(paste(summ$n_loci, summ$odds_ratio, summ$raf, summ$threshold), paste(exact$n_loci, exact$odds_ratio, exact$raf, exact$threshold))
summ$gain_ceiling <- exact$gain[j]
summ$fraction_of_ceiling <- summ$gain_median / summ$gain_ceiling

chk$mcid <- unname(mcid[as.character(chk$threshold)])
check <- by_key(chk, c("n_loci", "odds_ratio", "raf", "threshold"), function(d) data.frame(mc_gain_mean = mean(d$gain), mc_gain_sd = sd(d$gain)))
jc <- match(paste(check$n_loci, check$odds_ratio, check$raf, check$threshold), paste(exact$n_loci, exact$odds_ratio, exact$raf, exact$threshold))
check$exact_gain <- exact$gain[jc]; check$difference <- check$mc_gain_mean - check$exact_gain

# Frontiers: smallest k (median gain, and at least 80% of replicates) and smallest OR
fr <- do.call(rbind, lapply(split(summ, list(summ$odds_ratio, summ$raf, summ$sigma, summ$threshold), drop = TRUE), function(x) { x <- x[order(x$n_loci), ]
  data.frame(kind = "min_loci", odds_ratio = x$odds_ratio[1], raf = x$raf[1], sigma = x$sigma[1], n_loci = NA_real_, threshold = x$threshold[1],
             median_reaches = first_reaching(x$n_loci, x$gain_median, x$mcid[1]), reliable_reaches = { i <- which(x$share_reaching_mcid >= 0.8); if (length(i)) x$n_loci[min(i)] else NA_real_ }) }))
fo <- do.call(rbind, lapply(split(summ, list(summ$n_loci, summ$raf, summ$sigma, summ$threshold), drop = TRUE), function(x) { x <- x[order(x$odds_ratio), ]
  data.frame(kind = "min_odds_ratio", odds_ratio = NA_real_, raf = x$raf[1], sigma = x$sigma[1], n_loci = x$n_loci[1], threshold = x$threshold[1],
             median_reaches = first_reaching(x$odds_ratio, x$gain_median, x$mcid[1]), reliable_reaches = { i <- which(x$share_reaching_mcid >= 0.8); if (length(i)) x$odds_ratio[min(i)] else NA_real_ }) }))
frontier <- rbind(fr, fo); rownames(frontier) <- NULL
# The oracle (sigma = 0) frontier from the exact sweep, on the same grid
ex <- exact[exact$odds_ratio %in% unique(grid$odds_ratio) & exact$raf %in% unique(grid$raf), ]
ex$sigma <- 0
exfr <- rbind(do.call(rbind, lapply(split(ex, list(ex$odds_ratio, ex$raf, ex$threshold), drop = TRUE), function(x) { x <- x[order(x$n_loci), ]
  data.frame(kind = "min_loci", odds_ratio = x$odds_ratio[1], raf = x$raf[1], sigma = 0, n_loci = NA_real_, threshold = x$threshold[1], median_reaches = first_reaching(x$n_loci, x$gain, x$mcid[1]), reliable_reaches = first_reaching(x$n_loci, x$gain, x$mcid[1])) })),
  do.call(rbind, lapply(split(ex, list(ex$n_loci, ex$raf, ex$threshold), drop = TRUE), function(x) { x <- x[order(x$odds_ratio), ]
  data.frame(kind = "min_odds_ratio", odds_ratio = NA_real_, raf = x$raf[1], sigma = 0, n_loci = x$n_loci[1], threshold = x$threshold[1], median_reaches = first_reaching(x$odds_ratio, x$gain, x$mcid[1]), reliable_reaches = first_reaching(x$odds_ratio, x$gain, x$mcid[1])) })))
frontier <- rbind(exfr, frontier); rownames(frontier) <- NULL

write.csv(summ, "output/aim4_realised_sweep.csv", row.names = FALSE)
write.csv(frontier, "output/aim4_realised_frontier.csv", row.names = FALSE)
write.csv(check, "output/aim4_realised_check.csv", row.names = FALSE)

pl <- frontier[frontier$kind == "min_loci" & frontier$odds_ratio >= 1.4, ]
pl$loci <- ifelse(is.na(pl$median_reaches), 21, pl$median_reaches)
pl$precision <- factor(ifelse(pl$sigma == 0, "exact weights", paste0("SE ", pl$sigma)), levels = c("exact weights", paste0("SE ", c(0.05, 0.10, 0.20, 0.40))))
pl$panel <- paste0("threshold ", round(pl$threshold * 100), "%, frequency ", pl$raf)
ggsave("output/figures/aim4_realised_frontier.png", ggplot(pl, aes(odds_ratio, loci, colour = precision)) + geom_line() + geom_point() +
         facet_wrap(~ panel, nrow = 2) + scale_y_continuous(breaks = c(1, 5, 10, 15, 20, 21), labels = c(1, 5, 10, 15, 20, ">20")) +
         labs(x = "Odds ratio per risk allele", y = "Fewest loci that reach the MCID (median replicate)", colour = NULL),
       width = 12, height = 6, dpi = 150)

message("\nSimulation against the exact ceiling at sigma = 0 (difference in gain):"); show_df(check[, c("n_loci", "odds_ratio", "raf", "threshold", "mc_gain_mean", "exact_gain", "difference")])
message("\nFewest loci reaching the 5% MCID (median replicate; NA = none within 20), by precision, frequency 0.3:")
w <- frontier[frontier$kind == "min_loci" & frontier$threshold == 0.05 & frontier$raf == 0.3, ]
print(reshape(w[c("odds_ratio", "sigma", "median_reaches")], idvar = "odds_ratio", timevar = "sigma", direction = "wide"), row.names = FALSE)
