# Aim 4: the pathogen-specificity scenario (scope decision 1). Run from the project root, after scripts 11 and 23.
#
# As script 23 (1 to 20 loci, identical independent loci, base cohort, refit baseline, MCID +0.00458 at 5% and +0.00160 at 20%), but only a fraction f of cases arises through the
# pathway the loci act on (for example Aspergillus among all invasive mould disease); the rest are genotype-independent. f = 1, 0.8, 0.6 and 0.4. The grid is odds ratio
# 1.2 to 3.0 (10 values) x frequency 0.1, 0.3, 0.5 x 1 to 20 loci. Two scores: "specific" (weights from the target-pathogen effect) and "diluted" (weights fitted to the
# all-mould outcome); see R/aim4_pathogen.R. f = 1 is checked against the exact sweep of script 23. The odds ratio is the effect on the target-pathogen pathway.
# f is an assumption: the cohort has no pathogen-level outcome, so it is not estimated here.
#
# Writes   output/aim4_pathogen_sweep.csv, output/aim4_pathogen_frontier.csv, output/aim4_pathogen_check.csv, output/figures/aim4_pathogen_frontier.png

source("R/baseline.R"); source("R/bs_data.R"); source("R/decision_curve.R"); source("R/boscore_data.R"); source("R/refit_baseline.R")
source("R/genetic_layer.R"); source("R/genetic_helpers.R"); source("R/aim4_sweep.R"); source("R/aim4_pathogen.R")
suppressPackageStartupMessages(library(ggplot2))

at <- c(0.05, 0.20)
base <- select_cohort(load_base_cohort(), include_twins = FALSE)
q <- tabulate(profile_index(base$X), 2L^length(PREDICTORS)) / nrow(base$X)
lp <- profile_logit(load_refit_params())
mcid <- setNames(baseline_mcid(q, lp, at), as.character(at))
exact <- read.csv("output/aim4_sweep.csv")

grid <- expand.grid(n_loci = 1:20, odds_ratio = c(1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 2.0, 2.5, 3.0), raf = c(0.1, 0.3, 0.5), case_fraction = c(1, 0.8, 0.6, 0.4))
message("grid: ", nrow(grid), " architectures")

cl <- start_workers(files = c("R/baseline.R", "R/decision_curve.R", "R/genetic_layer.R", "R/aim4_sweep.R", "R/aim4_pathogen.R"))
sw <- parallel_rows(cl, split(seq_len(nrow(grid)), rep_len(seq_len(60), nrow(grid))), function(ch, ctx)
  do.call(rbind, lapply(ch, function(i) pathogen_cell(ctx$q, ctx$lp, ctx$g$n_loci[i], ctx$g$odds_ratio[i], ctx$g$raf[i], ctx$g$case_fraction[i], ctx$at))),
  list(g = grid, q = q, lp = lp, at = at))
parallel::stopCluster(cl)

sw$mcid <- unname(mcid[as.character(sw$threshold)])
sw$ratio_specific <- sw$gain_specific / sw$mcid
sw$ratio_diluted <- sw$gain_diluted / sw$mcid
sw$meets_specific <- sw$gain_specific >= sw$mcid
sw$meets_diluted <- sw$gain_diluted >= sw$mcid

# Check: with f = 1 the specific score is the ceiling of script 23 and the diluted weight is the true one
chk <- sw[sw$case_fraction == 1, ]
j <- match(paste(chk$n_loci, chk$odds_ratio, chk$raf, chk$threshold), paste(exact$n_loci, exact$odds_ratio, exact$raf, exact$threshold))
chk$exact_gain <- exact$gain[j]
chk$difference <- chk$gain_specific - chk$exact_gain
chk$beta_difference <- chk$beta_diluted - log(chk$odds_ratio)
check <- data.frame(n_cells = sum(!is.na(chk$exact_gain)), max_abs_gain_difference = max(abs(chk$difference), na.rm = TRUE),
                    max_abs_beta_difference = max(abs(chk$beta_difference)))

# Frontier: the fewest loci reaching the MCID, by score
frontier <- do.call(rbind, lapply(split(sw, list(sw$odds_ratio, sw$raf, sw$case_fraction, sw$threshold), drop = TRUE), function(x) { x <- x[order(x$n_loci), ]
  data.frame(odds_ratio = x$odds_ratio[1], raf = x$raf[1], case_fraction = x$case_fraction[1], threshold = x$threshold[1],
             min_loci_specific = first_reaching(x$n_loci, x$gain_specific, x$mcid[1]), min_loci_diluted = first_reaching(x$n_loci, x$gain_diluted, x$mcid[1])) }))
rownames(frontier) <- NULL

write.csv(sw, "output/aim4_pathogen_sweep.csv", row.names = FALSE)
write.csv(frontier, "output/aim4_pathogen_frontier.csv", row.names = FALSE)
write.csv(check, "output/aim4_pathogen_check.csv", row.names = FALSE)

pl <- frontier[frontier$odds_ratio >= 1.4, ]
pl <- rbind(transform(pl, score = "specific weights", loci = ifelse(is.na(min_loci_specific), 21, min_loci_specific)),
            transform(pl, score = "all-mould weights", loci = ifelse(is.na(min_loci_diluted), 21, min_loci_diluted)))
pl$panel <- paste0("threshold ", round(pl$threshold * 100), "%, frequency ", pl$raf)
pl$fraction <- factor(pl$case_fraction)
dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)
ggsave("output/figures/aim4_pathogen_frontier.png", ggplot(pl, aes(odds_ratio, loci, colour = fraction, linetype = score)) + geom_line() + geom_point() +
         facet_wrap(~ panel, nrow = 2) + scale_y_continuous(breaks = c(1, 5, 10, 15, 20, 21), labels = c(1, 5, 10, 15, 20, ">20")) +
         labs(x = "Odds ratio per risk allele (target-pathogen pathway)", y = "Fewest loci that reach the MCID", colour = "Target-pathogen share of cases", linetype = NULL),
       width = 12, height = 6, dpi = 150)

message("\nCheck against script 23 at f = 1:"); show_df(check)
message("\nFewest loci reaching the 5% MCID (NA = none within 20), frequency 0.3, specific weights:")
w <- frontier[frontier$threshold == 0.05 & frontier$raf == 0.3, ]
print(reshape(w[c("odds_ratio", "case_fraction", "min_loci_specific")], idvar = "odds_ratio", timevar = "case_fraction", direction = "wide"), row.names = FALSE)
message("\nSame, all-mould weights:")
print(reshape(w[c("odds_ratio", "case_fraction", "min_loci_diluted")], idvar = "odds_ratio", timevar = "case_fraction", direction = "wide"), row.names = FALSE)
