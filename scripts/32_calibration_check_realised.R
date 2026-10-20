# The realised-score architecture frontier under the three comparisons (addendum 24, continued from script 31). Run from the project root,
# after scripts 11, 16, 23, 25 and 31.
#
# As script 25 (1 to 20 loci, identical independent loci, base cohort, refit baseline, MCID +0.00458 at 5% and +0.00160 at 20%, weight error
# sigma = 0.05, 0.10, 0.20, 0.40, 12 replicate discovery studies of 2,000 genotypes per architecture and the same seeds), but with the truth,
# comparator and score built under each comparison of R/calibration_check.R:
#   current    global-intercept truth; comparator expit(eta_X)      (reproduces script 25)
#   marginal   the same truth; comparator m(X) = E[p | X]
#   profile    profile-preserving truth; comparator expit(eta_X); the score has an intercept per profile
# The oracle (sigma = 0) frontier under each comparison is computed exactly on the same grid.
#
# Writes   output/calibration_check_realised_sweep.csv     per comparison, architecture, precision and threshold: median gain, 5% and 95%
#                                                          of replicates, share of replicates reaching the MCID
#          output/calibration_check_realised_frontier.csv  fewest loci (and smallest odds ratio) reaching the MCID, by comparison
#          output/calibration_check_realised_summary.csv   by comparison, precision and threshold: share of architectures reaching the
#                                                          MCID (median replicate) and verdicts changed against "current"
#          output/calibration_check_realised_check.csv     simulation against the exact ceiling at sigma = 0, by comparison

source("R/baseline.R"); source("R/bs_data.R"); source("R/decision_curve.R"); source("R/boscore_data.R"); source("R/refit_baseline.R")
source("R/genetic_layer.R"); source("R/genetic_helpers.R"); source("R/aim4_sweep.R"); source("R/calibration_check.R")

at <- c(0.05, 0.20)
comparisons <- c("current", "marginal", "profile")
base <- select_cohort(load_base_cohort(), include_twins = FALSE)
q <- tabulate(profile_index(base$X), 2L^length(PREDICTORS)) / nrow(base$X)
lp <- profile_logit(load_refit_params())
mcid <- setNames(baseline_mcid(q, lp, at), as.character(at))
ref_sweep <- read.csv("output/aim4_realised_sweep.csv")

grid <- expand.grid(n_loci = 1:20, odds_ratio = c(1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 2.0, 2.5, 3.0), raf = c(0.1, 0.3, 0.5), sigma = c(0.05, 0.10, 0.20, 0.40))
grid$seed <- 7000L + seq_len(nrow(grid))
check_grid <- expand.grid(n_loci = c(2, 5, 10, 20), odds_ratio = c(1.5, 2.0, 3.0), raf = c(0.1, 0.3, 0.5), sigma = 0)
check_grid$seed <- 9000L + seq_len(nrow(check_grid))
exact_grid <- unique(grid[, c("n_loci", "odds_ratio", "raf")])
message("grid: ", nrow(grid), " architectures x precisions; exact oracle: ", nrow(exact_grid), " architectures; check: ", nrow(check_grid))

cl <- start_workers(cores = 12, files = c("R/baseline.R", "R/decision_curve.R", "R/genetic_layer.R", "R/aim4_sweep.R", "R/calibration_check.R"))
chunks <- function(n) split(seq_len(n), rep_len(seq_len(60), n))
t0 <- Sys.time()
run <- function(g) parallel_rows(cl, chunks(nrow(g)), function(ch, ctx)
  do.call(rbind, lapply(ch, function(i) calibration_gains_realised(ctx$q, ctx$lp, ctx$g$n_loci[i], ctx$g$odds_ratio[i], ctx$g$raf[i],
                                                                   ctx$g$sigma[i], ctx$at, seed = ctx$g$seed[i]))),
  list(g = g, q = q, lp = lp, at = at))
rep_gains <- run(grid)
chk <- run(check_grid)
exact <- parallel_rows(cl, chunks(nrow(exact_grid)), function(ch, ctx)
  do.call(rbind, lapply(ch, function(i) {
    pmf <- dosage_pmf(ctx$g$n_loci[i], ctx$g$raf[i])
    gr <- list(G = matrix(seq_along(pmf) - 1, ncol = 1, dimnames = list(NULL, "D")), prob = pmf)
    b <- c(D = log(ctx$g$odds_ratio[i]))
    cbind(ctx$g[i, ], calibration_gains(ctx$q, ctx$lp, gr, b, b, ctx$at), row.names = NULL)
  })), list(g = exact_grid, q = q, lp = lp, at = at))
parallel::stopCluster(cl)
message(sprintf("simulation: %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins"))))

rep_gains$mcid <- unname(mcid[as.character(rep_gains$threshold)])
summ <- by_key(rep_gains, c("comparison", "n_loci", "odds_ratio", "raf", "sigma", "threshold"), function(d) data.frame(
  gain_median = median(d$gain), gain_lo = quantile(d$gain, 0.05, names = FALSE), gain_hi = quantile(d$gain, 0.95, names = FALSE),
  share_reaching_mcid = mean(d$gain >= d$mcid), mcid = d$mcid[1]))

# The current comparison must reproduce script 25
cur <- summ[summ$comparison == "current", ]
j <- match(paste(cur$n_loci, cur$odds_ratio, cur$raf, cur$sigma, cur$threshold),
           paste(ref_sweep$n_loci, ref_sweep$odds_ratio, ref_sweep$raf, ref_sweep$sigma, ref_sweep$threshold))
stopifnot(!anyNA(j), nrow(cur) == nrow(ref_sweep),
          max(abs(cur$gain_median - ref_sweep$gain_median[j])) < 1e-9, max(abs(cur$share_reaching_mcid - ref_sweep$share_reaching_mcid[j])) < 1e-12)

# Simulation against the exact oracle at sigma = 0
chk$mcid <- unname(mcid[as.character(chk$threshold)])
check <- by_key(chk, c("comparison", "n_loci", "odds_ratio", "raf", "threshold"), function(d) data.frame(mc_gain_mean = mean(d$gain), mc_gain_sd = sd(d$gain)))
je <- match(paste(check$comparison, check$n_loci, check$odds_ratio, check$raf, check$threshold),
            paste(exact$comparison, exact$n_loci, exact$odds_ratio, exact$raf, exact$threshold))
check$exact_gain <- exact$gain_oracle[je]; check$difference <- check$mc_gain_mean - check$exact_gain

# Frontiers: smallest number of loci (median gain, and at least 80% of replicates) and smallest odds ratio
frontier_of <- function(s, cmp) {
  s <- s[s$comparison == cmp, ]
  fr <- do.call(rbind, lapply(split(s, list(s$odds_ratio, s$raf, s$sigma, s$threshold), drop = TRUE), function(x) { x <- x[order(x$n_loci), ]
    data.frame(comparison = cmp, kind = "min_loci", odds_ratio = x$odds_ratio[1], raf = x$raf[1], sigma = x$sigma[1], n_loci = NA_real_, threshold = x$threshold[1],
               median_reaches = first_reaching(x$n_loci, x$gain_median, x$mcid[1]),
               reliable_reaches = { i <- which(x$share_reaching_mcid >= 0.8); if (length(i)) x$n_loci[min(i)] else NA_real_ }) }))
  fo <- do.call(rbind, lapply(split(s, list(s$n_loci, s$raf, s$sigma, s$threshold), drop = TRUE), function(x) { x <- x[order(x$odds_ratio), ]
    data.frame(comparison = cmp, kind = "min_odds_ratio", odds_ratio = NA_real_, raf = x$raf[1], sigma = x$sigma[1], n_loci = x$n_loci[1], threshold = x$threshold[1],
               median_reaches = first_reaching(x$odds_ratio, x$gain_median, x$mcid[1]),
               reliable_reaches = { i <- which(x$share_reaching_mcid >= 0.8); if (length(i)) x$odds_ratio[min(i)] else NA_real_ }) }))
  rbind(fr, fo)
}
oracle_of <- function(cmp) {
  ex <- exact[exact$comparison == cmp, ]
  ex$mcid <- unname(mcid[as.character(ex$threshold)])
  rbind(do.call(rbind, lapply(split(ex, list(ex$odds_ratio, ex$raf, ex$threshold), drop = TRUE), function(x) { x <- x[order(x$n_loci), ]
    r <- first_reaching(x$n_loci, x$gain_oracle, x$mcid[1])
    data.frame(comparison = cmp, kind = "min_loci", odds_ratio = x$odds_ratio[1], raf = x$raf[1], sigma = 0, n_loci = NA_real_, threshold = x$threshold[1], median_reaches = r, reliable_reaches = r) })),
    do.call(rbind, lapply(split(ex, list(ex$n_loci, ex$raf, ex$threshold), drop = TRUE), function(x) { x <- x[order(x$odds_ratio), ]
    r <- first_reaching(x$odds_ratio, x$gain_oracle, x$mcid[1])
    data.frame(comparison = cmp, kind = "min_odds_ratio", odds_ratio = NA_real_, raf = x$raf[1], sigma = 0, n_loci = x$n_loci[1], threshold = x$threshold[1], median_reaches = r, reliable_reaches = r) })))
}
frontier <- do.call(rbind, lapply(comparisons, function(cmp) rbind(oracle_of(cmp), frontier_of(summ, cmp))))
rownames(frontier) <- NULL

# Verdicts: does the median replicate reach the MCID? Changes against "current", by comparison, precision and threshold
summ$meets <- summ$gain_median >= summ$mcid
cur_meets <- with(summ[summ$comparison == "current", ], setNames(meets, paste(n_loci, odds_ratio, raf, sigma, threshold)))
summary_tab <- by_key(summ, c("comparison", "sigma", "threshold"), function(x) {
  key <- paste(x$n_loci, x$odds_ratio, x$raf, x$sigma, x$threshold)
  data.frame(architectures = nrow(x), share_reaching_mcid = mean(x$meets), verdicts_changed = sum(x$meets != cur_meets[key]),
             lost = sum(!x$meets & cur_meets[key]), gained = sum(x$meets & !cur_meets[key]))
})

dir.create("output", showWarnings = FALSE)
write.csv(summ[, setdiff(names(summ), "meets")], "output/calibration_check_realised_sweep.csv", row.names = FALSE)
write.csv(frontier, "output/calibration_check_realised_frontier.csv", row.names = FALSE)
write.csv(summary_tab, "output/calibration_check_realised_summary.csv", row.names = FALSE)
write.csv(check, "output/calibration_check_realised_check.csv", row.names = FALSE)

message("\nSimulation against the exact oracle at sigma = 0: largest absolute difference in gain by comparison")
print(aggregate(difference ~ comparison, check, function(x) max(abs(x))), row.names = FALSE)
message("\nShare of architectures whose median replicate reaches the MCID, and verdicts changed against 'current':")
print(summary_tab, row.names = FALSE)
