# Exploratory (addendum 48, not pre-specified): the genetic-layer gains over a grid of risk thresholds from 2% to 30% (including 10%),
# under the profile-preserving truth, for the four- and five-term panels. Run from the project root, after scripts 11, 16, 22 and 31.
# Needs the patient-level cohort (data/restricted/ or boscore.xlsx) for the profile distribution, as script 31 does.
#
# The pre-specified thresholds were 5% (primary) and 20% (secondary). This script reports every threshold on the same footing; nothing here
# changes the primary or secondary analyses. The MCID at each threshold is the pre-specified rule (20% more true cases at no extra false
# positives, 0.2 x TP_clin / n), applied to the refit baseline at that threshold. The run stops unless it reproduces, at 5% and 20%, the
# published gains of script 31 (output/calibration_check_layer.csv) and, at 5%, 10% and 20%, the MCID of script 11.
#
# Writes   output/threshold_curve_profile.csv            per panel, reading, score and threshold: gain (median, 90% range), MCID, gain / MCID,
#                                                         share of draws reaching the MCID, extra true cases and flagged share
#          output/figures/fig_threshold_curve_profile.png gain against the MCID across thresholds (pooled-mean reading)
#
# Dry run without patient data (code check only; the numbers are NOT the cohort's): PROSE_Q_FILE=output/profile_probabilities.csv
# (the copula reconstruction) with PROSE_OUT_DIR set to a scratch folder and PROSE_NDRAW to a small number. The reproduction checks are skipped.

source("R/baseline.R"); source("R/decision_curve.R"); source("R/refit_baseline.R"); source("R/boscore_data.R")
source("R/evidence_mapping.R"); source("R/genetic_layer.R"); source("R/genetic_helpers.R"); source("R/aim4_sweep.R")
source("R/prs_benchmark.R"); source("R/calibration_check.R")
suppressPackageStartupMessages({library(ggplot2)})
source("R/theme_jama.R")

thresholds <- round(seq(0.02, 0.30, by = 0.01), 2)
q_file <- Sys.getenv("PROSE_Q_FILE", "")
dry_run <- nzchar(q_file)
out_dir <- Sys.getenv("PROSE_OUT_DIR", "output")
n_draw <- as.integer(Sys.getenv("PROSE_NDRAW", "2000"))
if (!dry_run && out_dir != "output") stop("PROSE_OUT_DIR is for dry runs only", call. = FALSE)
if (dry_run && out_dir == "output") stop("a dry run must not write to output/: set PROSE_OUT_DIR", call. = FALSE)

lp <- profile_logit(load_refit_params())
q <- if (dry_run) {
  message("DRY RUN: profile distribution from ", q_file, "; results are not the cohort's")
  read.csv(q_file)$q
} else {
  base <- select_cohort(load_base_cohort(), include_twins = FALSE)
  tabulate(profile_index(base$X), 2L^length(PREDICTORS)) / nrow(base$X)
}
stopifnot(length(q) == length(lp), abs(sum(q) - 1) < 1e-9)

# The MCID at each threshold: 20% of the refit baseline's true positives (exact under the profile-preserving truth)
mcid <- setNames(baseline_mcid(q, lp, thresholds), as.character(thresholds))
if (!dry_run) {
  ref <- read.csv("output/dca_baseline_refit_summary.csv")
  stopifnot(max(abs(mcid[as.character(ref$threshold)] - ref$mcid_b)) < 1e-9)
}

cl <- start_workers(cores = max(1L, min(10L, parallel::detectCores() - 1L)),
                    files = c("R/baseline.R", "R/decision_curve.R", "R/genetic_layer.R", "R/aim4_sweep.R", "R/prs_benchmark.R", "R/calibration_check.R"))
on.exit(parallel::stopCluster(cl), add = TRUE)

run_draws <- function(grid, effects, beta_score) {
  jobs <- do.call(rbind, lapply(names(effects), function(h) data.frame(truth = h, draw = seq_len(min(n_draw, nrow(effects[[h]]))))))
  ctx <- list(q = q, lp = lp, grid = grid, effects = effects, beta_score = beta_score, at = thresholds, jobs = jobs)
  parallel_rows(cl, split(seq_len(nrow(jobs)), rep_len(seq_len(60), nrow(jobs))), function(ix, ctx) {
    do.call(rbind, lapply(ix, function(i) {
      j <- ctx$jobs[i, ]
      cbind(truth = j$truth, draw = j$draw,
            calibration_gains(ctx$q, ctx$lp, ctx$grid, ctx$effects[[j$truth]][j$draw, ], ctx$beta_score, ctx$at, comparisons = "profile", detail = TRUE))
    }))
  }, ctx)
}

eff <- read.csv("output/genetic_layer_effects.csv")
terms4 <- c("PTX3", "CLEC7A", "TLR4", "CD209")
effects4 <- lapply(c(mean = "mean", predictive = "predictive"), function(h) as.matrix(eff[eff$truth == h, terms4]))
d4 <- run_draws(dosage_grid(genetic_term_frequencies()), effects4, pooled_estimates(load_term_draws()))
s5 <- clec7a_second_setup()
d5 <- run_draws(s5$grid, s5$effects, s5$beta_hat)

summarise_panel <- function(d, panel) {
  long <- rbind(cbind(d, score = "ceiling", gain = d$gain_oracle, extra_tp = d$extra_tp_oracle, flagged = d$flagged_oracle),
                cbind(d, score = "realised", gain = d$gain_score, extra_tp = d$extra_tp_score, flagged = d$flagged_score))
  out <- by_key(long, c("truth", "threshold", "score"), function(x) {
    m <- mcid[[as.character(x$threshold[1])]]
    data.frame(gain = median(x$gain), lo = quantile(x$gain, 0.05, names = FALSE), hi = quantile(x$gain, 0.95, names = FALSE),
               mcid = m, gain_over_mcid = median(x$gain) / m, p_meets_mcid = mean(x$gain >= m),
               extra_cases_per_1000 = 1000 * median(x$extra_tp), mcid_cases_per_1000 = 1000 * m,
               flagged = median(x$flagged), flagged_clinical = median(x$flagged_clinical), nb_clinical = median(x$nb_clinical))
  })
  cbind(panel = panel, out)
}
curve <- rbind(summarise_panel(d4, "four terms"), summarise_panel(d5, "five terms"))
curve <- curve[order(curve$panel != "four terms", curve$truth, curve$score, curve$threshold), ]

# At 5% and 20% the curve must reproduce script 31's profile-preserving results
if (!dry_run) {
  ref <- subset(read.csv("output/calibration_check_layer.csv"), comparison == "profile")
  m <- merge(curve, ref, by = c("panel", "truth", "score", "threshold"), suffixes = c("", "_ref"))
  stopifnot(nrow(m) == nrow(ref), max(abs(m$gain - m$gain_ref)) < 1e-9, max(abs(m$p_meets_mcid - m$p_meets_mcid_ref)) < 1e-9,
            max(abs(m$extra_cases_per_1000 - m$extra_cases_per_1000_ref)) < 1e-6)
}

dir.create(file.path(out_dir, "figures"), recursive = TRUE, showWarnings = FALSE)
write.csv(curve, file.path(out_dir, "threshold_curve_profile.csv"), row.names = FALSE)

fig <- subset(curve, truth == "mean")
fig$panel <- factor(ifelse(fig$panel == "four terms", "Four terms (prespecified panel)", "Five terms (CLEC7A A-G added)"),
                    levels = c("Four terms (prespecified panel)", "Five terms (CLEC7A A-G added)"))
fig$score <- ifelse(fig$score == "ceiling", "true effects in the score", "pooled estimates in the score")
mc <- unique(fig[, c("panel", "threshold", "mcid")])
ggsave(file.path(out_dir, "figures", "fig_threshold_curve_profile.png"),
  ggplot(fig, aes(threshold, gain, colour = score, fill = score)) +
    geom_vline(xintercept = c(0.05, 0.20), colour = jama_cols[["light"]], linewidth = 0.8) +
    geom_hline(yintercept = 0, colour = jama_cols[["grey"]], linewidth = 0.3) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, colour = NA) +
    geom_line(linewidth = 0.8) +
    geom_line(data = mc, aes(threshold, mcid), inherit.aes = FALSE, linetype = "dashed", colour = jama_cols[["ink"]], linewidth = 0.5) +
    facet_wrap(~ panel) +
    scale_x_continuous(labels = function(x) paste0(round(100 * x), "%"), breaks = c(0.02, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30)) +
    scale_colour_manual(values = c("pooled estimates in the score" = jama_cols[["steel"]], "true effects in the score" = jama_cols[["navy"]])) +
    scale_fill_manual(values = c("pooled estimates in the score" = jama_cols[["steel"]], "true effects in the score" = jama_cols[["navy"]])) +
    labs(x = "Risk threshold", y = "Net-benefit gain over the clinical score", colour = NULL, fill = NULL) +
    theme_jama(),
  width = 10, height = 4.5, dpi = 300)

print(subset(curve, truth == "mean" & threshold %in% c(0.05, 0.10, 0.15, 0.20, 0.25, 0.30),
             c("panel", "score", "threshold", "gain", "mcid", "gain_over_mcid", "p_meets_mcid", "extra_cases_per_1000", "flagged_clinical")), row.names = FALSE)
