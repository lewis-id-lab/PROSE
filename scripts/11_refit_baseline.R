# Refit the baseline model on the base cohort (addendum 6, D10) and re-derive its performance and MCID.
# Run from the project root after scripts/08_reconcile_boscore.R.
#
# Base cohort: 4,078 unique admissions in 1,945 patients (twin rows removed, D9). The refit is a maximum-likelihood
# logistic regression on the seven predictors, with patient-clustered standard errors and a patient-level bootstrap.
# Optimism is reported, not subtracted (a bootstrap refit in every resample; Harrell's method).
#
# Writes   output/baseline_refit_parameters.csv     the fixed baseline: coefficients, robust SE, bootstrap interval, published value
#          output/baseline_refit_performance.csv    aROC, calibration slope, calibration in the large: apparent, optimism, corrected
#          output/dca_baseline_refit_summary.csv    net benefit and MCID at 5%, 10%, 20% for the refit baseline

source("R/baseline.R")
source("R/bs_data.R")
source("R/decision_curve.R")
source("R/boscore_data.R")
source("R/outcome_risk.R")
source("R/refit_baseline.R")

pars <- load_baseline_params()
base <- select_cohort(load_base_cohort(), include_twins = FALSE)
y <- base$rows$case; n <- length(y)
at <- c(0.05, 0.10, 0.20)
nprof <- 2L^length(PREDICTORS)

fit <- fit_baseline(base)
risk <- profile_risk_from_coef(fit$coef)
q <- tabulate(profile_index(base$X), nprof) / n
point <- dca_summary(dca_exact(q, risk, risk, thresholds = at), at = at)

# Patient-level bootstrap: refit in each resample, recompute the profile weights, the curve and the MCID
groups <- split(seq_len(n), base$rows$patient)
B <- 2000L
set.seed(3319)
coefs <- matrix(NA_real_, B, length(fit$coef), dimnames = list(NULL, names(fit$coef)))
curves <- matrix(NA_real_, B, 3 * length(at))
for (b in seq_len(B)) {
  rows <- resample_group_rows(groups)
  cf <- tryCatch(suppressWarnings(stats::coef(stats::glm(y[rows] ~ base$X[rows, , drop = FALSE], family = stats::binomial))), error = function(e) NULL)
  if (is.null(cf) || anyNA(cf)) next
  names(cf) <- names(fit$coef); coefs[b, ] <- cf
  rb <- profile_risk_from_coef(cf)
  s <- dca_exact(tabulate(profile_index(base$X[rows, , drop = FALSE]), nprof) / length(rows), rb, rb, thresholds = at)
  curves[b, ] <- c(s$nb, s$mcid_b, s$flagged)
}
ok <- stats::complete.cases(coefs)
colnames(curves) <- c(paste0("nb_", at), paste0("mcid_", at), paste0("flagged_", at))
qci <- function(m, prefix) t(apply(curves[ok, paste0(prefix, "_", at), drop = FALSE], 2, stats::quantile, c(0.05, 0.95)))

params <- data.frame(term = names(fit$coef), published = c(pars$intercept, pars$coef[PREDICTORS]), refit = unname(fit$coef),
                     se_robust = unname(fit$se_robust), lo = apply(coefs[ok, ], 2, stats::quantile, 0.05),
                     hi = apply(coefs[ok, ], 2, stats::quantile, 0.95), row.names = NULL)
params$odds_ratio <- exp(params$refit)

# Optimism (observed outcomes), then the performance table
opt <- optimism_bootstrap(base, at, B = B, seed = 3320)
lp <- drop(fit$coef[1] + base$X %*% fit$coef[-1])
app_auc <- auc_rank(lp, y)
opt_auc <- mean(opt[, "auc_apparent"] - opt[, "auc_test"])
slope_corr <- mean(opt[, "slope_test"])
perf <- data.frame(
  metric = c("aROC", "calibration slope", "mean predicted risk", "observed prevalence",
             paste0("observed net benefit at ", at * 100, "%")),
  apparent = c(app_auc, 1, mean(fit$fitted), mean(y), observed_nb_rows(stats::plogis(lp), y, at)),
  optimism = c(opt_auc, 1 - slope_corr, 0, 0, vapply(at, function(t) mean(opt[, paste0("nb_apparent_", t)] - opt[, paste0("nb_test_", t)]), numeric(1)))
)
perf$corrected <- perf$apparent - perf$optimism

summary_tab <- cbind(point, nb_lo = qci(curves, "nb")[, 1], nb_hi = qci(curves, "nb")[, 2],
                     mcid_lo = qci(curves, "mcid")[, 1], mcid_hi = qci(curves, "mcid")[, 2],
                     flagged_lo = qci(curves, "flagged")[, 1], flagged_hi = qci(curves, "flagged")[, 2],
                     observed_nb_apparent = perf$apparent[5:7], observed_nb_corrected = perf$corrected[5:7])

dir.create("output", showWarnings = FALSE)
write.csv(params, "output/baseline_refit_parameters.csv", row.names = FALSE)
write.csv(perf, "output/baseline_refit_performance.csv", row.names = FALSE)
write.csv(summary_tab, "output/dca_baseline_refit_summary.csv", row.names = FALSE)

message("bootstrap resamples used: ", sum(ok), " of ", B)
message("\nRefit baseline (robust SE, patient-bootstrap 90% interval):")
print(cbind(params[1], round(params[-1], 3)), row.names = FALSE)
message("\nPerformance (apparent, optimism, corrected):")
print(cbind(perf[1], round(perf[-1], 4)), row.names = FALSE)
message("\nRefit baseline decision curve and MCID:")
print(round(summary_tab[, c("threshold", "nb", "tp_per_1000", "fp_per_1000", "flagged", "mcid_b", "mcid_lo", "mcid_hi")], 5), row.names = FALSE)
