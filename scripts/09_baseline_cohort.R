# Baseline decision curve and MCID on the base cohort (addendum 3, D9): the 4,078 unique admissions of
# boscore.xlsx, with the 49 non-event twins removed and patients resampled in the bootstrap. Run from the
# project root after scripts/08_reconcile_boscore.R.
#
# Reads    data/restricted/base_cohort_pseudonymised.csv (or boscore.xlsx if that is absent)
# Writes   output/dca_baseline_summary.csv        primary: base cohort, patient bootstrap
#          output/dca_baseline_sensitivity.csv    the 4,127 rows and the admission-level bootstrap, for comparison
# Only aggregates at the pre-specified thresholds are written; profile counts are not.

source("R/baseline.R")
source("R/bs_data.R")
source("R/decision_curve.R")
source("R/boscore_data.R")

pars <- load_baseline_params()
cohort <- load_base_cohort()
base <- select_cohort(cohort, include_twins = FALSE)
rows <- select_cohort(cohort, include_twins = TRUE)
risk <- plogis(profile_logit(pars, all_profiles()))
at <- c(0.05, 0.10, 0.20)
B <- 2000L

point <- function(x) {
  q <- tabulate(profile_index(x$X), 2L^length(PREDICTORS)) / nrow(x$X)
  dca_summary(dca_exact(q, risk, thresholds = at), at = at)
}
ci <- function(m, prefix) t(apply(m[, paste0(prefix, "_", at)], 2, quantile, c(0.05, 0.95)))
summarise <- function(x, m) {
  p <- point(x)
  cbind(p, nb_lo = ci(m, "nb")[, 1], nb_hi = ci(m, "nb")[, 2], mcid_lo = ci(m, "mcid")[, 1], mcid_hi = ci(m, "mcid")[, 2],
        flagged_lo = ci(m, "flagged")[, 1], flagged_hi = ci(m, "flagged")[, 2])
}

# Primary: base cohort, patients resampled
boot_base <- bootstrap_dca(base, risk, at, B = B, seed = 6612, by = "patient")
primary <- summarise(base, boot_base)

# Cross-check on the observed outcomes (not the model-implied risk)
lp <- as.numeric(pars$intercept + base$X %*% pars$coef[colnames(base$X)])
p <- plogis(lp); y <- base$rows$case; n <- length(y)
observed <- do.call(rbind, lapply(at, function(t) {
  f <- p >= t; tp <- sum(f & y == 1) / n; fp <- sum(f & y == 0) / n
  data.frame(observed_nb = tp - fp * t / (1 - t), observed_mcid = 0.2 * tp)
}))
primary <- cbind(cohort = "4,078 unique admissions (base cohort)", bootstrap = "patients", primary, observed)

# Sensitivity: the 4,127 rows with twins, patient and admission bootstrap; and the base cohort resampled by admission
sens <- rbind(
  cbind(cohort = "4,127 rows with twins", bootstrap = "patients", summarise(rows, bootstrap_dca(rows, risk, at, B = B, seed = 6613, by = "patient"))),
  cbind(cohort = "4,127 rows with twins", bootstrap = "admissions", summarise(rows, bootstrap_dca(rows, risk, at, B = B, seed = 6614, by = "admission"))),
  cbind(cohort = "4,078 unique admissions (base cohort)", bootstrap = "admissions", summarise(base, bootstrap_dca(base, risk, at, B = B, seed = 6615, by = "admission")))
)

dir.create("output", showWarnings = FALSE)
write.csv(primary, "output/dca_baseline_summary.csv", row.names = FALSE)
write.csv(sens, "output/dca_baseline_sensitivity.csv", row.names = FALSE)

message("Base cohort: ", nrow(base$X), " admissions, ", length(unique(base$rows$patient)), " patients, ", sum(y), " events")
print(round(primary[, c("threshold", "nb", "mcid_b", "mcid_lo", "mcid_hi", "flagged", "observed_nb", "observed_mcid")], 5), row.names = FALSE)
message("\nSensitivity (MCID with 90% interval):")
print(cbind(sens[, c("cohort", "bootstrap", "threshold")], round(sens[, c("mcid_b", "mcid_lo", "mcid_hi")], 5)), row.names = FALSE)
