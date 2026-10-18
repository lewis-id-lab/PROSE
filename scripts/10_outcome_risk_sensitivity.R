# How much does the MCID depend on the choice of outcome risk? Run from the project root after scripts/08.
#
# Base cohort: 4,078 unique admissions in 1,945 patients. The decision rule is the published equation's
# unless stated; "truth" is varied:
#   A published equation (adopted baseline), B refit on the cohort with the published rule, C refit rule and truth,
#   D the observed outcomes.
# Patients are resampled, and the refit is repeated inside every resample, so its uncertainty is included.
#
# Writes   output/mcid_outcome_risk_sensitivity.csv   scenario, threshold, point value, 90% interval
#          output/refit_coefficients.csv              published against refit, with a patient-bootstrap interval

source("R/baseline.R")
source("R/bs_data.R")
source("R/decision_curve.R")
source("R/boscore_data.R")
source("R/outcome_risk.R")

pars <- load_baseline_params()
base <- select_cohort(load_base_cohort(), include_twins = FALSE)
risk_pub <- plogis(profile_logit(pars, all_profiles()))
at <- c(0.05, 0.10, 0.20)
nprof <- 2L^length(PREDICTORS)

point_scenarios <- function(cohort, rows = seq_len(nrow(cohort$X)), coef = refit_coef(cohort)) {
  q <- tabulate(profile_index(cohort$X)[rows], nprof) / length(rows)
  list(table = rbind(mcid_scenarios(q, risk_pub, profile_risk_from_coef(coef), at),
                     observed_outcome_scenario(cohort, risk_pub, at, rows)[, c("scenario", "threshold", "nb", "tp", "fp", "flagged", "mcid_b")]),
       coef = coef)
}

est <- point_scenarios(base)

# Patient-level bootstrap with the refit redone in each resample
groups <- split(seq_len(nrow(base$X)), base$rows$patient)
B <- 2000L
set.seed(8421)
failed <- 0L
draws <- vector("list", B); coefs <- vector("list", B)
for (b in seq_len(B)) {
  rows <- resample_group_rows(groups)
  sub <- list(rows = base$rows[rows, , drop = FALSE], X = base$X[rows, , drop = FALSE])
  res <- tryCatch(suppressWarnings(point_scenarios(sub)), error = function(e) NULL)
  if (is.null(res) || any(!is.finite(res$coef))) { failed <- failed + 1L; next }
  draws[[b]] <- transform(res$table, rep = b); coefs[[b]] <- res$coef
}
draws <- do.call(rbind, draws); coefs <- do.call(rbind, coefs)
key <- paste(est$table$scenario, est$table$threshold)
lo <- tapply(draws$mcid_b, paste(draws$scenario, draws$threshold), quantile, 0.05)[key]
hi <- tapply(draws$mcid_b, paste(draws$scenario, draws$threshold), quantile, 0.95)[key]
nb_lo <- tapply(draws$nb, paste(draws$scenario, draws$threshold), quantile, 0.05)[key]
nb_hi <- tapply(draws$nb, paste(draws$scenario, draws$threshold), quantile, 0.95)[key]

out <- transform(est$table, tp_per_1000 = 1000 * tp, fp_per_1000 = 1000 * fp, mcid_lo = as.numeric(lo), mcid_hi = as.numeric(hi),
                 nb_lo = as.numeric(nb_lo), nb_hi = as.numeric(nb_hi))
out$mcid_relative_to_A <- out$mcid_b / out$mcid_b[match(paste("A: published truth", out$threshold), key)]
out <- out[, c("scenario", "threshold", "nb", "nb_lo", "nb_hi", "tp_per_1000", "fp_per_1000", "flagged", "mcid_b", "mcid_lo", "mcid_hi", "mcid_relative_to_A")]
coef_tab <- data.frame(term = names(est$coef), published = c(pars$intercept, pars$coef[PREDICTORS]), refit = unname(est$coef),
                       lo = apply(coefs, 2, quantile, 0.05), hi = apply(coefs, 2, quantile, 0.95), row.names = NULL)
dir.create("output", showWarnings = FALSE)
write.csv(out, "output/mcid_outcome_risk_sensitivity.csv", row.names = FALSE)
write.csv(coef_tab, "output/refit_coefficients.csv", row.names = FALSE)

message("bootstrap resamples used: ", B - failed, " of ", B)
print(cbind(out[, c("scenario", "threshold")], round(out[, c("mcid_b", "mcid_lo", "mcid_hi", "mcid_relative_to_A")], 5)), row.names = FALSE)
message("\nRefit coefficients against the published ones (patient-bootstrap 90% interval):")
print(cbind(coef_tab[1], round(coef_tab[-1], 2)), row.names = FALSE)
