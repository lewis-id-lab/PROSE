# Which outcome risk is "truth" for the baseline cohort? The published equation was fitted with the 49 twin rows in,
# and on the corrected cohort it under-predicts slightly (addendum 3, D9). These helpers score the published
# equation's decisions against alternative outcome risks, to see how much the MCID depends on the choice.
#
# Requires R/baseline.R (PREDICTORS, all_profiles, profile_index, profile_logit) and R/decision_curve.R (dca_exact).

# Logistic refit on the cohort's admissions, for verification and sensitivity only; the study model stays fixed.
refit_coef <- function(cohort) {
  fit <- stats::glm(cohort$rows$case ~ cohort$X, family = stats::binomial)
  stats::setNames(stats::coef(fit), c("(Intercept)", colnames(cohort$X)))
}

profile_risk_from_coef <- function(coef, profiles = all_profiles()) {
  as.numeric(stats::plogis(coef[["(Intercept)"]] + profiles %*% coef[colnames(profiles)]))
}

# Net benefit and MCID for one set of profile weights q under four readings of the outcome risk.
# The decision rule is the published equation's in A, B and D; it is the refit's in C.
#   A  truth = published equation        (the adopted baseline)
#   B  truth = refit                     ("truth equals the refit")
#   C  rule and truth both = refit       (a recalibrated baseline model)
#   D  truth = the observed outcomes     (non-parametric; needs the rows, so it is computed separately)
mcid_scenarios <- function(q, risk_pub, risk_refit, thresholds) {
  one <- function(name, rt, rp) cbind(scenario = name, dca_exact(q, rt, rp, thresholds)[, c("threshold", "nb", "tp", "fp", "flagged", "mcid_b")])
  rbind(one("A: published truth", risk_pub, risk_pub),
        one("B: refit truth, published rule", risk_refit, risk_pub),
        one("C: refit truth and rule", risk_refit, risk_refit))
}

# Scenario D on a set of rows: published predicted risk decides, the observed outcome is the truth
observed_outcome_scenario <- function(cohort, risk_pub, thresholds, rows = seq_len(nrow(cohort$X))) {
  p <- risk_pub[profile_index(cohort$X)][rows]; y <- cohort$rows$case[rows]; n <- length(rows)
  do.call(rbind, lapply(thresholds, function(t) {
    f <- p >= t; tp <- sum(f & y == 1) / n; fp <- sum(f & y == 0) / n
    data.frame(scenario = "D: observed outcomes", threshold = t, nb = tp - fp * t / (1 - t), tp = tp, fp = fp,
               flagged = mean(f), mcid_b = 0.2 * tp)
  }))
}
