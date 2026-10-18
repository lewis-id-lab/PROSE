# Empirical baseline cohort (addendum D8). The baseline cohort for the genetic layer is the observed
# distribution of the 128 clinical profiles in the 4,127 admissions of bs.xlsx, with each profile's
# outcome risk given by the fixed Stanzani equation. This script derives the baseline decision curve and
# MCID on that cohort and compares it with the fitted copula (kept as a sensitivity analysis).
# bs.xlsx is gitignored patient-level data. Only aggregates at the pre-specified thresholds and pairwise
# summaries are written; profile counts are not, because several profiles have fewer than 5 admissions.
# Run from the project root.

source("R/baseline.R")
source("R/copula.R")
source("R/decision_curve.R")
source("R/bs_data.R")

pars <- load_baseline_params()
bs <- read_bs_data("bs.xlsx")
n <- length(bs$y)
risk <- plogis(profile_logit(pars, all_profiles()))
at <- c(0.05, 0.10, 0.20)

cnt <- empirical_profile_counts(bs)
q_emp <- cnt / sum(cnt)
emp <- dca_summary(dca_exact(q_emp, risk, thresholds = at), at = at)

# Sampling uncertainty: bootstrap the admissions (drawing profile counts from a multinomial is the
# same thing). Admission clustering within patients is ignored until the patient link is available.
set.seed(3924)
B <- 2000
boot <- t(vapply(seq_len(B), function(b) {
  qb <- rmultinom(1, n, q_emp)[, 1] / n
  d <- dca_exact(qb, risk, thresholds = at)
  c(d$nb, d$mcid_b, d$flagged)
}, numeric(9)))
colnames(boot) <- c(paste0("nb_", at), paste0("mcid_", at), paste0("flagged_", at))
ci <- function(prefix) t(apply(boot[, paste0(prefix, "_", at)], 2, quantile, c(0.05, 0.95)))

# Cross-check on the observed outcomes (not the model-implied risk)
lp <- as.numeric(pars$intercept + bs$X %*% pars$coef[PREDICTORS])
p <- plogis(lp)
observed <- do.call(rbind, lapply(at, function(t) {
  f <- p >= t
  tp <- sum(f & bs$y == 1) / n; fp <- sum(f & bs$y == 0) / n
  data.frame(threshold = t, observed_nb = tp - fp * t / (1 - t), observed_mcid = 0.2 * tp)
}))

summary_tab <- cbind(emp, nb_lo = ci("nb")[, 1], nb_hi = ci("nb")[, 2],
                     mcid_lo = ci("mcid")[, 1], mcid_hi = ci("mcid")[, 2],
                     flagged_lo = ci("flagged")[, 1], flagged_hi = ci("flagged")[, 2],
                     observed[, c("observed_nb", "observed_mcid")])
dir.create("output", showWarnings = FALSE)
write.csv(summary_tab, "output/dca_empirical_summary.csv", row.names = FALSE)

# Comparison with the fitted copula and with simpler joint distributions
q_fit <- read.csv("output/profile_probabilities.csv")$q
q_ind <- independent_profile_probs(pars$marginal)
q_exch <- profile_probs_exact(exchangeable_corr(0.30), pars$marginal, abseps = 1e-6, seed = 1)
tv <- function(a, b) 0.5 * sum(abs(a - b))
cmp <- data.frame(
  metric = c("admissions", "profiles_observed_of_128", "profiles_with_fewer_than_5_admissions",
             "tv_distance_independence", "tv_distance_exchangeable_0.30", "tv_distance_copula_v2"),
  value = c(n, sum(cnt > 0), sum(cnt > 0 & cnt < 5), tv(q_ind, q_emp), tv(q_exch, q_emp), tv(q_fit, q_emp))
)
write.csv(cmp, "output/empirical_vs_copula_summary.csv", row.names = FALSE)

Rfit <- as.matrix(read.csv("output/fitted_correlation_matrix.csv", row.names = 1))
pairs <- t(combn(PREDICTORS, 2))
corr <- data.frame(predictor_1 = pairs[, 1], predictor_2 = pairs[, 2],
                   empirical_tetrachoric = apply(pairs, 1, function(pr) tetrachoric(bs$X[, pr[1]], bs$X[, pr[2]])),
                   copula_v2 = apply(pairs, 1, function(pr) Rfit[pr[1], pr[2]]))
corr$difference <- corr$copula_v2 - corr$empirical_tetrachoric
write.csv(corr, "output/empirical_vs_copula_correlations.csv", row.names = FALSE)

message("Empirical baseline (", sum(cnt > 0), " of 128 profiles observed):")
print(round(summary_tab[, c("threshold", "nb", "mcid_b", "mcid_lo", "mcid_hi", "flagged", "observed_nb", "observed_mcid")], 5),
      row.names = FALSE)
message("\nTotal variation distance from the empirical joint:")
print(cmp[grepl("^tv_", cmp$metric), ], row.names = FALSE)
