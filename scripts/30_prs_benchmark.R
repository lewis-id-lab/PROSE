# Polygenic score benchmark (addendum 23). Run from the project root, after scripts 08 and 11.
#
# How strong would a polygenic score have to be to reach the MCID on the refit baseline, and how large a discovery study would it
# need? The score is G ~ N(0, 1), independent of the clinical profile, with true effect beta (log OR per SD). Exact over the
# profile distribution of the base cohort (see R/prs_benchmark.R for the model and the accuracy formula).
#
#   1. Ceiling: gain against the OR per SD (true weights), and the OR per SD that reaches the MCID at 5% and 20%.
#   2. Estimated weights: for each true OR per SD, the accuracy rho (correlation of the estimated with the true score) needed.
#   3. Discovery size: the balanced case-control discovery study that gives that accuracy, for 1,000, 10,000 and 60,000 effective
#      independent loci; and the accuracy a discovery study the size of Fisher 2017 (2,609 recipients, 483 cases) would give.
#
# Writes   output/prs_ceiling.csv            per OR per SD and threshold: gain, share of the MCID, aROC of the score alone and combined
#          output/prs_benchmark_summary.csv  per threshold: the OR per SD needed (ceiling) and the aROCs at that point
#          output/prs_gain_by_accuracy.csv   gain against accuracy for selected true ORs per SD
#          output/prs_discovery_size.csv     per true OR per SD, threshold and M_e: accuracy needed and discovery cases needed

source("R/baseline.R")
source("R/bs_data.R")
source("R/decision_curve.R")
source("R/boscore_data.R")
source("R/refit_baseline.R")
source("R/evidence_mapping.R")
source("R/genetic_layer.R")
source("R/prs_benchmark.R")

at <- c(0.05, 0.20)
m_eff <- c(1e3, 1e4, 6e4)
fisher <- list(n = 2609, cases = 483)

base <- select_cohort(load_base_cohort(), include_twins = FALSE)
lp <- profile_logit(load_refit_params())
q <- tabulate(profile_index(base$X), 2L^length(PREDICTORS)) / nrow(base$X)
mcid_tab <- read.csv("output/dca_baseline_refit_summary.csv")
mcid <- setNames(mcid_tab$mcid_b[match(at, mcid_tab$threshold)], as.character(at))

# 1. Ceiling
or_grid <- seq(1.05, 6, by = 0.05)
ceiling_tab <- do.call(rbind, lapply(or_grid, function(o) {
  g <- prs_gain(q, lp, log(o), 1, at)
  data.frame(or_sd = o, g, frac_mcid = g$gain / mcid[as.character(g$threshold)])
}))
summary_tab <- do.call(rbind, lapply(at, function(t) {
  b <- uniroot(function(b) prs_gain(q, lp, b, 1, t)$gain - mcid[[as.character(t)]], c(0.05, 2.5), tol = 1e-6)$root
  g <- prs_gain(q, lp, b, 1, t)
  data.frame(threshold = t, mcid = mcid[[as.character(t)]], or_sd_needed = exp(b), auc_prs_alone = g$auc_prs,
             auc_baseline = g$auc_baseline, auc_combined = g$auc_combined)
}))

# 2. Accuracy needed, and 3. discovery size
beta_grid <- log(seq(1.5, 6, by = 0.1))
discovery <- do.call(rbind, lapply(beta_grid, function(b) {
  a <- prs_shift(q, lp, b)
  do.call(rbind, lapply(at, function(t) {
    rho <- prs_rho_needed(q, lp, b, mcid[[as.character(t)]], t, a = a)
    data.frame(or_sd = exp(b), threshold = t, rho_needed = rho, m_eff = m_eff,
               n_effective = prs_n_effective(rho, b, m_eff),
               cases_balanced = 2 * prs_n_effective(rho, b, m_eff),
               rho_fisher_sized = prs_accuracy(fisher$n, fisher$cases / fisher$n, b, m_eff))
  }))
}))

# Gain against accuracy, for the figure
by_acc <- do.call(rbind, lapply(c(2, 3, 4, 5), function(o) {
  a <- prs_shift(q, lp, log(o))
  do.call(rbind, lapply(seq(0, 1, by = 0.05), function(r) {
    g <- prs_gain(q, lp, log(o), r, at, a = a)
    data.frame(or_sd = o, rho = r, threshold = g$threshold, gain = g$gain, frac_mcid = g$gain / mcid[as.character(g$threshold)])
  }))
}))

dir.create("output", showWarnings = FALSE)
write.csv(ceiling_tab, "output/prs_ceiling.csv", row.names = FALSE)
write.csv(summary_tab, "output/prs_benchmark_summary.csv", row.names = FALSE)
write.csv(by_acc, "output/prs_gain_by_accuracy.csv", row.names = FALSE)
write.csv(discovery, "output/prs_discovery_size.csv", row.names = FALSE)

message("OR per SD needed with true weights:")
print(signif(summary_tab, 4), row.names = FALSE)
