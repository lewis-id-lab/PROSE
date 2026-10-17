# Track A3: decision curves for the baseline model, thresholds 1-30%. Every curve is an exact
# function of the profile distribution; the profile distribution itself comes from numerical
# integration of the fitted copula (script 01) and, for the ensemble, of each sampled matrix
# (script 02). No patient-level sampling is involved. Run from the project root.

source("R/baseline.R")
source("R/copula.R")
source("R/decision_curve.R")
suppressPackageStartupMessages({library(ggplot2)})

cores <- max(1L, min(15L, parallel::detectCores() - 1L))
pars <- load_baseline_params()
profiles <- all_profiles()
risk <- plogis(profile_logit(pars, profiles))
grid <- seq(0.01, 0.30, by = 0.001)
dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)

# 1. Fitted cohort
pt <- read.csv("output/profile_probabilities.csv")
stopifnot(all.equal(pt$risk, risk))
fit_curve <- dca_exact(pt$q, risk, thresholds = grid)
write.csv(fit_curve, "output/dca_fitted.csv", row.names = FALSE)

# 2. Ensemble of near-optimal matrices, each integrated exactly
members <- read.csv("output/correlation_ensemble_members.csv")
member_theta <- lapply(seq_len(nrow(members)), function(i) unlist(members[i, THETA_NAMES]))
member_curve <- function(th) {
  R <- pd_adjust(build_corr(th))$R
  q <- profile_probs_exact(R, pars$marginal, abseps = 1e-6, seed = 4401)
  cv <- dca_exact(q, plogis(profile_logit(pars)), thresholds = seq(0.01, 0.30, by = 0.005))
  m <- profile_metrics(q, pars)
  cv$prevalence_model <- m$prevalence
  cv$auc <- m$auc
  cv
}
curves <- run_parallel(member_theta, member_curve, cores = cores)
ens <- do.call(rbind, Map(function(cv, id) cbind(member = id, cv), curves, members$member))
write.csv(ens, "output/dca_ensemble.csv", row.names = FALSE)

band_of <- function(v) {
  tapply(ens[[v]], ens$threshold, function(x) quantile(x, c(0.05, 0.5, 0.95))) |>
    (\(l) data.frame(threshold = as.numeric(names(l)), lo = sapply(l, `[`, 1),
                     median = sapply(l, `[`, 2), hi = sapply(l, `[`, 3), row.names = NULL))()
}
bands <- lapply(c(nb = "nb", mcid_b = "mcid_b", tp = "tp", fp = "fp", flagged = "flagged"), band_of)
write.csv(do.call(rbind, Map(function(b, n) cbind(measure = n, b), bands, names(bands))),
          "output/dca_ensemble_band.csv", row.names = FALSE)

# 3. Summary at the pre-specified thresholds: fitted value and ensemble range
at <- c(0.05, 0.10, 0.20)
summary_tab <- do.call(rbind, lapply(at, function(t) {
  f <- dca_summary(fit_curve, t)
  e <- ens[abs(ens$threshold - t) < 1e-9, ]
  data.frame(f, nb_lo = quantile(e$nb, 0.05), nb_hi = quantile(e$nb, 0.95),
             mcid_lo = quantile(e$mcid_b, 0.05), mcid_hi = quantile(e$mcid_b, 0.95),
             flagged_lo = quantile(e$flagged, 0.05), flagged_hi = quantile(e$flagged, 0.95),
             row.names = NULL)
}))
write.csv(summary_tab, "output/dca_summary.csv", row.names = FALSE)

# 4. Figures
band_nb <- bands$nb[, c("threshold", "lo", "hi")]
band_mc <- bands$mcid_b[, c("threshold", "lo", "hi")]
ggsave("output/figures/dca_baseline.png", plot_dca(fit_curve, band_nb), width = 7, height = 4.5, dpi = 150)
ggsave("output/figures/dca_rates.png", plot_dca_rates(fit_curve), width = 8, height = 3.5, dpi = 150)
ggsave("output/figures/dca_mcid.png", plot_mcid(fit_curve, band_mc), width = 7, height = 4, dpi = 150)
