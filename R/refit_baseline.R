# The refit baseline model (addendum 6, D10). The clinical block is a maximum-likelihood logistic refit of the seven
# predictors on the base cohort (4,078 unique admissions in 1,945 patients), replacing the published equation as
# the fixed baseline. The published equation is kept as a sensitivity analysis.
#
# Requires R/baseline.R (PREDICTORS, all_profiles), R/bs_data.R (auc_rank), R/boscore_data.R (resample_group_rows).

design_matrix <- function(X) cbind(`(Intercept)` = 1, X)

# Cluster-robust (patient) variance of logistic-regression coefficients: bread %*% meat %*% bread, with the
# usual small-sample factor g/(g-1) * (n-1)/(n-k). small_sample = FALSE with one row per cluster gives HC0.
cluster_robust_vcov <- function(D, y, p, cluster, small_sample = TRUE) {
  bread <- solve(crossprod(D * (p * (1 - p)), D))
  meat <- crossprod(rowsum(D * (y - p), cluster))
  g <- length(unique(cluster)); n <- nrow(D); k <- ncol(D)
  adj <- if (small_sample) g / (g - 1) * (n - 1) / (n - k) else 1
  adj * bread %*% meat %*% bread
}

fit_baseline <- function(cohort) {
  y <- cohort$rows$case
  fit <- stats::glm(y ~ cohort$X, family = stats::binomial)
  cf <- stats::setNames(stats::coef(fit), c("(Intercept)", colnames(cohort$X)))
  p <- as.numeric(stats::fitted(fit))
  v <- cluster_robust_vcov(design_matrix(cohort$X), y, p, cohort$rows$patient)
  list(coef = cf, fitted = p, se_model = sqrt(diag(stats::vcov(fit))), se_robust = sqrt(diag(v)))
}

# A parameter list that works anywhere the published one does (profile_logit, net_benefit, ...)
refit_params <- function(coef) {
  list(intercept = coef[["(Intercept)"]], coef = coef[PREDICTORS])
}

load_refit_params <- function(path = "output/baseline_refit_parameters.csv") {
  d <- utils::read.csv(path, stringsAsFactors = FALSE)
  cf <- stats::setNames(d$refit, d$term)
  refit_params(cf)
}

# Observed net benefit of a rule: flag when the predicted risk reaches the threshold, outcomes as observed
observed_nb_rows <- function(p, y, thresholds) {
  n <- length(y)
  vapply(thresholds, function(t) {
    f <- p >= t
    sum(f & y == 1) / n - sum(f & y == 0) / n * t / (1 - t)
  }, numeric(1))
}

# Optimism of the refit by a patient-level bootstrap (Harrell): refit in each resample, then compare the
# performance on the resample with the performance on the original cohort. Net benefit uses the observed outcomes.
optimism_bootstrap <- function(cohort, thresholds, B = 2000L, seed) {
  X <- cohort$X; y <- cohort$rows$case
  groups <- split(seq_len(nrow(X)), cohort$rows$patient)
  set.seed(seed)
  out <- vector("list", B)
  for (b in seq_len(B)) {
    rows <- resample_group_rows(groups)
    Xb <- X[rows, , drop = FALSE]; yb <- y[rows]
    cf <- tryCatch(suppressWarnings(stats::coef(stats::glm(yb ~ Xb, family = stats::binomial))), error = function(e) NULL)
    if (is.null(cf) || anyNA(cf)) next
    lp_b <- drop(cf[1] + Xb %*% cf[-1]); lp_o <- drop(cf[1] + X %*% cf[-1])
    slope <- unname(suppressWarnings(stats::coef(stats::glm(y ~ lp_o, family = stats::binomial))[2]))
    out[[b]] <- c(auc_apparent = auc_rank(lp_b, yb), auc_test = auc_rank(lp_o, y), slope_test = slope,
                  stats::setNames(observed_nb_rows(stats::plogis(lp_b), yb, thresholds), paste0("nb_apparent_", thresholds)),
                  stats::setNames(observed_nb_rows(stats::plogis(lp_o), y, thresholds), paste0("nb_test_", thresholds)))
  }
  do.call(rbind, out)
}
