# Checks of the admission-level derivation data (bs.xlsx) against the published Stanzani 2019
# counts and coefficients. bs.xlsx is patient-level data and is gitignored; these functions only
# ever return aggregates.

# File column -> model predictor
BS_VARIABLES <- c(active_malig = "uncontrolled_malignancy", hr_chemo = "high_risk_chemo",
                  hd_steroids = "high_dose_steroids", lymph = "severe_lymphopenia", cmv = "cmv",
                  neutro = "prolonged_neutropenia", prior_IA = "prior_mould_disease")
BS_OUTCOME <- "eortc_prov_prob"

# Differences between the file and the published Table 2 that are known and unresolved
# (data minus published). Any other difference is a failure.
BS_KNOWN_COUNT_DIFFERENCES <- data.frame(
  predictor = c("high_risk_chemo", "prolonged_neutropenia"),
  group = c("non_cases", "non_cases"),
  diff = c(-1L, -1L)
)

read_bs_data <- function(path = "bs.xlsx") {
  if (!file.exists(path)) {
    stop("'", path, "' not found. It is gitignored patient-level data: copy it into the project root.",
         call. = FALSE)
  }
  d <- as.data.frame(readxl::read_excel(path))
  needed <- c("id", BS_OUTCOME, names(BS_VARIABLES))
  missing_cols <- setdiff(needed, names(d))
  if (length(missing_cols)) stop("Missing columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)

  coded <- d[, c(BS_OUTCOME, names(BS_VARIABLES))]
  if (anyNA(coded)) stop("Missing values in the outcome or predictors.", call. = FALSE)
  if (!all(unlist(coded) %in% 0:1)) stop("Outcome and predictors must be coded 0/1.", call. = FALSE)

  X <- as.matrix(d[, names(BS_VARIABLES)])
  colnames(X) <- unname(BS_VARIABLES)
  list(X = X[, PREDICTORS], y = d[[BS_OUTCOME]], id = d$id)
}

# Rank-based AUC (ties count one half)
auc_rank <- function(score, y) {
  r <- rank(score)
  n1 <- sum(y == 1); n0 <- sum(y == 0)
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

# Counts of each predictor in cases and non-cases, against the published counts
check_bs_counts <- function(bs, pars) {
  cases <- bs$y == 1
  out <- data.frame(
    predictor = PREDICTORS,
    cases_data = colSums(bs$X[cases, , drop = FALSE]),
    cases_published = as.numeric(pars$count_case[PREDICTORS]),
    non_cases_data = colSums(bs$X[!cases, , drop = FALSE]),
    non_cases_published = as.numeric(pars$count_non[PREDICTORS]),
    row.names = NULL
  )
  out$diff_cases <- out$cases_data - out$cases_published
  out$diff_non_cases <- out$non_cases_data - out$non_cases_published
  out
}

# Do the differences equal the documented ones, no more and no less?
bs_counts_as_expected <- function(counts) {
  observed <- rbind(
    data.frame(predictor = counts$predictor, group = "cases", diff = counts$diff_cases),
    data.frame(predictor = counts$predictor, group = "non_cases", diff = counts$diff_non_cases)
  )
  observed <- observed[observed$diff != 0, ]
  key <- function(d) paste(d$predictor, d$group, d$diff)
  setequal(key(observed), key(BS_KNOWN_COUNT_DIFFERENCES)) && nrow(observed) == nrow(BS_KNOWN_COUNT_DIFFERENCES)
}

# Refit the logistic model, for verification only. The model used in the study stays fixed at
# the published coefficients.
check_bs_refit <- function(bs, pars) {
  fit <- stats::glm(bs$y ~ bs$X, family = stats::binomial)
  co <- summary(fit)$coefficients
  data.frame(
    term = c("(Intercept)", PREDICTORS),
    published = c(pars$intercept, pars$coef[PREDICTORS]),
    refit = unname(co[, 1]),
    se = unname(co[, 2]),
    row.names = NULL
  ) |>
    transform(diff = refit - published)
}

# The published equation applied to the file
bs_apparent_performance <- function(bs, pars) {
  lp <- as.numeric(pars$intercept + bs$X %*% pars$coef[PREDICTORS])
  list(n = length(bs$y), cases = sum(bs$y), mean_predicted_risk = mean(stats::plogis(lp)),
       observed_prevalence = mean(bs$y), apparent_auc = auc_rank(lp, bs$y))
}

# Empirical distribution of the 2^7 clinical profiles in the analysis admissions. This replaces the
# fitted copula as the baseline cohort (addendum D8); the copula stays as a sensitivity analysis.
empirical_profile_counts <- function(bs) tabulate(profile_index(bs$X), 2L^length(PREDICTORS))

# Tetrachoric correlation of two binary variables: the latent normal correlation that reproduces the
# observed 2x2 table.
tetrachoric <- function(x, y) {
  p11 <- mean(x * y); ta <- stats::qnorm(1 - mean(x)); tb <- stats::qnorm(1 - mean(y))
  f <- function(r) mvtnorm::pmvnorm(lower = c(ta, tb), upper = c(Inf, Inf), corr = matrix(c(1, r, r, 1), 2))[1] - p11
  stats::uniroot(f, c(-0.95, 0.99), tol = 1e-8)$root
}
