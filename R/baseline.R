# Fixed Stanzani 2019 baseline model (J Infect 78:484-490).
# Coefficients are never refitted. Everything downstream works at the level of the
# 2^7 = 128 predictor profiles, because the clinical block is seven binary variables.

PREDICTORS <- c(
  "uncontrolled_malignancy", "high_risk_chemo", "high_dose_steroids",
  "severe_lymphopenia", "cmv", "prolonged_neutropenia", "prior_mould_disease"
)

load_baseline_params <- function(path = "data/baseline_model_parameters.csv") {
  p <- read.csv(path, stringsAsFactors = FALSE)
  icpt <- p[p$term == "(Intercept)", ]
  p <- p[match(PREDICTORS, p$term), ]
  stopifnot(!anyNA(p$term))
  list(
    intercept   = icpt$coefficient,
    coef        = setNames(p$coefficient, p$term),
    odds_ratio  = setNames(p$odds_ratio, p$term),
    marginal    = setNames(p$prev_overall, p$term),
    prev_case   = setNames(p$prev_in_IMD_cases, p$term),
    prev_non    = setNames(p$prev_in_non_cases, p$term),
    count_case  = setNames(p$n_IMD, p$term),
    count_non   = setNames(p$n_non_IMD, p$term),
    n_case      = icpt$n_IMD,
    n_non       = icpt$n_non_IMD
  )
}

load_validation_targets <- function(path = "data/baseline_validation_targets.csv") {
  t <- read.csv(path, stringsAsFactors = FALSE)
  val <- function(name) as.numeric(t$value[t$target == name])
  list(
    prevalence = val("IMD prevalence"),
    aroc       = val("aROC derivation"),
    aroc_boot  = val("aROC bootstrap internal validation"),
    # Published 95% CI 0.84-0.86 on the derivation aROC
    aroc_se    = (0.86 - 0.84) / (2 * 1.96)
  )
}

# All 128 predictor profiles; row i has index i, bit j <-> predictor j.
all_profiles <- function() {
  m <- as.matrix(expand.grid(rep(list(0:1), length(PREDICTORS)), KEEP.OUT.ATTRS = FALSE))
  dimnames(m) <- list(NULL, PREDICTORS)
  m
}

profile_index <- function(X) {
  as.integer(X %*% 2^(seq_len(ncol(X)) - 1)) + 1L
}

# Linear predictor of the fixed equation for every profile
profile_logit <- function(pars, profiles = all_profiles()) {
  as.numeric(pars$intercept + profiles %*% pars$coef[colnames(profiles)])
}

# Exact profile distribution if the seven predictors were independent
independent_profile_probs <- function(marginal, profiles = all_profiles()) {
  p <- marginal[colnames(profiles)]
  apply(profiles, 1, function(x) prod(ifelse(x == 1, p, 1 - p)))
}

# Expected discrimination from grouped weights, ties counted as one half
expected_auc <- function(lp, w_case, w_non) {
  key <- round(lp, 10)
  w1 <- as.numeric(rowsum(w_case, key))
  w0 <- as.numeric(rowsum(w_non, key))
  # rowsum sorts by key, so the vectors are in ascending score order
  below0 <- cumsum(w0) - w0
  sum(w1 * (below0 + 0.5 * w0)) / (sum(w1) * sum(w0))
}

# Everything the validation harness needs, from a profile distribution q
profile_metrics <- function(q, pars, profiles = all_profiles()) {
  lp <- profile_logit(pars, profiles)
  risk <- plogis(lp)
  w_case <- q * risk
  w_non  <- q * (1 - risk)
  prev <- sum(w_case)
  list(
    prevalence  = prev,
    case_prev   = setNames(as.numeric(crossprod(profiles, w_case)) / prev, colnames(profiles)),
    noncase_prev = setNames(as.numeric(crossprod(profiles, w_non)) / (1 - prev), colnames(profiles)),
    marginal    = setNames(as.numeric(crossprod(profiles, q)), colnames(profiles)),
    auc         = expected_auc(lp, w_case, w_non)
  )
}

# Decision-curve net benefit per admission at a threshold t, from a profile distribution.
# Flag = predicted risk >= t. Returns TP and FP rates per admission as well as NB.
net_benefit <- function(q, pars, threshold, profiles = all_profiles(), risk = NULL) {
  if (is.null(risk)) risk <- plogis(profile_logit(pars, profiles))
  flag <- risk >= threshold
  tp <- sum(q * risk * flag)
  fp <- sum(q * (1 - risk) * flag)
  prev <- sum(q * risk)
  list(
    nb = tp - fp * threshold / (1 - threshold),
    nb_all = prev - (1 - prev) * threshold / (1 - threshold),
    tp = tp, fp = fp, flagged = sum(q * flag)
  )
}
