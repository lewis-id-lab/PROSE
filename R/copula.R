# Gaussian-copula joint distribution of the seven binary clinical predictors.
# The joint is not published (see 02_current_status.md section 6), so the latent
# correlation matrix is a model assumption, fitted to the published targets.

# Clinical grouping used to structure the matrix. Group A is the marrow-suppression
# block; steroids and CMV are T-cell-immunosuppression variables; prior mould disease
# is treated as a history variable that is only weakly tied to current status.
GROUP_A <- c("uncontrolled_malignancy", "high_risk_chemo", "severe_lymphopenia", "prolonged_neutropenia")

# Version 1: eight parameters. Version 2 adds a separate chemotherapy-lymphopenia term so
# that malignancy's links to the rest of the block are not tied to a pair that needs to be
# strongly correlated for the matrix to stay positive definite.
THETA_NAMES_V1 <- c("rho_A", "d_malig_chemo", "d_chemo_neut", "d_neut_lymph",
                    "rho_A_steroid", "rho_A_cmv", "rho_steroid_cmv", "rho_prior")
THETA_NAMES <- c("rho_A", "d_malig_chemo", "d_chemo_neut", "d_neut_lymph", "d_chemo_lymph",
                 "rho_A_steroid", "rho_A_cmv", "rho_steroid_cmv", "rho_prior")

# Bounds encode the clinical constraints: all dependence is non-negative, and the named
# within-block pairs are at least as correlated as the baseline within-block pairs.
THETA_LOWER_V1 <- setNames(rep(0, 8), THETA_NAMES_V1)
THETA_UPPER_V1 <- c(rho_A = 0.8, d_malig_chemo = 0.5, d_chemo_neut = 0.5, d_neut_lymph = 0.5,
                    rho_A_steroid = 0.8, rho_A_cmv = 0.8, rho_steroid_cmv = 0.8, rho_prior = 0.5)

# Version 2 default bounds. The three block extras may reach 0.7, so a latent correlation
# of up to about 0.95 is permitted and positive definiteness, not the box, is what limits them.
THETA_LOWER <- setNames(rep(0, 9), THETA_NAMES)
THETA_UPPER <- c(rho_A = 0.8, d_malig_chemo = 0.5, d_chemo_neut = 0.7, d_neut_lymph = 0.7,
                 d_chemo_lymph = 0.7, rho_A_steroid = 0.8, rho_A_cmv = 0.8,
                 rho_steroid_cmv = 0.8, rho_prior = 0.5)
# Version 2 with the version 1 caps, to separate the effect of the new term from the caps
THETA_UPPER_CAPPED <- c(THETA_UPPER_V1[1:4], d_chemo_lymph = 0.5, THETA_UPPER_V1[5:8])

build_corr <- function(theta) {
  stopifnot(all(THETA_NAMES_V1 %in% names(theta)), all(names(theta) %in% THETA_NAMES))
  d_cl <- if ("d_chemo_lymph" %in% names(theta)) theta[["d_chemo_lymph"]] else 0
  R <- diag(length(PREDICTORS))
  dimnames(R) <- list(PREDICTORS, PREDICTORS)
  set_r <- function(a, b, r) { R[a, b] <<- r; R[b, a] <<- r }

  for (i in seq_along(GROUP_A)) for (j in seq_len(i - 1)) set_r(GROUP_A[i], GROUP_A[j], theta[["rho_A"]])
  set_r("uncontrolled_malignancy", "high_risk_chemo", theta[["rho_A"]] + theta[["d_malig_chemo"]])
  set_r("high_risk_chemo", "prolonged_neutropenia", theta[["rho_A"]] + theta[["d_chemo_neut"]])
  set_r("prolonged_neutropenia", "severe_lymphopenia", theta[["rho_A"]] + theta[["d_neut_lymph"]])
  set_r("high_risk_chemo", "severe_lymphopenia", theta[["rho_A"]] + d_cl)
  for (a in GROUP_A) {
    set_r(a, "high_dose_steroids", theta[["rho_A_steroid"]])
    set_r(a, "cmv", theta[["rho_A_cmv"]])
  }
  set_r("high_dose_steroids", "cmv", theta[["rho_steroid_cmv"]])
  for (o in setdiff(PREDICTORS, "prior_mould_disease")) set_r("prior_mould_disease", o, theta[["rho_prior"]])
  R
}

exchangeable_corr <- function(rho, d = length(PREDICTORS)) {
  R <- matrix(rho, d, d)
  diag(R) <- 1
  dimnames(R) <- list(PREDICTORS, PREDICTORS)
  R
}

# Continuous repair of a matrix that is not safely positive definite. Shifting the
# diagonal and rescaling keeps unit variances and changes the matrix continuously, so the
# fitting objective has no cliff at the feasibility boundary; a steep penalty keeps the
# optimiser from settling in the repaired region.
pd_adjust <- function(R, eps = 1e-3) {
  lmin <- min(eigen(R, symmetric = TRUE, only.values = TRUE)$values)
  if (lmin >= eps) return(list(R = R, penalty = 0, shift = 0))
  s <- eps - lmin
  Radj <- (R + s * diag(nrow(R))) / (1 + s)
  dimnames(Radj) <- dimnames(R)
  list(R = Radj, penalty = 1e6 * s^2, shift = s)
}

# Profile distribution under the copula, using fixed standard-normal draws Z0
# (common random numbers, so the fitting objective is smooth in R).
profile_probs_crn <- function(R, Z0, marginal) {
  L <- tryCatch(chol(R), error = function(e) NULL)
  if (is.null(L)) return(NULL)
  Z <- Z0 %*% L
  thr <- qnorm(1 - marginal[PREDICTORS])
  idx <- rep(1L, nrow(Z))
  for (j in seq_along(PREDICTORS)) idx <- idx + (Z[, j] > thr[j]) * 2L^(j - 1L)
  tabulate(idx, 2L^length(PREDICTORS)) / nrow(Z)
}

# Same distribution from fresh draws in chunks, for large final evaluations
profile_probs_mc <- function(R, marginal, n, seed, chunk = 1e6) {
  set.seed(seed)
  counts <- numeric(2L^length(PREDICTORS))
  done <- 0
  while (done < n) {
    m <- min(chunk, n - done)
    Z0 <- matrix(rnorm(m * length(PREDICTORS)), m)
    counts <- counts + profile_probs_crn(R, Z0, marginal) * m
    done <- done + m
  }
  counts / n
}

# Profile distribution by numerical integration of the 7-dimensional normal, with no
# patient-level sampling. mvtnorm uses a randomised quasi-Monte Carlo rule internally;
# at abseps = 1e-7 the profile probabilities repeat to about 1e-6 across seeds.
profile_probs_exact <- function(R, marginal, abseps = 1e-7, profiles = all_profiles(), seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  thr <- qnorm(1 - marginal[colnames(profiles)])
  algo <- mvtnorm::GenzBretz(abseps = abseps, releps = 0, maxpts = 2e6)
  q <- apply(profiles, 1, function(x) {
    as.numeric(mvtnorm::pmvnorm(lower = ifelse(x == 1, thr, -Inf), upper = ifelse(x == 1, Inf, thr),
                                corr = R, algorithm = algo))
  })
  q / sum(q)
}

# The 16 published targets, as z-scores against their published sampling error.
# Case and non-case prevalences use binomial SEs on the published group sizes.
target_table <- function(pars, targets) {
  c1 <- pars$prev_case; c0 <- pars$prev_non
  data.frame(
    target = c(paste0("case_prev:", PREDICTORS), paste0("noncase_prev:", PREDICTORS),
               "prevalence", "aroc"),
    published = c(c1, c0, targets$prevalence, targets$aroc),
    se = c(sqrt(c1 * (1 - c1) / pars$n_case), sqrt(c0 * (1 - c0) / pars$n_non),
           sqrt(targets$prevalence * (1 - targets$prevalence) / (pars$n_case + pars$n_non)),
           targets$aroc_se),
    row.names = NULL
  )
}

metrics_vector <- function(m) c(m$case_prev, m$noncase_prev, prevalence = m$prevalence, aroc = m$auc)

target_loss <- function(m, tt) sum(((metrics_vector(m) - tt$published) / tt$se)^2)

make_objective <- function(Z0, pars, targets, profiles = all_profiles(), theta_names = THETA_NAMES) {
  tt <- target_table(pars, targets)
  function(theta) {
    names(theta) <- theta_names
    pd <- pd_adjust(build_corr(theta))
    q <- profile_probs_crn(pd$R, Z0, pars$marginal)
    target_loss(profile_metrics(q, pars, profiles), tt) + pd$penalty
  }
}

# Parallel lapply on fresh R processes (PSOCK), because forking is unavailable in some IDEs.
# Workers start in the project root and source the model code.
run_parallel <- function(X, FUN, cores = 1L, root = getwd()) {
  if (cores <= 1L || length(X) <= 1L) return(lapply(X, FUN))
  cl <- parallel::makeCluster(min(cores, length(X)))
  on.exit(parallel::stopCluster(cl))
  parallel::clusterCall(cl, function(root) {
    setwd(root)
    for (f in c("R/baseline.R", "R/copula.R", "R/decision_curve.R")) if (file.exists(f)) source(f)
    invisible(NULL)
  }, root)
  parallel::parLapply(cl, X, FUN)
}

# Multi-start bounded fit, parallel over starts. L-BFGS-B is adequate because common random
# numbers make the objective smooth at the scale of the finite-difference step.
fit_structured_corr <- function(Z0, pars, targets, starts, profiles = all_profiles(),
                                lower = THETA_LOWER, upper = THETA_UPPER, cores = 1L) {
  nms <- names(lower)
  obj <- make_objective(Z0, pars, targets, profiles, theta_names = nms)
  one <- function(i) {
    f <- optim(unlist(starts[i, nms]), obj, method = "L-BFGS-B", lower = lower, upper = upper,
               control = list(ndeps = rep(2e-3, length(nms)), maxit = 100))
    c(setNames(f$par, nms), loss = f$value, convergence = f$convergence)
  }
  fits <- run_parallel(seq_len(nrow(starts)), one, cores = cores)
  out <- as.data.frame(do.call(rbind, fits))
  out[order(out$loss), ]
}
