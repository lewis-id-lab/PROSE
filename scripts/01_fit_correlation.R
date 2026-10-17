# Track A2: fit the structured latent correlation matrix of the seven clinical predictors
# to the published targets (D6: aROC target 0.85). Run from the project root.
#
# Three configurations are compared:
#   v1         eight parameters, block extras capped at 0.5 (the first fit)
#   v2_capped  adds a separate chemotherapy-lymphopenia term, same caps
#   v2         adds that term and lets the block extras reach 0.7
# FINAL names the configuration written to the unsuffixed output files.

source("R/baseline.R")
source("R/copula.R")

FINAL <- "v2"
cores <- max(1L, min(15L, parallel::detectCores() - 1L))

pars    <- load_baseline_params()
targets <- load_validation_targets()
profiles <- all_profiles()
tt <- target_table(pars, targets)
dir.create("output", showWarnings = FALSE)

configs <- list(
  v1        = list(lower = THETA_LOWER_V1, upper = THETA_UPPER_V1),
  v2_capped = list(lower = THETA_LOWER,    upper = THETA_UPPER_CAPPED),
  v2        = list(lower = THETA_LOWER,    upper = THETA_UPPER)
)

# 1. Multi-start search on fixed draws (common random numbers)
set.seed(3814)
Z_search <- matrix(rnorm(2e5 * length(PREDICTORS)), 2e5)
n_starts <- 32

make_starts <- function(cfg, seed) {
  set.seed(seed)
  nms <- names(cfg$lower)
  hi <- pmin(cfg$upper, 0.5)
  s <- as.data.frame(t(replicate(n_starts, runif(length(nms), cfg$lower, hi))))
  names(s) <- nms
  s
}

search <- Map(function(cfg, seed) {
  fit_structured_corr(Z_search, pars, targets, make_starts(cfg, seed), profiles,
                      lower = cfg$lower, upper = cfg$upper, cores = cores)
}, configs, c(9275, 9276, 9277))

multistart <- dplyr::bind_rows(lapply(names(search), function(n) {
  cbind(config = n, search[[n]])
}))
write.csv(multistart, "output/fit_multistart.csv", row.names = FALSE)

# 2. Polish each best start on a larger, independent set of draws (parallel over configurations)
set.seed(7461)
Z_polish <- matrix(rnorm(5e5 * length(PREDICTORS)), 5e5)
polish_one <- function(n) {
  cfg <- configs[[n]]
  start <- search[[n]][1, names(cfg$lower), drop = FALSE]
  fit_structured_corr(Z_polish, pars, targets, start, profiles, lower = cfg$lower, upper = cfg$upper)[1, ]
}
polished <- run_parallel(names(configs), polish_one, cores = length(configs))
names(polished) <- names(configs)
theta_of <- function(n) unlist(polished[[n]][names(configs[[n]]$lower)])

# 3. Final evaluation by numerical integration, with no patient-level sampling
adj <- lapply(names(configs), function(n) pd_adjust(build_corr(theta_of(n))))
names(adj) <- names(configs)

mats <- c(list(exchangeable_0.30 = exchangeable_corr(0.30)), lapply(adj, `[[`, "R"))
exact_one <- function(R) profile_probs_exact(R, pars$marginal, seed = 6127)
q_all <- run_parallel(mats, exact_one, cores = length(mats))
q_all$independent <- independent_profile_probs(pars$marginal, profiles)
m_all <- lapply(q_all, profile_metrics, pars = pars, profiles = profiles)

comparison <- do.call(rbind, lapply(names(configs), function(n) {
  th <- theta_of(n); cfg <- configs[[n]]
  z <- (metrics_vector(m_all[[n]]) - tt$published) / tt$se
  names(z) <- tt$target
  at_bound <- names(th)[abs(th - cfg$lower) < 0.01 | abs(th - cfg$upper[names(th)]) < 0.01]
  data.frame(config = n, n_params = length(th), chi2 = target_loss(m_all[[n]], tt),
             prevalence = m_all[[n]]$prevalence, auc = m_all[[n]]$auc,
             z_malignancy_case = z[["case_prev:uncontrolled_malignancy"]],
             max_abs_z = max(abs(z)),
             pd_shift = adj[[n]]$shift,
             at_bound = paste(at_bound, collapse = "; "))
}))
write.csv(comparison, "output/fit_comparison.csv", row.names = FALSE)

all_theta <- do.call(rbind, lapply(names(configs), function(n) {
  data.frame(config = n, parameter = names(theta_of(n)), value = unname(theta_of(n)))
}))
write.csv(all_theta, "output/fitted_theta_all.csv", row.names = FALSE)

# 4. Outputs for the chosen model
write.csv(data.frame(parameter = names(theta_of(FINAL)), value = unname(theta_of(FINAL))),
          "output/fitted_theta.csv", row.names = FALSE)
write.csv(adj[[FINAL]]$R, "output/fitted_correlation_matrix.csv")

q_fit <- q_all[[FINAL]]
prof_tab <- data.frame(profile = seq_len(nrow(profiles)), profiles, q = q_fit,
                       risk = plogis(profile_logit(pars, profiles)))
write.csv(prof_tab, "output/profile_probabilities.csv", row.names = FALSE)

report <- cbind(
  tt,
  independent = metrics_vector(m_all$independent),
  exchangeable_0.30 = metrics_vector(m_all$exchangeable_0.30),
  fitted_v1 = metrics_vector(m_all$v1),
  fitted = metrics_vector(m_all[[FINAL]])
)
for (col in c("independent", "exchangeable_0.30", "fitted_v1", "fitted")) {
  report[[paste0("z_", col)]] <- (report[[col]] - report$published) / report$se
}
names(report)[names(report) == "z_exchangeable_0.30"] <- "z_exchangeable"
write.csv(report, "output/reproduction_report.csv", row.names = FALSE)

loss <- c(independent = target_loss(m_all$independent, tt),
          exchangeable_0.30 = target_loss(m_all$exchangeable_0.30, tt),
          fitted_v1 = target_loss(m_all$v1, tt),
          fitted = target_loss(m_all[[FINAL]], tt))
write.csv(data.frame(model = names(loss), chi2 = loss), "output/reproduction_loss.csv", row.names = FALSE)
