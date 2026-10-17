# Track A2, identifiability: sample the near-optimal region of the structured correlation
# matrix with a random-walk Metropolis on exp(-chi2 / 2). The chi2 uses the published
# sampling errors and ignores admission clustering, so this is an ensemble of matrices
# that fit the targets about equally well, not a formal posterior. Run from project root.

source("R/baseline.R")
source("R/copula.R")

cores <- max(1L, min(15L, parallel::detectCores() - 1L))
pars    <- load_baseline_params()
targets <- load_validation_targets()
profiles <- all_profiles()

set.seed(7461)
Z <- matrix(rnorm(3e5 * length(PREDICTORS)), 3e5)
obj <- make_objective(Z, pars, targets, profiles)
theta0 <- setNames(read.csv("output/fitted_theta.csv")$value, THETA_NAMES)

n_chains <- 6; n_iter <- 3000; burn <- 1000; thin <- 5

run_chain <- function(i) {
  set.seed(1000 + i * 17)
  start <- pmin(pmax(theta0 + rnorm(length(theta0), 0, 0.05), THETA_LOWER), THETA_UPPER)
  step <- 0.04
  cur <- start; cur_loss <- obj(cur); acc <- 0
  draws <- matrix(NA_real_, n_iter, length(THETA_NAMES), dimnames = list(NULL, THETA_NAMES))
  loss <- numeric(n_iter)
  for (k in seq_len(n_iter)) {
    prop <- cur + rnorm(length(cur), 0, step)
    if (all(prop >= THETA_LOWER & prop <= THETA_UPPER)) {
      pl <- obj(prop)
      if (log(runif(1)) < (cur_loss - pl) / 2) { cur <- prop; cur_loss <- pl; acc <- acc + 1 }
    }
    draws[k, ] <- cur; loss[k] <- cur_loss
  }
  keep <- seq(burn + 1, n_iter, by = thin)
  list(draws = data.frame(chain = i, iter = keep, draws[keep, ], loss = loss[keep]), accept = acc / n_iter)
}

chains <- run_parallel(seq_len(n_chains), run_chain, cores = n_chains)
ens <- do.call(rbind, lapply(chains, `[[`, "draws"))
write.csv(ens, "output/correlation_ensemble_raw.csv", row.names = FALSE)
write.csv(data.frame(chain = seq_len(n_chains), accept = sapply(chains, `[[`, "accept")),
          "output/correlation_ensemble_acceptance.csv", row.names = FALSE)

# Thin to 150 members for the decision-curve evaluation (script 03)
set.seed(2208)
members <- ens[sample(nrow(ens), 150), ]
members$member <- seq_len(nrow(members))
write.csv(members, "output/correlation_ensemble_members.csv", row.names = FALSE)
