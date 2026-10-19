# Extensions of the genetic layer (addendum 12). Run from the project root, after scripts 14, 15 and 16.
#
#   A. Correlation ensemble. The 150 near-optimal copula matrices (script 02) give 150 profile distributions of the seven clinical
#      predictors; each simulated world draws one member, so the spread includes the unpublished correlation structure. The
#      MCID is the one of the baseline in that member's own distribution. This is a sensitivity analysis (D8, D9): the primary
#      cohort stays the empirical one.
#   B. Per-term attribution. Shapley decomposition of the gain into the four terms, for the ceiling and the realised score, on the
#      empirical cohort, with each term alone and the loss from leaving it out.
#   C. Sampling variability of a study of this size: 4,078-admission datasets built by resampling patients with their admissions,
#      assigning each patient a genotype and drawing outcomes from the true risk, then measuring the observed net-benefit gain of
#      fixed scores. Two kinds of run: effects held fixed (pure sampling variability) and effects drawn from the posteriors
#      (parameter and sampling variability together).
#
# Writes   output/correlation_ensemble_profile_probs.csv   the 150 profile distributions (cached)
#          output/genetic_ensemble_draws.csv, genetic_ensemble_summary.csv
#          output/genetic_attribution_draws.csv, genetic_attribution_summary.csv
#          output/genetic_sampling_draws.csv, genetic_sampling_summary.csv

source("R/baseline.R")
source("R/bs_data.R")
source("R/copula.R")
source("R/decision_curve.R")
source("R/boscore_data.R")
source("R/refit_baseline.R")
source("R/evidence_mapping.R")
source("R/genetic_layer.R")

at <- c(0.05, 0.20)
cores <- max(1L, min(15L, parallel::detectCores() - 1L))
n_chunks <- 30L
chunk_ids <- function(n, k = n_chunks) split(seq_len(n), rep_len(seq_len(k), n))

base <- select_cohort(load_base_cohort(), include_twins = FALSE)
pars_refit <- load_refit_params()
lp_profile <- profile_logit(pars_refit)
q_emp <- tabulate(profile_index(base$X), 2L^length(PREDICTORS)) / nrow(base$X)
mcid_tab <- read.csv("output/dca_baseline_refit_summary.csv")
mcid <- setNames(mcid_tab$mcid_b[match(at, mcid_tab$threshold)], as.character(at))

p <- genetic_term_frequencies()
grid <- dosage_grid(p)
cells <- genetic_cells(q_emp, lp_profile, grid)
draws <- load_term_draws()
beta_hat <- pooled_estimates(draws)
beta_hi <- vapply(draws, function(d) unname(stats::quantile(d$mu, 0.95)), numeric(1))
effects_file <- read.csv("output/genetic_layer_effects.csv")
term_names <- names(GENETIC_TERMS)

# Fresh worker processes (forking is unavailable in Positron) that source the model code
cl <- parallel::makeCluster(cores)
parallel::clusterCall(cl, function(root) {
  setwd(root)
  for (f in c("R/baseline.R", "R/decision_curve.R", "R/genetic_layer.R")) source(f)
  invisible(NULL)
}, getwd())
pl <- function(chunks, FUN, ctx) do.call(rbind, parallel::parLapply(cl, chunks, FUN, ctx))

# ---- A. Correlation ensemble ------------------------------------------------------------------------------------------------

members <- read.csv("output/correlation_ensemble_members.csv")
cache_path <- "output/correlation_ensemble_profile_probs.csv"
if (file.exists(cache_path)) {
  qm <- read.csv(cache_path)
  stopifnot(identical(qm$member, members$member))
} else {
  marginal <- load_baseline_params()$marginal
  thetas <- lapply(seq_len(nrow(members)), function(i) unlist(members[i, THETA_NAMES]))
  make_q <- function(marginal) function(th) profile_probs_exact(pd_adjust(build_corr(th))$R, marginal, abseps = 1e-6, seed = 4401)
  qs <- run_parallel(thetas, make_q(marginal), cores = cores)
  qm <- cbind(member = members$member, as.data.frame(do.call(rbind, qs)))
  names(qm)[-1] <- paste0("q", seq_len(2L^length(PREDICTORS)))
  write.csv(qm, cache_path, row.names = FALSE)
}
q_members <- as.matrix(qm[, -1]); stopifnot(all(abs(rowSums(q_members) - 1) < 1e-6))
tv <- apply(q_members, 1, function(q) 0.5 * sum(abs(q - q_emp)))
message("ensemble members: ", nrow(q_members), "; total-variation distance to the empirical profile distribution: median ",
        sprintf("%.3f", median(tv)), " (range ", sprintf("%.3f", min(tv)), " to ", sprintf("%.3f", max(tv)), ")")

n_ens <- 3000L
member_of_draw <- rep_len(seq_len(nrow(q_members)), n_ens)
ens_effects <- list(mean = draw_term_effects(draws, n_ens, "mean", seed = 4180L),
                    predictive = draw_term_effects(draws, n_ens, "predictive", seed = 7352L))
ens <- do.call(rbind, lapply(names(ens_effects), function(h) {
  ctx <- list(q = q_members, member = member_of_draw, lp = lp_profile, grid = grid, eff = ens_effects[[h]], bh = beta_hat, at = at)
  cbind(truth = h, pl(chunk_ids(n_ens), function(ch, ctx) {
    run_genetic_layer_ensemble(ctx$q, ctx$member[ch], ctx$lp, ctx$grid, ctx$eff[ch, , drop = FALSE], ctx$bh, ctx$at, draw_id = ch)
  }, ctx))
}))
ens$extra_cases_per_1000 <- 1000 * ens$delta_tp

summ_gain <- function(d) data.frame(
  gain = median(d$delta_nb), lo = quantile(d$delta_nb, 0.05, names = FALSE), hi = quantile(d$delta_nb, 0.95, names = FALSE),
  extra_cases_per_1000 = 1000 * median(d$delta_tp), p_meets_mcid = mean(d$meets_mcid), mcid = median(d$mcid), flagged = NA_real_)
by_key <- function(d, keys, f) {
  k <- unique(d[keys]); do.call(rbind, lapply(seq_len(nrow(k)), function(i) {
    s <- Reduce(`&`, lapply(keys, function(n) d[[n]] == k[[n]][i])); cbind(k[i, , drop = FALSE], f(d[s, ]))
  }))
}
ens_sum <- by_key(ens[ens$scenario != "baseline", ], c("truth", "scenario", "threshold"), summ_gain)[, -10]
# How much of the spread comes from the matrices: spread of each member's mean gain across members
member_gain <- aggregate(delta_nb ~ truth + scenario + threshold + member, data = ens[ens$scenario != "baseline", ], FUN = mean)
ens_sum$between_member_sd <- by_key(member_gain, c("truth", "scenario", "threshold"), function(d) data.frame(sd = stats::sd(d$delta_nb)))$sd
emp <- read.csv("output/genetic_layer_summary.csv")
ens_sum$empirical_gain <- emp$delta_nb[match(paste(ens_sum$truth, ens_sum$scenario, ens_sum$threshold), paste(emp$truth, emp$scenario, emp$threshold))]
ens_sum$empirical_p_meets_mcid <- emp$p_meets_mcid[match(paste(ens_sum$truth, ens_sum$scenario, ens_sum$threshold), paste(emp$truth, emp$scenario, emp$threshold))]

# ---- B. Per-term attribution --------------------------------------------------------------------------------------------------

n_attr <- 1000L
attr_all <- do.call(rbind, lapply(c("mean", "predictive"), function(h) {
  eff <- as.matrix(effects_file[effects_file$truth == h & effects_file$draw <= n_attr, term_names])
  ctx <- list(cells = cells, eff = eff, bh = beta_hat, at = at, terms = term_names)
  pl(chunk_ids(n_attr), function(ch, ctx) {
    do.call(rbind, lapply(ch, function(i) do.call(rbind, lapply(c("ceiling", "realised"), function(sc) {
      s <- shapley_gain(ctx$cells, ctx$eff[i, ], if (sc == "ceiling") ctx$eff[i, ] else ctx$bh, ctx$at)
      cbind(draw = i, score = sc, threshold = ctx$at, total = s$total,
            setNames(as.data.frame(s$phi), paste0("phi_", ctx$terms)),
            setNames(as.data.frame(s$alone), paste0("alone_", ctx$terms)),
            setNames(as.data.frame(s$leave_out), paste0("loo_", ctx$terms)))
    }))))
  }, ctx) |> (\(d) cbind(truth = h, d))()
}))
attr_sum <- by_key(attr_all, c("truth", "score", "threshold"), function(d) {
  out <- data.frame(total = median(d$total))
  for (t in term_names) {
    out[[paste0("phi_", t)]] <- median(d[[paste0("phi_", t)]])
    out[[paste0("share_", t)]] <- sum(d[[paste0("phi_", t)]]) / sum(d$total)
    out[[paste0("alone_", t)]] <- median(d[[paste0("alone_", t)]])
    out[[paste0("loo_", t)]] <- median(d[[paste0("loo_", t)]])
  }
  out$p_phi_negative_any <- mean(rowSums(d[paste0("phi_", term_names)] < -1e-9) > 0)
  out
})

# ---- C. Sampling variability of a study of this size ----------------------------------------------------------------------------

cg <- cohort_groups(base)
n_sim <- 2000L
runs <- list(
  fixed_medians = list(truth = matrix(beta_hat, n_sim, 4, byrow = TRUE, dimnames = list(NULL, term_names)), seed = 5290L),
  fixed_upper = list(truth = matrix(beta_hi, n_sim, 4, byrow = TRUE, dimnames = list(NULL, term_names)), seed = 6841L),
  draws_mean = list(truth = as.matrix(effects_file[effects_file$truth == "mean" & effects_file$draw <= n_sim, term_names]), seed = 2736L),
  draws_predictive = list(truth = as.matrix(effects_file[effects_file$truth == "predictive" & effects_file$draw <= n_sim, term_names]), seed = 9417L))
samp <- do.call(rbind, lapply(names(runs), function(r) {
  ctx <- list(cg = cg, lp = lp_profile, p = p, cells = cells, eff = runs[[r]]$truth, bh = beta_hat, at = at, mcid = mcid, seed = runs[[r]]$seed, terms = term_names)
  cbind(analysis = r, pl(chunk_ids(n_sim), function(ch, ctx) {
    set.seed(ctx$seed + ch[1])
    do.call(rbind, lapply(ch, function(i) {
      sh <- genetic_scores(ctx$cells, ctx$eff[i, ], ctx$bh)
      s <- simulate_study(ctx$cg, ctx$lp, ctx$p, sh, ctx$eff[i, ], ctx$bh, ctx$at, ctx$mcid)
      cbind(draw = i, s[s$scenario != "baseline", ])
    }))
  }, ctx))
}))
exact_ref <- do.call(rbind, lapply(c("fixed_medians", "fixed_upper"), function(r) {
  b <- runs[[r]]$truth[1, ]
  x <- genetic_nb(cells, b, beta_hat, at)
  g <- x[x$scenario != "baseline", ]
  g$delta_nb <- g$nb - x$nb[x$scenario == "baseline"][match(g$threshold, x$threshold[x$scenario == "baseline"])]
  data.frame(analysis = r, scenario = g$scenario, threshold = g$threshold, exact_gain = g$delta_nb)
}))
samp_sum <- by_key(samp, c("analysis", "scenario", "threshold"), function(d) data.frame(
  mean_gain = mean(d$delta_nb), median_gain = median(d$delta_nb), sd_gain = stats::sd(d$delta_nb),
  lo = quantile(d$delta_nb, 0.05, names = FALSE), hi = quantile(d$delta_nb, 0.95, names = FALSE),
  p_gain_positive = mean(d$delta_nb > 0), p_reaches_mcid = mean(d$delta_nb >= d$mcid), p_reaches_sample_mcid = mean(d$delta_nb >= d$mcid_sample),
  mcid = d$mcid[1], median_events = median(d$events)))
samp_sum$exact_gain <- exact_ref$exact_gain[match(paste(samp_sum$analysis, samp_sum$scenario, samp_sum$threshold), paste(exact_ref$analysis, exact_ref$scenario, exact_ref$threshold))]

parallel::stopCluster(cl)

# ---- Write and report -----------------------------------------------------------------------------------------------------------

dir.create("output", showWarnings = FALSE)
write.csv(ens[ens$scenario != "baseline", c("truth", "draw", "member", "scenario", "threshold", "nb", "delta_nb", "delta_tp", "mcid", "meets_mcid", "auc")], "output/genetic_ensemble_draws.csv", row.names = FALSE)
write.csv(ens_sum, "output/genetic_ensemble_summary.csv", row.names = FALSE)
write.csv(attr_all, "output/genetic_attribution_draws.csv", row.names = FALSE)
write.csv(attr_sum, "output/genetic_attribution_summary.csv", row.names = FALSE)
write.csv(samp[, c("analysis", "draw", "scenario", "threshold", "tp", "fp", "nb", "delta_nb", "delta_tp", "mcid_sample", "mcid", "events")], "output/genetic_sampling_draws.csv", row.names = FALSE)
write.csv(samp_sum, "output/genetic_sampling_summary.csv", row.names = FALSE)

show <- function(x) print(cbind(x[sapply(x, is.character)], round(x[sapply(x, is.numeric)], 5)), row.names = FALSE)
message("\nA. Correlation ensemble (gain over the baseline, median and 90% range; MCID in each member's own distribution):")
show(ens_sum)
message("\nB. Shapley attribution (median contribution to the gain; share = total contribution over total gain across draws):")
show(attr_sum[, c("truth", "score", "threshold", "total", paste0("phi_", term_names), paste0("share_", term_names))])
message("\nC. Sampling variability of a 4,078-admission study (observed gain of fixed scores):")
show(samp_sum)
