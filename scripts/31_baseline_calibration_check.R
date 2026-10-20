# Does adding genetics preserve the clinical baseline fairly? (addendum 24). Run from the project root, after scripts 16, 22, 23
# and 30.
#
# The genetic layer recentres one global intercept (D5), so within a clinical profile the true risk averaged over genotype is not
# the clinical risk, and the clinical score is miscalibrated profile by profile under the simulated truth. Part of the reported
# gain could be the correction of that induced miscalibration. Three comparisons (R/calibration_check.R):
#   current    global-intercept truth, clinical comparator expit(eta_X)              (reproduces the reported results)
#   marginal   the same truth, clinical comparator m(X) = E[p | X] under it           (the best clinical-only rule)
#   profile    profile-preserving truth (intercept per profile, conditional genetic effects unchanged), comparator expit(eta_X)
# for the four-term panel (primary), the five-term panel, the Aim 4 sweep and the polygenic score benchmark. The pre-specified MCID
# (+0.00458 at 5%, +0.00160 at 20%) is fixed throughout; under the profile-preserving truth the baseline's true positives, and so
# the MCID, are exactly those of the refit baseline.
#
# Writes   output/calibration_check_layer.csv         four- and five-term panels: gain, MCID test and calibration by comparison
#          output/calibration_check_artefact.csv      per panel, reading and threshold: NB(m) - NB(expit(eta)) over draws
#          output/calibration_check_aim4.csv          Aim 4 sweep: share reaching the MCID and verdicts changed, by comparison
#          output/calibration_check_aim4_sweep.csv    Aim 4 sweep, per architecture, comparison and threshold: gain and MCID test
#          output/calibration_check_aim4_frontier.csv fewest loci reaching the MCID, by comparison
#          output/calibration_check_prs.csv           polygenic score: OR per SD needed (ceiling), accuracy and cases needed

source("R/baseline.R"); source("R/bs_data.R"); source("R/decision_curve.R"); source("R/boscore_data.R"); source("R/refit_baseline.R")
source("R/evidence_mapping.R"); source("R/genetic_layer.R"); source("R/genetic_helpers.R"); source("R/aim4_sweep.R")
source("R/prs_benchmark.R"); source("R/calibration_check.R")

at <- c(0.05, 0.20)
base <- select_cohort(load_base_cohort(), include_twins = FALSE)
lp <- profile_logit(load_refit_params())
q <- tabulate(profile_index(base$X), 2L^length(PREDICTORS)) / nrow(base$X)
mcid_tab <- read.csv("output/dca_baseline_refit_summary.csv")
mcid <- setNames(mcid_tab$mcid_b[match(at, mcid_tab$threshold)], as.character(at))

cl <- start_workers(cores = 10, files = c("R/baseline.R", "R/decision_curve.R", "R/genetic_layer.R", "R/aim4_sweep.R", "R/prs_benchmark.R",
                              "R/calibration_check.R"))

# ---- 1 and 2. The four- and five-term panels over the draws of the true effects ----------------------------------------------------
run_draws <- function(grid, effects, beta_score) {
  jobs <- do.call(rbind, lapply(names(effects), function(h) data.frame(truth = h, draw = seq_len(nrow(effects[[h]])))))
  ctx <- list(q = q, lp = lp, grid = grid, effects = effects, beta_score = beta_score, at = at, jobs = jobs)
  parallel_rows(cl, split(seq_len(nrow(jobs)), rep_len(seq_len(60), nrow(jobs))), function(ix, ctx) {
    do.call(rbind, lapply(ix, function(i) {
      j <- ctx$jobs[i, ]
      cbind(truth = j$truth, draw = j$draw,
            calibration_gains(ctx$q, ctx$lp, ctx$grid, ctx$effects[[j$truth]][j$draw, ], ctx$beta_score, ctx$at, calibration = TRUE, detail = TRUE))
    }))
  }, ctx)
}

eff <- read.csv("output/genetic_layer_effects.csv")
terms4 <- c("PTX3", "CLEC7A", "TLR4", "CD209")
effects4 <- lapply(c(mean = "mean", predictive = "predictive"), function(h) as.matrix(eff[eff$truth == h, terms4]))
d4 <- run_draws(dosage_grid(genetic_term_frequencies()), effects4, pooled_estimates(load_term_draws()))

s5 <- clec7a_second_setup()
d5 <- run_draws(s5$grid, s5$effects, s5$beta_hat)

summarise_panel <- function(d, panel) {
  long <- rbind(cbind(d, score = "ceiling", gain = d$gain_oracle, extra_tp = d$extra_tp_oracle, flagged = d$flagged_oracle, auc = d$auc_oracle),
                cbind(d, score = "realised", gain = d$gain_score, extra_tp = d$extra_tp_score, flagged = d$flagged_score, auc = d$auc_score))
  out <- by_key(long, c("truth", "comparison", "threshold", "score"), function(x) data.frame(
    gain = median(x$gain), lo = quantile(x$gain, 0.05, names = FALSE), hi = quantile(x$gain, 0.95, names = FALSE),
    extra_cases_per_1000 = 1000 * median(x$extra_tp), flagged = median(x$flagged), flagged_clinical = median(x$flagged_clinical),
    auc = median(x$auc), auc_clinical = median(x$auc_clinical),
    p_meets_mcid = mean(x$gain >= mcid[[as.character(x$threshold[1])]]),
    nb_clinical = median(x$nb_clinical), flag_disagree = median(x$flag_disagree), flag_disagree_max = max(x$flag_disagree),
    clin_slope = median(x$clin_slope), clin_citl = median(x$clin_citl),
    clin_mean_abs_err = median(x$clin_mean_abs_err), clin_max_abs_err = median(x$clin_max_abs_err),
    clin_max_abs_err_p95 = quantile(x$clin_max_abs_err, 0.95, names = FALSE),
    score_slope = if (x$score[1] == "realised") median(x$score_slope) else 1))
  cbind(panel = panel, out)
}
layer <- rbind(summarise_panel(d4, "four terms"), summarise_panel(d5, "five terms"))

# The current comparison must reproduce the reported layers
ref4 <- read.csv("output/genetic_layer_summary.csv")
ref5 <- subset(read.csv("output/clec7a_second_layer_summary.csv"), second_parameter == "estimated")
chk <- function(ref, panel_name, gain_col) {
  cur <- layer[layer$panel == panel_name & layer$comparison == "current", c("truth", "score", "threshold", "gain")]
  names(cur)[4] <- "gain_chk"
  m <- merge(cur, ref, by.x = c("truth", "score", "threshold"), by.y = c("truth", "scenario", "threshold"))
  stopifnot(nrow(m) == nrow(cur), nrow(m) > 0)
  max(abs(m$gain_chk - m[[gain_col]]))
}
stopifnot(chk(ref4, "four terms", "delta_nb") < 1e-9, chk(ref5, "five terms", "gain") < 1e-9)

# ... and so must the extra true cases, the flagged share and the aROC
chk_cols <- function(ref, panel_name) {
  cur <- layer[layer$panel == panel_name & layer$comparison == "current", c("truth", "score", "threshold", "extra_cases_per_1000", "flagged", "auc")]
  names(cur)[4:6] <- paste0(names(cur)[4:6], "_chk")
  m <- merge(cur, ref, by.x = c("truth", "score", "threshold"), by.y = c("truth", "scenario", "threshold"))
  stopifnot(nrow(m) == nrow(cur), nrow(m) > 0)
  max(abs(m$extra_cases_per_1000_chk - m$extra_cases_per_1000), abs(m$flagged_chk - m$flagged), abs(m$auc_chk - m$auc))
}
stopifnot(chk_cols(ref4, "four terms") < 1e-8, chk_cols(ref5, "five terms") < 1e-8)

# The artefact: how much better the best clinical-only rule does than expit(eta) under the global truth
artefact <- do.call(rbind, lapply(list(`four terms` = d4, `five terms` = d5), function(d) {
  cur <- d[d$comparison == "current", ]; mar <- d[d$comparison == "marginal", ]
  stopifnot(identical(cur$draw, mar$draw), identical(cur$threshold, mar$threshold))
  x <- data.frame(truth = cur$truth, threshold = cur$threshold, artefact = mar$nb_clinical - cur$nb_clinical,
                  share_of_gain = ifelse(cur$gain_oracle > 0, (mar$nb_clinical - cur$nb_clinical) / cur$gain_oracle, NA))
  by_key(x, c("truth", "threshold"), function(y) data.frame(
    draws_with_artefact = mean(y$artefact > 1e-12), artefact_median = median(y$artefact),
    artefact_p95 = quantile(y$artefact, 0.95, names = FALSE), artefact_max = max(y$artefact),
    share_of_ceiling_gain_median = median(y$share_of_gain, na.rm = TRUE),
    share_of_ceiling_gain_p95 = quantile(y$share_of_gain, 0.95, na.rm = TRUE, names = FALSE)))
}))
artefact <- cbind(panel = rep(c("four terms", "five terms"), each = nrow(artefact) / 2), artefact)

# ---- 3. The Aim 4 sweep --------------------------------------------------------------------------------------------------------------
g4 <- expand.grid(n_loci = 1:20, odds_ratio = round(seq(1.1, 3.0, by = 0.1), 2), raf = round(seq(0.05, 0.50, by = 0.05), 2))
ctx <- list(q = q, lp = lp, at = at, g = g4)
aim <- parallel_rows(cl, split(seq_len(nrow(g4)), rep_len(seq_len(60), nrow(g4))), function(ch, ctx) {
  do.call(rbind, lapply(ch, function(i) {
    pmf <- dosage_pmf(ctx$g$n_loci[i], ctx$g$raf[i])
    grid <- list(G = matrix(seq_along(pmf) - 1, ncol = 1, dimnames = list(NULL, "D")), prob = pmf)
    b <- c(D = log(ctx$g$odds_ratio[i]))
    cbind(ctx$g[i, ], calibration_gains(ctx$q, ctx$lp, grid, b, b, ctx$at, calibration = TRUE), row.names = NULL)
  }))
}, ctx)
aim$mcid <- unname(mcid[as.character(aim$threshold)])
aim$meets <- aim$gain_oracle >= aim$mcid
ref_aim <- read.csv("output/aim4_sweep.csv")
m <- merge(subset(aim, comparison == "current"), ref_aim, by = c("n_loci", "odds_ratio", "raf", "threshold"))
stopifnot(nrow(m) == nrow(ref_aim), max(abs(m$gain_oracle - m$gain)) < 1e-9)

cur_meets <- with(subset(aim, comparison == "current"), setNames(meets, paste(n_loci, odds_ratio, raf, threshold)))
aim_summary <- by_key(aim, c("comparison", "threshold"), function(x) {
  key <- paste(x$n_loci, x$odds_ratio, x$raf, x$threshold)
  data.frame(share_reaching_mcid = mean(x$meets), verdicts_changed = sum(x$meets != cur_meets[key]),
             lost = sum(!x$meets & cur_meets[key]), gained = sum(x$meets & !cur_meets[key]),
             flag_disagree_max = max(x$flag_disagree), clin_slope_min = min(x$clin_slope, na.rm = TRUE),
             clin_max_abs_err_max = max(x$clin_max_abs_err))
})
frontier <- by_key(aim, c("comparison", "threshold", "odds_ratio", "raf"), function(x) {
  x <- x[order(x$n_loci), ]
  data.frame(min_loci = first_reaching(x$n_loci, x$gain_oracle, x$mcid[1]))
})

# ---- 4. The polygenic score benchmark --------------------------------------------------------------------------------------------------
comparisons <- c("current", "marginal", "profile")
z <- prs_nodes(400)
grid_z <- list(G = matrix(z, dimnames = list(NULL, "PRS")), prob = rep(1 / length(z), length(z)))
prs_ceiling <- do.call(rbind, lapply(comparisons, function(cmp) do.call(rbind, lapply(at, function(t) {
  f <- function(b) calibration_gains(q, lp, grid_z, c(PRS = b), c(PRS = b), t, comparisons = cmp)$gain_oracle - mcid[[as.character(t)]]
  b <- uniroot(f, c(0.05, 2.5), tol = 1e-6)$root
  cal <- calibration_gains(q, lp, grid_z, c(PRS = b), c(PRS = b), t, comparisons = cmp, calibration = TRUE)
  data.frame(comparison = cmp, threshold = t, or_sd_needed = exp(b), clin_slope = cal$clin_slope, clin_citl = cal$clin_citl,
             clin_max_abs_err = cal$clin_max_abs_err, flag_disagree = cal$flag_disagree)
}))))

# The quadrature of script 30 (400 x 64 nodes): cases needed scale with rho^2 / (1 - rho^2), which is steep near rho = 0.9, so a coarser
# grid moved them by up to 5%
s_nodes <- prs_nodes(400); u_nodes <- prs_nodes(64)
grid_su <- list(G = cbind(S = rep(s_nodes, times = length(u_nodes)), U = rep(u_nodes, each = length(s_nodes))),
                prob = rep(1 / (length(s_nodes) * length(u_nodes)), length(s_nodes) * length(u_nodes)))
tasks <- expand.grid(or_sd = c(2, 2.5, 3, 4, 5, 6), threshold = at, comparison = comparisons, stringsAsFactors = FALSE)
ctx <- list(q = q, lp = lp, grid = grid_su, mcid = mcid, tasks = tasks)
prs_rho <- parallel_rows(cl, as.list(seq_len(nrow(tasks))), function(i, ctx) {
  tk <- ctx$tasks[i, ]; b <- log(tk$or_sd)
  f <- function(rho) calibration_gains(ctx$q, ctx$lp, ctx$grid, c(S = b * rho, U = b * sqrt(1 - rho^2)), c(S = b * rho, U = 0),
                                       tk$threshold, comparisons = tk$comparison)$gain_score - ctx$mcid[[as.character(tk$threshold)]]
  rho <- if (f(1) < 0) NA_real_ else stats::uniroot(f, c(0, 1), tol = 1e-4)$root
  data.frame(tk, rho_needed = rho)
}, ctx)
# The current comparison must agree with script 30: the OR per SD needed with true weights (exactly), and the accuracy needed with estimated
# weights (to 0.01, a sanity bound). They do not agree exactly: script 30 takes the shift that holds prevalence from 400 nodes of the score
# but evaluates the truth on 400 x 64 nodes, and the 64-node part has variance 0.980, so its truth holds prevalence only to 0.2% to 1.0%
# (addendum 23). Here the shift is exact on the grid used. The accuracy needed differs by up to 0.004 (2% to 5% in cases needed).
ref_sum <- read.csv("output/prs_benchmark_summary.csv")
cur_ceil <- prs_ceiling[prs_ceiling$comparison == "current", ]
stopifnot(max(abs(cur_ceil$or_sd_needed - ref_sum$or_sd_needed[match(cur_ceil$threshold, ref_sum$threshold)])) < 1e-3)
ref_dz <- read.csv("output/prs_discovery_size.csv"); ref_dz <- ref_dz[ref_dz$m_eff == 1e3, ]
cur_rho <- prs_rho[prs_rho$comparison == "current", ]
j <- match(paste(cur_rho$or_sd, cur_rho$threshold), paste(round(ref_dz$or_sd, 1), ref_dz$threshold))
stopifnot(!anyNA(j), identical(is.na(cur_rho$rho_needed), is.na(ref_dz$rho_needed[j])),
          max(abs(cur_rho$rho_needed - ref_dz$rho_needed[j]), na.rm = TRUE) < 1e-2)

m_eff <- c(1e3, 1e4, 6e4)
prs <- merge(prs_rho, prs_ceiling[, c("comparison", "threshold", "or_sd_needed")], by = c("comparison", "threshold"))
prs <- do.call(rbind, lapply(m_eff, function(m) transform(prs, m_eff = m, cases_balanced = 2 * prs_n_effective(rho_needed, log(or_sd), m))))

parallel::stopCluster(cl)

dir.create("output", showWarnings = FALSE)
write.csv(layer, "output/calibration_check_layer.csv", row.names = FALSE)
write.csv(artefact, "output/calibration_check_artefact.csv", row.names = FALSE)
write.csv(aim_summary, "output/calibration_check_aim4.csv", row.names = FALSE)
write.csv(aim[, c("comparison", "n_loci", "odds_ratio", "raf", "threshold", "gain_oracle", "extra_tp_oracle", "mcid", "meets")],
          "output/calibration_check_aim4_sweep.csv", row.names = FALSE)
write.csv(frontier, "output/calibration_check_aim4_frontier.csv", row.names = FALSE)
write.csv(rbind(cbind(part = "ceiling", prs_ceiling[, 1:3], rho_needed = NA, or_sd = NA, m_eff = NA, cases_balanced = NA),
                cbind(part = "estimated weights", prs[, c("comparison", "threshold", "or_sd_needed", "rho_needed", "or_sd", "m_eff", "cases_balanced")])),
          "output/calibration_check_prs.csv", row.names = FALSE)

message("Four-term panel, ceiling, pooled-mean reading:")
print(subset(layer, panel == "four terms" & score == "ceiling" & truth == "mean")[, c("comparison", "threshold", "gain", "lo", "hi", "p_meets_mcid", "flag_disagree", "clin_slope", "clin_max_abs_err")], row.names = FALSE)
