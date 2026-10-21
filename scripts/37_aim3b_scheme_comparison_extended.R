# Script 36's stress test extended to the realised score and the five-term panel (addendum 28). Run from the project root, after scripts 22 and 36.
#
# Script 36 asked whether the global-intercept truth distorts the Aim 3b scenarios at effects as large as the one that exposed it (a per-copy effect 4.9
# times the pooled one), for the ceiling score of the four-term panel. This repeats it for
#   - the realised score, which carries the pooled estimates. In the placement grid (A) they are the event-weighted average of the stratum effects of the
#     truth (what a pooled estimate would be); in the six scenarios (B, C) they are the pooled estimates, times the scaling in the scaled version, as in script 21;
#   - the five-term panel (PTX3, CLEC7A, CLEC7A_AG, TLR4, CD209, with the joint CLEC7A genotypes of script 22). The scaling multiplies all five effects, so
#     the A-G term gets a log odds ratio near 3 per copy: a stress test, not a plausible effect.
# A and B use the pooled point estimates; C uses 200 posterior draws of the true effects (pooled-mean reading).
# The four-term ceiling rows must reproduce script 36 (the script stops if they do not).
#
# Writes   output/aim3b_scheme_comparison_ext_placement.csv, _ext_scenarios.csv, _ext_draws.csv

for (f in c("baseline", "bs_data", "decision_curve", "boscore_data", "refit_baseline", "evidence_mapping", "genetic_layer", "genetic_helpers", "aim3b",
            "calibration_check", "scheme_comparison")) source(file.path("R", paste0(f, ".R")))

base <- select_cohort(load_base_cohort(), include_twins = FALSE); lp <- profile_logit(load_refit_params())
q <- tabulate(profile_index(base$X), 2L^length(PREDICTORS)) / nrow(base$X)
strata <- neutropenia_strata(q, lp); w_ev <- strata$w; share <- c(N = strata$N$share, NN = strata$NN$share); th <- c(0.05, 0.20)
mc <- read.csv("output/dca_baseline_refit_summary.csv"); mcid5 <- mc$mcid_b[match(0.05, mc$threshold)]
n_draw <- 200L

terms4 <- c("PTX3", "CLEC7A", "TLR4", "CD209")
eff <- read.csv("output/genetic_layer_effects.csv")
s5 <- clec7a_second_setup(n_draw = n_draw)
mk_model <- function(name, grid, bhat, E) list(name = name, bhat = bhat, E = E,
  cells = lapply(list(N = strata$N, NN = strata$NN), function(s) trim_cells(genetic_cells(s$q, lp, grid))))
models <- list(`four terms` = mk_model("four terms", dosage_grid(genetic_term_frequencies()), pooled_estimates(load_term_draws()), as.matrix(eff[eff$truth == "mean" & eff$draw <= n_draw, terms4])),
               `five terms` = mk_model("five terms", s5$grid, s5$beta_hat, s5$effects$mean))
stopifnot(nrow(models[[1]]$E) == n_draw, nrow(models[[2]]$E) == n_draw)

# Rows (both scores, both thresholds) for one truth `tr` and one realised-score weight vector `bs`
cmp_rows <- function(md, label, tr, bs, extra = list()) {
  gg <- whole_gains(md$cells, share, tr, bs, "global", th); pp <- whole_gains(md$cells, share, tr, bs, "profile", th)
  do.call(rbind, lapply(c(ceiling = "ceil", realised = "real"), function(k) do.call(rbind, lapply(seq_along(th), function(i)
    data.frame(extra, model = md$name, score = names(which(c(ceiling = "ceil", realised = "real") == k)), label = label, threshold = th[i],
               global_vs_clinical = 1000 * gg[[paste0(k, "_clin")]][i], global_vs_marginal = 1000 * gg[[paste0(k, "_marg")]][i], profile = 1000 * pp[[paste0(k, "_clin")]][i],
               comparator_share = (gg[[paste0(k, "_clin")]][i] - gg[[paste0(k, "_marg")]][i]) / gg[[paste0(k, "_clin")]][i],
               profile_vs_global = pp[[paste0(k, "_clin")]][i] / gg[[paste0(k, "_clin")]][i] - 1,
               clinical_slope_N_global = gg$slope_N, clinical_slope_NN_global = gg$slope_NN)))))
}

# A. Placement and size
place <- list(`non-neutropenic only` = function(m, b) list(N = b * 0, NN = m * b), `neutropenic only` = function(m, b) list(N = m * b, NN = b * 0),
              `both strata` = function(m, b) list(N = m * b, NN = m * b))
A <- do.call(rbind, lapply(models, function(md) do.call(rbind, lapply(names(place), function(pl) do.call(rbind, lapply(c(1, 2, 3, 4.9), function(m) {
  tr <- place[[pl]](m, md$bhat); cmp_rows(md, pl, tr, w_ev * tr$N + (1 - w_ev) * tr$NN, list(multiplier = m)) }))))))

# B. The six scenarios, as in the paper and scaled so that the largest stratum effect is 4.9 times the pooled one
r_values <- list(`H_F only` = Inf, `H_F, ratio 4` = 4, homogeneous = 1, `Brunel direction, ratio 1/4` = 0.25, `Brunel only` = 0, `Brunel estimate` = c(N = -0.21, NN = 1.87))
mult <- function(tr, b) c(N = unname(tr$N[1] / b[1]), NN = unname(tr$NN[1] / b[1]))
k_of <- function(md, sc) 4.9 / max(abs(mult(scenario_effects(md$bhat, r_values[[sc]], w_ev), md$bhat)))
B <- do.call(rbind, lapply(models, function(md) do.call(rbind, lapply(names(r_values), function(sc) {
  tr <- scenario_effects(md$bhat, r_values[[sc]], w_ev); mu <- mult(tr, md$bhat); k <- k_of(md, sc)
  rbind(cmp_rows(md, sc, tr, md$bhat, list(version = "as in the paper", mult_N = mu[["N"]], mult_NN = mu[["NN"]])),
        cmp_rows(md, sc, list(N = k * tr$N, NN = k * tr$NN), k * md$bhat, list(version = "scaled to 4.9x", mult_N = k * mu[["N"]], mult_NN = k * mu[["NN"]])))
}))))

# C. 200 posterior draws on workers: median gain and the share of draws reaching the 5% MCID
tasks <- do.call(c, lapply(models, function(md) do.call(c, lapply(names(r_values), function(sc) lapply(c("as in the paper", "scaled to 4.9x"), function(ver)
  list(model = md$name, scenario = sc, version = ver, k = if (ver == "as in the paper") 1 else k_of(md, sc)))))))
cl <- start_workers(cores = 12, files = c("R/baseline.R", "R/decision_curve.R", "R/genetic_layer.R", "R/aim3b.R", "R/calibration_check.R", "R/scheme_comparison.R"))
t0 <- Sys.time()
C <- parallel_rows(cl, tasks, function(tk, ctx) {
  md <- ctx$models[[tk$model]]; k <- tk$k
  g <- vapply(seq_len(nrow(md$E)), function(i) {
    tr <- scenario_effects(md$E[i, ], ctx$r_values[[tk$scenario]], ctx$w_ev); tr <- list(N = k * tr$N, NN = k * tr$NN)
    gg <- whole_gains(md$cells, ctx$share, tr, k * md$bhat, "global", ctx$th); pp <- whole_gains(md$cells, ctx$share, tr, k * md$bhat, "profile", ctx$th)
    c(gc = gg$ceil_clin[1], pc = pp$ceil_clin[1], gr = gg$real_clin[1], pr = pp$real_clin[1])
  }, numeric(4))
  do.call(rbind, lapply(c(ceiling = "c", realised = "r"), function(s) data.frame(model = tk$model, score = if (s == "c") "ceiling" else "realised", scenario = tk$scenario, version = tk$version,
    median_global = 1000 * median(g[paste0("g", s), ]), median_profile = 1000 * median(g[paste0("p", s), ]),
    profile_vs_global = median(g[paste0("p", s), ]) / median(g[paste0("g", s), ]) - 1,
    share_reaching_mcid_global = mean(g[paste0("g", s), ] >= ctx$mcid5), share_reaching_mcid_profile = mean(g[paste0("p", s), ] >= ctx$mcid5))))
}, list(models = models, r_values = r_values, w_ev = w_ev, share = share, th = th, mcid5 = mcid5))
parallel::stopCluster(cl)
message(sprintf("draws: %.1f min on 12 workers", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
rownames(C) <- NULL

# The four-term ceiling rows must reproduce script 36
ref_a <- read.csv("output/aim3b_scheme_comparison_placement.csv"); ref_b <- read.csv("output/aim3b_scheme_comparison_scenarios.csv"); ref_c <- read.csv("output/aim3b_scheme_comparison_draws.csv")
a <- A[A$model == "four terms" & A$score == "ceiling", ]; b <- B[B$model == "four terms" & B$score == "ceiling", ]; cc <- C[C$model == "four terms" & C$score == "ceiling", ]
ja <- match(paste(ref_a$label, ref_a$multiplier, ref_a$threshold), paste(a$label, a$multiplier, a$threshold))
jb <- match(paste(ref_b$label, ref_b$version, ref_b$threshold), paste(b$label, b$version, b$threshold))
jc <- match(paste(ref_c$scenario, ref_c$version), paste(cc$scenario, cc$version))
stopifnot(!anyNA(c(ja, jb, jc)),
          max(abs(a$global_vs_clinical[ja] - ref_a$global_vs_clinical), abs(a$profile[ja] - ref_a$profile)) < 1e-9,
          max(abs(b$global_vs_marginal[jb] - ref_b$global_vs_marginal), abs(b$profile[jb] - ref_b$profile)) < 1e-9,
          max(abs(cc$median_global[jc] - ref_c$median_global), abs(cc$median_profile[jc] - ref_c$median_profile), abs(cc$share_reaching_mcid_profile[jc] - ref_c$share_reaching_mcid_profile)) < 1e-9)

dir.create("output", showWarnings = FALSE)
write.csv(A, "output/aim3b_scheme_comparison_ext_placement.csv", row.names = FALSE)
write.csv(B, "output/aim3b_scheme_comparison_ext_scenarios.csv", row.names = FALSE)
write.csv(C, "output/aim3b_scheme_comparison_ext_draws.csv", row.names = FALSE)
message("A. Placement and size, 5% (x1000):"); print(A[A$threshold == 0.05, c("model", "score", "label", "multiplier", "global_vs_clinical", "global_vs_marginal", "profile", "comparator_share", "profile_vs_global")], row.names = FALSE, digits = 3)
message("\nC. Draws (5%, x1000; MCID ", round(1000 * mcid5, 2), "):"); print(C, row.names = FALSE, digits = 3)
