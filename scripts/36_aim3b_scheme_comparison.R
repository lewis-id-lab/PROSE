# Does the global-intercept construction distort the Aim 3b scenarios when their effects are as large as the one that exposed it? (addendum 27)
# Run from the project root, after scripts 08, 11, 14 and 21 (it reads the pooled estimates and the committed draws of the true effects).
#
# Addendum 26 found that the global intercept flatters "all of the effect in non-neutropenic admissions" because that scenario puts a per-copy effect
# about 4.9 times the pooled one into the low-risk stratum. This asks the same question of the other scenarios: for each, the whole-cohort ceiling gain
# under the global truth against the clinical score and against the best clinical-only rule m(X), and under the profile-preserving truth.
#   A. placement x size: the pooled effect times m (1, 2, 3, 4.9) in the non-neutropenic stratum, the neutropenic stratum, or both
#   B. the six scenarios of script 21 as in the paper, and scaled so that the largest stratum effect is 4.9 times the pooled one
#   C. the same with 200 posterior draws (pooled-mean reading): median gain and the share of draws reaching the 5% MCID
# Pooled point estimates in A and B. The Brunel-estimate stratum fractions are those script 21 prints (N -0.21, NN 1.87).
#
# Writes   output/aim3b_scheme_comparison_placement.csv, aim3b_scheme_comparison_scenarios.csv, aim3b_scheme_comparison_draws.csv

for (f in c("baseline", "bs_data", "decision_curve", "boscore_data", "refit_baseline", "evidence_mapping", "genetic_layer", "aim3b", "calibration_check")) source(file.path("R", paste0(f, ".R")))

base <- select_cohort(load_base_cohort(), include_twins = FALSE); lp <- profile_logit(load_refit_params())
q <- tabulate(profile_index(base$X), 2L^length(PREDICTORS)) / nrow(base$X)
grid <- dosage_grid(genetic_term_frequencies()); strata <- neutropenia_strata(q, lp)
cells <- lapply(list(N = strata$N, NN = strata$NN), function(s) trim_cells(genetic_cells(s$q, lp, grid)))
bhat <- pooled_estimates(load_term_draws()); th <- c(0.05, 0.20); w_ev <- strata$w; share <- c(N = strata$N$share, NN = strata$NN$share)
mc <- read.csv("output/dca_baseline_refit_summary.csv"); mcid5 <- mc$mcid_b[match(0.05, mc$threshold)]

nbf <- function(w, r, pred) nb_from(tp_fp(w, r, pred, th), th)
# Mean of a cell-level risk within each profile, mapped back to the cells: the marginal risk m(X)
profile_mean <- function(x, cl) as.numeric((tapply(cl$w * x, cl$profile, sum) / tapply(cl$w, cl$profile, sum))[as.character(cl$profile)])

# Ceiling gain of one stratum over the clinical score and over m(X), and the clinical score's calibration slope within the stratum
stratum_gain <- function(cl, b, scheme) {
  lp_t <- cl$lp_base + drop(cl$G %*% b[colnames(cl$G)])
  r <- plogis(lp_t + cell_shift(lp_t, cl$w, cl$profile, plogis(cl$lp_base), cl$prevalence, scheme))
  clin <- plogis(cl$lp_base); m <- profile_mean(r, cl); self <- nbf(cl$w, r, r)
  first <- !duplicated(cl$profile)
  cal <- calibration_fit(profile_mean(r, cl)[first][order(cl$profile[first])], clin[first][order(cl$profile[first])], as.numeric(tapply(cl$w, cl$profile, sum)))
  list(clin = self - nbf(cl$w, r, clin), marg = self - nbf(cl$w, r, m), slope = unname(cal["slope"]))
}
whole <- function(tr, scheme) {
  g <- list(N = stratum_gain(cells$N, tr$N, scheme), NN = stratum_gain(cells$NN, tr$NN, scheme))
  list(clin = share[["N"]] * g$N$clin + share[["NN"]] * g$NN$clin, marg = share[["N"]] * g$N$marg + share[["NN"]] * g$NN$marg, slope_N = g$N$slope, slope_NN = g$NN$slope)
}
cmp_row <- function(label, tr, extra = list()) {
  gg <- whole(tr, "global"); pp <- whole(tr, "profile")
  do.call(rbind, lapply(seq_along(th), function(i) data.frame(extra, label = label, threshold = th[i],
    global_vs_clinical = 1000 * gg$clin[i], global_vs_marginal = 1000 * gg$marg[i], profile = 1000 * pp$clin[i],
    comparator_share = (gg$clin[i] - gg$marg[i]) / gg$clin[i], profile_vs_global = pp$clin[i] / gg$clin[i] - 1,
    clinical_slope_N_global = gg$slope_N, clinical_slope_NN_global = gg$slope_NN)))
}

# A. Placement and size
zero <- bhat * 0
place <- list(`non-neutropenic only` = function(m) list(N = zero, NN = m * bhat), `neutropenic only` = function(m) list(N = m * bhat, NN = zero),
              `both strata` = function(m) list(N = m * bhat, NN = m * bhat))
A <- do.call(rbind, lapply(names(place), function(pl) do.call(rbind, lapply(c(1, 2, 3, 4.9), function(m) cmp_row(pl, place[[pl]](m), list(multiplier = m))))))

# B. The six scenarios, as in the paper and scaled to a largest stratum effect of 4.9 times the pooled one
r_values <- list(`H_F only` = Inf, `H_F, ratio 4` = 4, homogeneous = 1, `Brunel direction, ratio 1/4` = 0.25, `Brunel only` = 0, `Brunel estimate` = c(N = -0.21, NN = 1.87))
mult <- function(tr) c(N = unname(tr$N[1] / bhat[1]), NN = unname(tr$NN[1] / bhat[1]))
B <- do.call(rbind, lapply(names(r_values), function(sc) {
  tr <- scenario_effects(bhat, r_values[[sc]], w_ev); mu <- mult(tr); k <- 4.9 / max(abs(mu))
  rbind(cmp_row(sc, tr, list(version = "as in the paper", mult_N = mu[["N"]], mult_NN = mu[["NN"]])),
        cmp_row(sc, list(N = k * tr$N, NN = k * tr$NN), list(version = "scaled to 4.9x", mult_N = k * mu[["N"]], mult_NN = k * mu[["NN"]])))
}))

# C. 200 posterior draws
eff <- read.csv("output/genetic_layer_effects.csv"); terms4 <- c("PTX3", "CLEC7A", "TLR4", "CD209")
E <- as.matrix(eff[eff$truth == "mean" & eff$draw <= 200, terms4]); stopifnot(nrow(E) == 200)
ks <- sapply(names(r_values), function(sc) 4.9 / max(abs(mult(scenario_effects(bhat, r_values[[sc]], w_ev)))))
C <- do.call(rbind, lapply(names(r_values), function(sc) do.call(rbind, lapply(c("as in the paper", "scaled to 4.9x"), function(ver) {
  k <- if (ver == "as in the paper") 1 else ks[[sc]]
  g <- vapply(seq_len(nrow(E)), function(i) { tr <- scenario_effects(E[i, ], r_values[[sc]], w_ev); tr <- list(N = k * tr$N, NN = k * tr$NN)
    c(global = whole(tr, "global")$clin[1], profile = whole(tr, "profile")$clin[1]) }, numeric(2))
  data.frame(scenario = sc, version = ver, median_global = 1000 * median(g["global", ]), median_profile = 1000 * median(g["profile", ]),
             profile_vs_global = median(g["profile", ]) / median(g["global", ]) - 1, share_reaching_mcid_global = mean(g["global", ] >= mcid5),
             share_reaching_mcid_profile = mean(g["profile", ] >= mcid5))
}))))

dir.create("output", showWarnings = FALSE)
write.csv(A, "output/aim3b_scheme_comparison_placement.csv", row.names = FALSE)
write.csv(B, "output/aim3b_scheme_comparison_scenarios.csv", row.names = FALSE)
write.csv(C, "output/aim3b_scheme_comparison_draws.csv", row.names = FALSE)
message("A. Placement and size (whole-cohort ceiling gain at 5%, x1000):"); print(A[A$threshold == 0.05, setdiff(names(A), "threshold")], row.names = FALSE, digits = 3)
message("\nB. The six scenarios (5%, x1000):"); print(B[B$threshold == 0.05, setdiff(names(B), "threshold")], row.names = FALSE, digits = 3)
message("\nC. Draws (5%, x1000; MCID ", round(1000 * mcid5, 2), "):"); print(C, row.names = FALSE, digits = 3)
