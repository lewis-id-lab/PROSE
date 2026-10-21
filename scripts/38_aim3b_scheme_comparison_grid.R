# Where does the global-intercept truth start to distort the Aim 3b gain? (addendum 29) Run from the project root, after scripts 22, 36 and 37.
#
# Scripts 36 and 37 compared the truths at a few effect sizes (1, 2, 3 and 4.9 times the pooled effect), scaling all terms together. This runs a fine grid of
# multipliers (1 to 6 in steps of 0.25, and 4.9) and scales either all terms together or ONE term at a time (the others stay at their pooled estimates), with the
# effect in the non-neutropenic stratum, the neutropenic stratum or both, for the four- and five-term panels and for the ceiling and the realised score
# (whose weights are the event-weighted average of the truth's stratum effects). Pooled point estimates; whole-cohort gains over the clinical score.
#
# "Starts to matter" is defined before the results were seen, three ways, each as the smallest multiplier (linear interpolation on the grid) at which
#   (a) the profile gain differs from the global one by at least 10%;
#   (b) the clinical score's induced miscalibration is worth at least 5% of the global gain (the gain against the clinical score exceeds the gain against m(X) by 5%);
#   (c) the verdict against the 5% MCID differs: the global gain reaches the MCID and the profile gain does not (the interval between the two multipliers).
# The per-copy odds ratio of the scaled term at each crossing is reported so that the crossing can be set against the odds ratios in the literature.
# The all-terms rows at 1, 2, 3 and 4.9 must reproduce script 37 (the script stops if they do not).
#
# Writes   output/aim3b_scheme_comparison_grid.csv (every row), output/aim3b_scheme_comparison_grid_crossings.csv (the crossings, per series)

for (f in c("baseline", "bs_data", "decision_curve", "boscore_data", "refit_baseline", "evidence_mapping", "genetic_layer", "genetic_helpers", "aim3b",
            "calibration_check", "scheme_comparison")) source(file.path("R", paste0(f, ".R")))

base <- select_cohort(load_base_cohort(), include_twins = FALSE); lp <- profile_logit(load_refit_params())
q <- tabulate(profile_index(base$X), 2L^length(PREDICTORS)) / nrow(base$X)
strata <- neutropenia_strata(q, lp); w_ev <- strata$w; share <- c(N = strata$N$share, NN = strata$NN$share); th <- c(0.05, 0.20)
mc <- read.csv("output/dca_baseline_refit_summary.csv"); mcid5 <- 1000 * mc$mcid_b[match(0.05, mc$threshold)]

mk_model <- function(name, grid, bhat) list(name = name, bhat = bhat,
  cells = lapply(list(N = strata$N, NN = strata$NN), function(s) trim_cells(genetic_cells(s$q, lp, grid))))
models <- list(`four terms` = mk_model("four terms", dosage_grid(genetic_term_frequencies()), pooled_estimates(load_term_draws())),
               `five terms` = mk_model("five terms", clec7a_second_setup(n_draw = 10L)$grid, clec7a_second_setup(n_draw = 10L)$beta_hat))
place <- list(`non-neutropenic only` = function(b) list(N = b * 0, NN = b), `neutropenic only` = function(b) list(N = b, NN = b * 0), `both strata` = function(b) list(N = b, NN = b))
grid_m <- sort(unique(c(seq(1, 6, by = 0.25), 4.9)))

tasks <- do.call(c, lapply(names(models), function(md) do.call(c, lapply(names(place), function(pl)
  lapply(c("all terms", names(models[[md]]$bhat)), function(tg) list(model = md, placement = pl, target = tg))))))
cl <- start_workers(cores = 12, files = c("R/baseline.R", "R/decision_curve.R", "R/genetic_layer.R", "R/aim3b.R", "R/calibration_check.R", "R/scheme_comparison.R"))
t0 <- Sys.time()
G <- parallel_rows(cl, tasks, function(tk, ctx) {
  md <- ctx$models[[tk$model]]
  do.call(rbind, lapply(ctx$grid_m, function(m) {
    b <- md$bhat; if (tk$target == "all terms") b <- m * b else b[tk$target] <- m * b[tk$target]
    tr <- ctx$place[[tk$placement]](b); bs <- ctx$w_ev * tr$N + (1 - ctx$w_ev) * tr$NN
    gg <- whole_gains(md$cells, ctx$share, tr, bs, "global", ctx$th); pp <- whole_gains(md$cells, ctx$share, tr, bs, "profile", ctx$th)
    do.call(rbind, lapply(c(ceiling = "ceil", realised = "real"), function(k) do.call(rbind, lapply(seq_along(ctx$th), function(i) data.frame(
      model = tk$model, placement = tk$placement, target = tk$target, multiplier = m, score = if (k == "ceil") "ceiling" else "realised", threshold = ctx$th[i],
      or_per_copy_scaled = if (tk$target == "all terms") NA_real_ else exp(b[[tk$target]]), largest_abs_log_or = max(abs(b)),
      global_vs_clinical = 1000 * gg[[paste0(k, "_clin")]][i], global_vs_marginal = 1000 * gg[[paste0(k, "_marg")]][i], profile = 1000 * pp[[paste0(k, "_clin")]][i],
      comparator_share = (gg[[paste0(k, "_clin")]][i] - gg[[paste0(k, "_marg")]][i]) / gg[[paste0(k, "_clin")]][i],
      profile_vs_global = pp[[paste0(k, "_clin")]][i] / gg[[paste0(k, "_clin")]][i] - 1, clinical_slope_N_global = gg$slope_N, clinical_slope_NN_global = gg$slope_NN)))))
  }))
}, list(models = models, place = place, grid_m = grid_m, w_ev = w_ev, share = share, th = th))
parallel::stopCluster(cl)
rownames(G) <- NULL
message(sprintf("grid: %.1f min on 12 workers, %d rows", as.numeric(difftime(Sys.time(), t0, units = "mins")), nrow(G)))

# Gate: the all-terms rows at 1, 2, 3 and 4.9 reproduce script 37
ref <- read.csv("output/aim3b_scheme_comparison_ext_placement.csv"); a <- G[G$target == "all terms", ]
j <- match(paste(ref$model, ref$score, ref$label, ref$multiplier, ref$threshold), paste(a$model, a$score, a$placement, a$multiplier, a$threshold))
stopifnot(!anyNA(j), max(abs(a$global_vs_clinical[j] - ref$global_vs_clinical), abs(a$profile[j] - ref$profile), abs(a$global_vs_marginal[j] - ref$global_vs_marginal)) < 1e-9)

# Crossings: the smallest multiplier at which a metric first reaches a level, by linear interpolation on the grid (NA if never within the grid)
first_cross <- function(m, v, level) { i <- which(v >= level)[1]; if (is.na(i)) return(NA_real_); if (i == 1) return(m[1]); m[i - 1] + (level - v[i - 1]) / (v[i] - v[i - 1]) * (m[i] - m[i - 1]) }
key <- c("model", "placement", "target", "score", "threshold")
C <- do.call(rbind, lapply(split(G, G[key], drop = TRUE), function(d) { d <- d[order(d$multiplier), ]
  m_d <- first_cross(d$multiplier, abs(d$profile_vs_global), 0.10); m_c <- first_cross(d$multiplier, d$comparator_share, 0.05)
  m_g <- if (d$threshold[1] == 0.05) first_cross(d$multiplier, d$global_vs_clinical, mcid5) else NA_real_
  m_p <- if (d$threshold[1] == 0.05) first_cross(d$multiplier, d$profile, mcid5) else NA_real_
  or_at <- function(mm) if (is.na(mm) || d$target[1] == "all terms") NA_real_ else exp(mm * models[[d$model[1]]]$bhat[[d$target[1]]])
  data.frame(d[1, key], m_diff_10pct = m_d, m_comparator_5pct = m_c, m_global_reaches_mcid = m_g, m_profile_reaches_mcid = m_p,
             or_at_diff_10pct = or_at(m_d), or_at_comparator_5pct = or_at(m_c), or_at_global_mcid = or_at(m_g),
             diff_at_1x = d$profile_vs_global[d$multiplier == 1], diff_at_4.9x = d$profile_vs_global[d$multiplier == 4.9]) }))
rownames(C) <- NULL

dir.create("output", showWarnings = FALSE)
write.csv(G, "output/aim3b_scheme_comparison_grid.csv", row.names = FALSE)
write.csv(C, "output/aim3b_scheme_comparison_grid_crossings.csv", row.names = FALSE)
message("Crossings, ceiling score, 5% (NA = not within 1 to 6 times):")
print(C[C$score == "ceiling" & C$threshold == 0.05, c("model", "placement", "target", "m_diff_10pct", "m_comparator_5pct", "m_global_reaches_mcid", "m_profile_reaches_mcid", "or_at_diff_10pct")], row.names = FALSE, digits = 3)
