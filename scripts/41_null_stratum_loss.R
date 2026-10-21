# Why the five-term realised score loses 3.6 per 1,000 neutropenic admissions when the truth has no effect there (addendum 33). Run from the project root, after script 22.
#
# In the scenario with all of the effect in non-neutropenic admissions the truth has NO genetic effect in the neutropenic stratum (N), so the clinical score is exactly calibrated
# there and any change the genetic score makes can only cost. The realised score (pooled estimates, recentred globally to the stratum's mean risk) loses 1.07 per 1,000 N admissions
# with four terms and 3.61 with five. Pooled point estimates, 5%, per 1,000 N admissions. This script breaks the loss down:
#   1. by term: the score with one term only, all terms, and all terms without the A-G term;
#   2. by the scale of the weights (a multiplier k) with the share of N's patients whose recentred shift crosses the threshold;
#   3. by profile, using the identity that, with no effect in the truth, the net-benefit change is the sum over profiles of (clinical risk - 5%) / 0.95 times the change in that
#      profile's flagged share (the script stops if the identity does not reproduce the stratum gain).
#
# Writes   output/null_stratum_loss_terms.csv, null_stratum_loss_scaling.csv, null_stratum_loss_profiles.csv

for (f in c("baseline", "bs_data", "decision_curve", "boscore_data", "refit_baseline", "evidence_mapping", "genetic_layer", "aim3b", "calibration_check", "scheme_comparison")) source(file.path("R", paste0(f, ".R")))

base <- select_cohort(load_base_cohort(), include_twins = FALSE); lp <- profile_logit(load_refit_params())
q <- tabulate(profile_index(base$X), 2L^length(PREDICTORS)) / nrow(base$X); strata <- neutropenia_strata(q, lp); th <- 0.05
s5 <- clec7a_second_setup(n_draw = 10L); b4 <- pooled_estimates(load_term_draws()); b5 <- s5$beta_hat
cN <- list(`four terms` = trim_cells(genetic_cells(strata$N$q, lp, dosage_grid(genetic_term_frequencies()))), `five terms` = trim_cells(genetic_cells(strata$N$q, lp, s5$grid)))
bs <- list(`four terms` = b4, `five terms` = b5)
gain_n <- function(cl, w) 1000 * stratum_gains(cl, w * 0, w, "global", th)$real_clin        # the truth has no effect: b = 0
delta <- function(cl, w) { l <- cl$lp_base + drop(cl$G %*% w[colnames(cl$G)]); (l + cell_shift(l, cl$w, cl$profile, plogis(cl$lp_base), cl$pred_mean, "global")) - cl$lp_base }
lg_dist <- function(cl, p) qlogis(plogis(cl$lp_base[cl$profile == p][1])) - qlogis(th)

# 1. By term
terms <- rbind(data.frame(panel = "four terms", score = "all four terms", gain = gain_n(cN[[1]], b4)),
               data.frame(panel = "five terms", score = "all five terms", gain = gain_n(cN[[2]], b5)),
               data.frame(panel = "five terms", score = "without the A-G term", gain = gain_n(cN[[2]], replace(b5, "CLEC7A_AG", 0))),
               do.call(rbind, lapply(names(b5), function(t) data.frame(panel = "five terms", score = paste(t, "alone"), gain = gain_n(cN[[2]], b5 * (names(b5) == t))))))
sd_shift <- function(cl, w) { d <- drop(cl$G %*% w[colnames(cl$G)]); m <- sum(cl$w * d) / sum(cl$w); c(mean = m, sd = sqrt(sum(cl$w * (d - m)^2) / sum(cl$w))) }
spread <- rbind(four_terms = sd_shift(cN[[1]], b4), five_terms = sd_shift(cN[[2]], b5)); terms$sd_of_genetic_shift <- NA_real_
terms$sd_of_genetic_shift[terms$score == "all four terms"] <- spread["four_terms", "sd"]; terms$sd_of_genetic_shift[terms$score == "all five terms"] <- spread["five_terms", "sd"]

# 2. By the scale of the weights, with the share of profile 36's patients (a large profile just above the threshold) that the score unflags and the share of profile 35's it flags
need36 <- lg_dist(cN[[2]], 36); need35 <- lg_dist(cN[[2]], 35)
frac <- function(cl, d, p, test) { i <- cl$profile == p; sum(cl$w[i][test(d[i])]) / sum(cl$w[i]) }
scaling <- do.call(rbind, lapply(c(0.25, 0.5, 0.75, 1, 1.25, 1.5, 2), function(k) do.call(rbind, lapply(names(cN), function(m) { d <- delta(cN[[m]], k * bs[[m]])
  data.frame(panel = m, k = k, gain = gain_n(cN[[m]], k * bs[[m]]), share_profile36_unflagged = frac(cN[[m]], d, 36, function(x) x < -need36), share_profile35_flagged = frac(cN[[m]], d, 35, function(x) x > -need35)) }))))

# 3. By profile
prof_tab <- function(cl, w) { s <- plogis(cl$lp_base + delta(cl, w))
  do.call(rbind, lapply(sort(unique(cl$profile)), function(p) { i <- cl$profile == p; c_ <- plogis(cl$lp_base[i][1]); wt <- sum(cl$w[i]); fl_c <- if (c_ >= th) wt else 0; fl_s <- sum(cl$w[i][s[i] >= th])
    data.frame(profile = p, weight_per_1000_N = 1000 * wt, clinical_risk = c_, logit_distance_from_threshold = qlogis(c_) - qlogis(th), flagged_by_clinical = fl_c / wt, flagged_by_score = fl_s / wt,
               contribution = 1000 * (c_ - th) / (1 - th) * (fl_s - fl_c)) })) }
prof <- do.call(rbind, lapply(names(cN), function(m) cbind(panel = m, prof_tab(cN[[m]], bs[[m]]))))
for (m in names(cN)) stopifnot(abs(sum(prof$contribution[prof$panel == m]) - terms$gain[terms$panel == m & grepl("^all", terms$score)]) < 1e-8)

dir.create("output", showWarnings = FALSE)
write.csv(terms, "output/null_stratum_loss_terms.csv", row.names = FALSE); write.csv(scaling, "output/null_stratum_loss_scaling.csv", row.names = FALSE); write.csv(prof, "output/null_stratum_loss_profiles.csv", row.names = FALSE)
message("1. By term (gain per 1,000 N admissions):"); print(terms, row.names = FALSE, digits = 3)
message("\n2. By the scale of the weights:"); print(scaling, row.names = FALSE, digits = 3)
message("\n3. The profiles that carry the five-term loss:"); print(head(prof[prof$panel == "five terms", ][order(prof$contribution[prof$panel == "five terms"]), ], 6), row.names = FALSE, digits = 3)
