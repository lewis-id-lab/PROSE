# The genetic layer (addendum 11): four recipient-genotype terms added to the refit clinical baseline, evaluated exactly on the
# base cohort's profile distribution.
#
# Model. Each term j is the dosage (0, 1, 2) of its risk haplotype class, drawn under Hardy-Weinberg with the gnomAD NFE class
# frequency p_j. The four genes sit on different chromosomes (PTX3 3q, CLEC7A 12p, TLR4 9q, CD209 19p), so dosages are
# independent of each other, and they are taken to be independent of the clinical profile (an assumption: nothing in the
# evidence links genotype to the clinical predictors). Genotype is a property of the patient and covariates are properties of
# the admission, so the same dosages apply to all of a patient's admissions; this does not change the expected net benefit
# per admission, only its sampling variability, which assign_patient_genotypes() and simulate_nb() reproduce.
#
#   true logit = baseline logit(profile) + sum_j beta_j * dosage_j + a,   a chosen so that prevalence is unchanged (D5).
#
# Scores compared at the same truth:
#   baseline  the refit clinical model alone;
#   ceiling   the true effects are in the score (the best a genetic score could do);
#   realised  the pooled estimates (posterior medians, what a literature-derived score would use) are in the score, recentred
#             so its mean predicted risk equals the baseline prevalence.
#
# Requires R/baseline.R, R/decision_curve.R, R/refit_baseline.R, R/evidence_mapping.R (block_haplotypes).

# Terms and the risk class each one's effect refers to (addenda 7 to 9). Panel fixed by D14: four terms.
GENETIC_TERMS <- list(
  PTX3   = list(prefix = "ptx3_h2", block = "PTX3_trio", risk = c("A", NA, NA)),
  CLEC7A = list(prefix = "clec7a",  block = "CLEC7A",    risk = c("C", NA)),
  TLR4   = list(prefix = "tlr4",    block = "TLR4",      risk = c("G", "T")),
  CD209  = list(prefix = "cd209",   block = "CD209",     risk = c("G", NA))
)

# Frequency of each term's risk class on a chromosome (gnomAD NFE rescaled)
genetic_term_frequencies <- function(terms = GENETIC_TERMS, path = "data/haplotype_frequencies_phased.csv") {
  vapply(terms, function(t) {
    h <- block_haplotypes(t$block, t$risk, path = path, source = "gnomAD_NFE")
    sum(h$freq[h$risk == 1])
  }, numeric(1))
}

# Posterior draws of (mu, tau) for each term, one data frame per term
load_term_draws <- function(terms = GENETIC_TERMS, dir = "output") {
  lapply(terms, function(t) utils::read.csv(file.path(dir, paste0(t$prefix, "_posterior_draws.csv"))))
}

# Pooled estimate each term would carry in a literature-derived score: the posterior median of mu
pooled_estimates <- function(draws) vapply(draws, function(d) stats::median(d$mu), numeric(1))

# True per-copy effects for n_draw simulated worlds (rows) by term (columns). Draw i uses row i of every term's draws, which are
# independent posterior samples. "mean" takes the pooled mean effect as the truth; "predictive" adds between-study
# heterogeneity, theta ~ N(mu, tau), as for a new population.
draw_term_effects <- function(draws, n_draw, heterogeneity = c("mean", "predictive"), seed) {
  heterogeneity <- match.arg(heterogeneity)
  stopifnot(all(vapply(draws, nrow, integer(1)) >= n_draw))
  set.seed(seed)
  out <- vapply(draws, function(d) {
    mu <- d$mu[seq_len(n_draw)]
    if (heterogeneity == "predictive") stats::rnorm(n_draw, mu, d$tau[seq_len(n_draw)]) else mu
  }, numeric(n_draw))
  matrix(out, nrow = n_draw, dimnames = list(NULL, names(draws)))
}

# All dosage combinations (3^k) with their Hardy-Weinberg probabilities
dosage_grid <- function(p) {
  g <- as.matrix(expand.grid(rep(list(0:2), length(p)), KEEP.OUT.ATTRS = FALSE))
  colnames(g) <- names(p)
  list(G = g, prob = apply(g, 1, function(x) prod(stats::dbinom(x, 2, p))))
}

# Cells = clinical profile x dosage combination. q is the profile distribution (length 128).
#   offset: added to the logit of the true risk only, not to the scores (the allo-HSCT stratum, where the true odds are higher than the
#   baseline score knows). `prevalence` is the true event rate, `pred_mean` the mean risk the baseline score predicts.
genetic_cells <- function(q, lp_profile, grid, offset = 0) {
  np <- length(q); ng <- nrow(grid$G)
  prof <- rep(seq_len(np), times = ng); geno <- rep(seq_len(ng), each = np)
  list(w = as.numeric(outer(q, grid$prob)), lp_base = lp_profile[prof], G = grid$G[geno, , drop = FALSE],
       profile = prof, offset = offset, prevalence = sum(q * stats::plogis(lp_profile + offset)),
       pred_mean = sum(q * stats::plogis(lp_profile)))
}

# The shift a that holds mean risk at `target` when a linear predictor changes (D5)
recentre_shift <- function(lp, w, target) {
  # extendInt widens the bracket when a large effect in a small stratum needs a shift beyond +-8
  stats::uniroot(function(a) sum(w * stats::plogis(lp + a)) - target, c(-8, 8), tol = 1e-12, extendInt = "yes")$root
}

# Cell-level logits of the true risk, of the ceiling score (true effects, no allo offset) and of the realised score, each recentred
# (the truth to the true prevalence, the scores to the baseline score's mean risk, D5), and the shifts used.
genetic_scores <- function(cells, beta_true, beta_score) {
  w <- cells$w
  off <- if (is.null(cells$offset)) 0 else cells$offset
  pm <- if (is.null(cells$pred_mean)) cells$prevalence else cells$pred_mean
  gl_true <- drop(cells$G %*% beta_true[colnames(cells$G)])
  lp_true <- cells$lp_base + off + gl_true
  a_true <- recentre_shift(lp_true, w, cells$prevalence)
  lp_ceil <- cells$lp_base + gl_true
  a_ceil <- recentre_shift(lp_ceil, w, pm)
  lp_real <- cells$lp_base + drop(cells$G %*% beta_score[colnames(cells$G)])
  a_real <- recentre_shift(lp_real, w, pm)
  list(lp_true = lp_true + a_true, lp_ceiling = lp_ceil + a_ceil, lp_real = lp_real + a_real,
       a_true = a_true, a_ceil = a_ceil, a_real = a_real)
}

# Net benefit, true positives, false positives, flagged share and expected aROC of the three scores under one set of true effects.
#   beta_true   true per-copy effects (named by term); beta_score the effects carried by the "realised" score
genetic_nb <- function(cells, beta_true, beta_score, thresholds) {
  w <- cells$w
  s <- genetic_scores(cells, beta_true, beta_score)
  r_true <- stats::plogis(s$lp_true)
  scores <- list(baseline = cells$lp_base, ceiling = s$lp_ceiling, realised = s$lp_real)
  do.call(rbind, lapply(names(scores), function(sc) {
    d <- dca_exact(w, r_true, stats::plogis(scores[[sc]]), thresholds = thresholds)
    data.frame(scenario = sc, threshold = d$threshold, nb = d$nb, tp = d$tp, fp = d$fp, flagged = d$flagged,
               auc = expected_auc(scores[[sc]], w * r_true, w * (1 - r_true)), stringsAsFactors = FALSE)
  }))
}

# One row per draw, scenario and threshold, with the gain over the baseline under the same truth and the MCID test.
#   mcid  named by threshold: the pre-specified gain (20% more true cases at no extra false positives) from the refit baseline
run_genetic_layer <- function(cells, effects, beta_score, thresholds, mcid) {
  res <- lapply(seq_len(nrow(effects)), function(i) {
    r <- genetic_nb(cells, effects[i, ], beta_score, thresholds)
    base <- r[r$scenario == "baseline", ]
    r$draw <- i
    r$delta_nb <- r$nb - base$nb[match(r$threshold, base$threshold)]
    r$mcid_truth <- 0.2 * base$tp[match(r$threshold, base$threshold)]
    r
  })
  out <- do.call(rbind, res)
  out$mcid <- unname(mcid[as.character(out$threshold)])
  out$meets_mcid <- out$delta_nb >= out$mcid
  out
}

# One dosage vector per patient, copied to all of the patient's admissions
assign_patient_genotypes <- function(patient, p, seed) {
  set.seed(seed)
  u <- unique(patient)
  g <- vapply(p, function(pj) stats::rbinom(length(u), 2, pj), numeric(length(u)))
  g <- matrix(g, nrow = length(u), dimnames = list(NULL, names(p)))
  g[match(patient, u), , drop = FALSE]
}

# The cohort repeated R times, each repetition with its own patients (for checking the exact engine against simulation)
replicate_cohort <- function(cohort, R) {
  idx <- rep(seq_len(nrow(cohort$X)), times = R)
  rep_id <- rep(seq_len(R), each = nrow(cohort$X))
  list(rows = transform(cohort$rows[idx, , drop = FALSE], patient = paste(patient, rep_id, sep = "_")),
       X = cohort$X[idx, , drop = FALSE])
}

# Net benefit by simulating patient genotypes and outcomes on a cohort's rows. Returns the baseline and ceiling scores' NB.
simulate_nb <- function(cohort, lp_profile, p, beta_true, prevalence, thresholds, seed) {
  G <- assign_patient_genotypes(cohort$rows$patient, p, seed)
  lp_b <- lp_profile[profile_index(cohort$X)]
  lp_t <- lp_b + drop(G %*% beta_true[colnames(G)])
  # The same a as the exact engine: the shift that holds the population prevalence
  grid <- dosage_grid(p)
  q <- tabulate(profile_index(cohort$X), length(lp_profile)) / nrow(cohort$X)
  cells <- genetic_cells(q, lp_profile, grid)
  lp_cell <- cells$lp_base + drop(cells$G %*% beta_true[colnames(cells$G)])
  lp_t <- lp_t + recentre_shift(lp_cell, cells$w, cells$prevalence)
  set.seed(seed + 1L)
  y <- stats::rbinom(length(lp_t), 1, stats::plogis(lp_t))
  nb <- function(score) vapply(thresholds, function(t) { f <- stats::plogis(score) >= t
    mean(f & y == 1) - mean(f & y == 0) * t / (1 - t) }, numeric(1))
  data.frame(threshold = thresholds, baseline = nb(lp_b), ceiling = nb(lp_t))
}

# ---- Per-term attribution --------------------------------------------------------------------------------------------------

# True and false positive rates of a rule (flag when the predicted risk reaches the threshold), one row per threshold
tp_fp <- function(w, r_true, r_pred, thresholds) {
  t(vapply(thresholds, function(t) { f <- r_pred >= t; c(tp = sum(w * r_true * f), fp = sum(w * (1 - r_true) * f)) }, numeric(2)))
}
nb_from <- function(x, thresholds) x[, "tp"] - x[, "fp"] * thresholds / (1 - thresholds)

# Shapley decomposition of the net-benefit gain over the baseline into the four terms. The value of a subset S of terms is the
# gain of a score that carries only the effects of S (recentred), under the full true model; the empty set is the baseline.
# The four Shapley contributions sum to the gain of the full score. Net benefit is not additive in the terms (flagging is a
# threshold on the sum), so this is a fair division of the total rather than a separable effect. `beta_effects` is what the score
# carries: the true effects for the ceiling, the pooled estimates for the realised score.
# Also returns each term alone and the loss from leaving it out of the full score.
shapley_gain <- function(cells, beta_true, beta_effects, thresholds) {
  terms <- colnames(cells$G); k <- length(terms); w <- cells$w
  off <- if (is.null(cells$offset)) 0 else cells$offset
  pm <- if (is.null(cells$pred_mean)) cells$prevalence else cells$pred_mean
  lp_t <- cells$lp_base + off + drop(cells$G %*% beta_true[terms])
  r_true <- stats::plogis(lp_t + recentre_shift(lp_t, w, cells$prevalence))
  gl <- sweep(cells$G, 2, beta_effects[terms], "*")
  subsets <- as.matrix(expand.grid(rep(list(0:1), k), KEEP.OUT.ATTRS = FALSE))
  v <- vapply(seq_len(nrow(subsets)), function(i) {
    lp <- cells$lp_base
    if (any(subsets[i, ] == 1)) { lp <- lp + drop(gl %*% subsets[i, ]); lp <- lp + recentre_shift(lp, w, pm) }
    nb_from(tp_fp(w, r_true, stats::plogis(lp), thresholds), thresholds)
  }, numeric(length(thresholds)))
  v <- matrix(v, nrow = length(thresholds)); v <- v - v[, 1]                 # gain over the empty set
  code <- function(m) 1L + sum(m * 2L^(seq_len(k) - 1L))
  size <- rowSums(subsets)
  phi <- vapply(seq_len(k), function(j) {
    without <- which(subsets[, j] == 0)
    contrib <- vapply(without, function(i) {
      m1 <- subsets[i, ]; m1[j] <- 1L
      factorial(size[i]) * factorial(k - size[i] - 1) / factorial(k) * (v[, code(m1)] - v[, i])
    }, numeric(length(thresholds)))
    rowSums(matrix(contrib, nrow = length(thresholds)))
  }, numeric(length(thresholds)))
  phi <- matrix(phi, nrow = length(thresholds), dimnames = list(thresholds, terms))
  full <- code(rep(1L, k))
  list(phi = phi, total = v[, full],
       alone = matrix(vapply(seq_len(k), function(j) v[, code(replace(integer(k), j, 1L))], numeric(length(thresholds))), nrow = length(thresholds), dimnames = list(thresholds, terms)),
       leave_out = matrix(vapply(seq_len(k), function(j) v[, full] - v[, code(replace(rep(1L, k), j, 0L))], numeric(length(thresholds))), nrow = length(thresholds), dimnames = list(thresholds, terms)))
}

# ---- Correlation ensemble ---------------------------------------------------------------------------------------------------

# MCID (20% more true cases at no extra false positives) of the baseline model in a profile distribution, per threshold
baseline_mcid <- function(q, lp_profile, thresholds, offset = 0) dca_exact(q, stats::plogis(lp_profile + offset), stats::plogis(lp_profile), thresholds = thresholds)$mcid_b

# The genetic layer over an ensemble of profile distributions (rows of q_members): draw i of the true effects is evaluated in the
# member member_of_draw[i], so the spread reflects both the effects and the unpublished clinical correlation structure. The MCID
# is the one of the baseline in that member's own distribution.
run_genetic_layer_ensemble <- function(q_members, member_of_draw, lp_profile, grid, effects, beta_score, thresholds, draw_id = seq_len(nrow(effects))) {
  cache <- new.env()
  member_cells <- function(m) {
    key <- as.character(m)
    if (is.null(cache[[key]])) {
      cc <- genetic_cells(q_members[m, ], lp_profile, grid)
      cc$mcid <- baseline_mcid(q_members[m, ], lp_profile, thresholds)
      cache[[key]] <- cc
    }
    cache[[key]]
  }
  out <- lapply(seq_len(nrow(effects)), function(i) {
    cc <- member_cells(member_of_draw[i])
    r <- genetic_nb(cc, effects[i, ], beta_score, thresholds)
    b <- match(r$threshold, r$threshold[r$scenario == "baseline"])
    base <- r[r$scenario == "baseline", ]
    r$draw <- draw_id[i]; r$member <- member_of_draw[i]
    r$delta_nb <- r$nb - base$nb[b]; r$delta_tp <- r$tp - base$tp[b]
    r$mcid <- cc$mcid[b]
    r
  })
  out <- do.call(rbind, out)
  out$meets_mcid <- out$delta_nb >= out$mcid
  out
}

# ---- Sampling variability of a study of this size -------------------------------------------------------------------------------

# Patient groups of a cohort, for resampling
cohort_groups <- function(cohort) list(prof = profile_index(cohort$X), groups = split(seq_len(nrow(cohort$X)), cohort$rows$patient))

# One simulated study of the cohort's size: patients resampled with their admissions (so the clinical mix varies), a fresh genotype
# per patient, outcomes drawn from the true risk, and the observed net benefit of each score. The scores are fixed, not fitted
# in the sample. `shifts` are the recentring shifts of the exact engine (genetic_scores); mcid is the population MCID per threshold.
simulate_study <- function(cg, lp_profile, p, shifts, beta_true, beta_score, thresholds, mcid, resample = TRUE) {
  by <- cg$groups
  idx <- if (resample) sample.int(length(by), replace = TRUE) else seq_along(by)
  rows <- unlist(by[idx], use.names = FALSE)
  pat <- rep(seq_along(idx), lengths(by[idx]))
  g <- matrix(vapply(p, function(pj) stats::rbinom(length(idx), 2, pj), numeric(length(idx))), nrow = length(idx), dimnames = list(NULL, names(p)))
  G <- g[pat, , drop = FALSE]
  lp_b <- lp_profile[cg$prof[rows]]
  lp_t <- lp_b + drop(G %*% beta_true[names(p)]) + shifts$a_true
  lp_r <- lp_b + drop(G %*% beta_score[names(p)]) + shifts$a_real
  y <- stats::rbinom(length(lp_t), 1, stats::plogis(lp_t))
  n <- length(y)
  one <- function(lp, name) {
    r <- stats::plogis(lp)
    d <- t(vapply(thresholds, function(t) { f <- r >= t; c(tp = sum(f & y == 1) / n, fp = sum(f & y == 0) / n) }, numeric(2)))
    data.frame(scenario = name, threshold = thresholds, tp = d[, "tp"], fp = d[, "fp"], nb = nb_from(d, thresholds), stringsAsFactors = FALSE)
  }
  r <- rbind(one(lp_b, "baseline"), one(lp_t, "ceiling"), one(lp_r, "realised"))
  b <- match(r$threshold, r$threshold[r$scenario == "baseline"])
  base <- r[r$scenario == "baseline", ]
  r$delta_nb <- r$nb - base$nb[b]
  r$delta_tp <- r$tp - base$tp[b]
  r$mcid_sample <- 0.2 * base$tp[b]
  r$mcid <- unname(mcid[as.character(r$threshold)])
  r$events <- sum(y)
  r
}

# ---- Ensemble members as the source of each study's clinical mix ------------------------------------------------------------------

# One simulated study whose clinical mix comes from a profile distribution q (a copula member, or the empirical one) instead of the
# real rows. Patients are resampled from `counts` (admissions per patient in the real cohort, so the clustering of genotype is kept);
# each admission's profile is drawn independently from q, which ignores any covariate correlation within a patient. A fresh genotype
# per patient (donor genotype for the allo-HSCT stratum), outcomes from the true risk, observed net benefit of fixed scores.
#   shifts  from genetic_scores() on the cells of the same q (a_true, a_ceil, a_real); offset as in genetic_cells()
simulate_study_q <- function(counts, q, lp_profile, p, shifts, beta_true, beta_score, thresholds, mcid, offset = 0, resample_counts = TRUE) {
  cnt <- if (resample_counts) sample(counts, length(counts), replace = TRUE) else counts
  n <- sum(cnt)
  prof <- sample.int(length(q), n, replace = TRUE, prob = q)
  pat <- rep(seq_along(cnt), cnt)
  g <- matrix(vapply(p, function(pj) stats::rbinom(length(cnt), 2, pj), numeric(length(cnt))), nrow = length(cnt), dimnames = list(NULL, names(p)))
  G <- g[pat, , drop = FALSE]
  lp_b <- lp_profile[prof]
  lp_t <- lp_b + offset + drop(G %*% beta_true[names(p)]) + shifts$a_true
  lp_c <- lp_b + drop(G %*% beta_true[names(p)]) + shifts$a_ceil
  lp_r <- lp_b + drop(G %*% beta_score[names(p)]) + shifts$a_real
  y <- stats::rbinom(n, 1, stats::plogis(lp_t))
  one <- function(lp, name) {
    r <- stats::plogis(lp)
    d <- t(vapply(thresholds, function(t) { f <- r >= t; c(tp = sum(f & y == 1) / n, fp = sum(f & y == 0) / n) }, numeric(2)))
    data.frame(scenario = name, threshold = thresholds, tp = d[, "tp"], fp = d[, "fp"], nb = nb_from(d, thresholds), stringsAsFactors = FALSE)
  }
  r <- rbind(one(lp_b, "baseline"), one(lp_c, "ceiling"), one(lp_r, "realised"))
  b <- match(r$threshold, r$threshold[r$scenario == "baseline"])
  base <- r[r$scenario == "baseline", ]
  r$delta_nb <- r$nb - base$nb[b]
  r$delta_tp <- r$tp - base$tp[b]
  r$mcid_sample <- 0.2 * base$tp[b]
  r$mcid <- unname(mcid[as.character(r$threshold)])
  r$events <- sum(y)
  r$admissions <- n
  r
}

# ---- Donor genotype in the allogeneic HSCT stratum (D1, secondary analysis) ---------------------------------------------------------

# Published Table 1 counts (Stanzani 2019): allogeneic HSCT is the reason for admission in 44 of 133 IMD admissions and 248 of 3,994 others.
ALLO_COUNTS <- c(case_allo = 44, case = 133, non_allo = 248, non = 3994)

# The allo-HSCT stratum from the published margins alone. No admission-level allo flag is in the base cohort, so allo status is taken
# to depend on the outcome only, given the clinical profile: P(allo | case) = 44/133 and P(allo | no case) = 248/3994. That gives the
# stratum's profile distribution, q_allo(x) proportional to q(x) [r(x) a1 + (1 - r(x)) a0], a true risk whose odds are the baseline's
# times a1/a0 (5.3) at every profile, and a share of admissions that matches Table 1 (7.1%).
allo_stratum <- function(q, lp_profile, counts = ALLO_COUNTS) {
  a1 <- counts[["case_allo"]] / counts[["case"]]; a0 <- counts[["non_allo"]] / counts[["non"]]
  r <- stats::plogis(lp_profile)
  w <- q * (r * a1 + (1 - r) * a0)
  offset <- log(a1 / a0)
  qa <- w / sum(w)
  list(share = sum(w), q = qa, offset = offset, a1 = a1, a0 = a0, event_rate = sum(qa * stats::plogis(lp_profile + offset)))
}

# Donor-genotype posterior draws of (mu, tau), one data frame per term; NULL where the evidence has no donor-genotype study (CD209)
load_donor_draws <- function(terms = GENETIC_TERMS, dir = "output") {
  lapply(terms, function(t) {
    f <- file.path(dir, paste0(t$prefix, "_donor_posterior_draws.csv"))
    if (file.exists(f)) utils::read.csv(f) else NULL
  })
}

# Donor effects in the same form as draw_term_effects(); a term without donor evidence has no estimable effect and is held at zero
draw_donor_effects <- function(donor_draws, n_draw, heterogeneity = c("mean", "predictive"), seed) {
  heterogeneity <- match.arg(heterogeneity)
  set.seed(seed)
  out <- vapply(donor_draws, function(d) {
    if (is.null(d)) return(numeric(n_draw))
    stopifnot(nrow(d) >= n_draw)
    mu <- d$mu[seq_len(n_draw)]
    if (heterogeneity == "predictive") stats::rnorm(n_draw, mu, d$tau[seq_len(n_draw)]) else mu
  }, numeric(n_draw))
  matrix(out, nrow = n_draw, dimnames = list(NULL, names(donor_draws)))
}

donor_pooled_estimates <- function(donor_draws) vapply(donor_draws, function(d) if (is.null(d)) 0 else stats::median(d$mu), numeric(1))

# The allo-HSCT stratum from the cohort's own flag (addendum 14): the empirical profile distribution of the flagged admissions, and a
# recalibration offset on the baseline logit chosen so that the stratum's predicted event rate equals its observed rate (the baseline
# score has no allo item). The profile distribution rests on a few hundred admissions, so many profiles are empty.
allo_stratum_data <- function(cohort, lp_profile) {
  a <- which(cohort$rows$allo_hsct)
  if (!length(a)) stop("the cohort carries no allo-HSCT flag", call. = FALSE)
  q <- tabulate(profile_index(cohort$X[a, , drop = FALSE]), length(lp_profile)) / length(a)
  events <- sum(cohort$rows$case[a])
  rate <- events / length(a)
  offset <- stats::uniroot(function(o) sum(q * stats::plogis(lp_profile + o)) - rate, c(-6, 8), tol = 1e-12)$root
  list(share = length(a) / nrow(cohort$X), q = q, offset = offset, event_rate = rate, n = length(a), events = events,
       profiles_present = sum(q > 0), predicted_rate_baseline = sum(q * stats::plogis(lp_profile)))
}

# ---- Terms whose genotypes are not independent binomials (the two CLEC7A parameters, addendum 16) --------------------------------------

# Genotype distribution of one gene with several haplotype classes: the two chromosomes are drawn from the haplotype frequencies
# (Hardy-Weinberg) and each term is the number of copies of its haplotype. `classes` is a named character vector, term -> haplotype name.
haplotype_gene_part <- function(hap, classes) {
  ij <- expand.grid(i = seq_len(nrow(hap)), j = seq_len(nrow(hap)))
  G <- vapply(classes, function(h) as.numeric(hap$haplotype[ij$i] == h) + as.numeric(hap$haplotype[ij$j] == h), numeric(nrow(ij)))
  G <- matrix(G, nrow = nrow(ij), dimnames = list(NULL, names(classes)))
  prob <- hap$freq[ij$i] * hap$freq[ij$j]
  key <- apply(G, 1, paste, collapse = "_")
  u <- !duplicated(key)
  list(G = G[u, , drop = FALSE], prob = as.numeric(tapply(prob, factor(key, levels = key[u]), sum)))
}

# Product of independent genes, each a list(G, prob): the joint dosage grid in the form dosage_grid() returns
joint_grid <- function(parts) {
  idx <- expand.grid(lapply(parts, function(p) seq_along(p$prob)), KEEP.OUT.ATTRS = FALSE)
  list(G = do.call(cbind, Map(function(p, i) p$G[i, , drop = FALSE], parts, idx)),
       prob = Reduce(`*`, Map(function(p, i) p$prob[i], parts, idx)))
}

# The allo-HSCT stratum with a smoothed profile distribution (addendum 16): a logistic model of allo status on the seven predictors (main effects) gives each
# of the 128 profiles a probability of being an allo admission, and the stratum's distribution is the cohort's profile distribution times that probability.
# Profiles that no allo admission happens to show get weight if the model says they are plausible. The offset recalibrates the baseline as in allo_stratum_data().
# The share equals the observed share (the model has an intercept). Allo status is almost a part of the high-risk chemotherapy item, so the fit is close to separated
# for profiles without it; the fitted probabilities there are tiny, which is the intended behaviour.
allo_stratum_smooth <- function(cohort, lp_profile, profiles = all_profiles()) {
  a <- cohort$rows$allo_hsct
  if (all(is.na(a))) stop("the cohort carries no allo-HSCT flag", call. = FALSE)
  fit <- suppressWarnings(stats::glm(as.integer(a) ~ cohort$X, family = stats::binomial))
  cf <- stats::setNames(stats::coef(fit), c("(Intercept)", colnames(cohort$X)))
  p_allo <- as.numeric(stats::plogis(cf[["(Intercept)"]] + profiles %*% cf[colnames(profiles)]))
  q_all <- tabulate(profile_index(cohort$X), nrow(profiles)) / nrow(cohort$X)
  w <- q_all * p_allo
  q <- w / sum(w)
  events <- sum(cohort$rows$case[which(a)]); n <- sum(a); rate <- events / n
  offset <- stats::uniroot(function(o) sum(q * stats::plogis(lp_profile + o)) - rate, c(-6, 8), tol = 1e-12)$root
  list(share = sum(w), q = q, offset = offset, event_rate = rate, n = n, events = events, profiles_present = sum(q > 1e-6),
       predicted_rate_baseline = sum(q * stats::plogis(lp_profile)), coefficients = cf, converged = fit$converged)
}

# Inputs for the genetic layer with two CLEC7A parameters (addendum 16): the joint dosage grid of PTX3, CLEC7A (Y238X copies), CLEC7A_AG (A-G copies),
# TLR4 and CD209; the true effects drawn for each reading (the four-term draws plus beta2 from its pooled posterior); and the pooled estimates carried by the realised score.
clec7a_second_setup <- function(n_draw = 2000L, seed = 6724L, dir = "output", hap_path = "data/haplotype_frequencies_phased.csv") {
  hap <- block_haplotypes("CLEC7A", c("C", NA), path = hap_path, source = "gnomAD_NFE")
  p <- genetic_term_frequencies(path = hap_path)
  terms <- c("PTX3", "CLEC7A", "CLEC7A_AG", "TLR4", "CD209")
  grid <- joint_grid(list(dosage_grid(p["PTX3"]), haplotype_gene_part(hap, c(CLEC7A = "C-G", CLEC7A_AG = "A-G")), dosage_grid(p["TLR4"]), dosage_grid(p["CD209"])))
  grid$G <- grid$G[, terms]
  b2 <- utils::read.csv(file.path(dir, "clec7a_second_posterior_draws.csv"))
  ef <- utils::read.csv(file.path(dir, "genetic_layer_effects.csv"))
  stopifnot(nrow(b2) >= n_draw)
  set.seed(seed)
  b2_draw <- list(mean = b2$mu[seq_len(n_draw)], predictive = stats::rnorm(n_draw, b2$mu[seq_len(n_draw)], b2$tau[seq_len(n_draw)]))
  effects <- lapply(c(mean = "mean", predictive = "predictive"), function(h) {
    m <- ef[ef$truth == h & ef$draw <= n_draw, c("PTX3", "CLEC7A", "TLR4", "CD209")]
    cbind(PTX3 = m$PTX3, CLEC7A = m$CLEC7A, CLEC7A_AG = b2_draw[[h]], TLR4 = m$TLR4, CD209 = m$CD209)
  })
  pooled <- load_term_draws(dir = dir)
  pe <- pooled_estimates(pooled)
  beta_hat <- c(PTX3 = pe[["PTX3"]], CLEC7A = pe[["CLEC7A"]], CLEC7A_AG = stats::median(b2$mu), TLR4 = pe[["TLR4"]], CD209 = pe[["CD209"]])
  beta_hi <- c(PTX3 = unname(stats::quantile(pooled$PTX3$mu, 0.95)), CLEC7A = unname(stats::quantile(pooled$CLEC7A$mu, 0.95)), CLEC7A_AG = unname(stats::quantile(b2$mu, 0.95)),
               TLR4 = unname(stats::quantile(pooled$TLR4$mu, 0.95)), CD209 = unname(stats::quantile(pooled$CD209$mu, 0.95)))
  list(grid = grid, terms = terms, effects = effects, beta_hat = beta_hat, beta_hi = beta_hi, hap = hap, p = p)
}
