# Aim 4: the pathogen-specificity scenario (scope decision 1; not in addenda 16 and 17). Requires R/baseline.R, R/decision_curve.R, R/genetic_layer.R, R/aim4_sweep.R
# and R/calibration_check.R (profile_shifts).
#
# Two constructions of the truth, as everywhere (genetic_scheme(), environment variable PROSE_SCHEME): "global" holds the event rate with one intercept (D5, the
# first version, addenda 35 to 39), "profile" (the primary, D27) gives the target-pathogen component an intercept per clinical profile so that the clinical
# score stays calibrated in every profile. The mixture below is then calibrated in every profile as well, because the genotype-insensitive part is the baseline risk.
#
# The sweep of script 23 treats every case of invasive mould disease as genotype-sensitive. Here only a fraction f of cases (the target pathogen, for example
# Aspergillus) arises through the pathway the loci act on; the other 1 - f arise independently of the genotype. The outcome is still "any proven or probable IMD".
#
#   true risk = (1 - f) * plogis(lp_base) + f * plogis(lp_base + log(OR) D + a)
#
# a is the shift that holds the target-pathogen component, and so the total event rate, at its baseline value (D5), so f = 1 reproduces script 23 exactly.
# Two scores are compared, both recentred to the baseline score's mean risk:
#   specific  carries the target-pathogen effect log(OR) (what a pathogen-specific discovery study would estimate), applied to every admission
#   diluted   carries the effect fitted to the all-mould outcome (what a discovery study on any IMD would estimate): the best single-weight score for this truth,
#             found by quasi-binomial regression of the true risk on the baseline logit (offset) and D
# The MCID is the baseline's, as in script 23 (20% more true cases at no extra false positives of the baseline under its own predicted risk).

pathogen_risk <- function(lp_base, gl, w, prevalence, f) {
  a <- recentre_shift(lp_base + gl, w, prevalence)
  (1 - f) * stats::plogis(lp_base) + f * stats::plogis(lp_base + gl + a)
}

# Effect of D on the all-mould outcome: the weighted fractional-logit fit with the baseline logit as an offset
diluted_beta <- function(lp_base, D, w, r_true) {
  if (all(D == D[1])) return(0)
  fit <- suppressWarnings(stats::glm(r_true ~ D + offset(lp_base), family = stats::quasibinomial(), weights = w, start = c(0, 0)))
  unname(stats::coef(fit)["D"])
}

pathogen_cell <- function(q, lp_profile, n_loci, odds_ratio, raf, f, thresholds, scheme = genetic_scheme()) {
  pmf <- dosage_pmf(n_loci, raf); ng <- length(pmf); np <- length(q)
  d <- seq_len(ng) - 1; beta <- log(odds_ratio)
  w <- as.numeric(outer(q, pmf)); lp_cell <- rep(lp_profile, ng); D <- rep(d, each = np)
  pred_mean <- sum(q * stats::plogis(lp_profile))
  clin <- stats::plogis(lp_cell)
  comp <- if (scheme == "global") stats::plogis(lp_cell + beta * D + recentre_shift(lp_cell + beta * D, w, pred_mean))
          else as.numeric(stats::plogis(outer(lp_profile + profile_shifts(lp_profile, beta * d, pmf), beta * d, "+")))
  r_true <- (1 - f) * clin + f * comp
  b_dil <- diluted_beta(lp_cell, D, w, r_true)
  pred <- function(b) {
    if (scheme == "global") { lp <- lp_cell + b * D; stats::plogis(lp + recentre_shift(lp, w, pred_mean)) }
    else as.numeric(stats::plogis(outer(lp_profile + profile_shifts(lp_profile, b * d, pmf), b * d, "+")))
  }
  one <- function(pr) { r <- tp_fp(w, r_true, pr, thresholds); list(r = r, nb = nb_from(r, thresholds)) }
  base <- one(clin); spec <- one(pred(beta)); dil <- one(pred(b_dil))
  data.frame(n_loci = n_loci, odds_ratio = odds_ratio, raf = raf, case_fraction = f, threshold = thresholds,
             beta_diluted = b_dil, or_diluted = exp(b_dil), baseline_nb = base$nb,
             gain_specific = spec$nb - base$nb, gain_diluted = dil$nb - base$nb,
             extra_cases_per_1000_specific = 1000 * (spec$r[, "tp"] - base$r[, "tp"]),
             extra_cases_per_1000_diluted = 1000 * (dil$r[, "tp"] - base$r[, "tp"]))
}

# ---- With estimation error in the weights (addendum 39; the counterpart of sweep_cell_realised in R/aim4_sweep.R) -------------------------------
# The score's weights are estimates with an independent error of standard deviation sigma on each locus (no winner's curse), around one of two centres:
#   specific  the target-pathogen effect log(OR)         (a discovery study that knows each case's pathogen)
#   all-mould the effect fitted to the all-mould outcome (a discovery study on any proven or probable mould disease)
# The truth is the pathogen mixture of pathogen_risk(); genotypes are simulated (N individuals crossed with the exact profile distribution), R replicate
# discovery studies, the same simulated genotypes and errors for both centres. Returns one row per replicate, score and threshold.
pathogen_cell_realised <- function(q, lp_profile, n_loci, odds_ratio, raf, f, sigma, thresholds, R = 12L, N = 2000L, seed = 1L, scheme = genetic_scheme()) {
  set.seed(seed)
  beta <- log(odds_ratio)
  pmf <- dosage_pmf(n_loci, raf); ng <- length(pmf); np <- length(q)
  d <- seq_len(ng) - 1
  w <- as.numeric(outer(q, pmf)); lp_cell <- rep(lp_profile, ng); D0 <- rep(d, each = np)
  pred_mean <- sum(q * stats::plogis(lp_profile))
  clin <- stats::plogis(lp_cell)
  comp <- if (scheme == "global") { a <- recentre_shift(lp_cell + beta * D0, w, pred_mean); stats::plogis(lp_cell + beta * D0 + a) }
          else { a <- profile_shifts(lp_profile, beta * d, pmf); stats::plogis(outer(lp_profile + a, beta * d, "+")) }
  r_cells <- (1 - f) * clin + f * as.numeric(comp)
  r_mat <- matrix(r_cells, np, ng)
  centres <- c(specific = beta, all_mould = diluted_beta(lp_cell, D0, w, r_cells))
  base_nb <- nb_from(tp_fp(w, r_cells, clin, thresholds), thresholds)
  wt <- rep(q, times = N) / N
  lpb <- rep(lp_profile, times = N)
  out <- lapply(seq_len(R), function(r) {
    G <- matrix(stats::rbinom(N * n_loci, 2, raf), N, n_loci)
    D <- rowSums(G)
    eps <- stats::rnorm(n_loci, 0, sigma)
    rt <- if (scheme == "global") (1 - f) * stats::plogis(lpb) + f * stats::plogis(lpb + beta * rep(D, each = np) + a)
          else as.numeric(r_mat[, D + 1L])
    do.call(rbind, lapply(names(centres), function(sc) {
      S <- drop(G %*% (centres[[sc]] + eps))
      pr <- if (scheme == "global") { lps <- lpb + rep(S, each = np); stats::plogis(lps + recentre_shift(lps, wt, pred_mean)) }
            else as.numeric(stats::plogis(outer(lp_profile + profile_shifts(lp_profile, S, rep(1 / N, N)), S, "+")))
      data.frame(replicate = r, score = sc, threshold = thresholds, gain = nb_from(tp_fp(wt, rt, pr, thresholds), thresholds) - base_nb)
    }))
  })
  cbind(n_loci = n_loci, odds_ratio = odds_ratio, raf = raf, case_fraction = f, sigma = sigma, do.call(rbind, out))
}

# ---- Stress tests of the realised score: locus-specific precision, the winner's curse, and incomplete tagging (addendum 42) -------------------------
# Profile-preserving truth only. The k true loci each carry the target-pathogen log odds ratio; a score is built from estimated weights as in
# pathogen_cell_realised(), with three departures that can each be switched on:
#   locus_specific  the loci differ in frequency (log-normal spread of 0.35 around raf, bounded to 0.02 to 0.50) and the error of each weight is inversely
#                   proportional to sqrt(2 p (1 - p)), equal to sigma at the frequency raf; the true dosage distribution is the exact convolution.
#   winners_curse   the discovery study screens candidates_per_locus * k loci, the k true ones and the rest null (same frequencies), and the score takes the k
#                   with the largest z, with their estimated weights: true loci come out with inflated weights, and null loci can be selected.
#   tag_r           the scored genotype of each true locus is a tag correlated r with the causal one (each chromosome copies the causal allele with
#                   probability r, otherwise draws independently); its marginal per-allele effect, and so the estimated weight, is about r times the causal one.
#                   Linkage among the causal loci themselves is not modelled.
# Returns one row per replicate, score and threshold.
pathogen_cell_stress <- function(q, lp_profile, n_loci, odds_ratio, raf, f, sigma, thresholds, locus_specific = FALSE, winners_curse = FALSE, tag_r = 1,
                                 candidates_per_locus = 5L, R = 12L, N = 2000L, seed = 1L) {
  if (genetic_scheme() != "profile") stop("pathogen_cell_stress() builds the profile-preserving truth; set PROSE_SCHEME = \"profile\"", call. = FALSE)
  set.seed(seed)
  k <- n_loci; beta <- log(odds_ratio); np <- length(q)
  spread <- function(n) pmin(pmax(raf * exp(stats::rnorm(n, 0, 0.35)), 0.02), 0.5)
  pj <- if (locus_specific) spread(k) else rep(raf, k)
  pmf <- dosage_pmf(k, pj); ng <- length(pmf); d <- seq_len(ng) - 1
  w <- as.numeric(outer(q, pmf)); lp_cell <- rep(lp_profile, ng); D0 <- rep(d, each = np)
  clin <- stats::plogis(lp_cell)
  comp <- stats::plogis(outer(lp_profile + profile_shifts(lp_profile, beta * d, pmf), beta * d, "+"))
  r_cells <- (1 - f) * clin + f * as.numeric(comp)
  r_mat <- matrix(r_cells, np, ng)
  centres <- c(specific = beta, all_mould = diluted_beta(lp_cell, D0, w, r_cells))
  base_nb <- nb_from(tp_fp(w, r_cells, clin, thresholds), thresholds)
  M <- if (winners_curse) as.integer(candidates_per_locus * k) else k
  pM <- c(pj, if (M > k) { if (locus_specific) spread(M - k) else rep(raf, M - k) })
  se <- sigma * sqrt(raf * (1 - raf) / (pM * (1 - pM)))
  true_idx <- seq_len(k)
  wt <- rep(q, times = N) / N
  out <- lapply(seq_len(R), function(r) {
    P <- matrix(rep(pM, each = N), N, M)
    c1 <- matrix(stats::rbinom(N * M, 1L, P), N, M); c2 <- matrix(stats::rbinom(N * M, 1L, P), N, M)
    if (tag_r < 1) {
      swap <- function(cc) { cc[, true_idx] <- ifelse(matrix(stats::runif(N * k), N, k) < tag_r, cc[, true_idx], matrix(stats::rbinom(N * k, 1L, P[, true_idx]), N, k)); cc }
      t1 <- swap(c1); t2 <- swap(c2)
    } else { t1 <- c1; t2 <- c2 }
    Gt <- t1 + t2
    Dtrue <- rowSums((c1 + c2)[, true_idx, drop = FALSE])
    rt <- as.numeric(r_mat[, Dtrue + 1L])
    z <- stats::rnorm(M)
    do.call(rbind, lapply(names(centres), function(sc) {
      mu <- numeric(M); mu[true_idx] <- tag_r * centres[[sc]]
      est <- mu + z * se
      sel <- if (winners_curse) order(est / se, decreasing = TRUE)[seq_len(k)] else true_idx
      S <- drop(Gt[, sel, drop = FALSE] %*% est[sel])
      pr <- as.numeric(stats::plogis(outer(lp_profile + profile_shifts(lp_profile, S, rep(1 / N, N)), S, "+")))
      data.frame(replicate = r, score = sc, threshold = thresholds, gain = nb_from(tp_fp(wt, rt, pr, thresholds), thresholds) - base_nb,
                 true_selected = sum(sel <= k))
    }))
  })
  cbind(n_loci = n_loci, odds_ratio = odds_ratio, raf = raf, case_fraction = f, sigma = sigma, do.call(rbind, out))
}

# ---- Linkage among the causal loci (addendum 43) ------------------------------------------------------------------------------------------------
# The k causal loci are grouped into blocks of block_size adjacent loci. Within a block, on each chromosome, the alleles are correlated rho (identical frequency raf):
# a shared draw u ~ Bernoulli(raf) is copied to a locus with probability sqrt(rho) and otherwise the locus draws independently, so the pairwise allele correlation is
# exactly rho. Blocks are independent, the two chromosomes are independent, and the dosage distribution of the true genetic score is exact: per chromosome and block,
# X ~ raf * Binomial(m, w + (1 - w) raf) + (1 - raf) * Binomial(m, (1 - w) raf) with w = sqrt(rho), and the dosage is the convolution of two such chromosomes and then of
# the blocks. The variance of the dosage is 2 k raf (1 - raf) (1 + (m - 1) rho) for blocks of m.
dosage_pmf_linked <- function(n_loci, p, block_size, rho) {
  w <- sqrt(rho)
  blocks <- split(seq_len(n_loci), ceiling(seq_len(n_loci) / block_size))
  one <- function(m) { ch <- p * stats::dbinom(0:m, m, w + (1 - w) * p) + (1 - p) * stats::dbinom(0:m, m, (1 - w) * p); convolve_pmf(ch, ch) }
  Reduce(convolve_pmf, lapply(blocks, function(ix) one(length(ix))))
}

# The realised score when the causal loci are linked. The true effect of each locus is the causal (conditional) one. Two ways of estimating the weights:
#   marginal  each locus estimated on its own, as in a single-SNP discovery study: the estimate carries the effect of its linked neighbours, centred on
#             (1 + (m - 1) rho) times the causal effect, with errors correlated rho within a block (covariance sigma^2 times the genotype correlation matrix)
#   joint     all loci of a block estimated together: centred on the causal effect, with error covariance sigma^2 times the inverse of that matrix (the variance of
#             each weight is inflated by its variance inflation factor)
# sigma is the standard error of a single locus's marginal estimate. gain_exact is the gain with the causal weights known (the ceiling under this linkage).
# Profile-preserving truth only. Returns one row per replicate, estimation, score and threshold.
pathogen_cell_linked <- function(q, lp_profile, n_loci, odds_ratio, raf, f, sigma, thresholds, block_size = 4L, rho = 0.3, R = 12L, N = 2000L, seed = 1L) {
  if (genetic_scheme() != "profile") stop("pathogen_cell_linked() builds the profile-preserving truth; set PROSE_SCHEME = \"profile\"", call. = FALSE)
  set.seed(seed)
  k <- n_loci; beta <- log(odds_ratio); np <- length(q)
  pmf <- dosage_pmf_linked(k, raf, block_size, rho); ng <- length(pmf); d <- seq_len(ng) - 1
  w <- as.numeric(outer(q, pmf)); lp_cell <- rep(lp_profile, ng); D0 <- rep(d, each = np)
  clin <- stats::plogis(lp_cell)
  comp <- stats::plogis(outer(lp_profile + profile_shifts(lp_profile, beta * d, pmf), beta * d, "+"))
  r_cells <- (1 - f) * clin + f * as.numeric(comp)
  r_mat <- matrix(r_cells, np, ng)
  centres <- c(specific = beta, all_mould = diluted_beta(lp_cell, D0, w, r_cells))
  base_nb <- nb_from(tp_fp(w, r_cells, clin, thresholds), thresholds)
  exact_gain <- lapply(centres, function(b) { pr <- as.numeric(stats::plogis(outer(lp_profile + profile_shifts(lp_profile, b * d, pmf), b * d, "+"))); nb_from(tp_fp(w, r_cells, pr, thresholds), thresholds) - base_nb })
  blocks <- split(seq_len(k), ceiling(seq_len(k) / block_size)); wsh <- sqrt(rho)
  sig <- lapply(blocks, function(ix) { m <- length(ix); (1 - rho) * diag(m) + rho * matrix(1, m, m) })
  chol_marg <- lapply(sig, chol); chol_joint <- lapply(sig, function(s) chol(solve(s)))
  mult <- unlist(lapply(blocks, function(ix) rep(1 + (length(ix) - 1) * rho, length(ix))))
  wt <- rep(q, times = N) / N
  out <- lapply(seq_len(R), function(r) {
    G <- matrix(0L, N, k)
    for (ix in blocks) for (ch in 1:2) {
      u <- stats::rbinom(N, 1L, raf)
      copy <- matrix(stats::runif(N * length(ix)) < wsh, N, length(ix))
      G[, ix] <- G[, ix] + ifelse(copy, u, matrix(stats::rbinom(N * length(ix), 1L, raf), N, length(ix)))
    }
    rt <- as.numeric(r_mat[, rowSums(G) + 1L])
    z <- stats::rnorm(k)
    err <- list(marginal = unlist(Map(function(R_, ix) sigma * drop(t(R_) %*% z[ix]), chol_marg, blocks)),
                joint = unlist(Map(function(R_, ix) sigma * drop(t(R_) %*% z[ix]), chol_joint, blocks)))
    do.call(rbind, lapply(names(err), function(est) do.call(rbind, lapply(names(centres), function(sc) {
      wv <- (if (est == "marginal") centres[[sc]] * mult else centres[[sc]]) + err[[est]]
      S <- drop(G %*% wv)
      pr <- as.numeric(stats::plogis(outer(lp_profile + profile_shifts(lp_profile, S, rep(1 / N, N)), S, "+")))
      data.frame(replicate = r, estimation = est, score = sc, threshold = thresholds, gain = nb_from(tp_fp(wt, rt, pr, thresholds), thresholds) - base_nb,
                 gain_exact = exact_gain[[sc]])
    }))))
  })
  cbind(n_loci = n_loci, odds_ratio = odds_ratio, raf = raf, case_fraction = f, sigma = sigma, block_size = block_size, rho = rho, do.call(rbind, out))
}
