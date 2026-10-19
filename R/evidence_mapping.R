# Evidence-to-parameter mapping for the genetic layer (addendum 7). The target parameter for PTX3 is the log odds ratio
# per copy of the h2 haplotype class, recipient genotype, on the odds-ratio scale.
#
# h2 DEFINITION. PTX3 has three loci in strong linkage (rs2305619, rs3816527, rs1840680; forward strand). In 1000G EUR
# the common haplotypes are G-A-G (54%), A-C-A (40%) and A-A-A (6%). A chromosome is h2 if it carries the A allele at
# rs2305619; h1 is G-A-G. Because rs2305619 and rs1840680 are the same class (r2 = 0.996), the h2 dosage is also the
# number of A alleles at rs1840680. rs3816527 C marks only the A-C-A sub-haplotype (about 87% of h2), so it is not
# interchangeable and is mapped to h2 through the haplotype structure.
#
# Each study contrast (dominant, recessive, allelic, counts) is converted to this one scale by finding the per-copy
# h2 effect beta that reproduces the study's contrast under a haplotype model: additive on the log-odds scale per h2
# copy, Hardy-Weinberg haplotype pairs, and a case-control (rare outcome) sampling.
#
# Requires nothing outside base R and stats.

PTX3_LOCI <- c("rs2305619", "rs3816527", "rs1840680")

# The haplotype class whose per-copy effect is estimated is marked 1 in a column named `risk` (or `h2` for PTX3)
risk_of <- function(hap) if ("risk" %in% names(hap)) hap$risk else hap$h2

complement_base <- function(x) chartr("ACGT", "TGCA", x)

# Read one table from the evidence bundle without unpacking it
read_bundle_table <- function(name, zip = "IMD_genetics_analysis_bundle (1).zip") {
  # read.csv opens and closes an unopened connection itself, so it is not closed again here
  utils::read.csv(unz(zip, paste0("IMD_genetics_bundle/data/", name)), stringsAsFactors = FALSE, na.strings = c("", "NA"))
}

# PTX3 haplotypes and their frequencies; h2 = A at rs2305619. Kept for the earlier code and tests; block_haplotypes() below
# is the general form.
ptx3_haplotypes <- function(path = "data/haplotype_frequencies_phased.csv", source = c("gnomAD_NFE", "1000G_EUR")) {
  source <- match.arg(source)
  d <- utils::read.csv(path, stringsAsFactors = FALSE)
  d <- d[d$block == "PTX3_trio", ]
  col <- if (source == "gnomAD_NFE") "freq_gnomAD_NFE_rescaled" else "freq_1000G_EUR"
  h <- data.frame(haplotype = d$haplotype, a1 = d$allele_1, a2 = d$allele_2, a3 = d$allele_3, freq = d[[col]],
                  stringsAsFactors = FALSE)
  h <- h[h$freq > 0, , drop = FALSE]
  h$freq <- h$freq / sum(h$freq)
  h$h2 <- as.integer(h$a1 == "A")
  rownames(h) <- NULL
  h
}

# Fixed log-odds effects carried by a haplotype that is not the class being estimated (the second CLEC7A parameter, addendum 16): a column named
# `offset` in the haplotype table, zero when absent
offset_of <- function(hap) if ("offset" %in% names(hap)) hap$offset else rep(0, nrow(hap))

# All ordered pairs of haplotypes (Hardy-Weinberg): probability, genotype at one locus (alleles sorted), risk-class copies, fixed offset of the pair
genotype_pairs <- function(hap, locus) {
  al <- hap[[paste0("a", locus)]]
  ij <- expand.grid(i = seq_len(nrow(hap)), j = seq_len(nrow(hap)))
  off <- offset_of(hap)
  data.frame(prob = hap$freq[ij$i] * hap$freq[ij$j],
             genotype = vapply(seq_len(nrow(ij)), function(r) paste(sort(c(al[ij$i[r]], al[ij$j[r]])), collapse = ""), ""),
             k = risk_of(hap)[ij$i] + risk_of(hap)[ij$j], o = off[ij$i] + off[ij$j], stringsAsFactors = FALSE)
}

# Log OR of one genotype group against another, under risk-class effect beta (and any fixed offsets). Cases are weighted by exp(effect of the
# pair); the null is no effect at all.
implied_log_or_genotype <- function(beta, pairs, numerator, reference) {
  o <- if ("o" %in% names(pairs)) pairs$o else rep(0, nrow(pairs))
  f <- function(b, offsets = TRUE) {
    n <- pairs$genotype %in% numerator; r <- pairs$genotype %in% reference
    oo <- if (offsets) o else rep(0, length(o))
    log(sum(pairs$prob[n] * exp(b * pairs$k[n] + oo[n])) / sum(pairs$prob[r] * exp(b * pairs$k[r] + oo[r])))
  }
  f(beta) - f(0, offsets = FALSE)
}

# Log OR per copy of an allele, chromosome level (exact for an allelic test under multiplicative risk)
implied_log_or_allele <- function(beta, hap, locus, allele) {
  has <- hap[[paste0("a", locus)]] == allele
  off <- offset_of(hap)
  f <- function(b, offsets = TRUE) { w <- hap$freq * exp(b * risk_of(hap) + if (offsets) off else 0); log(sum(w[has]) / sum(w[!has])) }
  f(beta) - f(0, offsets = FALSE)
}

# Solve implied(beta) = observed for beta, if implied() is strictly monotone over the search range
solve_beta <- function(observed, implied, range = c(-6, 6)) {
  grid <- seq(range[1], range[2], length.out = 41)
  v <- vapply(grid, implied, numeric(1))
  if (any(!is.finite(v))) return(list(beta = NA_real_, slope = NA_real_, flag = "a genotype group has zero frequency under the haplotype model"))
  d <- diff(v)
  if (!(all(d > 0) || all(d < 0))) return(list(beta = NA_real_, slope = NA_real_, flag = "implied effect not monotone in beta"))
  if (observed < min(v) || observed > max(v)) return(list(beta = NA_real_, slope = NA_real_, flag = "observed effect outside the range the model can produce"))
  b <- stats::uniroot(function(x) implied(x) - observed, range, tol = 1e-10)$root
  h <- 1e-4
  list(beta = b, slope = (implied(b + h) - implied(b - h)) / (2 * h), flag = "")
}

# Genotype tokens in free text: "AG", "A/G", "G / A" -> sorted two-letter strings
parse_genotype_tokens <- function(txt) {
  m <- gregexpr("[ACGT]\\s*/\\s*[ACGT]|(?<![A-Za-z])[ACGT]{2}(?![A-Za-z])", txt, perl = TRUE)
  tok <- regmatches(txt, m)[[1]]
  unique(vapply(tok, function(x) paste(sort(strsplit(gsub("[^ACGT]", "", x), "")[[1]]), collapse = ""), ""))
}

# Which genotype group is compared with which? Read the comparison text; if it cannot be read, use the model label.
parse_contrast <- function(comparison, genetic_model, effect_allele, other_allele, flip_text = FALSE) {
  # The alleles are given on the forward strand. If the study wrote its text on the opposite strand, flip_text complements
  # the genotypes read from the comparison text.
  e <- effect_allele; o <- other_allele
  tok <- function(txt) { t <- parse_genotype_tokens(txt); if (flip_text) unique(vapply(t, function(g) paste(sort(strsplit(complement_base(g), "")[[1]]), collapse = ""), "")) else t }
  valid <- if (!is.na(e) && !is.na(o)) c(paste0(e, e), paste(sort(c(e, o)), collapse = ""), paste0(o, o)) else character()
  parts <- if (!is.na(comparison)) strsplit(comparison, "(?i)\\s+(vs\\.?|versus)\\s+", perl = TRUE)[[1]] else character()
  if (length(parts) == 2 && length(valid)) {
    L <- tok(parts[1]); R <- tok(parts[2])
    if (length(L) && length(R) && all(c(L, R) %in% valid) && !length(intersect(L, R))) {
      left_ref <- grepl("ref", parts[1], ignore.case = TRUE) && !grepl("ref", parts[2], ignore.case = TRUE)
      return(list(type = "genotype", numerator = if (left_ref) R else L, reference = if (left_ref) L else R, method = "text", flag = ""))
    }
  }
  if (length(parts) > 2) return(list(type = NA, numerator = NULL, reference = NULL, method = "unreadable", flag = "more than two groups"))
  if (!length(valid)) return(list(type = NA, numerator = NULL, reference = NULL, method = "unreadable", flag = "effect allele not stated"))
  het <- paste(sort(c(e, o)), collapse = ""); ee <- paste0(e, e); oo <- paste0(o, o)
  switch(as.character(genetic_model),
         dominant = , carrier = list(type = "genotype", numerator = c(ee, het), reference = oo, method = "model_label", flag = "contrast taken from the model label"),
         recessive = list(type = "genotype", numerator = ee, reference = c(het, oo), method = "model_label", flag = "contrast taken from the model label"),
         allelic = list(type = "allele", allele = e, method = "model_label", flag = "contrast taken from the model label"),
         additive_per_allele = list(type = "allele", allele = e, method = "model_label", flag = "contrast taken from the model label"),
         list(type = NA, numerator = NULL, reference = NULL, method = "unreadable", flag = paste("model", genetic_model, "has no single contrast")))
}

# The implied log OR function for a parsed contrast at one locus
contrast_function <- function(contrast, hap, locus) {
  if (identical(contrast$type, "genotype")) {
    pairs <- genotype_pairs(hap, locus)
    return(function(b) implied_log_or_genotype(b, pairs, contrast$numerator, contrast$reference))
  }
  if (identical(contrast$type, "allele")) return(function(b) implied_log_or_allele(b, hap, locus, contrast$allele))
  NULL
}

# From an observed log OR with its standard error to the per-copy h2 effect, by the delta method
map_to_h2 <- function(log_or, se, contrast, hap, locus) {
  fn <- contrast_function(contrast, hap, locus)
  if (is.null(fn)) return(data.frame(beta = NA_real_, se_beta = NA_real_, flag = contrast$flag))
  s <- solve_beta(log_or, fn)
  data.frame(beta = s$beta, se_beta = if (is.na(s$beta)) NA_real_ else se / abs(s$slope),
             flag = paste(c(contrast$flag, s$flag)[nzchar(c(contrast$flag, s$flag))], collapse = "; "))
}

# Log OR and standard error from a published estimate and its 95% interval
log_or_from_estimate <- function(est, lo, hi) {
  est <- as.numeric(est); lo <- as.numeric(lo); hi <- as.numeric(hi)
  if (anyNA(c(est, lo, hi)) || est <= 0 || lo <= 0 || hi <= lo) return(list(y = NA_real_, se = NA_real_, flag = "estimate or interval missing or not positive"))
  if (est < lo - 1e-9 || est > hi + 1e-9 || est >= hi - 1e-9 || est <= lo + 1e-9) return(list(y = NA_real_, se = NA_real_, flag = "estimate not strictly inside its interval"))
  list(y = log(est), se = (log(hi) - log(lo)) / (2 * stats::qnorm(0.975)), flag = "")
}

# Additive-trend log OR per effect allele from a 2x3 table (cases and controls by effect-allele copies 2, 1, 0)
log_or_from_counts_2x3 <- function(cases, controls) {
  d <- data.frame(g = c(2, 1, 0), cases = as.numeric(cases), controls = as.numeric(controls))
  d <- d[d$cases + d$controls > 0, ]
  if (nrow(d) < 2 || sum(d$cases) == 0 || sum(d$controls) == 0) return(list(y = NA_real_, se = NA_real_, flag = "table cannot give a trend"))
  fit <- suppressWarnings(stats::glm(cbind(cases, controls) ~ g, family = stats::binomial, data = d))
  s <- summary(fit)$coefficients
  list(y = unname(s["g", 1]), se = unname(s["g", 2]), flag = "")
}

# Log OR and standard error from a 2x2 table: events and non-events in the group against the rest (0.5 added if any cell is zero)
log_or_from_counts_2x2 <- function(events_group, events_rest, nonevents_group, nonevents_rest) {
  v <- c(events_group, events_rest, nonevents_group, nonevents_rest)
  if (anyNA(v)) return(list(y = NA_real_, se = NA_real_, flag = "incomplete 2x2 table"))
  if (any(v == 0)) v <- v + 0.5
  list(y = log((v[1] / v[3]) / (v[2] / v[4])), se = sqrt(sum(1 / v)), flag = "")
}

# Fixed-effect inverse-variance combination of estimates from independent cohorts of one study
combine_fixed <- function(y, se) {
  w <- 1 / se^2
  list(y = sum(w * y) / sum(w), se = sqrt(1 / sum(w)), q = sum(w * (y - sum(w * y) / sum(w))^2), k = length(y))
}

# Posterior of mu given tau for the normal-normal model, with a N(0, prior_sd^2) prior on mu
mu_given_tau <- function(y, se, tau, prior_sd) {
  w <- 1 / (se^2 + tau^2)
  prec <- sum(w) + 1 / prior_sd^2
  list(mean = sum(w * y) / prec, sd = sqrt(1 / prec))
}

# Bayesian random-effects pooling on a grid for tau, with mu integrated out analytically.
#   y_i ~ N(theta_i, se_i^2), theta_i ~ N(mu, tau^2), mu ~ N(0, prior_sd_mu^2), tau ~ half-normal(prior_scale_tau)
bayes_meta_grid <- function(y, se, prior_sd_mu = 1, prior_scale_tau = 0.5, n_draw = 4000L, seed = 1L, tau_max = 3, n_grid = 1201L) {
  stopifnot(length(y) == length(se), length(y) >= 1, all(se > 0))
  tau <- seq(0, tau_max, length.out = n_grid)
  v <- prior_sd_mu^2
  logpost <- vapply(tau, function(t) {
    w <- 1 / (se^2 + t^2)
    a <- sum(w); b <- sum(w * y); q <- sum(w * y^2)
    # marginal likelihood after integrating mu (matrix determinant lemma and Sherman-Morrison)
    -0.5 * sum(log(se^2 + t^2)) - 0.5 * log(1 + v * a) - 0.5 * (q - v * b^2 / (1 + v * a)) +
      stats::dnorm(t, 0, prior_scale_tau, log = TRUE)
  }, numeric(1))
  p <- exp(logpost - max(logpost)); p <- p / sum(p)
  set.seed(seed)
  tdraw <- sample(tau, n_draw, replace = TRUE, prob = p) + stats::runif(n_draw, -0.5, 0.5) * (tau[2] - tau[1])
  tdraw <- pmax(tdraw, 0)
  mu_par <- lapply(tdraw, function(t) mu_given_tau(y, se, t, prior_sd_mu))
  mu <- vapply(mu_par, function(m) stats::rnorm(1, m$mean, m$sd), numeric(1))
  pred <- stats::rnorm(n_draw, mu, tdraw)
  q <- function(x) stats::quantile(x, c(0.05, 0.5, 0.95), names = FALSE)
  list(draws = data.frame(mu = mu, tau = tdraw), mu = q(mu), tau = q(tdraw), prediction = q(pred),
       prob_mu_below_zero = mean(mu < 0), k = length(y))
}

# Haplotypes and frequencies of any block in the phased table, with the risk class marked.
#   locus_alleles: the allele at each locus that defines the risk class, e.g. c(A = "A") for PTX3 (rs2305619 A) or
#   c("C", NA) for the CLEC7A Y238X allele; NA means the locus is not part of the rule.
block_haplotypes <- function(block, risk_alleles, path = "data/haplotype_frequencies_phased.csv", source = c("gnomAD_NFE", "1000G_EUR")) {
  source <- match.arg(source)
  d <- utils::read.csv(path, stringsAsFactors = FALSE)
  d <- d[d$block == block, , drop = FALSE]
  if (!nrow(d)) stop("block not found: ", block, call. = FALSE)
  col <- if (source == "gnomAD_NFE") "freq_gnomAD_NFE_rescaled" else "freq_1000G_EUR"
  n_loci <- sum(grepl("^allele_[0-9]+$", names(d)))
  h <- data.frame(haplotype = d$haplotype, freq = d[[col]], stringsAsFactors = FALSE)
  for (j in seq_len(n_loci)) h[[paste0("a", j)]] <- d[[paste0("allele_", j)]]
  h <- h[h$freq > 0, , drop = FALSE]
  h$freq <- h$freq / sum(h$freq)
  inrule <- Reduce(`&`, lapply(seq_len(n_loci), function(j) is.na(risk_alleles[j]) | h[[paste0("a", j)]] == risk_alleles[j]))
  h$risk <- as.integer(inrule)
  rownames(h) <- NULL
  h
}
