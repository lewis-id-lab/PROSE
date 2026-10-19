# Run the evidence-to-parameter mapping for one genetic term and write its outputs (addenda 8 and 9).
# Requires R/evidence_mapping.R, R/ptx3_evidence.R and R/source_verification.R.

# Haplotype table for a single SNP, in the form block_haplotypes() returns, from the panel allele frequencies.
#   risk_allele: forward-strand allele whose per-copy effect is estimated.
single_snp_haplotypes <- function(rsid, risk_allele, path = "data/panel_allele_frequencies.csv", source = c("gnomAD_NFE", "1000G_EUR")) {
  source <- match.arg(source)
  d <- utils::read.csv(path, stringsAsFactors = FALSE)
  d <- d[d$rsid == rsid, , drop = FALSE]
  if (nrow(d) != 1) stop("SNP not found in ", path, ": ", rsid, call. = FALSE)
  f <- unlist(jsonlite::fromJSON(d[[if (source == "gnomAD_NFE") "gnomad_nfe_freqs" else "tg_eur_freqs"]]))
  f <- f[f > 0]
  if (!risk_allele %in% names(f)) stop(risk_allele, " is not an observed allele of ", rsid, call. = FALSE)
  h <- data.frame(haplotype = names(f), freq = as.numeric(f) / sum(f), a1 = names(f), stringsAsFactors = FALSE)
  h$risk <- as.integer(h$a1 == risk_allele)
  rownames(h) <- NULL
  h
}

# Mapping for one term.
#   hap, loci, snp_order, alt_snp: as in scripts 12 and 14. `alt_snp` (or NULL) is the cross-check SNP.
#   hr_fallback_min: if the odds-ratio-only primary analysis has fewer studies than this, hazard-ratio studies are added to
#     the primary analysis and the scenario is labelled as such (decision for CD209 and IFNG, addendum 9). NULL = never.
run_term <- function(prefix, hap, loci, snp_order, alt_snp = NULL, seed, vr, st, applied_ids = character(),
                     scenarios_extra = list(), hr_fallback_min = NULL, counts_model = "carrier", out_dir = "output") {
  re <- ptx3_row_effects(vr, st, hap, loci = loci, counts_model = counts_model, corrected_ids = applied_ids)
  scenarios <- c(list(
    "primary: recipient, OR scale, whole cohort, without Herrero-Sanchez" = ptx3_study_effects(re, snp_order = snp_order),
    "with Herrero-Sanchez 2018 (S035)" = ptx3_study_effects(re, snp_order = snp_order, exclude_studies = character()),
    "with hazard-ratio studies" = ptx3_study_effects(re, snp_order = snp_order, allow_hr = TRUE),
    "European-ancestry cohorts only" = ptx3_study_effects(re, snp_order = snp_order, allow_non_european = FALSE)),
    if (!is.null(alt_snp)) list("second SNP as the tag, mapped to the risk class" = ptx3_study_effects(re, snp_order = alt_snp)),
    list("donor genotype (secondary analysis)" = ptx3_study_effects(re, snp_order = snp_order, genotype = "donor", allow_hr = TRUE, exclude_studies = character())),
    lapply(scenarios_extra, function(f) f(re)))
  if (!is.null(hr_fallback_min) && nrow(scenarios[[1]]) < hr_fallback_min) {
    label <- sprintf("primary (fallback): hazard ratios included, odds-ratio evidence has %d of %d studies needed", nrow(scenarios[[1]]), hr_fallback_min)
    scenarios[[1]] <- scenarios[["with hazard-ratio studies"]]
    names(scenarios)[1] <- label
  }
  pooled <- Map(pool_study_effects, scenarios, seq_along(scenarios) + seed)
  keep <- !vapply(pooled, is.null, logical(1))
  tab <- if (any(keep)) do.call(rbind, Map(function(p, nm) cbind(scenario = nm, p$summary), pooled[keep], names(scenarios)[keep])) else data.frame()
  primary <- scenarios[[1]]
  study_tab <- if (nrow(primary)) {
    x <- primary[, c("study_id", "first_author", "year", "rsid", "cohort", "tier", "model", "measure", "method", "y", "se", "beta", "se_beta", "flag", "ancestry_as_reported")]
    x$or_per_copy <- exp(x$beta); x$ancestry_as_reported <- substr(x$ancestry_as_reported, 1, 40); x
  } else data.frame()
  cross <- data.frame()
  if (!is.null(alt_snp)) {
    tag <- ptx3_study_effects(re, snp_order = snp_order, exclude_studies = character())
    alt <- ptx3_study_effects(re, snp_order = alt_snp, exclude_studies = character())
    if (nrow(tag) && nrow(alt)) {
      cross <- merge(tag[, c("study_id", "rsid", "beta", "se_beta")], alt[, c("study_id", "rsid", "beta", "se_beta")], by = "study_id", suffixes = c("_first", "_second"))
      if (nrow(cross)) cross$difference <- cross$beta_second - cross$beta_first
    }
  }
  dir.create(out_dir, showWarnings = FALSE)
  write.csv(re[, setdiff(names(re), "ancestry_as_reported")], file.path(out_dir, paste0(prefix, "_row_effects.csv")), row.names = FALSE)
  write.csv(study_tab, file.path(out_dir, paste0(prefix, "_study_effects.csv")), row.names = FALSE)
  write.csv(tab, file.path(out_dir, paste0(prefix, "_pooled.csv")), row.names = FALSE)
  if (nrow(cross)) write.csv(cross, file.path(out_dir, paste0(prefix, "_snp_crosscheck.csv")), row.names = FALSE)
  if (!is.null(pooled[[1]])) write.csv(pooled[[1]]$draws, file.path(out_dir, paste0(prefix, "_posterior_draws.csv")), row.names = FALSE)
  # Donor genotype (the D1 secondary analysis): the pooled draws are kept for the allo-HSCT stratum analysis
  donor_draws <- pooled[["donor genotype (secondary analysis)"]]
  if (!is.null(donor_draws)) write.csv(donor_draws$draws, file.path(out_dir, paste0(prefix, "_donor_posterior_draws.csv")), row.names = FALSE)
  list(hap = hap, re = re, study = study_tab, pooled = tab, cross = cross, scenarios = scenarios)
}

# Console summary of one term's result
show_term <- function(x, label) {
  message("\n=====", label, "=====")
  message("risk class frequency (gnomAD NFE): ", sprintf("%.3f", sum(x$hap$freq[x$hap$risk == 1])), " | usable rows: ", nrow(x$re), " | mapped: ", sum(!is.na(x$re$beta)))
  if (nrow(x$study)) print(cbind(x$study[, c("study_id", "rsid", "tier", "model", "measure")], round(x$study[, c("y", "se", "beta", "se_beta", "or_per_copy")], 3)), row.names = FALSE)
  if (nrow(x$pooled)) print(cbind(x$pooled[, c("scenario", "k")], round(x$pooled[, c("mu_lo", "mu", "mu_hi", "or", "or_lo", "or_hi", "tau", "prob_below_zero")], 3)), row.names = FALSE)
  if (nrow(x$cross)) { message("cross-check against the second SNP:"); print(cbind(study = x$cross$study_id, round(x$cross[, c("beta_first", "beta_second", "difference")], 3)), row.names = FALSE) }
}
