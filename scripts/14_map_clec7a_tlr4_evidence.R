# Evidence-to-parameter mapping for CLEC7A and TLR4 (addendum 8). Run from the project root, after scripts/13.
#
# CLEC7A: the term is the Y238X haplotype. Y238X (rs16910526, forward-strand C allele) sits entirely on the G background of
#   rs7309123, so the three haplotypes are C-G (Y238X, 7.7% in gnomAD NFE), A-G (34.9%) and A-C (57.4%). The risk class is
#   C-G. rs16910526 rows map to its per-copy effect directly. rs7309123 rows would identify a second parameter, the effect of
#   G on the non-Y238X background, but rs7309123 is a C/G palindrome and no row records a control allele frequency, so its
#   strand cannot be resolved and its rows are not used. That second parameter is held at zero.
# TLR4: the term is the G-T haplotype of rs4986790 (G) and rs4986791 (T), 6.1% in gnomAD NFE. Both SNPs are used through
#   carriers against non-carriers with 0.5 added to a zero cell (decision D3c), because no cohort has minor homozygotes.
#
# Writes   output/{clec7a,tlr4}_row_effects.csv, _study_effects.csv, _pooled.csv, _posterior_draws.csv, _snp_crosscheck.csv
# The extraction rows are round 1 and mostly unverified, so every number is provisional.

source("R/evidence_mapping.R")
source("R/ptx3_evidence.R")
source("R/source_verification.R")

vr <- read_bundle_table("extraction_variants_round1.csv")
st <- read_bundle_table("extraction_studies_round1.csv")
corrections <- if (file.exists("data/extraction_corrections.csv")) read.csv("data/extraction_corrections.csv", stringsAsFactors = FALSE) else NULL
applied <- apply_corrections(vr, corrections)
vr <- applied$data

run_term <- function(prefix, block, risk_alleles, loci, snp_order, alt_snp, scenarios_extra = list(), seed) {
  hap <- block_haplotypes(block, risk_alleles = risk_alleles, source = "gnomAD_NFE")
  re <- ptx3_row_effects(vr, st, hap, loci = loci, counts_model = "carrier", corrected_ids = applied$applied)
  scenarios <- c(list(
    "primary: recipient, OR scale, whole cohort, without Herrero-Sanchez" = ptx3_study_effects(re, snp_order = snp_order),
    "with Herrero-Sanchez 2018 (S035)" = ptx3_study_effects(re, snp_order = snp_order, exclude_studies = character()),
    "with hazard-ratio studies" = ptx3_study_effects(re, snp_order = snp_order, allow_hr = TRUE),
    "European-ancestry cohorts only" = ptx3_study_effects(re, snp_order = snp_order, allow_non_european = FALSE)),
    if (!is.null(alt_snp)) list("second SNP as the tag, mapped to the risk class" = ptx3_study_effects(re, snp_order = alt_snp)),
    list("donor genotype (secondary analysis)" = ptx3_study_effects(re, snp_order = snp_order, genotype = "donor", allow_hr = TRUE, exclude_studies = character())),
    lapply(scenarios_extra, function(f) f(re)))
  pooled <- Map(pool_study_effects, scenarios, seq_along(scenarios) + seed)
  keep <- !vapply(pooled, is.null, logical(1))
  tab <- do.call(rbind, Map(function(p, nm) cbind(scenario = nm, p$summary), pooled[keep], names(scenarios)[keep]))
  primary <- scenarios[[1]]
  study_tab <- if (nrow(primary)) {
    x <- primary[, c("study_id", "first_author", "year", "rsid", "cohort", "tier", "model", "measure", "method", "y", "se", "beta", "se_beta", "flag", "ancestry_as_reported")]
    x$or_per_copy <- exp(x$beta); x$ancestry_as_reported <- substr(x$ancestry_as_reported, 1, 40); x
  } else data.frame()
  cross <- data.frame()
  if (!is.null(alt_snp)) {
    tag <- ptx3_study_effects(re, snp_order = snp_order, exclude_studies = character())
    alt <- ptx3_study_effects(re, snp_order = alt_snp, exclude_studies = character())
    cross <- merge(tag[, c("study_id", "rsid", "beta", "se_beta")], alt[, c("study_id", "rsid", "beta", "se_beta")], by = "study_id", suffixes = c("_first", "_second"))
    if (nrow(cross)) cross$difference <- cross$beta_second - cross$beta_first
  }
  dir.create("output", showWarnings = FALSE)
  write.csv(re[, setdiff(names(re), "ancestry_as_reported")], paste0("output/", prefix, "_row_effects.csv"), row.names = FALSE)
  write.csv(study_tab, paste0("output/", prefix, "_study_effects.csv"), row.names = FALSE)
  write.csv(tab, paste0("output/", prefix, "_pooled.csv"), row.names = FALSE)
  if (nrow(cross)) write.csv(cross, paste0("output/", prefix, "_snp_crosscheck.csv"), row.names = FALSE)
  if (!is.null(pooled[[1]])) write.csv(pooled[[1]]$draws, paste0("output/", prefix, "_posterior_draws.csv"), row.names = FALSE)
  list(hap = hap, re = re, study = study_tab, pooled = tab, cross = cross, scenarios = scenarios)
}

clec <- run_term("clec7a", "CLEC7A", c("C", NA), c("rs16910526", "rs7309123"), snp_order = "rs16910526", alt_snp = NULL, seed = 7241L)
tlr4 <- run_term("tlr4", "TLR4", c("G", "T"), c("rs4986790", "rs4986791"), snp_order = c("rs4986790", "rs4986791"), alt_snp = "rs4986791", seed = 3086L,
                 scenarios_extra = list("with de Boer 2011 (S061, genotype source unstated, treated as recipient)" = function(re) {
                   re$genotype[re$study_id == "S061"] <- "recipient"
                   ptx3_study_effects(re, snp_order = c("rs4986790", "rs4986791"))
                 }))

show <- function(x, label) {
  message("\n=====", label, "=====")
  message("risk class frequency (gnomAD NFE): ", sprintf("%.3f", sum(x$hap$freq[x$hap$risk == 1])), " | usable rows: ", nrow(x$re), " | mapped: ", sum(!is.na(x$re$beta)))
  if (nrow(x$study)) print(cbind(x$study[, c("study_id", "rsid", "tier", "model")], round(x$study[, c("y", "se", "beta", "se_beta", "or_per_copy")], 3)), row.names = FALSE)
  print(cbind(x$pooled[, c("scenario", "k")], round(x$pooled[, c("mu_lo", "mu", "mu_hi", "or", "or_lo", "or_hi", "tau", "prob_below_zero")], 3)), row.names = FALSE)
  if (nrow(x$cross)) { message("cross-check against the second SNP:"); print(cbind(x$cross[, c("study_id")], round(x$cross[, c("beta_first", "beta_second", "difference")], 3)), row.names = FALSE) }
}
show(clec, "CLEC7A (risk class: Y238X, haplotype C-G)")
show(tlr4, "TLR4 (risk class: G-T haplotype)")

r73 <- clec$re[clec$re$rsid == "rs7309123", ]
message("\nCLEC7A rs7309123: ", nrow(r73), " usable rows from ", length(unique(r73$study_id)), " studies are not used (strand unresolved, no control allele frequency).")
