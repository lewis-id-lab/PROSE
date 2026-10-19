# PTX3 pilot of the evidence-to-parameter mapping (addendum 7). Run from the project root.
#
# Reads    the extraction tables inside IMD_genetics_analysis_bundle (1).zip, and data/haplotype_frequencies_phased.csv
# Writes   output/ptx3_row_effects.csv           every usable PTX3 row, mapped to a per-copy h2 effect where possible
#          output/ptx3_study_effects.csv         the primary analysis: one effect per study
#          output/ptx3_h2_pooled.csv             pooled per-copy h2 effect, primary analysis and sensitivity analyses
#          output/ptx3_snp_crosscheck.csv        the same study mapped through rs2305619 and through rs3816527
#          output/ptx3_h2_posterior_draws.csv    draws of (mu, tau) from the primary analysis
#
# The extraction rows are round-1. Corrections found by checking rows against their sources (scripts/13) are applied first;
# everything else is unverified, so every number here is provisional.

source("R/evidence_mapping.R")
source("R/ptx3_evidence.R")
source("R/source_verification.R")

vr <- read_bundle_table("extraction_variants_round1.csv")
st <- read_bundle_table("extraction_studies_round1.csv")
corrections <- if (file.exists("data/extraction_corrections.csv")) read.csv("data/extraction_corrections.csv", stringsAsFactors = FALSE) else NULL
applied <- apply_corrections(vr, corrections)
vr <- applied$data
hap <- block_haplotypes("PTX3_trio", risk_alleles = c("A", NA, NA), source = "gnomAD_NFE")
re <- ptx3_row_effects(vr, st, hap, loci = PTX3_LOCI, counts_model = "trend", corrected_ids = applied$applied)

scenarios <- list(
  "primary: recipient, OR scale, whole cohort, without Herrero-Sanchez" = ptx3_study_effects(re),
  "with Herrero-Sanchez 2018 (S035), detection endpoint" = ptx3_study_effects(re, exclude_studies = character()),
  "with Herrero-Sanchez 2018 (S035), invasive aspergillosis endpoint" = ptx3_study_effects(use_s035_ia_endpoint(re), exclude_studies = character()),
  "without Brunel 2018 (S029)" = ptx3_study_effects(re, exclude_studies = c("S035", "S029")),
  "European-ancestry cohorts only" = ptx3_study_effects(re, allow_non_european = FALSE),
  "rs3816527 as the tag SNP, mapped to h2" = ptx3_study_effects(re, snp_order = "rs3816527"),
  "with hazard-ratio studies" = ptx3_study_effects(re, allow_hr = TRUE)
)

# Donor genotype is the pre-specified secondary analysis (D1, allo-HSCT only); one study, so it is described, not pooled
donor <- ptx3_study_effects(re, genotype = "donor", allow_hr = TRUE, exclude_studies = character())

pool <- function(effects, seed) {
  if (nrow(effects) < 1) return(NULL)
  b <- bayes_meta_grid(effects$beta, effects$se_beta, prior_sd_mu = 1, prior_scale_tau = 0.5, seed = seed)
  list(summary = data.frame(k = nrow(effects), studies = paste(effects$study_id, collapse = " "),
                            mu_lo = b$mu[1], mu = b$mu[2], mu_hi = b$mu[3],
                            or = exp(b$mu[2]), or_lo = exp(b$mu[1]), or_hi = exp(b$mu[3]),
                            tau_lo = b$tau[1], tau = b$tau[2], tau_hi = b$tau[3],
                            pred_lo = b$prediction[1], pred_hi = b$prediction[3], prob_below_zero = b$prob_mu_below_zero),
       draws = b$draws)
}
pooled <- Map(pool, scenarios, seq_along(scenarios) + 5517L)
pooled_tab <- do.call(rbind, Map(function(p, nm) cbind(scenario = nm, p$summary), pooled[!vapply(pooled, is.null, logical(1))], names(scenarios)[!vapply(pooled, is.null, logical(1))]))

primary <- scenarios[[1]]
study_tab <- primary[, c("study_id", "first_author", "year", "rsid", "cohort", "tier", "model", "measure", "method", "y", "se", "beta", "se_beta", "flag", "ancestry_as_reported")]
study_tab$or_per_h2_copy <- exp(study_tab$beta)
study_tab$ancestry_as_reported <- substr(study_tab$ancestry_as_reported, 1, 40)

# Cross-check: the same study through the h2 tag SNP and through rs3816527 (should agree if h2 explains both)
tag <- ptx3_study_effects(re, exclude_studies = character())
alt <- ptx3_study_effects(re, snp_order = "rs3816527", exclude_studies = character())
cross <- merge(tag[, c("study_id", "rsid", "beta", "se_beta")], alt[, c("study_id", "rsid", "beta", "se_beta")], by = "study_id", suffixes = c("_tag", "_rs3816527"))
cross$difference <- cross$beta_rs3816527 - cross$beta_tag

dir.create("output", showWarnings = FALSE)
write.csv(re[, setdiff(names(re), c("ancestry_as_reported"))], "output/ptx3_row_effects.csv", row.names = FALSE)
write.csv(study_tab, "output/ptx3_study_effects.csv", row.names = FALSE)
write.csv(pooled_tab, "output/ptx3_h2_pooled.csv", row.names = FALSE)
write.csv(cross, "output/ptx3_snp_crosscheck.csv", row.names = FALSE)
if (nrow(donor)) {
  donor$or_per_h2_copy <- exp(donor$beta)
  write.csv(donor[, c("study_id", "first_author", "year", "rsid", "tier", "model", "measure", "method", "y", "se", "beta", "se_beta", "or_per_h2_copy", "flag")],
            "output/ptx3_donor_effects.csv", row.names = FALSE)
}
write.csv(pooled[[1]]$draws, "output/ptx3_h2_posterior_draws.csv", row.names = FALSE)
# Donor draws for the allo-HSCT stratum analysis (one study: its estimate with the prior's shrinkage)
donor_pooled <- pool(donor, 5617L)
if (!is.null(donor_pooled)) write.csv(donor_pooled$draws, "output/ptx3_h2_donor_posterior_draws.csv", row.names = FALSE)

message("haplotype frequencies used (gnomAD NFE): h2 = ", sprintf("%.3f", sum(hap$freq[hap$risk == 1])), " | rows corrected: ", length(applied$applied))
message("usable PTX3 rows: ", nrow(re), " | mapped to a per-copy h2 effect: ", sum(!is.na(re$beta)))
message("\nPrimary analysis, one effect per study:")
print(cbind(study_tab[, c("study_id", "rsid", "tier", "model")], round(study_tab[, c("y", "se", "beta", "se_beta", "or_per_h2_copy")], 3)), row.names = FALSE)
message("\nPooled per-copy h2 effect (log OR median and 90% interval; OR; tau; prediction interval):")
print(cbind(pooled_tab[, c("scenario", "k")], round(pooled_tab[, c("mu_lo", "mu", "mu_hi", "or", "or_lo", "or_hi", "tau", "pred_lo", "pred_hi", "prob_below_zero")], 3)), row.names = FALSE)
message("\nCross-check, rs2305619/rs1840680 against rs3816527:")
print(cbind(cross[, c("study_id", "rsid_tag")], round(cross[, c("beta_tag", "beta_rs3816527", "difference")], 3)), row.names = FALSE)
