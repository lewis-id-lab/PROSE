# Sensitivity of the pooled donor PTX3 effect to Cunha 2014 (S079). Run from the project root, after script 12 (addendum 37).
#
# S079 reports donor PTX3 results (the study that first described the h2/h2 haplotype) but its extraction rows carry no rsid and are flagged
# `rsid_unresolved`, so the pipeline leaves them out of the donor pooling. Here, in memory only (the extraction is not changed), its rows are
# given the rsids of the SNPs they report (+281A>G = rs2305619, +734A>C = rs3816527), a donor label and a forward-strand flag, and the donor
# effect is pooled with and without them. Two S079 results are tried: the discovery hazard ratio (V1698, +281 GG, the study that picked the
# SNP from 22 candidates, so optimistic) and the confirmation odds ratio (V1703, a multicentre case-control study).
#
# Writes   output/ptx3_donor_s079_sensitivity.csv (pooled donor effect per h2 copy for each set) and output/ptx3_donor_s079_effects.csv (the study effects)

source("R/baseline.R"); source("R/bs_data.R"); source("R/decision_curve.R"); source("R/boscore_data.R"); source("R/refit_baseline.R")
source("R/evidence_mapping.R"); source("R/ptx3_evidence.R"); source("R/source_verification.R"); source("R/genetic_layer.R"); source("R/genetic_helpers.R")

vr <- read_bundle_table("extraction_variants_round1.csv"); st <- read_bundle_table("extraction_studies_round1.csv")
vr <- apply_corrections(vr, read.csv("data/extraction_corrections.csv", stringsAsFactors = FALSE))$data
vr <- drop_duplicate_reports(vr)
vr <- add_extraction_rows(vr)
hap <- block_haplotypes("PTX3_trio", risk_alleles = c("A", NA, NA), source = "gnomAD_NFE")

s79 <- c(V1698 = "rs2305619", V1700 = "rs3816527", V1703 = "rs2305619", V1704 = "rs3816527")
i <- match(names(s79), vr$variant_row_id)
vr$rsid[i] <- unname(s79); vr$orientation_flag[i] <- "fwd_ok"; vr$subgroup[i] <- "donor genotype; whole cohort"
re <- ptx3_row_effects(vr, st, hap, loci = PTX3_LOCI, counts_model = "trend")

pool <- function(effects, seed) {
  b <- bayes_meta_grid(effects$beta, effects$se_beta, prior_sd_mu = 1, prior_scale_tau = 0.5, seed = seed)
  data.frame(k = nrow(effects), studies = paste(effects$study_id, collapse = " "), or = exp(b$mu[2]), or_lo = exp(b$mu[1]), or_hi = exp(b$mu[3]), tau = b$tau[2])
}
donor_with <- function(keep79) ptx3_study_effects(re[re$study_id != "S079" | re$variant_row_id %in% keep79, ], genotype = "donor", allow_hr = TRUE, exclude_studies = character())
sets <- list("S080 and the Leiden letter (the pipeline)" = character(), "plus S079 discovery hazard ratio (V1698)" = "V1698", "plus S079 confirmation odds ratio (V1703)" = "V1703")
eff <- Map(donor_with, sets)
tab <- do.call(rbind, Map(function(e, nm) cbind(set = nm, pool(e, 5617L)), eff, names(sets)))
rownames(tab) <- NULL
effects <- do.call(rbind, Map(function(e, nm) cbind(set = nm, e[, c("study_id", "rsid", "model", "measure", "beta", "se_beta")]), eff, names(sets)))
effects$or_per_h2_copy <- exp(effects$beta)

write.csv(tab, "output/ptx3_donor_s079_sensitivity.csv", row.names = FALSE)
write.csv(effects, "output/ptx3_donor_s079_effects.csv", row.names = FALSE)
message("Pooled donor effect per h2 copy (OR, 90% interval):"); show_df(tab)
message("\nStudy effects:"); show_df(effects)
