# Evidence-to-parameter mapping for CD209, IL10 and IFNG (addendum 9). Run from the project root, after scripts/13.
#
# CD209: rs4804800 and rs7248637 give four haplotypes (gnomAD NFE: G-A 7.8%, A-A 2.6%, G-G 4.7%, A-G 84.9%). The risk class is
#   the chromosomes carrying G at rs4804800 (G-A and G-G, 12.5%), decided with the user on 19 October 2026. rs4804800 rows map
#   to it directly; rs7248637 rows (A allele, r2 = 0.41 with rs4804800) map through the haplotype structure and are used where
#   a study has no rs4804800 row. The rs4804800 and rs7248637 evidence is compared in the SNP cross-check.
# IL10 rs1800896 and IFNG rs2069705: optional single-SNP terms. The risk allele is the forward-strand allele that most studies
#   report as their effect allele (after complementing strand-flipped studies): IL10 C (reported as G), IFNG G (reported as C).
#   A study reporting the other allele is mapped to the same class with the sign handled by the haplotype model.
# Hazard-ratio rule (decided 19 October 2026, per term): odds-ratio evidence is the primary analysis unless it has fewer than
#   three studies, in which case hazard-ratio studies are added and the scenario is labelled as a fallback.
#
# Writes output/{cd209,il10,ifng}_row_effects.csv, _study_effects.csv, _pooled.csv, _posterior_draws.csv (and _snp_crosscheck.csv for CD209).
# The extraction rows carry the round-1 values plus the recorded corrections, so every number is provisional.

source("R/evidence_mapping.R")
source("R/ptx3_evidence.R")
source("R/source_verification.R")
source("R/term_mapping.R")

vr <- read_bundle_table("extraction_variants_round1.csv")
st <- read_bundle_table("extraction_studies_round1.csv")
corrections <- if (file.exists("data/extraction_corrections.csv")) read.csv("data/extraction_corrections.csv", stringsAsFactors = FALSE) else NULL
applied <- apply_corrections(vr, corrections)
vr <- applied$data

cd209 <- run_term("cd209", block_haplotypes("CD209", c("G", NA)), c("rs4804800", "rs7248637"), snp_order = c("rs4804800", "rs7248637"),
                  alt_snp = "rs7248637", seed = 5918L, vr = vr, st = st, applied_ids = applied$applied, hr_fallback_min = 3)
il10 <- run_term("il10", single_snp_haplotypes("rs1800896", "C"), "rs1800896", snp_order = "rs1800896",
                 seed = 2473L, vr = vr, st = st, applied_ids = applied$applied, hr_fallback_min = 3,
                 scenarios_extra = list("with de Boer 2011 (S061, genotype source unstated, treated as recipient)" = function(re) {
                   re$genotype[re$study_id == "S061"] <- "recipient"
                   ptx3_study_effects(re, snp_order = "rs1800896")
                 }))
ifng <- run_term("ifng", single_snp_haplotypes("rs2069705", "G"), "rs2069705", snp_order = "rs2069705",
                 seed = 8630L, vr = vr, st = st, applied_ids = applied$applied, hr_fallback_min = 3)

show_term(cd209, "CD209 (risk class: G at rs4804800)")
show_term(il10, "IL10 rs1800896 (risk allele: forward-strand C)")
show_term(ifng, "IFNG rs2069705 (risk allele: forward-strand G)")
