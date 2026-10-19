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

source("R/term_mapping.R")

clec <- run_term("clec7a", block_haplotypes("CLEC7A", c("C", NA)), c("rs16910526", "rs7309123"), snp_order = "rs16910526",
                 seed = 7241L, vr = vr, st = st, applied_ids = applied$applied)
tlr4 <- run_term("tlr4", block_haplotypes("TLR4", c("G", "T")), c("rs4986790", "rs4986791"), snp_order = c("rs4986790", "rs4986791"),
                 alt_snp = "rs4986791", seed = 3086L, vr = vr, st = st, applied_ids = applied$applied,
                 scenarios_extra = list("with de Boer 2011 (S061, genotype source unstated, treated as recipient)" = function(re) {
                   re$genotype[re$study_id == "S061"] <- "recipient"
                   ptx3_study_effects(re, snp_order = c("rs4986790", "rs4986791"))
                 }))

show_term(clec, "CLEC7A (risk class: Y238X, haplotype C-G)")
show_term(tlr4, "TLR4 (risk class: G-T haplotype)")

r73 <- clec$re[clec$re$rsid == "rs7309123", ]
message("\nCLEC7A rs7309123: ", nrow(r73), " usable rows from ", length(unique(r73$study_id)), " studies are not used (strand unresolved, no control allele frequency).")
