# Verify PTX3 extraction rows against their source papers, starting with Herrero-Sanchez 2018 (S035) and Brunel 2018 (S029),
# and record the corrections. Run from the project root, before scripts/12 and 14.
#
# Sources read on 18 October 2026: Herrero-Sanchez 2018 from PMC (PMC6246897), Brunel 2018 from the Haematologica
# site (doi 10.3324/haematol.2018.195453). These are read-throughs by the assistant, not the human round-2 check, which is
# still pending. The counts in data/source_checks_S035.csv were transcribed from Tables 3, 4 and 5 of the S035 paper.
#
# Writes   data/extraction_corrections.csv        one row per corrected field, with the evidence
#          output/s035_direction_check.csv        the direction of each printed odds ratio, from the genotype counts
#          output/s029_genotype_frequency_check.csv  Brunel's homozygote groups against Hardy-Weinberg expectations
#          output/ptx3_verification_log.csv       what was checked, what it showed, and what follows

source("R/evidence_mapping.R")
source("R/source_verification.R")

vr <- read_bundle_table("extraction_variants_round1.csv")

# 1. Herrero-Sanchez 2018: which group does each printed odds ratio belong to?
chk <- read.csv("data/source_checks_S035.csv", stringsAsFactors = FALSE)
dir <- do.call(rbind, lapply(seq_len(nrow(chk)), function(i) with(chk[i, ], printed_or_direction(reference_events, reference_nonevents, other_events, other_nonevents, printed_or))))
s035_check <- cbind(chk[, c("table", "snp", "model", "reference_label", "printed_or", "printed_ci_lo", "printed_ci_hi")], dir)
write.csv(s035_check, "output/s035_direction_check.csv", row.names = FALSE)

# 2. The corrections: in every S035 dominant and recessive PTX3 row the "(ref)" label is on the wrong group
rows <- vr[vr$study_id == "S035" & vr$rsid %in% PTX3_LOCI & vr$genetic_model %in% c("dominant", "recessive") &
             grepl("\\(ref\\)", vr$comparison_as_reported), ]
corrections <- data.frame(
  variant_row_id = rows$variant_row_id, study_id = rows$study_id, field = "comparison_as_reported",
  original_value = rows$comparison_as_reported, corrected_value = clean_reversed_comparison(rows$comparison_as_reported),
  evidence = paste0("The printed odds ratio belongs to the reference-labelled group against the other group, not the reverse. ",
                    "Checked from the genotype counts in Tables 3, 4 and 5 for all three SNPs (output/s035_direction_check.csv); ",
                    "the paper's text agrees (the AA genotype of rs2305619 and rs1840680 and the CC genotype of rs3816527 are called the risk genotypes). ",
                    "Corrected text keeps only the genotype contrast, with the first group as the numerator."),
  source_location = rows$page_or_table,
  checked_on = "2026-10-18", checked_by = "assistant read of PMC6246897 (not the human round-2 check)",
  stringsAsFactors = FALSE
)
write.csv(corrections, "data/extraction_corrections.csv", row.names = FALSE)

# 3. Brunel 2018: do the extracted homozygote groups look like the minor-allele homozygotes the paper describes?
# The paper reports minor allele frequencies of 0.48 (rs2305619) and 0.44 (rs3816527) in the whole cohort.
b <- vr[vr$study_id == "S029" & vr$rsid %in% c("rs2305619", "rs3816527") & vr$usable_tier == "counts_2x2" & !is.na(vr$usable_tier), ]
num <- function(x) suppressWarnings(as.numeric(x))
b$group_total <- num(b$n_cases_hom_effect) + num(b$n_ctrl_hom_effect)
b$n_total <- b$group_total + num(b$n_cases_het) + num(b$n_ctrl_het)
maf <- c(rs2305619 = 0.48, rs3816527 = 0.44)
hw <- do.call(rbind, lapply(seq_len(nrow(b)), function(i) cbind(data.frame(rsid = b$rsid[i], analysis = substr(b$subgroup[i], 1, 40), group_total = b$group_total[i], n_total = b$n_total[i]),
                                                              hwe_group_check(b$group_total[i], b$n_total[i], maf[[b$rsid[i]]]))))
write.csv(hw, "output/s029_genotype_frequency_check.csv", row.names = FALSE)

# 4. The log
log <- data.frame(
  study = c(rep("S035", 7), rep("S029", 6)),
  item = c("extracted odds ratios and intervals, Table 3, all three SNPs and models",
           "reference-group label on the odds ratios, Tables 3 to 8",
           "endpoint of the extracted primary rows (Table 3)",
           "invasive aspergillosis endpoint (Table 5)",
           "genotype source and covariates",
           "rs1840680 recessive interval in Table 3",
           "CLEC7A results (supplementary tables)",
           "design, size and outcome",
           "minor allele frequencies and linkage",
           "per-SNP genotype counts",
           "combined exposure and hazard ratio in the non-neutropenic stratum",
           "extracted rs3816527 homozygote group against Hardy-Weinberg expectations",
           "genotype source"),
  source_says = c("Table 3 prints the same ORs and CIs as the extraction (rs2305619 dominant 1.83, 0.89 to 3.76; recessive 1.55; rs3816527 dominant 1.74, recessive 2.02; rs1840680 dominant 1.94, recessive 1.76)",
                  "The reference group carries 1.00 but the printed OR is for that group against the others. The group counts show it for every SNP and model checked, and the text calls AA (rs2305619, rs1840680) and CC (rs3816527) the risk genotypes",
                  "Detection of Candida and/or Aspergillus, which includes colonisation (Table 3); 198 patients",
                  "Invasive aspergillosis in the whole cohort: 12 patients with EORTC/MSG invasive disease against 186 without; adjusted for sex, age and malignancy",
                  "Samples taken at diagnosis, before HSCT in transplant patients; adjusted for sex, age, haematological malignancy",
                  "The paper prints 1.76 (0.92 to 1.76): the upper limit equals the estimate, a printing error in the source",
                  "Not read: they are in the supplementary file, which was not retrieved",
                  "185 patients and 212 chemotherapy courses (2007 to 2017); 26 probable or proven invasive mould infections (16 probable, 10 proven), all moulds; 6-month cumulative incidence 14%",
                  "Minor allele frequencies 0.44 (rs3816527) and 0.48 (rs2305619); linkage R = 0.78; both in Hardy-Weinberg equilibrium. These fit A and C as the minor alleles",
                  "The main text prints none: per-SNP results are in Online Supplementary Figure S2, an image",
                  "Exposure is two copies of the minor allele of rs2305619 and/or rs3816527; adjusted HR 5.06 (1.68 to 15.2) without pre-existing neutropenia",
                  "With a minor allele frequency of 0.44, CC homozygotes should be about 19% of courses; the extracted group is about 30%, which fits AA homozygotes (about 31%) and not CC (z about 3.8). For rs2305619 the minor and major homozygote shares (23% and 27%) are too close to tell apart",
                  "Not stated in the text. Samples from allo-HSCT relapse patients were used only if taken after engraftment, when blood DNA may be the donor's"),
  verdict = c("confirmed", "corrected", "consistent with D3b", "available as a sensitivity analysis", "confirmed", "source error, kept as printed",
              "not verified", "confirmed", "supports the inferred alleles, weakly for rs2305619", "not verifiable from the text",
              "confirmed, but it cannot be mapped to one SNP", "unresolved: the group may not be the minor homozygotes", "unresolved"),
  consequence = c("none", "the 30 dominant and recessive rows are corrected; the S035 effect is positive (A allele raises risk), in line with the other studies, not the opposite",
                  "S035 stays out of the primary analysis", "a separate sensitivity analysis uses the Table 5 rows",
                  "S035 is a recipient-genotype study", "that row is excluded by the interval check",
                  "S035 CLEC7A rows stay unverified", "none", "Brunel's inferred effect allele (A at rs2305619) is plausible but not proven",
                  "Brunel's per-SNP counts cannot be checked against the text", "none for the per-SNP analysis",
                  "Brunel's rs3816527 counts may be the wrong homozygote group; a sensitivity analysis drops Brunel",
                  "Brunel is treated as a recipient study, flagged"),
  stringsAsFactors = FALSE
)
log$read_how <- ifelse(log$study == "S035", "full article text from PMC, tables read directly",
                       "summary of the article page returned by the fetch tool; figures and quotes not verbatim, to be re-read")
write.csv(log, "output/ptx3_verification_log.csv", row.names = FALSE)

message("S035 rows corrected: ", nrow(corrections))
message("\nDirection of the printed odds ratios (decisive = the crude OR is far enough from 1 to tell the readings apart):")
print(table(reading = s035_check$reading, decisive = s035_check$decisive))
message("\nS035 detail:")
print(cbind(s035_check[, c("table", "snp", "model", "printed_or")], round(s035_check[, c("crude_ref_vs_other", "distance_if_reference_numerator", "distance_if_other_numerator")], 2), decisive = s035_check$decisive), row.names = FALSE)
message("\nBrunel homozygote groups against Hardy-Weinberg:")
print(cbind(hw[, c("rsid", "analysis")], round(hw[, c("observed", "expected_minor_homozygotes", "expected_major_homozygotes", "z_minor_homozygotes", "z_major_homozygotes")], 2)), row.names = FALSE)
