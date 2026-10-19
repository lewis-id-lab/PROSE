source(file.path(proj_root, "R", "evidence_mapping.R"))
source(file.path(proj_root, "R", "ptx3_evidence.R"))
source(file.path(proj_root, "R", "source_verification.R"))

test_that("crude_or follows the formula and adds 0.5 only when a cell is zero", {
  expect_equal(crude_or(23, 18, 66, 91), (23 / 18) / (66 / 91))
  expect_equal(crude_or(0, 10, 5, 20), (0.5 / 10.5) / (5.5 / 20.5))
})

test_that("a printed odds ratio is attributed to the reference-labelled group when the counts say so", {
  # Herrero-Sanchez 2018 Table 3, rs2305619 dominant: AA labelled reference, printed OR 1.83
  d <- printed_or_direction(23, 18, 66, 91, 1.83)
  expect_equal(d$reading, "reference-labelled group is the numerator")
  expect_true(d$decisive)
  expect_lt(d$distance_if_reference_numerator, 0.1); expect_gt(d$distance_if_other_numerator, 1)
  # the opposite case: a printed value close to the reciprocal
  r <- printed_or_direction(23, 18, 66, 91, 0.57)
  expect_equal(r$reading, "other group is the numerator")
  # near 1 the two readings cannot be separated, and the check says so
  expect_false(printed_or_direction(10, 90, 11, 89, 1.0)$decisive)
})

test_that("Hardy-Weinberg group check gives the expected shares and z-scores", {
  h <- hwe_group_check(60, 200, 0.44)
  expect_equal(h$observed, 0.30)
  expect_equal(h$expected_minor_homozygotes, 0.44^2); expect_equal(h$expected_major_homozygotes, 0.56^2)
  expect_equal(h$z_minor_homozygotes, (0.30 - 0.1936) / sqrt(0.1936 * 0.8064 / 200))
  expect_gt(h$z_minor_homozygotes, 3); expect_lt(abs(h$z_major_homozygotes), 1)
})

test_that("corrections are applied, each one checks what it replaces, and failures are loud", {
  vr <- data.frame(variant_row_id = c("V1", "V2"), comparison_as_reported = c("AA (ref) vs AG+GG", "CC (ref) vs CA+AA"), stringsAsFactors = FALSE)
  cr <- data.frame(variant_row_id = "V1", field = "comparison_as_reported", original_value = "AA (ref) vs AG+GG", corrected_value = "AA vs AG+GG", stringsAsFactors = FALSE)
  out <- apply_corrections(vr, cr)
  expect_equal(out$data$comparison_as_reported, c("AA vs AG+GG", "CC (ref) vs CA+AA")); expect_equal(out$applied, "V1")
  expect_equal(vr$comparison_as_reported[1], "AA (ref) vs AG+GG")                          # the input is not modified
  stale <- transform(cr, original_value = "something else")
  expect_error(apply_corrections(vr, stale), "expected 'something else'")
  expect_error(apply_corrections(vr, transform(cr, variant_row_id = "V9")), "not found")
  expect_error(apply_corrections(vr, transform(cr, field = "nope")), "no field")
  none <- apply_corrections(vr, NULL)
  expect_identical(none$data, vr); expect_length(none$applied, 0)
})

test_that("the reversed comparison text is cleaned to a plain genotype contrast the parser reads", {
  expect_equal(clean_reversed_comparison("AA (ref) vs AG+GG combined, IA vs no IA"), "AA vs AG+GG combined")
  expect_equal(clean_reversed_comparison("AA+AG combined (ref) vs GG"), "AA+AG combined vs GG")
  expect_equal(clean_reversed_comparison("CC (ref) vs CA+AA combined, invasive candidiasis vs no candidiasis"), "CC vs CA+AA combined")
  ct <- parse_contrast(clean_reversed_comparison("AA (ref) vs AG+GG combined, IA vs no IA"), "dominant", "A", "G")
  expect_equal(ct$numerator, "AA"); expect_setequal(ct$reference, c("AG", "GG")); expect_equal(ct$method, "text")
  # before the correction the same text could not be read at all, because it contains a second "vs"
  expect_equal(parse_contrast("AA (ref) vs AG+GG combined, IA vs no IA", "dominant", "A", "G")$flag, "more than two groups")
})

bundle <- file.path(proj_root, "IMD_genetics_analysis_bundle (1).zip")
hap_path <- file.path(proj_root, "data", "haplotype_frequencies_phased.csv")
checks_path <- file.path(proj_root, "data", "source_checks_S035.csv")
corr_path <- file.path(proj_root, "data", "extraction_corrections.csv")

test_that("the transcribed S035 counts add up and every printed odds ratio belongs to the reference-labelled group", {
  skip_if_not(file.exists(checks_path), "the source checks are not present")
  ck <- read.csv(checks_path, stringsAsFactors = FALSE)
  expect_equal(nrow(ck), 18L)
  # the groups of a dominant and a recessive row of the same SNP and table overlap on one genotype, so totals must match
  tot <- function(d) d$reference_events + d$reference_nonevents + d$other_events + d$other_nonevents
  by <- split(ck, paste(ck$table, ck$snp))
  for (nm in names(by)) expect_equal(length(unique(tot(by[[nm]]))), 1L, info = nm)       # the same patients in both models
  expect_equal(unique(tot(ck[ck$table == "Table 3", ])), 198); expect_equal(unique(tot(ck[ck$table == "Table 5", ])), 198)
  dir <- do.call(rbind, lapply(seq_len(nrow(ck)), function(i) with(ck[i, ], printed_or_direction(reference_events, reference_nonevents, other_events, other_nonevents, printed_or))))
  expect_true(all(dir$reading == "reference-labelled group is the numerator"))
  expect_gte(sum(dir$decisive), 16L)
  # printed odds ratios sit inside their own intervals, except the one source typo
  inside <- ck$printed_or > ck$printed_ci_lo & ck$printed_or < ck$printed_ci_hi
  expect_equal(ck$snp[!inside], "rs1840680"); expect_equal(ck$model[!inside], "recessive"); expect_equal(ck$table[!inside], "Table 3")
})

test_that("every recorded correction matches the real extraction and carries its evidence", {
  skip_if_not(file.exists(bundle) && file.exists(corr_path), "the evidence bundle or the corrections file is not present")
  vr <- read_bundle_table("extraction_variants_round1.csv", bundle)
  cr <- read.csv(corr_path, stringsAsFactors = FALSE)
  expect_equal(nrow(cr), 43L); expect_equal(sum(cr$study_id == "S035"), 30L); expect_equal(sum(cr$study_id == "S073"), 12L); expect_equal(sum(cr$study_id == "S071"), 1L)
  expect_setequal(cr$study_id, c("S035", "S073", "S071"))
  expect_equal(cr$corrected_value[cr$variant_row_id == "V0261"], "G/G vs C/C+C/G")
  expect_true(all(nzchar(cr$evidence)) && all(nzchar(cr$checked_by)) && all(nzchar(cr$checked_on)))
  expect_false(anyDuplicated(paste(cr$variant_row_id, cr$field)) > 0)
  out <- apply_corrections(vr, cr)                                                       # stops if any original value has changed
  expect_setequal(out$applied, cr$variant_row_id)
  expect_false(any(grepl("\\(ref\\)", out$data$comparison_as_reported[out$data$variant_row_id %in% cr$variant_row_id])))
  # genotypic three-level rows are not corrected: they cannot be mapped to one contrast
  expect_false(any(vr$genetic_model[vr$variant_row_id %in% cr$variant_row_id] == "genotypic_3level"))
  # Fischer 2016 (S073): the corrected rows put the compared group first, so the risk genotype is the numerator
  s73 <- out$data[out$data$variant_row_id %in% cr$variant_row_id[cr$study_id == "S073"], ]
  expect_true(all(s73$comparison_as_reported %in% c("A/C vs A/A", "G/G vs C/C + C/G", "G/G + C/G vs C/C")))
  expect_equal(sort(unique(s73$rsid)), c("rs16910526", "rs7309123"))
  # Sainz 2012 (S071): the recessive rs7309123 row now reads G/G against the rest, like its duplicate V0322
  v <- out$data[out$data$variant_row_id %in% c("V0261", "V0322"), ]
  expect_equal(unique(v$comparison_as_reported), "G/G vs C/C+C/G")
})

test_that("after the S073 correction Fischer 2016 reads as risk for rs7309123 G and as no excess for Y238X carriers", {
  skip_if_not(file.exists(bundle) && file.exists(corr_path) && file.exists(hap_path), "the evidence bundle, corrections or haplotype table is not present")
  vr <- apply_corrections(read_bundle_table("extraction_variants_round1.csv", bundle), read.csv(corr_path, stringsAsFactors = FALSE))$data
  st <- read_bundle_table("extraction_studies_round1.csv", bundle)
  hc <- block_haplotypes("CLEC7A", c("C", NA), hap_path)
  rc <- ptx3_row_effects(vr, st, hc, loci = c("rs16910526", "rs7309123"), counts_model = "carrier")
  y238x <- ptx3_study_effects(rc, snp_order = "rs16910526")
  expect_lt(y238x$beta[y238x$study_id == "S073"], 0)                                     # carriers OR 0.7
  expect_equal(exp(y238x$beta[y238x$study_id == "S073"]), 0.7, tolerance = 0.01)
  g <- rc[rc$study_id == "S073" & rc$rsid == "rs7309123" & !is.na(rc$beta), ]
  expect_true(all(g$y > 0))                                                              # G/G, or G carriers, against the rest: OR above 1
})
