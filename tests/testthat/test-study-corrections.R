source(file.path(proj_root, "R", "source_verification.R"))

st <- data.frame(study_id = c("S1", "S2", "S3"), n_total = c(100, NA, 50), n_cases = c(NA, NA, 10), note = c("a", "b", "c"), stringsAsFactors = FALSE)
cr <- function(...) {
  d <- data.frame(study_id = "S1", field = "n_total", original_value = "100.0", corrected_value = "90", status = "corrected", stringsAsFactors = FALSE)
  for (n in names(list(...))) d[[n]] <- list(...)[[n]]
  d
}

test_that("a correction replaces the value, comparing numbers as numbers", {
  r <- apply_study_corrections(st, cr())
  expect_equal(as.numeric(r$data$n_total[1]), 90)
  expect_equal(nrow(r$applied), 1)
})

test_that("a blank expected value must match a blank cell", {
  ok <- apply_study_corrections(st, cr(study_id = "S2", field = "n_cases", original_value = NA_character_, corrected_value = "12"))
  expect_equal(as.numeric(ok$data$n_cases[2]), 12)
  expect_error(apply_study_corrections(st, cr(study_id = "S3", field = "n_cases", original_value = NA_character_, corrected_value = "12")), "expected")
})

test_that("an unexpected value, unknown study or unknown field stops the run", {
  expect_error(apply_study_corrections(st, cr(original_value = "101")), "expected")
  expect_error(apply_study_corrections(st, cr(study_id = "S9")), "not found")
  expect_error(apply_study_corrections(st, cr(field = "nope")), "no field")
})

test_that("needs_decision rows are reported and not applied", {
  r <- apply_study_corrections(st, cr(status = "needs_decision"))
  expect_equal(r$data$n_total[1], 100)
  expect_equal(nrow(r$pending), 1)
})

test_that("the committed study corrections apply cleanly to the bundle's study table", {
  skip_if_not(file.exists(file.path(proj_root, "IMD_genetics_analysis_bundle (1).zip")))
  source(file.path(proj_root, "R", "evidence_mapping.R"))
  old <- setwd(proj_root); on.exit(setwd(old))
  r <- apply_study_corrections(read_bundle_table("extraction_studies_round1.csv"), read.csv("data/study_corrections.csv", stringsAsFactors = FALSE))
  expect_gt(nrow(r$applied), 0)
  expect_equal(as.numeric(r$data$n_total[r$data$study_id == "S074"]), 172)
})

test_that("rows marked dup_of_<study> are dropped and other duplicate groups are kept", {
  vr <- data.frame(variant_row_id = paste0("V", 1:5), study_id = c("S002", "S002", "S079", "S073", "S073"),
                   duplicate_group = c("dup_of_S079", "dup_of_S079", NA, "D001", ""), stringsAsFactors = FALSE)
  out <- drop_duplicate_reports(vr)
  expect_equal(out$variant_row_id, c("V3", "V4", "V5"))
  expect_equal(nrow(drop_duplicate_reports(vr[3:5, ])), 3L)
  expect_error(drop_duplicate_reports(vr[c("variant_row_id", "study_id")]), "duplicate_group")
})

test_that("on the real extraction the 17 S002 rows are dropped after the corrections and nothing else is", {
  bundle <- file.path(proj_root, "IMD_genetics_analysis_bundle (1).zip")
  skip_if_not(file.exists(bundle))
  source(file.path(proj_root, "R", "evidence_mapping.R"))
  vr <- read_bundle_table("extraction_variants_round1.csv", bundle)
  after <- drop_duplicate_reports(apply_corrections(vr, read.csv(file.path(proj_root, "data", "extraction_corrections.csv"), stringsAsFactors = FALSE))$data)
  expect_equal(nrow(vr) - nrow(after), 17L)
  expect_false("S002" %in% after$study_id)
  expect_true(all(vr$variant_row_id[vr$study_id == "S079"] %in% after$variant_row_id))
})

test_that("add_extraction_rows appends new rows, keeps numeric columns numeric and refuses clashes", {
  vr <- data.frame(variant_row_id = c("V1", "V2"), n = c(1, 2), note = c("a", "b"), stringsAsFactors = FALSE)
  f <- tempfile(fileext = ".csv")
  write.csv(data.frame(variant_row_id = "L1", n = 5, note = "c"), f, row.names = FALSE)
  out <- add_extraction_rows(vr, f)
  expect_equal(out$variant_row_id, c("V1", "V2", "L1")); expect_true(is.numeric(out$n)); expect_equal(out$n[3], 5)
  write.csv(data.frame(variant_row_id = "V2", n = 5, note = "c"), f, row.names = FALSE)
  expect_error(add_extraction_rows(vr, f), "new and unique")
  write.csv(data.frame(variant_row_id = "L1", n = 5, other = "c"), f, row.names = FALSE)
  expect_error(add_extraction_rows(vr, f), "columns")
})

test_that("the Leiden letter's rows are entered correctly: counts total as printed and reproduce the printed crude odds ratios", {
  bundle <- file.path(proj_root, "IMD_genetics_analysis_bundle (1).zip")
  skip_if_not(file.exists(bundle))
  source(file.path(proj_root, "R", "evidence_mapping.R"))
  vr <- read_bundle_table("extraction_variants_round1.csv", bundle)
  out <- add_extraction_rows(vr, file.path(proj_root, "data", "extraction_added_rows.csv"))
  expect_equal(nrow(out) - nrow(vr), 3L)
  l <- out[out$variant_row_id %in% c("L0001", "L0002"), ]
  cases <- as.matrix(l[c("n_cases_hom_effect", "n_cases_het", "n_cases_hom_other")]); ctrls <- as.matrix(l[c("n_ctrl_hom_effect", "n_ctrl_het", "n_ctrl_hom_other")])
  expect_equal(unname(rowSums(cases)), c(42, 41)); expect_equal(unname(rowSums(ctrls)), c(54, 59))
  # the letter prints crude odds ratios of 0.88 (+281A>G GG) and 0.90 (+734A>C AA)
  crude <- vapply(1:2, function(i) crude_or(cases[i, 1], sum(cases[i, 2:3]), ctrls[i, 1], sum(ctrls[i, 2:3])), numeric(1))
  expect_equal(round(crude, 2), c(0.88, 0.90))
  expect_true(all(out$study_id[out$variant_row_id %in% c("L0001", "L0002", "L0003")] == "S002"))
})
