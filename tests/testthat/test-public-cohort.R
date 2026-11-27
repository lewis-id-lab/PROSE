source(file.path(proj_root, "R", "public_cohort.R"))
source(file.path(proj_root, "R", "boscore_data.R"))

fake_X <- function(n = 600, seed = 4410) {
  set.seed(seed)
  X <- sapply(PREDICTORS, function(p) rbinom(n, 1, if (p == "prior_mould_disease") 0.04 else 0.3))
  matrix(as.numeric(X), ncol = length(PREDICTORS), dimnames = list(NULL, PREDICTORS))
}

test_that("no published cell is below the minimum and suppressed cells carry no count", {
  X <- fake_X()
  tab <- public_profile_table(X)
  n <- tabulate(profile_index(X), 128L)
  expect_equal(nrow(tab), 128L)
  expect_true(all(is.na(tab$admissions[tab$suppressed])))
  shown <- tab$admissions[!tab$suppressed]
  expect_false(any(shown > 0 & shown < PUBLIC_MIN_CELL))
  expect_equal(tab$suppressed, n > 0 & n < PUBLIC_MIN_CELL)
  expect_equal(tab$admissions[!tab$suppressed], n[!tab$suppressed])
})

test_that("the table has no outcome, patient or date column", {
  tab <- public_profile_table(fake_X())
  expect_identical(names(tab), c("profile", PREDICTORS, "admissions", "suppressed"))
})

test_that("check_public_table refuses a cell below the minimum and a lone suppressed cell", {
  X <- fake_X()
  tab <- public_profile_table(X); meta <- public_profile_meta(X, tab)
  bad <- tab; i <- which(bad$suppressed)[1]; bad$suppressed[i] <- FALSE; bad$admissions[i] <- 2L
  expect_error(check_public_table(bad, meta), "below the minimum")
  one <- tab; s <- which(one$suppressed)
  one$suppressed[s[-1]] <- FALSE; one$admissions[s[-1]] <- 0L
  expect_error(check_public_table(one, meta), "exactly one suppressed cell")
})

test_that("the stand-in cohort keeps the total and every published cell exactly", {
  X <- fake_X()
  dir <- tempfile(); dir.create(dir)
  p <- file.path(dir, "t.csv"); m <- file.path(dir, "m.csv")
  out <- write_public_profile_table(X, p, m)
  co <- cohort_from_public(p, m)
  expect_true(is_public_cohort(co))
  expect_equal(nrow(co$X), nrow(X))
  n_true <- tabulate(profile_index(X), 128L); n_pub <- tabulate(profile_index(co$X), 128L)
  keep <- !out$table$suppressed
  expect_equal(n_pub[keep], n_true[keep])
  expect_true(all(n_pub[out$table$suppressed] >= 1 & n_pub[out$table$suppressed] < PUBLIC_MIN_CELL))
  expect_true(all(is.na(co$rows$case)))
})

test_that("pseudo-patients follow the published admissions-per-patient distribution exactly", {
  app <- data.frame(admissions = c("1", "2", "3", "8+"), patients = c(50L, 20L, 10L, 3L))
  sizes <- pseudo_patient_sizes(app, 50 + 40 + 30 + 27)
  expect_equal(sum(sizes), 147); expect_equal(length(sizes), 83L)
  expect_equal(sort(tail(sort(sizes), 3)), c(9, 9, 9))
  pp <- pseudo_patients(147, app)
  expect_equal(sort(as.integer(table(pp))), sort(sizes))
  expect_error(pseudo_patient_sizes(app, 100), "does not add up")
})

test_that("the committed admissions-per-patient table matches the base cohort", {
  app <- read.csv(file.path(proj_root, "output", "admissions_per_patient.csv"), stringsAsFactors = FALSE)
  sizes <- pseudo_patient_sizes(app, 4078)
  expect_equal(length(sizes), 1945L)
  expect_true(all(tail(sort(sizes), 14) >= 8))
})

test_that("load_base_cohort falls back to the public table and select_cohort refuses twins on it", {
  X <- fake_X()
  dir <- tempfile(); dir.create(dir)
  p <- file.path(dir, "t.csv"); m <- file.path(dir, "m.csv")
  write_public_profile_table(X, p, m)
  co <- suppressMessages(load_base_cohort(restricted = file.path(dir, "none.csv"), boscore = file.path(dir, "none.xlsx"), public = p, public_meta = m))
  expect_true(is_public_cohort(co))
  expect_true(is_public_cohort(select_cohort(co, include_twins = FALSE)))
  expect_error(select_cohort(co, include_twins = TRUE), "twins")
  expect_error(require_patient_data(co, "Refitting"), "patient-level")
})

test_that("the committed public table, when present, follows the rule", {
  p <- file.path(proj_root, "data", "public", "profile_admissions.csv")
  skip_if_not(file.exists(p), "public profile table not yet exported")
  pub <- read_public_profile_table(p, file.path(proj_root, "data", "public", "profile_admissions_meta.csv"))
  expect_equal(pub$meta$value[pub$meta$key == "total_admissions"], 4078)
})
