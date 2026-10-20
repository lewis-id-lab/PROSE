source(file.path(proj_root, "R", "bs_data.R"))
source(file.path(proj_root, "R", "decision_curve.R"))
source(file.path(proj_root, "R", "boscore_data.R"))

# Synthetic stand-in for boscore.xlsx: the same column names, 15 rows with a missing BOSCORE item (excluded),
# and 8 event admissions that also appear as a non-event twin.
make_fake_boscore_raw <- function(seed = 3071, n_patients = 60, n_twins = 8) {
  set.seed(seed)
  n_adm <- sample(1:4, n_patients, replace = TRUE)
  patient <- rep(sprintf("P%03d", seq_len(n_patients)), n_adm)
  n <- length(patient)
  X <- sapply(PREDICTORS, function(p) rbinom(n, 1, pars$marginal[[p]]))
  case <- rbinom(n, 1, 0.06)
  case[sample(which(case == 0), 10)] <- 1                      # make sure there are enough events to twin
  cols <- setNames(BOSCORE_PREDICTORS, names(BOSCORE_PREDICTORS))
  raw <- data.frame(id = seq_len(n), IDRICO = sprintf("A%05d", seq_len(n)), IDIND = patient,
                    DRICEM = format(as.Date("2010-01-01") + sample(0:3000, n), "%m/%d/%Y"),
                    stringsAsFactors = FALSE)
  raw$yearDRICEM <- as.character(as.integer(format(as.Date(raw$DRICEM, "%m/%d/%Y"), "%Y")))
  raw$EORTC_PROV_PROB <- as.character(case)
  for (cn in names(cols)) raw[[cn]] <- as.character(X[, cols[[cn]]])
  incomplete <- sample(n, 15)
  for (i in incomplete) {
    cn <- sample(names(cols), 1)              # choose the column once; indexing with sample() inline evaluates it twice
    raw[[cn]][i] <- NA_character_
  }
  raw$pred <- ifelse(seq_len(n) %in% incomplete, NA_character_, "0.04")
  # twins: copy an included event row, make the copy a non-event, same admission number and patient
  src <- sample(setdiff(which(case == 1), incomplete), n_twins)
  twin <- raw[src, ]; twin$id <- n + seq_len(n_twins); twin$EORTC_PROV_PROB <- "0"
  rbind(raw, twin)
}
fake <- make_fake_boscore_raw()
bc <- parse_boscore(fake)

test_that("the BOSCORE columns map onto the right model predictors", {
  expect_equal(colnames(bc$X), PREDICTORS)
  for (cn in names(BOSCORE_PREDICTORS)) {
    expect_equal(unname(bc$X[, BOSCORE_PREDICTORS[[cn]]]), suppressWarnings(as.numeric(fake[[cn]])), info = cn)
  }
  expect_setequal(unname(BOSCORE_PREDICTORS), PREDICTORS)
  expect_equal(BOSCORE_PREDICTORS[["BOSCORE_8"]], "uncontrolled_malignancy")
  expect_equal(BOSCORE_PREDICTORS[["BOSCORE_9"]], "high_risk_chemo")
  expect_equal(BOSCORE_PREDICTORS[["BOSCORE_4"]], "prior_mould_disease")
})

test_that("the inclusion rule is complete data on the seven items, and matches where pred is present", {
  inc <- boscore_included(bc)
  expect_equal(sum(!inc), 15L)
  expect_true(all(is.na(bc$rows$pred) == !inc))
  expect_true(all(complete.cases(bc$X[inc, ])))
})

test_that("twin pairs are found and the pattern check separates real twins from anything else", {
  pairs <- boscore_twin_pairs(bc)
  expect_equal(nrow(pairs), 8L)
  expect_true(all(unlist(check_twins(pairs))))

  bad <- pairs; bad$n_cases[1] <- 2
  expect_false(check_twins(bad)$each_pair_has_one_event)
  bad <- pairs; bad$same_patient[2] <- FALSE
  expect_false(check_twins(bad)$same_patient)
  bad <- pairs; bad$predictors_identical[3] <- FALSE
  expect_false(check_twins(bad)$predictors_identical)

  # A genuine duplicate with a different profile is reported as not identical
  f2 <- fake; i <- which(duplicated(f2$IDRICO) | duplicated(f2$IDRICO, fromLast = TRUE))[1]
  f2$BOSCORE_8[i] <- ifelse(f2$BOSCORE_8[i] == "1", "0", "1")
  expect_false(check_twins(boscore_twin_pairs(parse_boscore(f2)))$predictors_identical)
})

test_that("the base cohort keeps the event row of every pair and drops only the twin", {
  cohort <- boscore_cohort(bc)
  inc_rows <- sum(boscore_included(bc))
  expect_equal(nrow(cohort$rows), inc_rows)
  expect_equal(sum(cohort$rows$twin_duplicate), 8L)
  expect_true(all(cohort$rows$case[cohort$rows$twin_duplicate] == 0))

  base <- select_cohort(cohort)
  expect_equal(nrow(base$X), inc_rows - 8L)
  expect_false(anyDuplicated(base$rows$admission) > 0)
  expect_equal(sum(base$rows$case), sum(cohort$rows$case))                   # no event is lost
  expect_equal(length(unique(base$rows$patient)), length(unique(cohort$rows$patient)))
  expect_equal(nrow(select_cohort(cohort, include_twins = TRUE)$X), inc_rows)
})

test_that("row order does not decide which row of a pair is kept", {
  rev_bc <- parse_boscore(fake[rev(seq_len(nrow(fake))), ])
  cohort <- boscore_cohort(rev_bc)
  expect_equal(sum(cohort$rows$twin_duplicate), 8L)
  expect_true(all(cohort$rows$case[cohort$rows$twin_duplicate] == 0))
  expect_equal(sum(select_cohort(cohort)$rows$case), sum(select_cohort(boscore_cohort(bc))$rows$case))
})

test_that("cohort_summary agrees with a hand calculation", {
  base <- select_cohort(boscore_cohort(bc))
  s <- cohort_summary(base, pars)
  expect_equal(s$admissions, nrow(base$X))
  expect_equal(s$events, sum(base$rows$case))
  expect_equal(s$prevalence, mean(base$rows$case))
  expect_equal(s$patients, length(unique(base$rows$patient)))
  lp <- as.numeric(pars$intercept + base$X %*% pars$coef[PREDICTORS])
  expect_equal(s$mean_predicted_risk, mean(plogis(lp)))
})

test_that("admissions per patient pools the top category", {
  base <- select_cohort(boscore_cohort(bc))
  app <- admissions_per_patient(base, top = 3L)
  expect_equal(sum(app$patients), length(unique(base$rows$patient)))
  expect_true("3+" %in% app$admissions)
  expect_false(any(c("4", "5") %in% app$admissions))
})

test_that("resampling groups draws whole patients, with every admission of a drawn patient", {
  base <- select_cohort(boscore_cohort(bc))
  groups <- split(seq_len(nrow(base$X)), base$rows$patient)
  set.seed(1801)
  rows <- resample_group_rows(groups)
  times <- table(factor(rows, levels = seq_len(nrow(base$X))))
  # every row of a patient appears equally often, so whole patients were drawn
  for (g in groups) expect_equal(length(unique(as.integer(times[g]))), 1L)
  # and the draw is a real resample: some patients are drawn more than once, some not at all
  mult <- vapply(groups, function(g) as.integer(times[g[1]]), integer(1))
  expect_gt(max(mult), 1L)
  expect_true(any(mult == 0L))
})

risk_all <- plogis(profile_logit(pars, all_profiles()))

test_that("bootstrap_dca is reproducible and returns the documented columns", {
  base <- select_cohort(boscore_cohort(bc))
  a <- bootstrap_dca(base, risk_all, c(0.05, 0.20), B = 40, seed = 77)
  b <- bootstrap_dca(base, risk_all, c(0.05, 0.20), B = 40, seed = 77)
  expect_identical(a, b)
  expect_equal(dim(a), c(40L, 6L))
  expect_equal(colnames(a), c("nb_0.05", "nb_0.2", "mcid_0.05", "mcid_0.2", "flagged_0.05", "flagged_0.2"))
  expect_true(all(is.finite(a)))
})

test_that("resampling patients gives wider intervals than resampling admissions when admissions cluster", {
  # 80 patients with six identical admissions each, profiles differing between patients
  set.seed(2954)
  P <- 80; k <- 6
  prof <- all_profiles()[sample(128, P, replace = TRUE, prob = independent_profile_probs(pars$marginal)), ]
  X <- prof[rep(seq_len(P), each = k), ]
  cohort <- list(rows = data.frame(patient = rep(sprintf("P%02d", seq_len(P)), each = k), case = 0), X = X)
  width <- function(m) diff(quantile(m[, "mcid_0.05"], c(0.05, 0.95)))
  w_pat <- width(bootstrap_dca(cohort, risk_all, 0.05, B = 400, seed = 5, by = "patient"))
  w_adm <- width(bootstrap_dca(cohort, risk_all, 0.05, B = 400, seed = 5, by = "admission"))
  expect_gt(w_pat / w_adm, 1.5)
})

test_that("reconciliation passes for a consistent bs file and fails for each kind of tampering", {
  inc <- boscore_included(bc)
  bs <- list(X = bc$X[inc, ], y = bc$rows$case[inc], id = bc$rows$id[inc])
  expect_true(all(unlist(reconcile_boscore_bs(bc, bs))))
  expect_true(all(boscore_predictor_agreement(bc, bs)$agreement == 1))

  t1 <- bs; t1$X[5, "cmv"] <- 1 - t1$X[5, "cmv"]
  expect_false(reconcile_boscore_bs(bc, t1)$predictors_identical)
  expect_lt(min(boscore_predictor_agreement(bc, t1)$agreement), 1)

  t2 <- bs; t2$y[3] <- 1 - t2$y[3]
  expect_false(reconcile_boscore_bs(bc, t2)$outcome_identical)

  t3 <- bs; t3$id[1] <- 999999
  expect_false(reconcile_boscore_bs(bc, t3)$bs_ids_all_in_boscore)

  t4 <- list(X = bs$X[-1, ], y = bs$y[-1], id = bs$id[-1])
  expect_false(reconcile_boscore_bs(bc, t4)$included_rows_are_the_bs_rows)

  f2 <- fake; excluded <- which(is.na(f2$pred))[1]; f2$pred[excluded] <- "0.05"
  expect_false(reconcile_boscore_bs(parse_boscore(f2), bs)$pred_missing_exactly_when_excluded)
})

test_that("the pseudonymised cohort carries no identifier, preserves structure and round-trips", {
  cohort <- boscore_cohort(bc)
  ps <- pseudonymise_cohort(cohort, seed = 11)
  expect_false(any(c("id", "IDRICO", "IDIND", "patient", "admission", "date") %in% names(ps)))
  everything <- unlist(lapply(ps, as.character))
  expect_false(any(grepl("^A[0-9]{5}$|^P[0-9]{3}$", everything) & everything %in% c(fake$IDRICO, fake$IDIND)))
  expect_true(all(grepl("^P[0-9]+$", ps$patient_code)) && all(grepl("^A[0-9]+$", ps$admission_code)))
  one_to_one <- function(a, b) length(unique(a)) == length(unique(b)) && length(unique(a)) == nrow(unique(data.frame(a, b)))
  expect_true(one_to_one(cohort$rows$patient, ps$patient_code))
  expect_true(one_to_one(admission_key(cohort$rows), ps$admission_code))
  expect_equal(as.vector(table(table(ps$patient_code))), as.vector(table(table(cohort$rows$patient))))
  expect_identical(pseudonymise_cohort(cohort, seed = 11), ps)

  path <- tempfile(fileext = ".csv"); utils::write.csv(ps, path, row.names = FALSE)
  back <- cohort_from_pseudonymised(path)
  expect_equal(tabulate(profile_index(back$X), 128L), tabulate(profile_index(cohort$X), 128L))
  expect_equal(sum(back$rows$twin_duplicate), sum(cohort$rows$twin_duplicate))
  expect_equal(sum(select_cohort(back)$rows$case), sum(select_cohort(cohort)$rows$case))
  expect_equal(cohort_summary(select_cohort(back), pars)$patients, cohort_summary(select_cohort(cohort), pars)$patients)
  # the loader prefers the pseudonymised file when it exists
  loaded <- load_base_cohort(restricted = path, boscore = file.path(tempdir(), "none.xlsx"))
  expect_equal(nrow(loaded$X), nrow(cohort$X))
})

test_that("read_boscore fails clearly when the file or a needed column is missing", {
  expect_error(read_boscore(file.path(tempdir(), "no-such-file.xlsx")), "gitignored")
})

boscore_path <- file.path(proj_root, "boscore.xlsx")
bs_path <- file.path(proj_root, "bs.xlsx")

test_that("the real boscore.xlsx reconciles with bs.xlsx exactly", {
  skip_if_not(file.exists(boscore_path) && file.exists(bs_path), "boscore.xlsx and bs.xlsx are gitignored patient-level data and are not present")
  real <- read_boscore(boscore_path); bs <- read_bs_data(bs_path)
  expect_false(any(c("Connome", "Dnas") %in% names(real$rows)))                 # names and birth dates are never loaded
  expect_equal(nrow(real$rows), BOSCORE_EXPECTED$rows)
  expect_equal(sum(boscore_included(real)), BOSCORE_EXPECTED$included)
  expect_true(all(unlist(reconcile_boscore_bs(real, bs))))
  expect_true(all(boscore_predictor_agreement(real, bs)$agreement == 1))
})

test_that("the real file has 49 twin pairs, and removing them gives 4,078 admissions with all 133 events", {
  skip_if_not(file.exists(boscore_path), "boscore.xlsx is gitignored patient-level data and is not present")
  real <- read_boscore(boscore_path)
  pairs <- boscore_twin_pairs(real)
  expect_equal(nrow(pairs), BOSCORE_EXPECTED$twin_pairs)
  expect_true(all(unlist(check_twins(pairs))))
  base <- select_cohort(boscore_cohort(real))
  s <- cohort_summary(base, pars)
  expect_equal(s$admissions, BOSCORE_EXPECTED$unique_admissions)
  expect_equal(s$events, BOSCORE_EXPECTED$events)
  expect_equal(s$patients, BOSCORE_EXPECTED$patients)
  expect_equal(s$prevalence, 133 / 4078)                                          # the paper's stated 3.3%
  expect_equal(round(100 * s$prevalence, 1), 3.3)
})

test_that("the pseudonymised base cohort matches the identified file", {
  restricted <- file.path(proj_root, "data", "restricted", "base_cohort_pseudonymised.csv")
  skip_if_not(file.exists(boscore_path) && file.exists(restricted), "the identified file or the pseudonymised cohort is not present")
  from_xlsx <- boscore_cohort(read_boscore(boscore_path)); from_csv <- cohort_from_pseudonymised(restricted)
  expect_equal(tabulate(profile_index(select_cohort(from_xlsx)$X), 128L), tabulate(profile_index(select_cohort(from_csv)$X), 128L))
  expect_equal(nrow(from_xlsx$X), nrow(from_csv$X))
  expect_equal(as.vector(table(table(select_cohort(from_xlsx)$rows$patient))), as.vector(table(table(select_cohort(from_csv)$rows$patient))))
})

test_that("the allo-HSCT flag is the composite item with a transplant date present, and NA without the columns", {
  raw <- data.frame(recode2SCORE9 = c("2", "2", "0", "2", NA, "0"), DTMO = c("44500", NA, "44501", "", "44502", NA), stringsAsFactors = FALSE)
  expect_equal(allo_hsct_flag(raw), c(TRUE, FALSE, FALSE, FALSE, NA, FALSE))
  expect_true(all(is.na(allo_hsct_flag(data.frame(x = 1:3)))))
  expect_true(all(is.na(bc$rows$allo_hsct)))                    # the fixture has neither column
  expect_true(all(c("DTMO", "recode2SCORE9") %in% BOSCORE_KEEP))
})

test_that("the flag travels through the cohort and the pseudonymised file, without a date, and an older file without it still loads", {
  f <- fake
  f$recode2SCORE9 <- ifelse(seq_len(nrow(f)) %% 5 == 0, "2", "0")
  f$DTMO <- ifelse(seq_len(nrow(f)) %% 10 == 0, NA_character_, "44500")
  b <- parse_boscore(f)
  expect_equal(sum(b$rows$allo_hsct), sum(f$recode2SCORE9 == "2" & !is.na(f$DTMO)))
  cohort <- boscore_cohort(b)
  ps <- pseudonymise_cohort(cohort, seed = 11)
  expect_true("allo_hsct" %in% names(ps)); expect_false(any(c("DTMO", "recode2SCORE9") %in% names(ps)))
  # adding the flag does not change the codes
  expect_identical(ps[setdiff(names(ps), "allo_hsct")], pseudonymise_cohort(boscore_cohort(bc), seed = 11)[setdiff(names(ps), "allo_hsct")])
  path <- tempfile(fileext = ".csv"); utils::write.csv(ps, path, row.names = FALSE)
  back <- cohort_from_pseudonymised(path)
  expect_equal(back$rows$allo_hsct, cohort$rows$allo_hsct)
  expect_equal(allo_summary(select_cohort(back)), allo_summary(select_cohort(cohort)))
  s <- allo_summary(select_cohort(cohort))
  expect_equal(s$allo_admissions, sum(select_cohort(cohort)$rows$allo_hsct)); expect_lte(s$allo_events, s$allo_admissions)
  ps$allo_hsct <- NULL; old <- tempfile(fileext = ".csv"); utils::write.csv(ps, old, row.names = FALSE)
  expect_true(all(is.na(cohort_from_pseudonymised(old)$rows$allo_hsct)))
  expect_true(is.na(allo_summary(cohort_from_pseudonymised(old))$allo_admissions))
})
