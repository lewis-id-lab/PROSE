# The base dataset: boscore.xlsx (the parent of bs.xlsx), with the twin admissions corrected.
#
# boscore.xlsx has 4,694 rows (admissions 2007 to 2016). The 4,127 rows with complete data on the seven
# BOSCORE_* items are the analysis set of the paper (bs.xlsx). Among them, 49 admission numbers appear on two
# rows that share a patient and a predictor profile, with exactly one event row and one non-event twin. The
# base cohort keeps each admission once (the event row), giving 4,078 admissions. See phase3_decisions.md,
# addendum 3.
#
# boscore.xlsx is gitignored patient-level data with direct identifiers. read_boscore() reads only the columns
# it needs, so names and birth dates are never loaded. Everything returned is a plain list; the aggregate
# helpers below return no identifiers.
#
# Requires R/baseline.R (PREDICTORS, profile_index, ...), R/bs_data.R (auc_rank) and R/decision_curve.R (dca_exact).

# BOSCORE_* column -> model predictor (exact match to the bs.xlsx variables on all 4,127 rows)
BOSCORE_PREDICTORS <- c(BOSCORE_8 = "uncontrolled_malignancy", BOSCORE_9 = "high_risk_chemo",
                        BOSCORE_6 = "high_dose_steroids", BOSCORE_11 = "severe_lymphopenia", BOSCORE_15 = "cmv",
                        BOSCORE_10 = "prolonged_neutropenia", BOSCORE_4 = "prior_mould_disease")
BOSCORE_KEEP <- c("id", "IDRICO", "IDIND", "DRICEM", "yearDRICEM", "EORTC_PROV_PROB", names(BOSCORE_PREDICTORS), "pred")

# What the review of 17 to 18 October 2026 found; the reconciliation script checks the file against these.
BOSCORE_EXPECTED <- list(rows = 4694L, included = 4127L, twin_pairs = 49L, unique_admissions = 4078L,
                         events = 133L, patients = 1945L)

# Admission dates are Excel serial numbers or m/d/Y text
parse_admission_date <- function(x) {
  n <- suppressWarnings(as.numeric(x))
  d <- as.Date(x, format = "%m/%d/%Y")
  d[!is.na(n)] <- as.Date(n[!is.na(n)], origin = "1899-12-30")
  d
}

read_boscore <- function(path = "boscore.xlsx") {
  if (!file.exists(path)) {
    stop("'", path, "' not found. It is gitignored patient-level data: copy it into the project root.", call. = FALSE)
  }
  nm <- names(readxl::read_excel(path, n_max = 0, .name_repair = "minimal"))
  missing_cols <- setdiff(BOSCORE_KEEP, nm)
  if (length(missing_cols)) stop("Missing columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  if (anyDuplicated(nm[nm %in% BOSCORE_KEEP])) stop("A needed column name is duplicated in the sheet.", call. = FALSE)
  # Everything not needed (names, birth dates, free text) is skipped at read time
  raw <- suppressWarnings(readxl::read_excel(path, col_types = ifelse(nm %in% BOSCORE_KEEP, "text", "skip"),
                                             .name_repair = "minimal"))
  parse_boscore(raw)
}

parse_boscore <- function(raw) {
  num <- function(x) suppressWarnings(as.numeric(x))
  X <- vapply(names(BOSCORE_PREDICTORS), function(cn) num(raw[[cn]]), numeric(nrow(raw)))
  X <- matrix(X, nrow = nrow(raw), dimnames = list(NULL, unname(BOSCORE_PREDICTORS)))[, PREDICTORS, drop = FALSE]
  date <- parse_admission_date(raw$DRICEM)
  year <- num(raw$yearDRICEM)
  year[is.na(year)] <- as.integer(format(date, "%Y"))[is.na(year)]
  list(rows = data.frame(id = num(raw$id), admission = as.character(raw$IDRICO), patient = as.character(raw$IDIND),
                         date = date, year = year, case = num(raw$EORTC_PROV_PROB), pred = num(raw$pred),
                         stringsAsFactors = FALSE),
       X = X)
}

# Inclusion rule of the paper's analysis set: complete data on the seven BOSCORE items
boscore_included <- function(bc) stats::complete.cases(bc$X)

admission_key <- function(rows) ifelse(is.na(rows$admission), paste0("row-", rows$id), rows$admission)

# Admission numbers that appear on more than one included row
boscore_twin_pairs <- function(bc, rows = boscore_included(bc)) {
  d <- bc$rows[rows, , drop = FALSE]; X <- bc$X[rows, , drop = FALSE]
  key <- admission_key(d)
  dup <- unique(key[duplicated(key)])
  empty <- data.frame(admission = character(), n_rows = integer(), n_cases = numeric(),
                      same_patient = logical(), predictors_identical = logical())
  if (!length(dup)) return(empty)
  prof <- apply(X, 1, paste, collapse = "")
  do.call(rbind, lapply(dup, function(k) {
    r <- which(key == k)
    data.frame(admission = k, n_rows = length(r), n_cases = sum(d$case[r]),
               same_patient = length(unique(d$patient[r])) == 1,
               predictors_identical = length(unique(prof[r])) == 1)
  }))
}

# The pattern that makes the second row of a pair a twin: two rows, one event, same patient and profile
check_twins <- function(pairs) {
  list(all_pairs_have_two_rows = all(pairs$n_rows == 2),
       each_pair_has_one_event = all(pairs$n_cases == 1),
       same_patient = all(pairs$same_patient),
       predictors_identical = all(pairs$predictors_identical))
}

# The included rows, with the non-event twin of each pair flagged (the event row is kept)
boscore_cohort <- function(bc) {
  inc <- which(boscore_included(bc))
  d <- bc$rows[inc, , drop = FALSE]
  key <- admission_key(d)
  ord <- order(key, -d$case, d$id)
  keep <- ord[!duplicated(key[ord])]
  d$twin_duplicate <- !(seq_along(inc) %in% keep)
  list(rows = d, X = bc$X[inc, , drop = FALSE])
}

# The base cohort drops the twins; include_twins = TRUE gives the 4,127 rows of the paper
select_cohort <- function(cohort, include_twins = FALSE) {
  s <- if (include_twins) rep(TRUE, nrow(cohort$rows)) else !cohort$rows$twin_duplicate
  list(rows = cohort$rows[s, , drop = FALSE], X = cohort$X[s, , drop = FALSE])
}

cohort_summary <- function(cohort, pars) {
  lp <- as.numeric(pars$intercept + cohort$X %*% pars$coef[colnames(cohort$X)])
  list(admissions = nrow(cohort$X), patients = length(unique(cohort$rows$patient)),
       events = sum(cohort$rows$case), prevalence = mean(cohort$rows$case),
       mean_predicted_risk = mean(stats::plogis(lp)), apparent_auc = auc_rank(lp, cohort$rows$case))
}

# Patients by number of admissions, with the top category pooled so no cell is tiny
admissions_per_patient <- function(cohort, top = 8L) {
  k <- as.integer(table(cohort$rows$patient))
  k <- pmin(k, top)
  tb <- table(k)
  data.frame(admissions = ifelse(as.integer(names(tb)) >= top, paste0(top, "+"), names(tb)), patients = as.integer(tb))
}

# Rows drawn by resampling whole groups (patients) with replacement
resample_group_rows <- function(rows_by_group) {
  unlist(rows_by_group[sample.int(length(rows_by_group), replace = TRUE)], use.names = FALSE)
}

# Bootstrap of net benefit, MCID and flagged share. by = "patient" resamples patients with all their
# admissions (the primary method); by = "admission" resamples rows independently, for comparison.
bootstrap_dca <- function(cohort, risk, thresholds, B = 2000L, seed, by = c("patient", "admission")) {
  by <- match.arg(by)
  idx <- profile_index(cohort$X)
  n <- nrow(cohort$X)
  groups <- if (by == "patient") split(seq_len(n), cohort$rows$patient) else as.list(seq_len(n))
  set.seed(seed)
  res <- vapply(seq_len(B), function(b) {
    rows <- resample_group_rows(groups)
    s <- dca_exact(tabulate(idx[rows], 2L^length(PREDICTORS)) / length(rows), risk, thresholds = thresholds)
    c(s$nb, s$mcid_b, s$flagged)
  }, numeric(3 * length(thresholds)))
  res <- t(res)
  colnames(res) <- c(paste0("nb_", thresholds), paste0("mcid_", thresholds), paste0("flagged_", thresholds))
  res
}

# Compare the parent file with bs.xlsx (bs is the result of read_bs_data())
reconcile_boscore_bs <- function(bc, bs) {
  inc <- boscore_included(bc)
  m <- match(bs$id, bc$rows$id)
  found <- !anyNA(m)
  list(
    bs_ids_all_in_boscore = found,
    included_rows_are_the_bs_rows = setequal(bc$rows$id[inc], bs$id),
    predictors_identical = found && all(bc$X[m, PREDICTORS, drop = FALSE] == bs$X[, PREDICTORS, drop = FALSE]),
    outcome_identical = found && all(bc$rows$case[m] == bs$y),
    pred_missing_exactly_when_excluded = all(is.na(bc$rows$pred) == !inc)
  )
}

# Per predictor: the BOSCORE column that supplies it and how often it equals the bs.xlsx variable
boscore_predictor_agreement <- function(bc, bs) {
  m <- match(bs$id, bc$rows$id)
  data.frame(predictor = PREDICTORS,
             boscore_column = names(BOSCORE_PREDICTORS)[match(PREDICTORS, BOSCORE_PREDICTORS)],
             agreement = vapply(PREDICTORS, function(p) mean(bc$X[m, p] == bs$X[, p]), numeric(1)),
             row.names = NULL)
}

# A pseudonymised copy of the included rows (twins flagged) that downstream analysis can use without the
# identified file. Codes are a seeded random permutation; the year replaces the admission date.
pseudonymise_cohort <- function(cohort, seed) {
  set.seed(seed)
  code <- function(x, prefix) {
    u <- sort(unique(x))
    sprintf(paste0(prefix, "%0", nchar(as.character(length(u))), "d"), sample.int(length(u)))[match(x, u)]
  }
  d <- cohort$rows
  data.frame(patient_code = code(d$patient, "P"), admission_code = code(admission_key(d), "A"),
             admission_year = d$year, cohort$X, case = d$case, twin_duplicate = d$twin_duplicate,
             check.names = FALSE, row.names = NULL)
}

cohort_from_pseudonymised <- function(path) {
  d <- utils::read.csv(path, stringsAsFactors = FALSE)
  list(rows = data.frame(id = NA_real_, admission = d$admission_code, patient = d$patient_code,
                         date = as.Date(NA), year = d$admission_year, case = d$case, pred = NA_real_,
                         twin_duplicate = as.logical(d$twin_duplicate), stringsAsFactors = FALSE),
       X = as.matrix(d[, PREDICTORS]))
}

# Prefer the pseudonymised copy; otherwise build the cohort from the identified file
load_base_cohort <- function(restricted = "data/restricted/base_cohort_pseudonymised.csv", boscore = "boscore.xlsx") {
  if (file.exists(restricted)) return(cohort_from_pseudonymised(restricted))
  boscore_cohort(read_boscore(boscore))
}
