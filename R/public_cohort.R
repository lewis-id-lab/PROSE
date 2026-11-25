# The public, de-identified view of the base cohort: admissions per clinical profile, with small cells suppressed.
#
# The simulation uses the cohort only through the distribution of its 4,078 admissions over the 128 profiles of
# the seven binary predictors (the outcome is simulated from the fixed equation, never read from the cohort).
# This file holds that distribution in a form that can be published:
#
#   data/public/profile_admissions.csv       one row per profile: the seven predictors, admissions, suppressed
#   data/public/profile_admissions_meta.csv  total admissions, the minimum cell size and the number of suppressed cells
#
# A profile with 1 to 4 admissions is published as suppressed (admissions = NA); an empty profile is published as 0.
# No outcome, patient, date, ward or transplant information is written. Event counts are not published because
# 55 of the 133 events sit in profiles with fewer than five events (decision of 25 November 2026).
#
# When neither data/restricted/ nor boscore.xlsx is present, load_base_cohort() rebuilds a stand-in cohort from
# the public table: suppressed cells share the suppressed admissions equally (integers, the remainder going to the
# lowest profile indices), so the total is exact and every published cell is exact. Outcomes, patients and the
# allo-HSCT flag are NA in that cohort, and scripts that need them stop through require_patient_data().
#
# Requires R/baseline.R (PREDICTORS, all_profiles, profile_index).

PUBLIC_MIN_CELL <- 5L
PUBLIC_PROFILE_PATH <- "data/public/profile_admissions.csv"
PUBLIC_META_PATH <- "data/public/profile_admissions_meta.csv"

# The publishable table for a cohort (a list with X, as from select_cohort())
public_profile_table <- function(X, min_cell = PUBLIC_MIN_CELL) {
  stopifnot(identical(colnames(X), PREDICTORS))
  n <- tabulate(profile_index(X), 2L^length(PREDICTORS))
  suppressed <- n > 0L & n < min_cell
  prof <- all_profiles()
  storage.mode(prof) <- "integer"
  data.frame(profile = seq_along(n), prof, admissions = ifelse(suppressed, NA_integer_, n), suppressed = suppressed,
             check.names = FALSE, row.names = NULL)
}

public_profile_meta <- function(X, tab, min_cell = PUBLIC_MIN_CELL) {
  data.frame(key = c("total_admissions", "min_cell", "suppressed_cells", "profiles"),
             value = c(nrow(X), min_cell, sum(tab$suppressed), nrow(tab)))
}

# Refuses a table that would publish a cell below the minimum (the check run before anything is written)
check_public_table <- function(tab, meta) {
  m <- setNames(meta$value, meta$key)
  shown <- tab$admissions[!tab$suppressed]
  problems <- c(
    if (any(is.na(shown))) "a cell that is not suppressed has no count",
    if (any(shown > 0 & shown < m[["min_cell"]])) "a published cell is below the minimum cell size",
    if (any(!is.na(tab$admissions[tab$suppressed]))) "a suppressed cell carries a count",
    if (sum(tab$suppressed) == 1L) "exactly one suppressed cell: its count follows from the total",
    if (sum(shown) > m[["total_admissions"]]) "published cells exceed the total",
    if (any(tab$suppressed) && m[["total_admissions"]] - sum(shown) < sum(tab$suppressed)) "suppressed total is smaller than the number of suppressed cells",
    if (!identical(names(tab), c("profile", PREDICTORS, "admissions", "suppressed"))) "unexpected columns",
    if (nrow(tab) != 2L^length(PREDICTORS)) "not one row per profile")
  if (length(problems)) stop("Public profile table refused: ", paste(problems, collapse = "; "), call. = FALSE)
  invisible(TRUE)
}

write_public_profile_table <- function(X, path = PUBLIC_PROFILE_PATH, meta_path = PUBLIC_META_PATH, min_cell = PUBLIC_MIN_CELL) {
  tab <- public_profile_table(X, min_cell)
  meta <- public_profile_meta(X, tab, min_cell)
  check_public_table(tab, meta)
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  utils::write.csv(tab, path, row.names = FALSE, na = "")
  utils::write.csv(meta, meta_path, row.names = FALSE)
  invisible(list(table = tab, meta = meta))
}

read_public_profile_table <- function(path = PUBLIC_PROFILE_PATH, meta_path = PUBLIC_META_PATH) {
  tab <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  tab$suppressed <- as.logical(tab$suppressed)
  meta <- utils::read.csv(meta_path, stringsAsFactors = FALSE)
  check_public_table(tab, meta)
  list(table = tab, meta = meta)
}

# Admissions per profile with the suppressed admissions shared equally over the suppressed cells
impute_profile_counts <- function(tab, total) {
  n <- tab$admissions
  s <- which(tab$suppressed)
  rest <- total - sum(n[-s], na.rm = TRUE)
  if (length(s)) {
    each <- rest %/% length(s)
    n[s] <- each + (seq_along(s) <= rest %% length(s))
  }
  as.integer(n)
}

# A stand-in for load_base_cohort(): one row per admission, predictors only
cohort_from_public <- function(path = PUBLIC_PROFILE_PATH, meta_path = PUBLIC_META_PATH) {
  pub <- read_public_profile_table(path, meta_path)
  total <- pub$meta$value[pub$meta$key == "total_admissions"]
  n <- impute_profile_counts(pub$table, total)
  X <- all_profiles()[rep(seq_along(n), n), , drop = FALSE]
  storage.mode(X) <- "double"
  rows <- data.frame(id = NA_real_, admission = sprintf("PUB%04d", seq_len(nrow(X))), patient = NA_character_,
                     date = as.Date(NA), year = NA_real_, case = NA_real_, pred = NA_real_, twin_duplicate = FALSE,
                     allo_hsct = NA, stringsAsFactors = FALSE)
  structure(list(rows = rows, X = X), public = TRUE)
}

is_public_cohort <- function(cohort) isTRUE(attr(cohort, "public"))

# Stops a script that needs outcomes, patients, twins or the allo-HSCT flag when only the public table is present
require_patient_data <- function(cohort, what) {
  if (is_public_cohort(cohort)) {
    stop(what, " needs the patient-level cohort (data/restricted/ or boscore.xlsx); only the public profile table is present. ",
         "The committed outputs of this step are used instead (see reproduce.R).", call. = FALSE)
  }
  invisible(TRUE)
}
