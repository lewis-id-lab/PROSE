# Pseudonymisation of the full BOSCORE database export (BOSCORE_recoded_2019.xlsx).
#
# The export holds direct identifiers (name, birth date, hospital patient and admission numbers) and
# calendar dates. The extract keeps the structure that matters for analysis (which rows belong to the
# same patient, admission and episode) and every clinical variable, and removes or generalises the rest:
#   names and birth dates      dropped
#   patient / admission / episode / ward identifiers    replaced by random codes
#   admission date             generalised to the year
#   every other date           converted to days from the admission date
#   age                        top-coded at 90
# Pseudonymised data are still personal data: whoever holds the original file can re-link them.

BOSCORE_REQUIRED <- c("Reparto", "IDRICO", "IDIND", "Connome", "Dnas", "AGEAMM", "DRICEM", "DDIMEM")
BOSCORE_DATES_POSIX <- c("DTMO", "DEPISOD", "Expr1064", "TabMicetiIsolatiDISOL", "TabBatteriIsolatiDISOL")
BOSCORE_DATES_CHAR  <- c("DMOR", "DLFUP", "DSTATUS")              # dd/mm/yyyy text
BOSCORE_DROPPED <- c("Connome", "Dnas", "EMODES", "IDMICE", "IDBACT", "DDIMEM")
BOSCORE_REPLACED <- c("IDIND", "IDRICO", "IDEPISOD", "Reparto", "AGEAMM", "DRICEM")

pseudonymise_boscore <- function(d, seed) {
  d <- as.data.frame(d)
  stopifnot(all(BOSCORE_REQUIRED %in% names(d)))
  set.seed(seed)

  make_code <- function(x, prefix) {
    u <- sort(unique(x[!is.na(x)]))
    width <- nchar(as.character(length(u)))
    codes <- sprintf(paste0(prefix, "%0", width, "d"), sample.int(length(u)))
    codes[match(x, u)]
  }
  adm <- as.Date(d$DRICEM)
  offset <- function(x) {
    dt <- if (inherits(x, "POSIXct")) as.Date(x) else as.Date(x, format = "%d/%m/%Y")
    as.integer(dt - adm)
  }

  out <- data.frame(
    patient_code = make_code(d$IDIND, "P"),
    admission_code = make_code(d$IDRICO, "A"),
    episode_code = if ("IDEPISOD" %in% names(d)) make_code(d$IDEPISOD, "E") else NA_character_,
    ward_code = make_code(d$Reparto, "W"),
    admission_year = as.integer(format(adm, "%Y")),
    age_years = pmin(d$AGEAMM, 90),
    stringsAsFactors = FALSE
  )
  date_cols <- intersect(c(BOSCORE_DATES_POSIX, BOSCORE_DATES_CHAR), names(d))
  keep <- setdiff(names(d), c(BOSCORE_DROPPED, BOSCORE_REPLACED, date_cols))
  out <- cbind(out, d[, keep, drop = FALSE])
  for (v in date_cols) out[[paste0("days_from_admission_", v)]] <- offset(d[[v]])

  manifest <- data.frame(column = names(d), stringsAsFactors = FALSE)
  manifest$treatment <- vapply(manifest$column, function(v) {
    if (v %in% c("Connome", "Dnas")) "dropped: direct identifier"
    else if (v == "EMODES") "dropped: free-text diagnosis (coded version kept)"
    else if (v %in% c("IDMICE", "IDBACT")) "dropped: database record key"
    else if (v == "DDIMEM") "dropped: length of stay is kept in GGEM"
    else if (v %in% c("IDIND", "IDRICO", "IDEPISOD")) "replaced by a random code"
    else if (v == "Reparto") "replaced by a random ward code"
    else if (v == "AGEAMM") "kept as age_years, top-coded at 90"
    else if (v == "DRICEM") "generalised to admission_year"
    else if (v %in% date_cols) paste0("converted to days_from_admission_", v)
    else "kept unchanged"
  }, "")
  list(data = out, manifest = manifest)
}

# Does the extract still contain anything it should not? All checks must be TRUE.
check_pseudonymised <- function(x) {
  # A date uses the same separator twice (dd/mm/yyyy, dd-mm-yyyy, dd.mm.yyyy) or is ISO. Requiring the
  # same separator avoids flagging category labels such as the temperature range "38.6-40°C".
  date_like <- "(^|[^0-9])[0-9]{1,2}([/.-])[0-9]{1,2}\\2[0-9]{2,4}($|[^0-9])|[0-9]{4}-[0-9]{2}-[0-9]{2}"
  chr <- names(x)[vapply(x, is.character, logical(1))]
  list(
    no_date_columns = !any(vapply(x, function(c) inherits(c, c("Date", "POSIXct")), logical(1))),
    no_identifier_columns = !any(c("Connome", "Dnas", "IDIND", "IDRICO", "IDEPISOD", "Reparto", "DRICEM",
                                   "DDIMEM", "AGEAMM", "IDMICE", "IDBACT") %in% names(x)),
    no_date_strings = !any(vapply(chr, function(v) any(grepl(date_like, x[[v]], perl = TRUE), na.rm = TRUE),
                                  logical(1))),
    codes_well_formed = all(grepl("^P[0-9]+$", x$patient_code)) && all(grepl("^A[0-9]+$", x$admission_code)) &&
      all(grepl("^W[0-9]+$", x$ward_code)) &&
      all(grepl("^E[0-9]+$", x$episode_code[!is.na(x$episode_code) & nzchar(x$episode_code)])),
    age_top_coded = max(x$age_years, na.rm = TRUE) <= 90
  )
}

# Is the structure preserved, and is every mapping one to one?
verify_pseudonymisation <- function(orig, ps) {
  orig <- as.data.frame(orig)
  one_to_one <- function(a, b) {
    ok <- !is.na(a)
    n <- length(unique(a[ok]))
    n == length(unique(b[ok])) && n == nrow(unique(data.frame(a = a[ok], b = b[ok])))
  }
  adm_per_patient <- function(patient, admission) {
    u <- unique(data.frame(p = patient, a = admission)); as.vector(table(table(u$p)))
  }
  rows_per_adm <- function(admission) as.vector(table(table(admission)))
  list(
    same_rows = nrow(orig) == nrow(ps),
    patients_one_to_one = one_to_one(orig$IDIND, ps$patient_code),
    admissions_one_to_one = one_to_one(orig$IDRICO, ps$admission_code),
    episodes_one_to_one = !("IDEPISOD" %in% names(orig)) || one_to_one(orig$IDEPISOD, ps$episode_code),
    admissions_per_patient_preserved = identical(adm_per_patient(orig$IDIND, orig$IDRICO),
                                                 adm_per_patient(ps$patient_code, ps$admission_code)),
    rows_per_admission_preserved = identical(rows_per_adm(orig$IDRICO), rows_per_adm(ps$admission_code)),
    no_values_lost = all(vapply(intersect(c(BOSCORE_DATES_POSIX, BOSCORE_DATES_CHAR), names(orig)), function(v) {
      sum(!is.na(orig[[v]])) == sum(!is.na(ps[[paste0("days_from_admission_", v)]]))
    }, logical(1)))
  )
}
