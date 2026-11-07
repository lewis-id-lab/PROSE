# Checks of extracted rows against their source papers, and the corrections that follow from them.
#
# A correction is a recorded change to one field of one extraction row, with the evidence for it. The extraction file
# itself (inside the evidence bundle) is never edited; corrections are applied on top of it when the mapping runs, and
# each one checks that the value it replaces is still the one it expected.

# Crude odds ratio of group A against group B (events and non-events), 0.5 added to every cell if any is zero
crude_or <- function(a_events, a_nonevents, b_events, b_nonevents) {
  v <- c(a_events, a_nonevents, b_events, b_nonevents)
  if (any(v == 0)) v <- v + 0.5
  (v[1] / v[2]) / (v[3] / v[4])
}

# Which group does a printed odds ratio belong to? Compares the printed value with the crude odds ratio of the
# reference-labelled group against the other group, and with its reciprocal. `decisive` is FALSE when the crude odds
# ratio is too close to 1 for the two readings to be told apart (outside 0.82 to 1.22 is decisive).
printed_or_direction <- function(ref_events, ref_nonevents, other_events, other_nonevents, printed_or) {
  crude <- crude_or(ref_events, ref_nonevents, other_events, other_nonevents)
  d_ref <- abs(log(printed_or) - log(crude))
  d_other <- abs(log(printed_or) - log(1 / crude))
  data.frame(crude_ref_vs_other = crude,
             reading = ifelse(d_ref < d_other, "reference-labelled group is the numerator", "other group is the numerator"),
             decisive = abs(log(crude)) > 0.2,
             distance_if_reference_numerator = d_ref, distance_if_other_numerator = d_other)
}

# Observed share of a group against Hardy-Weinberg expectations for the minor-allele and major-allele homozygotes
hwe_group_check <- function(group_total, n_total, minor_allele_freq) {
  obs <- group_total / n_total
  z <- function(p) (obs - p) / sqrt(p * (1 - p) / n_total)
  data.frame(observed = obs, expected_minor_homozygotes = minor_allele_freq^2, expected_major_homozygotes = (1 - minor_allele_freq)^2,
             z_minor_homozygotes = z(minor_allele_freq^2), z_major_homozygotes = z((1 - minor_allele_freq)^2))
}

# Apply corrections to the extraction. Each correction names a row and a field, the value expected there and the new value.
apply_corrections <- function(vr, corrections) {
  if (is.null(corrections) || !nrow(corrections)) return(list(data = vr, applied = character()))
  for (i in seq_len(nrow(corrections))) {
    cr <- corrections[i, ]
    j <- which(vr$variant_row_id == cr$variant_row_id)
    if (length(j) != 1) stop("correction for ", cr$variant_row_id, ": row not found or not unique", call. = FALSE)
    if (!(cr$field %in% names(vr))) stop("correction for ", cr$variant_row_id, ": no field ", cr$field, call. = FALSE)
    current <- vr[[cr$field]][j]
    if (!identical(as.character(current), as.character(cr$original_value))) {
      stop("correction for ", cr$variant_row_id, " expected '", cr$original_value, "' in ", cr$field, " but found '", current, "'", call. = FALSE)
    }
    vr[[cr$field]][j] <- cr$corrected_value
  }
  list(data = vr, applied = unique(corrections$variant_row_id))
}

# The comparison text with the reversed reference marker removed, and only the genotype contrast kept (the text before
# the first comma). "AA (ref) vs AG+GG combined, IA vs no IA" becomes "AA vs AG+GG combined".
clean_reversed_comparison <- function(text) {
  x <- sub(",.*$", "", text)
  trimws(gsub("\\s*\\(ref\\)", "", x))
}

# Apply study-level corrections (data/study_corrections.csv) to the study table. Each correction names a study and a field, the value expected
# there and the new value. Only rows with status "corrected" are applied; "needs_decision" rows are returned untouched in `pending`. An
# expected value that is blank in the correction must be blank in the table, and numbers are compared as numbers.
apply_study_corrections <- function(st, corrections) {
  if (is.null(corrections) || !nrow(corrections)) return(list(data = st, applied = data.frame(), pending = data.frame()))
  same <- function(a, b) {
    blank <- function(x) is.na(x) || !nzchar(trimws(as.character(x)))
    if (blank(a) || blank(b)) return(blank(a) && blank(b))
    na <- suppressWarnings(as.numeric(a)); nb <- suppressWarnings(as.numeric(b))
    if (!is.na(na) && !is.na(nb)) na == nb else identical(trimws(as.character(a)), trimws(as.character(b)))
  }
  todo <- corrections[corrections$status == "corrected", , drop = FALSE]
  for (i in seq_len(nrow(todo))) {
    cr <- todo[i, ]
    j <- which(st$study_id == cr$study_id)
    if (length(j) != 1) stop("correction for ", cr$study_id, ": study not found or not unique", call. = FALSE)
    if (!(cr$field %in% names(st))) stop("correction for ", cr$study_id, ": no field ", cr$field, call. = FALSE)
    if (!same(st[[cr$field]][j], cr$original_value)) {
      stop("correction for ", cr$study_id, " expected '", cr$original_value, "' in ", cr$field, " but found '", st[[cr$field]][j], "'", call. = FALSE)
    }
    st[[cr$field]] <- as.character(st[[cr$field]])
    st[[cr$field]][j] <- as.character(cr$corrected_value)
  }
  list(data = st, applied = todo[c("study_id", "field", "original_value", "corrected_value")],
       pending = corrections[corrections$status != "corrected", c("study_id", "field", "status"), drop = FALSE])
}

# Drop rows that duplicate another study's rows (duplicate_group "dup_of_<study>", set in data/extraction_corrections.csv), so the same
# result is not pooled twice. Run after apply_corrections(). Rows with other duplicate_group values (D001 and so on) are kept.
drop_duplicate_reports <- function(vr) {
  if (!("duplicate_group" %in% names(vr))) stop("no duplicate_group column", call. = FALSE)
  vr[is.na(vr$duplicate_group) | !grepl("^dup_of_", vr$duplicate_group), , drop = FALSE]
}

# Add extraction rows that are not in the round-1 bundle (data/extraction_added_rows.csv: rows read from a source PDF by hand, with their own
# variant_row_id, here L0001 and so on). The file must have exactly the bundle's columns, and its ids must be new. Numeric columns of the extraction
# stay numeric. Called after drop_duplicate_reports() in scripts 12, 14, 15, 21, 22 and 26 (addendum 37).
add_extraction_rows <- function(vr, path = "data/extraction_added_rows.csv") {
  add <- utils::read.csv(path, stringsAsFactors = FALSE, na.strings = c("", "NA"), colClasses = "character")
  if (!setequal(names(add), names(vr))) stop("the added rows must have the extraction's columns: ", paste(setdiff(union(names(add), names(vr)), intersect(names(add), names(vr))), collapse = ", "), call. = FALSE)
  if (anyDuplicated(add$variant_row_id) || any(add$variant_row_id %in% vr$variant_row_id)) stop("added variant_row_id values must be new and unique", call. = FALSE)
  add <- add[names(vr)]
  for (n in names(vr)) add[[n]] <- if (is.numeric(vr[[n]])) suppressWarnings(as.numeric(add[[n]])) else if (is.logical(vr[[n]])) as.logical(add[[n]]) else as.character(add[[n]])
  rbind(vr, add)
}
