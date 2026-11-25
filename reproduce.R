# Rerun the analysis from public inputs only and compare every regenerated output with the committed one.
#
#   Rscript reproduce.R                 every step that runs without patient data (several hours on one machine)
#   Rscript reproduce.R 16 23 27        only these scripts, in the order given
#   Rscript reproduce.R --list          show the plan and stop
#
# Inputs, all in this repository:
#   - the extraction bundle `IMD_genetics_analysis_bundle (1).zip` and the dated corrections in data/
#   - 1000 Genomes and gnomAD inputs in data/
#   - the published Stanzani 2019 equation and validation targets (data/baseline_*.csv)
#   - the de-identified cohort table data/public/profile_admissions.csv (admissions per clinical profile, cells of
#     1 to 4 admissions suppressed; R/public_cohort.R)
#   - the committed aggregate outputs of the restricted steps (below), which later steps read
#
# Restricted steps need the patient-level cohort and are not rerun here. Their committed outputs (aggregates only) are
# what the public steps read: the refit baseline equation (output/baseline_refit_parameters.csv) and the MCID
# (output/dca_baseline_refit_summary.csv) come from scripts 09 to 11.
#
# Without the patient files, load_base_cohort() builds the cohort from the public table. Suppressed profiles share their
# admissions equally, so results that depend on the profile mix can differ slightly from the committed ones; the
# report says by how much. Columns that need observed outcomes (for example the observed events of script 21) come out NA.
#
# Writes output/reproduce_report.csv: one row per output file the run rewrote, with the largest absolute and relative
# difference from the committed version. The run overwrites files in output/; `git checkout output` restores them.

RESTRICTED <- c("05", "06", "07", "08", "09", "10", "11", "19", "20")
# Script 13 wrote the first version of data/extraction_corrections.csv; later hand-checked corrections were added to that
# file, so it is an input here and script 13 is not rerun (it would overwrite the file)
NOT_RERUN <- c("13")
# Scripts that take PROSE_SCHEME; the manuscript reports the profile-preserving truth
PROFILE <- c("17", "18", "21", "22", "24", "26", "27", "34", "43", "47", "48", "49", "50", "51")
PLAN <- c("01", "02", "03", "04", "28",                         # correlation structure, decision curves, haplotypes
          "12", "14", "15", "45", "46",                         # evidence to parameters
          "16", "17", "18", "21", "22", "24", "26",             # genetic layer and Aim 3b
          "23", "25", "27", "29", "30",                         # Aim 4, figures, NRI, polygenic benchmark
          "31", "32", "33", "34", "35",                         # calibration checks
          "36", "37", "38", "39", "40", "41", "42",             # scheme comparison
          "43", "47", "48", "49", "50", "51",                   # Aim 4 pathogen specificity and stress tests
          "44", "52", "55")                                     # PRISMA flow, manuscript figures, threshold curve

script_path <- function(n) {
  f <- list.files("scripts", pattern = paste0("^", n, "_.*\\.R$"), full.names = TRUE)
  if (length(f) != 1L) stop("No unique script ", n, call. = FALSE)
  f
}

snapshot_outputs <- function() {
  files <- list.files("output", pattern = "\\.csv$", full.names = TRUE)
  setNames(lapply(files, function(f) list(mtime = file.mtime(f), data = utils::read.csv(f, check.names = FALSE))), files)
}

compare_csv <- function(old, new) {
  if (!identical(dim(old), dim(new)) || !identical(names(old), names(new))) return(c(abs = NA, rel = NA, same_shape = 0, na_changed = NA))
  num <- vapply(old, is.numeric, logical(1)) & vapply(new, is.numeric, logical(1))
  if (!any(num)) return(c(abs = 0, rel = 0, same_shape = 1, na_changed = 0))
  a <- as.matrix(old[num]); b <- as.matrix(new[num])
  d <- abs(a - b)
  rel <- d / pmax(abs(a), 1e-12)
  na_changed <- sum(is.na(a) != is.na(b))
  c(abs = suppressWarnings(max(d, na.rm = TRUE)), rel = suppressWarnings(max(rel, na.rm = TRUE)), same_shape = 1, na_changed = na_changed)
}

run_script <- function(n) {
  scheme <- if (n %in% PROFILE) "profile" else "global"
  old <- Sys.getenv("PROSE_SCHEME", unset = NA)
  Sys.setenv(PROSE_SCHEME = scheme)
  on.exit(if (is.na(old)) Sys.unsetenv("PROSE_SCHEME") else Sys.setenv(PROSE_SCHEME = old))
  t0 <- Sys.time()
  env <- new.env(parent = globalenv())
  ok <- tryCatch({ sys.source(script_path(n), envir = env); "ok" }, error = function(e) paste("error:", conditionMessage(e)))
  data.frame(script = n, scheme = scheme, status = ok, minutes = round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2))
}

main <- function(args = commandArgs(trailingOnly = TRUE)) {
  plan <- if (length(args) && !identical(args, "--list")) sprintf("%02d", as.integer(args)) else PLAN
  if (length(intersect(plan, NOT_RERUN))) stop("Script 13 would overwrite data/extraction_corrections.csv; it is not rerun.", call. = FALSE)
  bad <- intersect(plan, RESTRICTED)
  if (length(bad)) stop("Scripts ", paste(bad, collapse = ", "), " need the patient-level cohort and are not part of the public rerun.", call. = FALSE)
  message("Plan: ", paste(plan, collapse = " "), "\nRestricted, not rerun: ", paste(RESTRICTED, collapse = " "))
  if (identical(args, "--list")) return(invisible(plan))
  if (!file.exists("data/public/profile_admissions.csv")) stop("data/public/profile_admissions.csv is missing.", call. = FALSE)

  before <- snapshot_outputs()
  runs <- do.call(rbind, lapply(plan, function(n) { message("== script ", n); r <- run_script(n); message("   ", r$status); r }))

  changed <- Filter(function(f) !file.exists(f) || !(f %in% names(before)) || file.mtime(f) > before[[f]]$mtime,
                    list.files("output", pattern = "\\.csv$", full.names = TRUE))
  report <- do.call(rbind, lapply(changed, function(f) {
    if (!(f %in% names(before))) return(data.frame(file = f, status = "new", max_abs_diff = NA, max_rel_diff = NA, na_changed = NA))
    cmp <- compare_csv(before[[f]]$data, utils::read.csv(f, check.names = FALSE))
    status <- if (cmp[["same_shape"]] == 0) "shape differs" else if (cmp[["abs"]] <= 1e-8 && cmp[["na_changed"]] == 0) "reproduced" else "differs"
    data.frame(file = f, status = status, max_abs_diff = cmp[["abs"]], max_rel_diff = cmp[["rel"]], na_changed = cmp[["na_changed"]])
  }))
  utils::write.csv(runs, "output/reproduce_runs.csv", row.names = FALSE)
  if (!is.null(report)) utils::write.csv(report, "output/reproduce_report.csv", row.names = FALSE)
  message("\nScripts: ", sum(runs$status == "ok"), " of ", nrow(runs), " ran. Outputs rewritten: ", NROW(report),
          if (!is.null(report)) paste0(" (", paste(names(table(report$status)), table(report$status), sep = " ", collapse = ", "), ")"),
          ". See output/reproduce_report.csv.")
  invisible(list(runs = runs, report = report))
}

if (sys.nframe() == 0L) main()
