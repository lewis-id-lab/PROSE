# Reconcile boscore.xlsx (the parent file) with bs.xlsx and build the base cohort. Run from the project root.
#
# Reads    boscore.xlsx, bs.xlsx                                     (gitignored; patient-level)
# Writes   output/boscore_reconciliation.csv                         check, value, expected, pass
#          output/boscore_predictor_mapping.csv                      BOSCORE column behind each predictor
#          output/boscore_cohort_summary.csv                         rows (4,127) against unique admissions (4,078)
#          output/admissions_per_patient.csv                         top-coded at 8+
#          output/allo_hsct_summary.csv                              allo-HSCT admissions and events (counts only)
#          data/restricted/base_cohort_pseudonymised.csv             (gitignored) the cohort with codes, no identifiers
#
# The script stops if any check fails, so a changed file cannot slip through.

source("R/baseline.R")
source("R/bs_data.R")
source("R/decision_curve.R")
source("R/boscore_data.R")

pars <- load_baseline_params()
bc <- read_boscore("boscore.xlsx")
bs <- read_bs_data("bs.xlsx")
inc <- boscore_included(bc)

# 1. Predictor mapping, outcome and inclusion rule against bs.xlsx
rec <- reconcile_boscore_bs(bc, bs)
mapping <- boscore_predictor_agreement(bc, bs)

# 2. The twin admissions
pairs <- boscore_twin_pairs(bc, inc)
twins <- check_twins(pairs)

# 3. Cohorts
cohort <- boscore_cohort(bc)
rows_cohort <- select_cohort(cohort, include_twins = TRUE)
base_cohort <- select_cohort(cohort, include_twins = FALSE)
s_rows <- cohort_summary(rows_cohort, pars)
s_base <- cohort_summary(base_cohort, pars)

E <- BOSCORE_EXPECTED
checks <- data.frame(
  check = c("rows in boscore.xlsx", "ids run 1 to the number of rows", "rows complete on the seven BOSCORE items",
            "admission numbers on two rows", "unique admissions after removing twins", "IMD events in the base cohort",
            "patients in the base cohort", "allo-HSCT admissions among the 4,127 rows", "events among them",
            names(rec), paste0("twin_pairs_", names(twins))),
  value = c(nrow(bc$rows), all(sort(bc$rows$id) == seq_len(nrow(bc$rows))), sum(inc), nrow(pairs), s_base$admissions, s_base$events,
            s_base$patients, allo_summary(rows_cohort)$allo_admissions, allo_summary(rows_cohort)$allo_events, unlist(rec), unlist(twins)),
  expected = c(E$rows, TRUE, E$included, E$twin_pairs, E$unique_admissions, E$events, E$patients, E$allo_rows, E$allo_events, rep(TRUE, length(rec) + length(twins))),
  stringsAsFactors = FALSE
)
checks$pass <- as.numeric(checks$value) == as.numeric(checks$expected)

dir.create("output", showWarnings = FALSE)
write.csv(checks, "output/boscore_reconciliation.csv", row.names = FALSE)
write.csv(mapping, "output/boscore_predictor_mapping.csv", row.names = FALSE)
summary_tab <- rbind(
  data.frame(cohort = "4,127 rows (with the 49 non-event twins)", as.data.frame(s_rows)),
  data.frame(cohort = "4,078 unique admissions (base cohort)", as.data.frame(s_base))
)
write.csv(summary_tab, "output/boscore_cohort_summary.csv", row.names = FALSE)
write.csv(admissions_per_patient(base_cohort), "output/admissions_per_patient.csv", row.names = FALSE)
allo_tab <- rbind(data.frame(cohort = "4,127 rows", as.data.frame(allo_summary(rows_cohort))), data.frame(cohort = "4,078 unique admissions (base cohort)", as.data.frame(allo_summary(base_cohort))))
write.csv(allo_tab, "output/allo_hsct_summary.csv", row.names = FALSE)

if (!all(checks$pass)) {
  print(checks[!checks$pass, ], row.names = FALSE)
  stop("One or more reconciliation checks failed; the pseudonymised cohort was not written.", call. = FALSE)
}

# 4. The pseudonymised cohort (all 4,127 included rows, twins flagged), for analysis without the identified file
ps <- pseudonymise_cohort(cohort, seed = 4318)
dir.create("data/restricted", recursive = TRUE, showWarnings = FALSE)
write.csv(ps, "data/restricted/base_cohort_pseudonymised.csv", row.names = FALSE)

message("All ", nrow(checks), " checks passed.")
message("Predictor agreement with bs.xlsx: ", paste(sprintf("%s %.0f%%", mapping$boscore_column, 100 * mapping$agreement), collapse = ", "))
print(summary_tab, row.names = FALSE)
message("\nWrote data/restricted/base_cohort_pseudonymised.csv (", nrow(ps), " rows, twins flagged)")
