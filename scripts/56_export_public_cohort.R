# Export the de-identified profile table that lets the simulation be rerun without patient data (R/public_cohort.R).
# Run from the project root on the machine that holds the patient-level data.
#
# Reads    data/restricted/base_cohort_pseudonymised.csv (or boscore.xlsx if that is absent)
# Writes   data/public/profile_admissions.csv       admissions per clinical profile; profiles with 1 to 4 admissions suppressed
#          data/public/profile_admissions_meta.csv  total admissions, minimum cell size, number of suppressed cells
#
# Rule (25 November 2026): no published cell below five admissions; outcomes, patients, dates and the allo-HSCT flag are not
# exported. check_public_table() refuses to write a table that breaks the rule.

source("R/baseline.R")
source("R/bs_data.R")
source("R/decision_curve.R")
source("R/boscore_data.R")
source("R/public_cohort.R")

full <- load_base_cohort()
if (is_public_cohort(full)) stop("This script exports from the patient-level cohort; only the public table was found.", call. = FALSE)
base <- select_cohort(full, include_twins = FALSE)
stopifnot(nrow(base$X) == BOSCORE_EXPECTED$unique_admissions)

out <- write_public_profile_table(base$X)
tab <- out$table
message("wrote ", PUBLIC_PROFILE_PATH, ": ", sum(tab$admissions > 0, na.rm = TRUE), " published profiles (",
        sum(tab$admissions, na.rm = TRUE), " admissions), ", sum(tab$suppressed), " suppressed (",
        nrow(base$X) - sum(tab$admissions, na.rm = TRUE), " admissions), ", sum(tab$admissions == 0, na.rm = TRUE), " empty")

# How far the imputed stand-in moves the profile distribution
q_true <- tabulate(profile_index(base$X), 2L^length(PREDICTORS)) / nrow(base$X)
q_pub <- tabulate(profile_index(cohort_from_public()$X), 2L^length(PREDICTORS)) / nrow(base$X)
message(sprintf("total variation distance between the true and the public profile distribution: %.5f", sum(abs(q_true - q_pub)) / 2))
