# Check the admission-level derivation data (bs.xlsx) against the published Stanzani 2019 counts
# and coefficients. bs.xlsx is gitignored patient-level data: copy it into the project root.
# Only aggregates are written to output/. Run from the project root.

source("R/baseline.R")
source("R/bs_data.R")

pars <- load_baseline_params()
targets <- load_validation_targets()
bs <- read_bs_data("bs.xlsx")

n_admissions <- length(bs$y)
checks <- list()

# 1. Structure
checks$rows <- n_admissions == pars$n_case + pars$n_non
checks$outcome <- sum(bs$y) == pars$n_case && sum(bs$y == 0) == pars$n_non
checks$id_is_admission_number <- !anyDuplicated(bs$id)   # unique per row, so patients cannot be recovered

# 2. Counts against Table 2
counts <- check_bs_counts(bs, pars)
checks$counts_as_documented <- bs_counts_as_expected(counts)

# 3. Refit against the published coefficients (verification only; the study model stays fixed)
refit <- check_bs_refit(bs, pars)
checks$coefficients_within_0.01 <- all(abs(refit$diff) < 0.01)

# 4. The published equation applied to the file
perf <- bs_apparent_performance(bs, pars)
checks$apparent_auc_near_published <- abs(perf$apparent_auc - targets$aroc) < 0.005
checks$prevalence_consistent <- abs(perf$mean_predicted_risk - perf$observed_prevalence) < 0.002

dir.create("output", showWarnings = FALSE)
write.csv(counts, "output/bs_data_check_counts.csv", row.names = FALSE)
write.csv(refit, "output/bs_data_check_coefficients.csv", row.names = FALSE)

message("\nRows: ", n_admissions, " | cases: ", perf$cases, " | observed prevalence: ",
        sprintf("%.2f%%", 100 * perf$observed_prevalence),
        " | mean predicted risk: ", sprintf("%.2f%%", 100 * perf$mean_predicted_risk),
        " | apparent aROC: ", sprintf("%.4f", perf$apparent_auc))
message("\nCounts that differ from Table 2 (data minus published):")
differ <- counts[counts$diff_cases != 0 | counts$diff_non_cases != 0, ]
if (nrow(differ)) print(differ, row.names = FALSE) else message("  none")
message("\nLargest coefficient difference from the published equation: ",
        sprintf("%.4f", max(abs(refit$diff))))
message("\nChecks:")
for (nm in names(checks)) message(sprintf("  %-32s %s", nm, if (checks[[nm]]) "pass" else "FAIL"))
if (!all(unlist(checks))) stop("One or more checks failed.", call. = FALSE)
message("\nAll checks passed.")
