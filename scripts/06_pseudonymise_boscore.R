# Build a pseudonymised extract of the full BOSCORE export for analysis. Run from the project root.
#
# Reads    BOSCORE_recoded_2019.xlsx                       (gitignored; direct identifiers)
# Writes   data/restricted/boscore_pseudonymised.csv       (gitignored; still personal data)
#          data/boscore_pseudonymisation_manifest.csv      (tracked; column names and treatments only)
#
# The codes are a seeded random permutation of the sorted original identifiers, so the mapping can be
# regenerated from the original file and no key file needs to exist. Keep the original file where it is
# governed, and keep the extract out of version control and out of email.

source("R/pseudonymise.R")

SEED <- 6283
src <- "BOSCORE_recoded_2019.xlsx"
if (!file.exists(src)) stop("'", src, "' not found. It is gitignored patient-level data: copy it into the project root.", call. = FALSE)

orig <- readxl::read_excel(src, .name_repair = "minimal")
names(orig) <- make.unique(names(orig))

ps <- pseudonymise_boscore(orig, seed = SEED)

checks <- c(check_pseudonymised(ps$data), verify_pseudonymisation(orig, ps$data))
for (nm in names(checks)) message(sprintf("  %-34s %s", nm, if (isTRUE(checks[[nm]])) "pass" else "FAIL"))
if (!all(unlist(checks))) stop("The extract failed a check and was not written.", call. = FALSE)

dir.create("data/restricted", recursive = TRUE, showWarnings = FALSE)
write.csv(ps$data, "data/restricted/boscore_pseudonymised.csv", row.names = FALSE, na = "")
write.csv(ps$manifest, "data/boscore_pseudonymisation_manifest.csv", row.names = FALSE)

message("\nWrote ", nrow(ps$data), " rows and ", ncol(ps$data), " columns to data/restricted/boscore_pseudonymised.csv")
message("Columns dropped: ", sum(startsWith(ps$manifest$treatment, "dropped")),
        " | replaced or generalised: ", sum(grepl("^(replaced|kept as|generalised)", ps$manifest$treatment)),
        " | dates converted to day offsets: ", sum(startsWith(ps$manifest$treatment, "converted")))
