# Apply the student's study-level corrections to the bundle's study table. Run from the project root.
#
# Reads    data/study_corrections.csv and extraction_studies_round1.csv from the evidence bundle
# Writes   data/extraction_studies_corrected.csv (the bundle's study table with the corrections applied; needs_decision rows are left as they were)

source("R/evidence_mapping.R"); source("R/source_verification.R")

st <- read_bundle_table("extraction_studies_round1.csv")
res <- apply_study_corrections(st, read.csv("data/study_corrections.csv", stringsAsFactors = FALSE))
write.csv(res$data, "data/extraction_studies_corrected.csv", row.names = FALSE)
message(nrow(res$applied), " corrections applied; ", nrow(res$pending), " left pending")
print(res$pending, row.names = FALSE)
