# PRISMA 2020 flow diagram (databases only). Run from the project root.
#
# Counts from 02_current_status.md, section 2: PubMed 267 records, Embase 785, 831 unique after title deduplication (so 221 removed), 81 included.
# Anything marked PLACEHOLDER is not recorded in the repository and must be filled in from Covidence before the figure is used.
# Placeholders print as "XX". Fill them in `counts` below, as numbers, then rerun.
#
# Writes   output/figures/prisma_flow.html, output/figures/prisma_flow.png, output/prisma_counts.csv

library(PRISMA2020)

PH <- "XX"

counts <- list(
  database_specific_results = "PubMed, 267; Embase, 785",
  database_results          = 267 + 785,
  register_results          = 0,
  duplicates                = 267 + 785 - 831,
  excluded_automatic        = 0,
  excluded_other            = 0,
  records_screened          = 831,
  records_excluded          = PH,   # PLACEHOLDER: excluded at title and abstract (Covidence)
  dbr_sought_reports        = PH,   # PLACEHOLDER: sent to full-text review
  dbr_notretrieved_reports  = PH,   # PLACEHOLDER: full text not retrieved (12 of the 81 included records had no PDF supplied; not confirmed as the same thing)
  dbr_assessed              = PH,   # PLACEHOLDER: full texts assessed for eligibility
  dbr_excluded              = PH,   # PLACEHOLDER: excluded at full text, with reasons ("Reason, n; Reason, n")
  new_studies               = 81,   # included records; the number of distinct studies is lower if duplicate reports are confirmed (status log, section 4)
  new_reports               = 81
)

tpl <- read.csv(system.file("extdata", "PRISMA.csv", package = "PRISMA2020"))
for (k in names(counts)) tpl$n[tpl$data == k] <- counts[[k]]

data <- PRISMA_data(tpl)
flow <- PRISMA_flowdiagram(data, interactive = FALSE, previous = FALSE, other = FALSE, fontsize = 12, side_boxes = TRUE)

dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)
write.csv(tpl[!is.na(tpl$data), c("data", "n")], "output/prisma_counts.csv", row.names = FALSE)
PRISMA_save(flow, filename = "output/figures/prisma_flow.html", filetype = "HTML", overwrite = TRUE)
PRISMA_save(flow, filename = "output/figures/prisma_flow.png", filetype = "PNG", overwrite = TRUE)
flow
