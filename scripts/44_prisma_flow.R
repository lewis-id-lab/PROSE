# PRISMA 2020 flow diagram (databases only). Run from the project root.
#
# Counts from prisma_summary.xlsx (7 November 2026; gitignored with the other spreadsheets, so they are typed in here with their derivation):
#   PubMed 267 and Embase 785 records (1,052); 221 duplicates removed (264 and 586 unique within database, 850 to 831 on title matching);
#   831 screened, 750 excluded at title and abstract (reasons not recorded: the automated screening decision file was never assembled);
#   81 reports sought; 18 not retrieved (12 with no full text, 6 whose link pointed to the wrong document); 63 assessed (56 verified, 7 re-extracted);
#   2 excluded with reasons (an eQTL study with no clinical outcome; a paediatric AML study with no SNP analysis of the fungal endpoint); 61 included.
# Still open: the PubMed search date is unrecorded (Embase was 8 October 2026), and S010 is counted as included here although addendum 37 excludes it
# from the evidence base (it would be 60 included and 3 excluded at full text). The diagram has no box for either.
#
# Writes   output/figures/prisma_flow.html, output/figures/prisma_flow.png, output/prisma_counts.csv

library(PRISMA2020)

counts <- list(
  database_specific_results = "PubMed, 267; Embase, 785",
  database_results          = 267 + 785,
  register_results          = 0,
  duplicates                = 221,
  excluded_automatic        = 0,
  excluded_other            = 0,
  records_screened          = 831,
  records_excluded          = 750,
  dbr_sought_reports        = 81,
  dbr_notretrieved_reports  = 18,
  dbr_assessed              = 63,
  dbr_excluded              = "No invasive mould disease outcome (eQTL study), 1; No SNP analysis for the fungal endpoint (paediatric AML), 1",
  new_studies               = 61,
  new_reports               = 61
)
stopifnot(counts$database_results - counts$duplicates == counts$records_screened,
          counts$records_screened - counts$records_excluded == counts$dbr_sought_reports,
          counts$dbr_sought_reports - counts$dbr_notretrieved_reports == counts$dbr_assessed,
          counts$dbr_assessed - 2 == counts$new_studies)

tpl <- read.csv(system.file("extdata", "PRISMA.csv", package = "PRISMA2020"))
for (k in names(counts)) tpl$n[tpl$data == k] <- counts[[k]]

data <- PRISMA_data(tpl)
flow <- PRISMA_flowdiagram(data, interactive = FALSE, previous = FALSE, other = FALSE, fontsize = 12, side_boxes = TRUE)

dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)
write.csv(tpl[!is.na(tpl$data), c("data", "n")], "output/prisma_counts.csv", row.names = FALSE)
PRISMA_save(flow, filename = "output/figures/prisma_flow.html", filetype = "HTML", overwrite = TRUE)
PRISMA_save(flow, filename = "output/figures/prisma_flow.png", filetype = "PNG", overwrite = TRUE)
flow
