# PRISMA 2020 flow diagram (databases only). Run from the project root.
#
# Counts from prisma_summary.xlsx (7 November 2026; gitignored with the other spreadsheets, so they are typed in here with their derivation):
#   PubMed 267 and Embase 785 records (1,052); 221 duplicates removed (264 and 586 unique within database, 850 to 831 on title matching);
#   831 screened, 750 excluded at title and abstract (reasons not recorded: the automated screening decision file was never assembled);
#   81 reports sought; 18 not retrieved (12 with no full text, 6 whose link pointed to the wrong document); 63 assessed (56 verified, 7 re-extracted);
#   3 excluded with reasons (an eQTL study with no clinical outcome; a paediatric AML study with no SNP analysis of the fungal endpoint; S010, Ozdemir
#   2011, excluded on 6 November 2026, addendum 37: leukaemia patients compared with healthy controls, no infection-by-genotype analysis, abstract only,
#   which the spreadsheet still counts as included); 60 included.
# Search dates (not on the diagram, which has no box for them): PubMed 7 July 2026, Embase 8 October 2026.
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
  dbr_excluded              = "No invasive mould disease outcome (eQTL study), 1; No SNP analysis for the fungal endpoint (paediatric AML), 1; No infection-by-genotype analysis (leukaemia patients against healthy controls), 1",
  new_studies               = 60,
  new_reports               = 60
)
stopifnot(counts$database_results - counts$duplicates == counts$records_screened,
          counts$records_screened - counts$records_excluded == counts$dbr_sought_reports,
          counts$dbr_sought_reports - counts$dbr_notretrieved_reports == counts$dbr_assessed,
          counts$dbr_assessed - 3 == counts$new_studies)

tpl <- read.csv(system.file("extdata", "PRISMA.csv", package = "PRISMA2020"))
for (k in names(counts)) tpl$n[tpl$data == k] <- counts[[k]]

data <- PRISMA_data(tpl)
# Blue-grey colours to match the other figures (R/theme_jama.R). Graphviz colour names: the package does not quote hex codes, so they fail.
# The side boxes keep the package's light steel blue
flow <- PRISMA_flowdiagram(data, interactive = FALSE, previous = FALSE, other = FALSE, fontsize = 12, side_boxes = TRUE,
                           title_colour = "SlateGray3", greybox_colour = "Gainsboro", main_colour = "SteelBlue4", arrow_colour = "SteelBlue4")

dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)
write.csv(tpl[!is.na(tpl$data), c("data", "n")], "output/prisma_counts.csv", row.names = FALSE)
PRISMA_save(flow, filename = "output/figures/prisma_flow.html", filetype = "HTML", overwrite = TRUE)
# The PNG is rendered from the SVG at 300 dpi (7 inches wide; the package's own PNG export is a screen-resolution image)
source("R/figure_export.R")
svg_to_png(prisma_side_labels(DiagrammeRsvg::export_svg(flow)), "output/figures/prisma_flow.png", width_in = 7, dpi = 300)
flow
