# PRISMA 2020 flow diagram (databases only). Run from the project root.
#
# Counts from the Covidence export and the thesis (19 November 2026). They are typed in here with their derivation, because spreadsheets are gitignored:
#   PubMed 267 and Embase 785 records (1,052); 221 duplicates removed before screening;
#   622 excluded by the automated (AI-assisted) title-and-abstract screen, run by the supervisor alone; 209 imported into Covidence;
#   1 excluded on title and abstract in Covidence; 208 full texts assessed by two reviewers; 127 excluded at full text (reasons: 59 review articles, 24 wrong
#   patient population, 15 conference proceedings; the other 29 have NO REASON STATED YET and must be filled from Covidence); 81 included.
# Extraction stage (the box added below the package's flow): 18 of the 81 had no usable source document (12 with no full text, 6 linked to the wrong document),
#   63 were extracted, 3 were excluded when read for extraction (the third, S010, Ozdemir 2011, on 6 November 2026, addendum 37), 60 form the evidence base.
# The 440 records from citation chasing (6 October 2026) did not enter screening and are not in the diagram.
# Search dates (not on the diagram, which has no box for them): PubMed 7 July 2026, Embase 8 October 2026.
#
# Writes   output/figures/prisma_flow.html, output/figures/prisma_flow.png, output/prisma_counts.csv

library(PRISMA2020)

counts <- list(
  database_specific_results = "PubMed, 267; Embase, 785",
  database_results          = 267 + 785,
  register_results          = 0,
  duplicates                = 221,
  excluded_automatic        = 622,
  excluded_other            = 0,
  records_screened          = 209,
  records_excluded          = 1,
  dbr_sought_reports        = 208,
  dbr_notretrieved_reports  = 0,
  dbr_assessed              = 208,
  dbr_excluded              = "Review article, 59; Wrong patient population, 24; Conference proceedings, 15; Reason not yet recorded, 29",
  new_studies               = 81,
  new_reports               = 81
)
extraction <- c(not_matched = 18, extracted = 63, excluded = 3, evidence_base = 60)
stopifnot(counts$database_results - counts$duplicates - counts$excluded_automatic == counts$records_screened,
          counts$records_screened - counts$records_excluded == counts$dbr_sought_reports,
          counts$dbr_sought_reports - counts$dbr_notretrieved_reports == counts$dbr_assessed,
          59 + 24 + 15 + 29 == counts$dbr_assessed - counts$new_studies,
          counts$new_studies - extraction[["not_matched"]] == extraction[["extracted"]],
          extraction[["extracted"]] - extraction[["excluded"]] == extraction[["evidence_base"]])

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
extraction_lines <- c("Data extraction",
                      sprintf("Studies sent for extraction (n = %d)", counts$new_studies),
                      sprintf("No usable source document (n = %d)", extraction[["not_matched"]]),
                      sprintf("Excluded at extraction (n = %d)", extraction[["excluded"]]),
                      sprintf("Studies in the evidence base (n = %d)", extraction[["evidence_base"]]))
svg_to_png(prisma_side_labels(prisma_extraction_box(DiagrammeRsvg::export_svg(flow), extraction_lines)), "output/figures/prisma_flow.png", width_in = 7, dpi = 300)
flow
