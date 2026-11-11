# Figures of the manuscript and the website

PNG files at 300 dpi in the blue-grey style of `R/theme_jama.R`, made by `scripts/53_export_figures.R` (run from the project root; it redraws the figures that sit
inside the pages and copies the ones made by `scripts/27`, `34`, `44` and `52`).

- `manuscript/`: the figures of the Draft manuscript page, numbered in the order of the page (`Figure_<n>_<id>.png`). The numbers are those of the rendered page.
- `website/`: the figures drawn inside the other pages of the site (baseline, decision curves, evidence, mapping), named `<page>_Figure_<n>_<id>.png`.
- `manifest.csv`: each file's source, pixel size, the resolution stamped in the file and the print size in inches that follows (pixels divided by 300).

Each file carries its resolution in the PNG header (pHYs chunk), so a viewer or a journal's checker reads 300 dpi at the print size in the manifest. The captions are in
`manuscript.qmd` and the page sources: the figures carry no title or caption of their own. Tables are not included.
