# Figures of the journal submission (submission/) as SVG, PDF and PNG (300 dpi), numbered as in submission/manuscript.qmd and submission/supplement.qmd.
# Run from the project root, after the figure scripts (27 with PROSE_SCHEME=profile, 34, 44, 52), whose outputs it redraws.
#
# The figures are drawn by the same code as the PNGs in output/figures/: each figure script is run here with ggsave() (and, for the PRISMA diagram,
# svg_to_png()) replaced by a function that writes the submission copies instead, so a submission figure cannot drift from the website figure. The
# scripts' other writes (CSV summaries, the PRISMA HTML) are switched off, so nothing in output/ is touched. The NRI grid is drawn inside the page
# manuscript.qmd and not by a script, so its code is repeated below (keep the two in step).
#
# Writes   submission/figures/Figure_<n>.{svg,pdf,png}, submission/figures/Figure_S<n>.{svg,pdf,png},
#          submission/figures/PROSE_figures_main.pdf (the main figures, one per page), submission/figures/manifest.csv

suppressPackageStartupMessages({library(ggplot2); library(dplyr); library(tidyr)})
source("R/theme_jama.R")
source("R/figure_export.R")
stopifnot(requireNamespace("svglite", quietly = TRUE), requireNamespace("rsvg", quietly = TRUE))

out <- "submission/figures"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
for (f in list.files(out, pattern = "^(Figure_|PROSE_figures)", full.names = TRUE)) file.remove(f)

# PNG written by a figure script -> submission figure (order of first citation in the submission text)
targets <- c(
  "fig_headline_profile.png"        = "Figure_1",
  "fig_icon_array_profile.png"      = "Figure_2",
  "fig_aim3b_gain_profile.png"      = "Figure_3",
  "aim4_sweep_5pct_profile.png"     = "Figure_4",
  "fig_frontier_moves_profile.png"  = "Figure_5",
  "prisma_flow.png"                 = "Figure_S1",
  "fig_forest.png"                  = "Figure_S2",
  "fig_attribution_profile.png"     = "Figure_S3",
  "fig_aim3b_power_profile.png"     = "Figure_S4",
  "aim4_sweep_20pct_profile.png"    = "Figure_S5",
  "nri_grid"                        = "Figure_S6")

manifest <- list()
record <- function(name, source, width, height) {
  manifest[[length(manifest) + 1L]] <<- data.frame(figure = name, source = source, width_in = round(width, 2), height_in = round(height, 2))
}

save_submission <- function(plot, name, width, height, source) {
  base <- file.path(out, name)
  ggplot2::ggsave(paste0(base, ".svg"), plot, width = width, height = height, units = "in", device = svglite::svglite)
  ggplot2::ggsave(paste0(base, ".pdf"), plot, width = width, height = height, units = "in", device = grDevices::cairo_pdf)
  # The PNG is rasterised from the SVG: the R png device can measure text differently from the vector devices (axis titles then overlap the tick labels)
  rsvg::rsvg_png(paste0(base, ".svg"), paste0(base, ".png"), width = round(width * 300))
  png_set_dpi(paste0(base, ".png"), 300)
  record(name, source, width, height)
}

# Run a figure script with its file writes redirected: only the figures in `targets` are written, to submission/figures/
run_figure_script <- function(script) {
  env <- new.env(parent = globalenv())
  env$ggsave <- function(filename, plot = ggplot2::last_plot(), width, height, dpi = 300, ...) {
    key <- basename(filename)
    if (key %in% names(targets)) save_submission(plot, targets[[key]], width, height, paste0(script, " (", key, ")"))
    invisible(filename)
  }
  env$write.csv <- function(...) invisible(NULL)
  env$PRISMA_save <- function(...) invisible(NULL)
  env$svg_to_png <- function(svg, file, width_in, dpi = 300) {
    key <- basename(file)
    if (!key %in% names(targets)) return(invisible(file))
    base <- file.path(out, targets[[key]])
    writeLines(svg, paste0(base, ".svg"))
    rsvg::rsvg_pdf(charToRaw(svg), paste0(base, ".pdf"))
    rsvg::rsvg_png(charToRaw(svg), paste0(base, ".png"), width = round(width_in * dpi))
    png_set_dpi(paste0(base, ".png"), dpi)
    i <- png_info(paste0(base, ".png"))
    record(targets[[key]], paste0(script, " (", key, ")"), width_in, i[["height_px"]] / dpi)
    invisible(file)
  }
  sys.source(script, envir = env)
}

old_scheme <- Sys.getenv("PROSE_SCHEME", NA)
Sys.setenv(PROSE_SCHEME = "profile")
for (s in c("scripts/27_writeup_figures.R", "scripts/34_calibration_check_figures.R", "scripts/44_prisma_flow.R", "scripts/52_manuscript_figures.R")) {
  message("running ", s)
  run_figure_script(s)
}
if (is.na(old_scheme)) Sys.unsetenv("PROSE_SCHEME") else Sys.setenv(PROSE_SCHEME = old_scheme)

# NRI grid: the chunk fig-nri-grid of manuscript.qmd (7 x 5.5 inches, the page's size)
rd <- function(f) read.csv(file.path("output", f))
mcid5 <- rd("dca_baseline_refit_summary.csv") |> filter(threshold == 0.05) |> pull(mcid_b)
panels <- c(nri_cont = "Continuous NRI", nri_cat3 = "NRI, cuts at 5% and 20%", nri_cat5 = "NRI at 5%", dnb_05 = "Net-benefit gain at 5%")
nri_grid <- rd("calibration_check_nri_grid.csv") |>
  filter(comparison == "profile") |>
  pivot_longer(all_of(names(panels)), names_to = "metric") |>
  mutate(metric = factor(panels[metric], panels), or = exp(log_or),
         structure = recode(structure, "PTX3 h2 haplotype" = "h2 haplotype (correct)",
                            "PTX3 three SNPs, independent" = "three SNPs, independent (earlier analysis)")) |>
  ggplot(aes(or, value, colour = structure)) +
  geom_hline(data = data.frame(metric = factor(panels, panels), y = c(0.225, 0.225, 0.225, mcid5)),
             aes(yintercept = y), linetype = "dashed", colour = jama_cols[["ink"]], linewidth = 0.4) +
  geom_vline(xintercept = 1.35, colour = jama_cols[["grey"]], linewidth = 0.4) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~metric, scales = "free_y") +
  scale_x_log10(breaks = c(1, 1.5, 2, 3, 4, 5)) +
  scale_colour_manual(values = c("h2 haplotype (correct)" = jama_cols[["navy"]], "three SNPs, independent (earlier analysis)" = jama_cols[["steel"]])) +
  labs(x = "Odds ratio per copy (log scale)", y = NULL, colour = NULL) +
  theme_jama(base_size = 12)
save_submission(nri_grid, targets[["nri_grid"]], 7, 5.5, "manuscript.qmd, chunk fig-nri-grid (code repeated here)")

manifest <- bind_rows(manifest) |> arrange(grepl("_S", figure), as.integer(sub("^Figure_S?", "", figure)))
missing <- setdiff(unname(targets), manifest$figure)
if (length(missing)) stop("not drawn: ", paste(missing, collapse = ", "), call. = FALSE)
write.csv(manifest, file.path(out, "manifest.csv"), row.names = FALSE)

# The main figures in one PDF, one per page, for reviewers (poppler's pdfunite, if installed)
main <- file.path(out, paste0(manifest$figure[!grepl("_S", manifest$figure)], ".pdf"))
if (nzchar(Sys.which("pdfunite"))) {
  system2("pdfunite", c(shQuote(main), shQuote(file.path(out, "PROSE_figures_main.pdf"))))
} else message("pdfunite not found: PROSE_figures_main.pdf not written")
print(manifest)
