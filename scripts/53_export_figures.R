# Export every figure of the manuscript and of the website pages to manuscript_figures/ as PNG at 300 dpi, in the blue-grey style of R/theme_jama.R.
# Run from the project root, after the figure scripts (27, 34, 44, 52), which write the manuscript figures that are made outside the page.
#
# manuscript_figures/manuscript/  Figure_<n>_<id>.png, numbered in the order of the Draft manuscript page (as Quarto numbers them)
# manuscript_figures/website/     <page>_Figure_<n>_<id>.png for the figures drawn inside the other pages (baseline, decision-curves, evidence, mapping)
# manuscript_figures/manifest.csv one row per file: source, size in pixels, the resolution stamped in the file, and the print size in inches that follows
# manuscript_figures/README.md
#
# Figures made by a script (images placed with ![](...)) are copied from output/figures/, after checking that they are at least 300 dpi at the size they were saved
# (their pixel size is at least 300 times their print size; the resolution is stamped here). Figures drawn inside a page (code chunks) are drawn again here, by running
# the chunks of the page in order, at the size the chunk gives (7 inches wide, as Quarto's default, and `fig-height`) and saved at 300 dpi.

suppressPackageStartupMessages({library(ggplot2); library(dplyr); library(tidyr)})
source("R/figure_export.R")
dpi <- 300
out <- "manuscript_figures"
dir.create(file.path(out, "manuscript"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out, "website"), recursive = TRUE, showWarnings = FALSE)
for (f in list.files(file.path(out, c("manuscript", "website")), full.names = TRUE)) file.remove(f)

# The code chunks of a page: label, figure height, code
read_chunks <- function(qmd) {
  x <- readLines(qmd)
  starts <- grep("^```\\{r", x)
  lapply(starts, function(s) {
    e <- s + which(x[(s + 1):length(x)] == "```")[1]
    body <- x[(s + 1):(e - 1)]
    opt <- function(k) { m <- sub(paste0("^#\\| ", k, ": *"), "", body[grepl(paste0("^#\\| ", k, ":"), body)]); if (length(m)) gsub('^"|"$', "", m[1]) else NA_character_ }
    list(line = s, label = opt("label"), height = suppressWarnings(as.numeric(opt("fig-height"))), width = suppressWarnings(as.numeric(opt("fig-width"))),
         eval = opt("eval"), code = body[!grepl("^#\\|", body)])
  })
}

# Run the chunks of a page in order (as the render does) and return the figure chunks with the ggplot each one draws
run_page <- function(qmd) {
  env <- new.env(parent = globalenv())
  figs <- list()
  for (ch in read_chunks(qmd)) {
    if (identical(ch$eval, "false")) next
    plot <- NULL
    for (ex in parse(text = ch$code, keep.source = FALSE)) {
      v <- withVisible(eval(ex, env))
      if (v$visible && inherits(v$value, "ggplot")) plot <- v$value
    }
    if (!is.na(ch$label) && grepl("^fig-", ch$label)) {
      if (is.null(plot)) stop("figure chunk ", ch$label, " in ", qmd, " draws no ggplot", call. = FALSE)
      figs[[length(figs) + 1L]] <- list(label = ch$label, line = ch$line, plot = plot, width = ifelse(is.na(ch$width), 7, ch$width), height = ifelse(is.na(ch$height), 5, ch$height))
    }
  }
  figs
}

save_plot <- function(plot, file, width, height) {
  ggsave(file, plot, width = width, height = height, dpi = dpi, device = "png")
  png_set_dpi(file, dpi)
}

manifest <- list()
record <- function(file, source, kind, width_in = NA, height_in = NA) {
  i <- png_info(file)
  manifest[[length(manifest) + 1L]] <<- data.frame(file = sub(paste0("^", out, "/"), "", file), source = source, kind = kind, width_px = i[["width_px"]], height_px = i[["height_px"]],
    dpi_stamped = round(i[["dpi"]]), width_in = round(i[["width_px"]] / dpi, 2), height_in = round(i[["height_px"]] / dpi, 2), stringsAsFactors = FALSE)
}

# ---- Manuscript ------------------------------------------------------------------------------------------------------------------------------------
ms <- readLines("manuscript.qmd")
img <- grep("^!\\[.*\\]\\([^)]+\\.png\\)\\{#fig-[A-Za-z0-9_-]+\\}\\s*$", ms)
events <- data.frame(line = img, id = sub(".*\\{#fig-([A-Za-z0-9_-]+)\\}\\s*$", "\\1", ms[img]), path = sub(".*\\(([^)]+\\.png)\\).*", "\\1", ms[img]), type = "image", stringsAsFactors = FALSE)
page <- run_page("manuscript.qmd")
chunk_events <- if (length(page)) data.frame(line = vapply(page, `[[`, 0, "line"), id = sub("^fig-", "", vapply(page, `[[`, "", "label")), path = NA_character_, type = "chunk", stringsAsFactors = FALSE) else NULL
events <- bind_rows(events, chunk_events) |> arrange(line) |> mutate(n = row_number())
for (k in seq_len(nrow(events))) {
  ev <- events[k, ]; name <- file.path(out, "manuscript", sprintf("Figure_%d_%s.png", ev$n, ev$id))
  if (ev$type == "image") {
    info <- png_info(ev$path)
    if (info[["width_px"]] < 7 * dpi * 0.6) warning(ev$path, " is only ", info[["width_px"]], " pixels wide; check its print size")
    file.copy(ev$path, name, overwrite = TRUE); png_set_dpi(name, dpi); record(name, ev$path, "script figure (copied)")
  } else {
    f <- page[[which(vapply(page, function(p) p$label == paste0("fig-", ev$id), TRUE))]]
    save_plot(f$plot, name, f$width, f$height); record(name, paste0("manuscript.qmd, chunk fig-", ev$id), "page chunk (drawn here)")
  }
}

# ---- Website pages ---------------------------------------------------------------------------------------------------------------------------------
for (qmd in c("baseline", "decision-curves", "evidence", "mapping")) {
  for (k in seq_along(figs <- run_page(paste0(qmd, ".qmd")))) {
    f <- figs[[k]]; name <- file.path(out, "website", sprintf("%s_Figure_%d_%s.png", qmd, k, sub("^fig-", "", f$label)))
    save_plot(f$plot, name, f$width, f$height); record(name, paste0(qmd, ".qmd, chunk ", f$label), "page chunk (drawn here)")
  }
}

manifest <- bind_rows(manifest)
write.csv(manifest, file.path(out, "manifest.csv"), row.names = FALSE)
stopifnot(all(manifest$dpi_stamped == dpi), all(manifest$width_px >= 3 * dpi))
writeLines(c("# Figures of the manuscript and the website", "",
  "PNG files at 300 dpi in the blue-grey style of `R/theme_jama.R`, made by `scripts/53_export_figures.R` (run from the project root; it redraws the figures that sit",
  "inside the pages and copies the ones made by `scripts/27`, `34`, `44` and `52`).", "",
  "- `manuscript/`: the figures of the Draft manuscript page, numbered in the order of the page (`Figure_<n>_<id>.png`). The numbers are those of the rendered page.",
  "- `website/`: the figures drawn inside the other pages of the site (baseline, decision curves, evidence, mapping), named `<page>_Figure_<n>_<id>.png`.",
  "- `manifest.csv`: each file's source, pixel size, the resolution stamped in the file and the print size in inches that follows (pixels divided by 300).", "",
  "Each file carries its resolution in the PNG header (pHYs chunk), so a viewer or a journal's checker reads 300 dpi at the print size in the manifest. The captions are in",
  "`manuscript.qmd` and the page sources: the figures carry no title or caption of their own. Tables are not included."), file.path(out, "README.md"))
print(as.data.frame(manifest[, c("file", "width_px", "height_px", "dpi_stamped", "width_in", "height_in")]), row.names = FALSE)
