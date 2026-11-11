# A blue-grey figure style in the manner of JAMA: sans-serif, thin axis lines, no panel background, no title inside the figure (it belongs in the caption).
# Sourced by _setup.R (the Quarto pages) and by the figure scripts. Requires ggplot2.

jama_cols <- c(navy = "#1B365D", blue = "#3C6E9E", steel = "#7FA3C7", pale = "#C9D8E6", grey = "#6B7785", light = "#E9EDF1", ink = "#2B2B2B")

# n colours for unordered groups (most contrast first) or for ordered groups (dark to pale blue)
jama_qual <- function(n) unname(jama_cols[c("navy", "steel", "grey", "blue", "pale")][seq_len(n)])
jama_seq <- function(n) grDevices::colorRampPalette(jama_cols[c("navy", "blue", "steel", "pale")])(n)

scale_colour_jama <- function(..., ordered = FALSE) ggplot2::discrete_scale("colour", palette = if (ordered) jama_seq else jama_qual, ...)
scale_fill_jama <- function(..., ordered = FALSE) ggplot2::discrete_scale("fill", palette = if (ordered) jama_seq else jama_qual, ...)

# grid: light reference lines across the value axis ("y" for vertical plots, "x" for horizontal ones, "none")
theme_jama <- function(base_size = 11, base_family = "", grid = c("y", "x", "none")) {
  grid <- match.arg(grid)
  line <- ggplot2::element_line(colour = jama_cols[["light"]], linewidth = 0.3)
  ggplot2::theme_classic(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      text = ggplot2::element_text(colour = jama_cols[["ink"]]),
      axis.text = ggplot2::element_text(colour = jama_cols[["ink"]]),
      axis.line = ggplot2::element_line(colour = jama_cols[["ink"]], linewidth = 0.35),
      axis.ticks = ggplot2::element_line(colour = jama_cols[["ink"]], linewidth = 0.35),
      panel.grid.major.y = if (grid == "y") line else ggplot2::element_blank(),
      panel.grid.major.x = if (grid == "x") line else ggplot2::element_blank(),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold", hjust = 0, colour = jama_cols[["navy"]], size = ggplot2::rel(1)),
      legend.position = "top", legend.justification = "left", legend.title = ggplot2::element_text(size = ggplot2::rel(0.9)),
      legend.key.size = ggplot2::unit(0.9, "lines"), legend.margin = ggplot2::margin(0, 0, 0, 0),
      plot.caption = ggplot2::element_text(colour = jama_cols[["grey"]], hjust = 0),
      plot.margin = ggplot2::margin(8, 12, 8, 8))
}

# Table style for gt tables in the same scheme: heavy navy rules above and below the table, a thin rule under the column labels on a pale blue-grey band,
# light horizontal rules between rows, no vertical lines, sans-serif text. Requires gt.
gt_jama_style <- function(data, font_size = 13) {
  gt::tab_options(data,
    table.font.names = c("Helvetica", "Arial", "sans-serif"), table.font.size = gt::px(font_size), table.font.color = jama_cols[["ink"]],
    table.background.color = "white", table.width = NULL,
    table.border.top.style = "solid", table.border.top.color = jama_cols[["navy"]], table.border.top.width = gt::px(2),
    table.border.bottom.style = "solid", table.border.bottom.color = jama_cols[["navy"]], table.border.bottom.width = gt::px(2),
    heading.align = "left", heading.border.bottom.style = "none",
    column_labels.font.weight = "bold", column_labels.background.color = jama_cols[["light"]],
    column_labels.border.top.style = "none", column_labels.border.bottom.style = "solid",
    column_labels.border.bottom.color = jama_cols[["navy"]], column_labels.border.bottom.width = gt::px(1),
    column_labels.padding = gt::px(6),
    row_group.font.weight = "bold", row_group.background.color = "white", row_group.border.top.color = jama_cols[["pale"]], row_group.border.bottom.color = jama_cols[["pale"]],
    table_body.hlines.style = "solid", table_body.hlines.color = jama_cols[["light"]], table_body.hlines.width = gt::px(1),
    table_body.vlines.style = "none", table_body.border.bottom.color = jama_cols[["navy"]], table_body.border.bottom.width = gt::px(1),
    data_row.padding = gt::px(4), source_notes.font.size = gt::px(font_size - 2))
}
