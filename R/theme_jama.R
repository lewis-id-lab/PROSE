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
