# Figures of the write-up under the profile-preserving truth (addendum 24): the headline gain against the MCID, and the Aim 4 sweep.
# Run from the project root, after scripts 31 and 32. The figures of scripts 23 and 27 stay as the global-intercept versions.
#
# Writes   output/figures/fig_headline_profile.png, aim4_sweep_5pct_profile.png, aim4_sweep_20pct_profile.png

suppressPackageStartupMessages({library(ggplot2); library(dplyr)})
source("R/theme_jama.R")
dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)
thr_lab <- function(x) factor(paste0(round(100 * x), "% threshold"), levels = c("5% threshold", "20% threshold"))

mcid <- read.csv("output/dca_baseline_refit_summary.csv") |> filter(threshold %in% c(0.05, 0.20)) |> transmute(threshold = thr_lab(threshold), mcid = mcid_b)
# The all-ancestry estimate of the second CLEC7A parameter under the profile truth comes from script 22 (PROSE_SCHEME = "profile")
all_anc <- read.csv("output/clec7a_second_layer_summary_profile.csv") |>
  filter(second_parameter == "estimated from all ancestries", truth == "mean", threshold %in% c(0.05, 0.20)) |>
  transmute(panel = "five terms, all ancestries", comparison = "profile", truth, threshold, score = scenario, gain, lo, hi)
lay <- read.csv("output/calibration_check_layer.csv") |>
  filter(comparison == "profile", truth == "mean", threshold %in% c(0.05, 0.20)) |>
  select(panel, comparison, truth, threshold, score, gain, lo, hi) |>
  bind_rows(all_anc) |>
  mutate(panel = factor(recode(panel, "four terms" = "four terms", "five terms" = "five terms,\nEuropean A-G effect", "five terms, all ancestries" = "five terms,\nall-ancestry A-G effect"),
                        levels = c("four terms", "five terms,\nEuropean A-G effect", "five terms,\nall-ancestry A-G effect")),
         score = recode(score, ceiling = "true effects in the score", realised = "pooled estimates in the score"), threshold = thr_lab(threshold))
ggsave("output/figures/fig_headline_profile.png",
  ggplot(lay, aes(panel, gain, colour = score)) +
    geom_hline(yintercept = 0, colour = jama_cols[["grey"]], linewidth = 0.3) +
    geom_hline(data = mcid, aes(yintercept = mcid), linetype = "dashed", colour = jama_cols[["ink"]], linewidth = 0.4) +
    geom_text(data = mcid, aes(x = Inf, y = mcid, label = "MCID"), inherit.aes = FALSE, hjust = 1.1, vjust = -0.5, size = 3.2, colour = jama_cols[["ink"]]) +
    geom_pointrange(aes(ymin = lo, ymax = hi), position = position_dodge(width = 0.45), linewidth = 0.7, size = 0.55) +
    facet_wrap(~ threshold, scales = "free_y") +
    scale_colour_manual(values = c("pooled estimates in the score" = jama_cols[["steel"]], "true effects in the score" = jama_cols[["navy"]])) +
    labs(x = NULL, y = "Net-benefit gain over the clinical score", colour = NULL) +
    theme_jama(), width = 10, height = 4.5, dpi = 300)

sw <- read.csv("output/calibration_check_aim4_sweep.csv") |> filter(comparison == "profile") |> mutate(ratio = gain_oracle / mcid)
# Bands of the gain relative to the MCID; the darkest band reaches it. `mark`: an architecture to mark (the four-term panel's equivalent, 5% only).
plot_sweep <- function(d, th, file, mark = NULL) {
  x <- d[abs(d$threshold - th) < 1e-9, ]
  x$band <- cut(x$ratio, c(-Inf, 0.25, 0.5, 1, Inf), labels = c("under a quarter of the MCID", "a quarter to half", "half to the MCID", "reaches the MCID"), right = FALSE)
  x$panel <- paste0("Risk-allele frequency ", formatC(x$raf, format = "f", digits = 2))
  p <- ggplot(x, aes(n_loci, odds_ratio)) +
    geom_tile(aes(fill = band), colour = "white", linewidth = 0.15) +
    facet_wrap(~ panel, nrow = 2) +
    scale_fill_manual(values = unname(jama_cols[c("light", "pale", "steel", "navy")]), name = NULL, drop = FALSE) +
    scale_x_continuous(breaks = c(1, 5, 10, 15, 20), expand = c(0, 0)) + scale_y_continuous(breaks = seq(1.5, 3, 0.5), expand = c(0, 0)) +
    labs(x = "Number of independent risk loci", y = "Odds ratio per risk allele") + theme_jama(grid = "none")
  if (!is.null(mark)) { mark$panel <- paste0("Risk-allele frequency ", formatC(mark$raf, format = "f", digits = 2))
    p <- p + geom_point(data = mark, shape = 21, size = 2.4, stroke = 0.9, fill = "white", colour = jama_cols[["ink"]]) +
      geom_text(data = mark, aes(label = label), hjust = 0, vjust = 0.5, nudge_x = 1, size = 2.9, lineheight = 0.9, colour = jama_cols[["ink"]]) }
  ggsave(file, p, width = 11, height = 5.5, dpi = 300)
}
# The four-term panel's ceiling gain at 5% (0.00056) equals that of 7 loci of odds ratio 1.3 at frequency 0.3
plot_sweep(sw, 0.05, "output/figures/aim4_sweep_5pct_profile.png", mark = data.frame(n_loci = 7, odds_ratio = 1.3, raf = 0.3, label = "four-term panel"))
plot_sweep(sw, 0.20, "output/figures/aim4_sweep_20pct_profile.png")
