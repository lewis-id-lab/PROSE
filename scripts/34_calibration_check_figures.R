# Figures of the write-up under the profile-preserving truth (addendum 24): the headline gain against the MCID, and the Aim 4 sweep.
# Run from the project root, after scripts 31 and 32. The figures of scripts 23 and 27 stay as the global-intercept versions.
#
# Writes   output/figures/fig_headline_profile.png, aim4_sweep_5pct_profile.png, aim4_sweep_20pct_profile.png

suppressPackageStartupMessages({library(ggplot2); library(dplyr)})
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
    geom_hline(data = mcid, aes(yintercept = mcid), linetype = "dashed") +
    geom_pointrange(aes(ymin = lo, ymax = hi), position = position_dodge(width = 0.5)) +
    facet_wrap(~ threshold, scales = "free_y") +
    labs(x = NULL, y = "Net-benefit gain over the clinical score", colour = NULL,
         caption = "Profile-preserving truth. Dashed line: the minimum clinically important difference. Points: median over draws of the true effects; bars: 90% range.") +
    theme(legend.position = "bottom"), width = 10, height = 4.5, dpi = 150)

sw <- read.csv("output/calibration_check_aim4_sweep.csv") |> filter(comparison == "profile") |> mutate(ratio = gain_oracle / mcid)
plot_sweep <- function(d, th, file) {
  x <- d[abs(d$threshold - th) < 1e-9, ]
  p <- ggplot(x, aes(n_loci, odds_ratio)) +
    geom_raster(aes(fill = ratio)) +
    geom_contour(aes(z = ratio), breaks = 1, colour = "white", linewidth = 0.5) +
    facet_wrap(~ paste0("risk-allele frequency ", raf), nrow = 2) +
    scale_fill_viridis_c(trans = "log10", name = "gain / MCID") +
    labs(x = "Number of loci", y = "Odds ratio per risk allele")
  ggsave(file, p, width = 11, height = 5.5, dpi = 150)
}
plot_sweep(sw, 0.05, "output/figures/aim4_sweep_5pct_profile.png")
plot_sweep(sw, 0.20, "output/figures/aim4_sweep_20pct_profile.png")
