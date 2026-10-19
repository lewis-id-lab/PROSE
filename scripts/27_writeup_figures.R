# Figures for the write-up (addendum 21): the headline gain against the MCID, the Aim 3b gain by effect-modification scenario, the power of the Aim 3b interaction test, and the Shapley
# attribution. Run from the project root, after scripts 16, 21, 22 and 24. The Aim 4 figures are made by scripts 23 and 25.
#
# Writes   output/figures/fig_headline.png, fig_aim3b_gain.png, fig_aim3b_power.png, fig_attribution.png

suppressPackageStartupMessages({library(ggplot2); library(dplyr); library(tidyr)})
dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)
pct <- function(x) paste0(round(100 * x), "%")
thr_lab <- function(x) factor(paste0(round(100 * x), "% threshold"), levels = c("5% threshold", "20% threshold"))
mcid <- read.csv("output/dca_baseline_refit_summary.csv") |> filter(threshold %in% c(0.05, 0.20)) |> transmute(threshold = thr_lab(threshold), mcid = mcid_b)

# 1. Headline: expected gain (median and 90% range of the true effects) against the MCID
l2 <- read.csv("output/clec7a_second_layer_summary.csv") |>
  filter(truth == "mean", scenario %in% c("ceiling", "realised"), threshold %in% c(0.05, 0.20)) |>
  mutate(panel = recode(second_parameter, `held at zero` = "four terms", estimated = "five terms,\nEuropean A-G effect", `estimated from all ancestries` = "five terms,\nall-ancestry A-G effect"),
         panel = factor(panel, levels = c("four terms", "five terms,\nEuropean A-G effect", "five terms,\nall-ancestry A-G effect")),
         score = recode(scenario, ceiling = "true effects in the score", realised = "pooled estimates in the score"), threshold = thr_lab(threshold))
ggsave("output/figures/fig_headline.png",
  ggplot(l2, aes(panel, gain, colour = score)) +
    geom_hline(data = mcid, aes(yintercept = mcid), linetype = "dashed") +
    geom_pointrange(aes(ymin = lo, ymax = hi), position = position_dodge(width = 0.5)) +
    facet_wrap(~ threshold, scales = "free_y") +
    labs(x = NULL, y = "Net-benefit gain over the clinical score", colour = NULL, caption = "Dashed line: the minimum clinically important difference. Points: median over draws of the true effects; bars: 90% range.") +
    theme(legend.position = "bottom"), width = 10, height = 4.5, dpi = 150)

# 2. Aim 3b: whole-cohort gain by scenario
sc_order <- c("H_F only", "H_F, ratio 4", "homogeneous", "Brunel direction, ratio 1/4", "Brunel only", "Brunel estimate")
sc_lab <- c("all in neutropenic", "ratio 4 (Feys)", "homogeneous", "ratio 1/4 (Brunel direction)", "all in non-neutropenic", "Brunel's stratum effects")
a3 <- bind_rows(read.csv("output/aim3b_layer_summary.csv") |> mutate(model = "four terms"), read.csv("output/aim3b_layer_summary_second_clec7a.csv") |> mutate(model = "five terms")) |>
  filter(truth == "mean", level == "whole", threshold %in% c(0.05, 0.20)) |>
  mutate(scenario = factor(sc_lab[match(scenario, sc_order)], levels = rev(sc_lab)), threshold = thr_lab(threshold),
         score = recode(score, ceiling = "true stratum effects (oracle)", common = "pooled effects, same in both strata", assumed_feys = "assumed Feys (ratio 4)", assumed_brunel = "assumed Brunel (ratio 1/4)"),
         model = factor(model, levels = c("four terms", "five terms")))
mc2 <- tidyr::crossing(mcid, model = factor(c("four terms", "five terms"), levels = c("four terms", "five terms")))
ggsave("output/figures/fig_aim3b_gain.png",
  ggplot(a3, aes(gain, scenario, colour = score)) +
    geom_vline(xintercept = 0, colour = "grey70") +
    geom_vline(data = mc2, aes(xintercept = mcid), linetype = "dashed") +
    geom_pointrange(aes(xmin = lo, xmax = hi), position = position_dodge(width = 0.7), size = 0.3) +
    facet_grid(model ~ threshold, scales = "free_x") +
    labs(x = "Whole-cohort net-benefit gain over the clinical score", y = "True effect modification by neutropenia", colour = NULL,
         caption = "Dashed line: the minimum clinically important difference. Points: median; bars: 90% range over draws of the true effects.") +
    theme(legend.position = "bottom") + guides(colour = guide_legend(nrow = 2)), width = 11, height = 8, dpi = 150)

# 3. Aim 3b: power of the interaction test in a study of the cohort's size
pw <- bind_rows(read.csv("output/aim3b_interaction_summary.csv") |> mutate(model = "four terms"), read.csv("output/aim3b_interaction_summary_second_clec7a.csv") |> mutate(model = "five terms")) |>
  mutate(rate = ifelse(true_difference_in_score_units > 0, p_sig_positive, ifelse(true_difference_in_score_units < 0, p_sig_negative, p_sig_positive + p_sig_negative)),
         scenario = factor(sc_lab[match(scenario, sc_order)], levels = sc_lab), model = factor(model, levels = c("four terms", "five terms")),
         effect_level = recode(effect_level, medians = "effects at the pooled medians", upper = "effects at the 95th percentile"))
ggsave("output/figures/fig_aim3b_power.png",
  ggplot(pw, aes(rate, scenario, colour = model)) +
    geom_vline(xintercept = 0.05, linetype = "dotted", colour = "grey50") +
    geom_point(position = position_dodge(width = 0.5), size = 2.2) +
    facet_wrap(~ effect_level) + scale_x_continuous(labels = pct, limits = c(0, 1)) +
    labs(x = "Studies with P < 0.05 in the true direction (false positives for the homogeneous row)", y = NULL, colour = NULL,
         caption = "1,000 simulated studies of 4,078 admissions per point; dotted line: 5%.") + theme(legend.position = "bottom"), width = 10, height = 4.5, dpi = 150)

# 4. Shapley attribution, four terms and five terms
at4 <- read.csv("output/genetic_attribution_summary.csv") |> filter(truth == "mean", threshold == 0.05) |> mutate(model = "four terms")
at5 <- read.csv("output/genetic_attribution_second_clec7a_summary.csv") |> filter(truth == "mean", threshold == 0.05) |> mutate(model = "five terms")
shares <- bind_rows(
  at4 |> transmute(model, score, `CLEC7A (Y238X)` = share_CLEC7A, `CLEC7A (A-G)` = 0, CD209 = share_CD209, PTX3 = share_PTX3, TLR4 = share_TLR4),
  at5 |> transmute(model, score, `CLEC7A (Y238X)` = share_CLEC7A, `CLEC7A (A-G)` = share_CLEC7A_AG, CD209 = share_CD209, PTX3 = share_PTX3, TLR4 = share_TLR4)) |>
  pivot_longer(-c(model, score), names_to = "term", values_to = "share") |>
  mutate(term = factor(term, levels = c("CLEC7A (A-G)", "CLEC7A (Y238X)", "CD209", "PTX3", "TLR4")), model = factor(model, levels = c("four terms", "five terms")),
         score = recode(score, ceiling = "true effects in the score", realised = "pooled estimates in the score"))
ggsave("output/figures/fig_attribution.png",
  ggplot(shares, aes(model, share, fill = term)) + geom_col(width = 0.7) + facet_wrap(~ score) + scale_y_continuous(labels = pct) +
    labs(x = NULL, y = "Share of the Shapley gain (5% threshold)", fill = NULL, caption = "Pooled-mean reading of the true effects.") + theme(legend.position = "bottom"), width = 9, height = 4.5, dpi = 150)
message("wrote four figures to output/figures/")
