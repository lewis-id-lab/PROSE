# Figures for the write-up (addendum 21): the headline gain against the MCID, the Aim 3b gain by effect-modification scenario, the power of the Aim 3b interaction test, and the Shapley
# attribution. Run from the project root, after scripts 16, 21, 22 and 24. The Aim 4 figures are made by scripts 23 and 25.
#
# Writes   output/figures/fig_headline.png, fig_aim3b_gain.png, fig_aim3b_power.png, fig_attribution.png

# Run with Sys.setenv(PROSE_SCHEME = "profile") to draw the Aim 3b and attribution figures from the profile-preserving results (addendum 25);
# they are then read from and written to files ending in _profile. The headline figure under that truth is made by script 34.
sfx <- if (Sys.getenv("PROSE_SCHEME", "global") == "profile") "_profile" else ""
suppressPackageStartupMessages({library(ggplot2); library(dplyr); library(tidyr)})
source("R/theme_jama.R")
dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)
pct <- function(x) paste0(round(100 * x), "%")
thr_lab <- function(x) factor(paste0(round(100 * x), "% threshold"), levels = c("5% threshold", "20% threshold"))
mcid <- read.csv("output/dca_baseline_refit_summary.csv") |> filter(threshold %in% c(0.05, 0.20)) |> transmute(threshold = thr_lab(threshold), mcid = mcid_b)

# 1. Headline: expected gain (median and 90% range of the true effects) against the MCID
if (sfx == "") {
l2 <- read.csv("output/clec7a_second_layer_summary.csv") |>
  filter(truth == "mean", scenario %in% c("ceiling", "realised"), threshold %in% c(0.05, 0.20)) |>
  mutate(panel = recode(second_parameter, `held at zero` = "four terms", estimated = "five terms,\nEuropean A-G effect", `estimated from all ancestries` = "five terms,\nall-ancestry A-G effect"),
         panel = factor(panel, levels = c("four terms", "five terms,\nEuropean A-G effect", "five terms,\nall-ancestry A-G effect")),
         score = recode(scenario, ceiling = "true effects in the score", realised = "pooled estimates in the score"), threshold = thr_lab(threshold))
ggsave("output/figures/fig_headline.png",
  ggplot(l2, aes(panel, gain, colour = score)) +
    geom_hline(yintercept = 0, colour = jama_cols[["grey"]], linewidth = 0.3) +
    geom_hline(data = mcid, aes(yintercept = mcid), linetype = "dashed", colour = jama_cols[["ink"]], linewidth = 0.4) +
    geom_text(data = mcid, aes(x = Inf, y = mcid, label = "MCID"), inherit.aes = FALSE, hjust = 1.1, vjust = -0.5, size = 3.2, colour = jama_cols[["ink"]]) +
    geom_pointrange(aes(ymin = lo, ymax = hi), position = position_dodge(width = 0.45), linewidth = 0.7, size = 0.55) +
    facet_wrap(~ threshold, scales = "free_y") +
    scale_colour_manual(values = c("pooled estimates in the score" = jama_cols[["steel"]], "true effects in the score" = jama_cols[["navy"]])) +
    labs(x = NULL, y = "Net-benefit gain over the clinical score", colour = NULL) +
    theme_jama(), width = 10, height = 4.5, dpi = 300)
}

# 2. Aim 3b: whole-cohort gain by scenario
sc_order <- c("H_F only", "H_F, ratio 4", "homogeneous", "Brunel direction, ratio 1/4", "Brunel only", "Brunel estimate")
sc_lab <- c("all in neutropenic", "ratio 4 (Feys)", "homogeneous", "ratio 1/4 (Brunel direction)", "all in non-neutropenic", "Brunel's stratum effects")
a3 <- bind_rows(read.csv(paste0("output/aim3b_layer_summary", sfx, ".csv")) |> mutate(model = "four terms"), read.csv(paste0("output/aim3b_layer_summary_second_clec7a", sfx, ".csv")) |> mutate(model = "five terms")) |>
  filter(truth == "mean", level == "whole", threshold %in% c(0.05, 0.20)) |>
  mutate(scenario = factor(sc_lab[match(scenario, sc_order)], levels = rev(sc_lab)), threshold = thr_lab(threshold),
         score = recode(score, ceiling = "true stratum effects (oracle)", common = "pooled effects, same in both strata", assumed_feys = "assumed Feys (ratio 4)", assumed_brunel = "assumed Brunel (ratio 1/4)"),
         score = factor(score, levels = c("true stratum effects (oracle)", "pooled effects, same in both strata", "assumed Feys (ratio 4)", "assumed Brunel (ratio 1/4)")),
         model = factor(model, levels = c("four terms", "five terms")))
mc2 <- tidyr::crossing(mcid, model = factor(c("four terms", "five terms"), levels = c("four terms", "five terms")))
ggsave(paste0("output/figures/fig_aim3b_gain", sfx, ".png"),
  ggplot(a3, aes(gain, scenario, colour = score)) +
    geom_vline(xintercept = 0, colour = jama_cols[["grey"]], linewidth = 0.3) +
    geom_vline(data = mc2, aes(xintercept = mcid), linetype = "dashed", colour = jama_cols[["ink"]], linewidth = 0.4) +
    geom_pointrange(aes(xmin = lo, xmax = hi), position = position_dodge2(width = 0.7, reverse = TRUE), size = 0.35, linewidth = 0.6) +
    facet_grid(model ~ threshold, scales = "free_x") +
    scale_colour_manual(values = c("true stratum effects (oracle)" = jama_cols[["navy"]], "pooled effects, same in both strata" = jama_cols[["blue"]],
                                   "assumed Feys (ratio 4)" = jama_cols[["steel"]], "assumed Brunel (ratio 1/4)" = jama_cols[["grey"]])) +
    labs(x = "Whole-cohort net-benefit gain over the clinical score", y = "True effect modification by neutropenia", colour = NULL) +
    theme_jama(grid = "x") + theme(panel.spacing = unit(1.1, "lines")) + guides(colour = guide_legend(nrow = 1)), width = 11, height = 7.5, dpi = 300)

# 3. Aim 3b: power of the interaction test in a study of the cohort's size
pw <- bind_rows(read.csv(paste0("output/aim3b_interaction_summary", sfx, ".csv")) |> mutate(model = "four terms"), read.csv(paste0("output/aim3b_interaction_summary_second_clec7a", sfx, ".csv")) |> mutate(model = "five terms")) |>
  mutate(rate = ifelse(true_difference_in_score_units > 0, p_sig_positive, ifelse(true_difference_in_score_units < 0, p_sig_negative, p_sig_positive + p_sig_negative)),
         scenario = factor(sc_lab[match(scenario, sc_order)], levels = sc_lab), model = factor(model, levels = c("four terms", "five terms")),
         effect_level = recode(effect_level, medians = "effects at the pooled medians", upper = "effects at the 95th percentile"))
ggsave(paste0("output/figures/fig_aim3b_power", sfx, ".png"),
  ggplot(pw, aes(rate, scenario, colour = model)) +
    geom_vline(xintercept = 0.05, linetype = "dotted", colour = jama_cols[["grey"]]) +
    geom_point(position = position_dodge(width = 0.5), size = 2.4) +
    facet_wrap(~ effect_level) + scale_x_continuous(labels = pct, limits = c(0, 1)) +
    scale_colour_manual(values = c("four terms" = jama_cols[["steel"]], "five terms" = jama_cols[["navy"]])) +
    labs(x = "Studies with P < 0.05 in the true direction (false positives for the homogeneous row)", y = NULL, colour = NULL,
         caption = "1,000 simulated studies of 4,078 admissions per point; dotted line: 5%.") + theme_jama(grid = "x"), width = 10, height = 4.5, dpi = 300)

# 4. Shapley attribution, four terms and five terms
at4 <- read.csv(paste0("output/genetic_attribution_summary", sfx, ".csv")) |> filter(truth == "mean", threshold == 0.05) |> mutate(model = "four terms")
at5 <- read.csv(paste0("output/genetic_attribution_second_clec7a_summary", sfx, ".csv")) |> filter(truth == "mean", threshold == 0.05) |> mutate(model = "five terms")
shares <- bind_rows(
  at4 |> transmute(model, score, `CLEC7A (Y238X)` = share_CLEC7A, `CLEC7A (A-G)` = 0, CD209 = share_CD209, PTX3 = share_PTX3, TLR4 = share_TLR4),
  at5 |> transmute(model, score, `CLEC7A (Y238X)` = share_CLEC7A, `CLEC7A (A-G)` = share_CLEC7A_AG, CD209 = share_CD209, PTX3 = share_PTX3, TLR4 = share_TLR4)) |>
  pivot_longer(-c(model, score), names_to = "term", values_to = "share") |>
  mutate(term = factor(term, levels = c("CLEC7A (A-G)", "CLEC7A (Y238X)", "CD209", "PTX3", "TLR4")), model = factor(model, levels = c("four terms", "five terms")),
         score = recode(score, ceiling = "true effects in the score", realised = "pooled estimates in the score"))
ggsave(paste0("output/figures/fig_attribution", sfx, ".png"),
  ggplot(shares, aes(share, model, fill = term)) + geom_col(width = 0.62, colour = "white", linewidth = 0.4) + facet_wrap(~ score, ncol = 1) +
    geom_text(aes(label = ifelse(share >= 0.06, pct(share), ""), colour = term), position = position_stack(vjust = 0.5), size = 3.1, show.legend = FALSE) +
    scale_x_continuous(labels = pct, expand = expansion(mult = c(0, 0.02))) + scale_y_discrete(limits = rev) +
    scale_fill_manual(values = c("CLEC7A (A-G)" = jama_cols[["navy"]], "CLEC7A (Y238X)" = jama_cols[["blue"]], "CD209" = jama_cols[["steel"]], "PTX3" = jama_cols[["pale"]], "TLR4" = jama_cols[["light"]])) +
    scale_colour_manual(values = c("CLEC7A (A-G)" = "white", "CLEC7A (Y238X)" = "white", "CD209" = "white", "PTX3" = jama_cols[["ink"]], "TLR4" = jama_cols[["ink"]])) +
    labs(y = NULL, x = "Share of the net-benefit gain at the 5% threshold (Shapley)", fill = NULL) + guides(fill = guide_legend(reverse = TRUE)) +
    theme_jama(grid = "none"), width = 9, height = 4.2, dpi = 300)
message("wrote four figures to output/figures/")
