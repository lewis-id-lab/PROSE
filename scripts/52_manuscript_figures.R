# Two further figures for the write-up, in the blue-grey style of R/theme_jama.R. Run from the project root, after scripts 12, 14, 15, 11 and 31.
#
# 1. Forest plot of the study effects and the pooled effect for the four terms of the panel (primary analysis: recipient genotype, odds-ratio scale, whole cohort),
#    per copy of each term's risk haplotype class. Study intervals are 90% (estimate +- 1.645 SE on the log scale) to match the 90% credible interval of the pooled effect.
# 2. Icon array of the primary result at the 5% threshold, per 10,000 admissions so that every icon is a whole case: cases the clinical score finds, cases the
#    four-term genetic score adds (true effects in the score, profile-preserving truth), further cases needed to reach the MCID, and cases still missed.
#
# 3. "How the bar moves": the banded architecture figure at frequency 0.3 (5% threshold, profile-preserving truth) with the frontiers of three departures drawn on it:
#    a scored variant that tags the causal one at r = 0.9 and at r = 0.7 (script 51; estimated weights, SE 0.10), and only 80% of cases due to the target pathogen
#    (script 43; exact all-mould weights). Run after scripts 31, 43 and 51.
#
# Writes   output/figures/fig_forest.png, output/figures/fig_icon_array_profile.png, output/fig_icon_array_counts.csv, output/figures/fig_frontier_moves_profile.png

suppressPackageStartupMessages({library(ggplot2); library(dplyr)})
source("R/theme_jama.R")
dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)

# ---- 1. Forest plot --------------------------------------------------------------------------------------------------------------------------------
terms <- list(
  list(study = "ptx3_study_effects.csv", pooled = "ptx3_h2_pooled.csv", label = "PTX3\nh2 haplotype"),
  list(study = "clec7a_study_effects.csv", pooled = "clec7a_pooled.csv", label = "CLEC7A\nY238X haplotype"),
  list(study = "tlr4_study_effects.csv", pooled = "tlr4_pooled.csv", label = "TLR4\nG-T haplotype"),
  list(study = "cd209_study_effects.csv", pooled = "cd209_pooled.csv", label = "CD209\nrs4804800 G"))
z90 <- stats::qnorm(0.95)
forest <- bind_rows(lapply(terms, function(t) {
  s <- read.csv(file.path("output", t$study)); p <- read.csv(file.path("output", t$pooled))[1, ]
  s <- s[order(s$year, s$first_author), ]
  bind_rows(data.frame(term = t$label, row = paste(sub(",.*$", "", s$first_author), s$year), kind = "study", or = exp(s$beta), lo = exp(s$beta - z90 * s$se_beta),
                       hi = exp(s$beta + z90 * s$se_beta), precision = 1 / s$se_beta),
            data.frame(term = t$label, row = paste0("Pooled (", p$k, " studies)"), kind = "pooled", or = p$or, lo = p$or_lo, hi = p$or_hi, precision = NA_real_))
}))
xlim <- c(0.1, 10)
forest <- forest |>
  mutate(term = factor(term, levels = vapply(terms, function(t) t$label, "")), key = factor(paste(term, row), levels = rev(unique(paste(term, row)))),
         lo_plot = pmax(lo, xlim[1]), hi_plot = pmin(hi, xlim[2]), text = sprintf("%.2f (%.2f to %.2f)", or, lo, hi))
labels <- setNames(forest$row, forest$key)
header <- data.frame(term = factor(levels(forest$term)[1], levels = levels(forest$term)), label = "Odds ratio (90% interval)")
ggsave("output/figures/fig_forest.png",
  ggplot(forest, aes(or, key)) +
    geom_vline(xintercept = 1, colour = jama_cols[["grey"]], linewidth = 0.35) +
    geom_segment(aes(x = lo_plot, xend = hi_plot, yend = key, colour = kind, linewidth = kind)) +
    geom_segment(data = filter(forest, lo < xlim[1]), aes(x = lo_plot * 1.25, xend = lo_plot, yend = key, colour = kind), linewidth = 0.5, arrow = arrow(length = unit(0.07, "inches"), type = "closed")) +
    geom_point(data = filter(forest, kind == "study"), aes(size = precision), shape = 15, colour = jama_cols[["blue"]]) +
    geom_point(data = filter(forest, kind == "pooled"), shape = 23, size = 3.4, fill = jama_cols[["navy"]], colour = jama_cols[["navy"]]) +
    geom_text(aes(x = xlim[2] * 1.25, label = text, fontface = ifelse(kind == "pooled", "bold", "plain")), hjust = 0, size = 3.1, colour = jama_cols[["ink"]]) +
    geom_text(data = header, aes(x = xlim[2] * 1.25, y = Inf, label = label), hjust = 0, vjust = -0.9, size = 3.1, fontface = "bold", colour = jama_cols[["navy"]]) +
    facet_grid(term ~ ., scales = "free_y", space = "free_y", switch = "y") +
    scale_x_log10(breaks = c(0.1, 0.25, 0.5, 1, 2, 4, 10), labels = c("0.1", "0.25", "0.5", "1", "2", "4", "10")) +
    scale_y_discrete(labels = labels) +
    scale_colour_manual(values = c(study = jama_cols[["steel"]], pooled = jama_cols[["navy"]]), guide = "none") +
    scale_linewidth_manual(values = c(study = 0.6, pooled = 1), guide = "none") + scale_size_continuous(range = c(1.6, 3.6), guide = "none") +
    coord_cartesian(xlim = xlim, clip = "off") +
    labs(x = "Odds ratio per copy of the risk haplotype (log scale)", y = NULL) +
    theme_jama(grid = "none") +
    theme(strip.placement = "outside", strip.text.y.left = element_text(angle = 0, hjust = 0, vjust = 1, face = "bold", colour = jama_cols[["navy"]]),
          axis.line.y = element_blank(), axis.ticks.y = element_blank(), panel.spacing = unit(0.9, "lines"), plot.margin = margin(22, 150, 8, 8)),
  width = 8.5, height = 5.6, dpi = 300)

# ---- 2. Icon array ---------------------------------------------------------------------------------------------------------------------------------
per <- 10000
base <- read.csv("output/dca_baseline_refit_summary.csv") |> filter(threshold == 0.05)
lay <- read.csv("output/calibration_check_layer.csv") |> filter(comparison == "profile", truth == "mean", threshold == 0.05, panel == "four terms", score == "ceiling")
# The event rate follows from the baseline's net benefit of treating everyone: nb_all = prevalence - (1 - prevalence) t / (1 - t)
prev <- (base$nb_all + 0.05 / 0.95) / (1 + 0.05 / 0.95)
n <- c(cases = prev * per, found = base$tp_per_1000 * per / 1000, added = lay$extra_cases_per_1000 * per / 1000, mcid = base$mcid_b * per,
       extra_flagged = (lay$flagged - lay$flagged_clinical) * per)
r <- round(n)
counts <- c("Found by the clinical score" = r[["found"]], "Added by the genetic score" = r[["added"]],
            "Further cases needed to reach the MCID" = r[["mcid"]] - r[["added"]], "Still missed" = r[["cases"]] - r[["found"]] - r[["mcid"]])
stopifnot(all(counts > 0), sum(counts) == r[["cases"]])
write.csv(data.frame(quantity = c(names(n), "extra_false_positives"), per_10000 = c(unname(n), n[["extra_flagged"]] - n[["added"]])), "output/fig_icon_array_counts.csv", row.names = FALSE)
ncol <- 33
icons <- data.frame(i = seq_len(sum(counts)) - 1, group = factor(rep(names(counts), counts), levels = names(counts)))
icons$x <- icons$i %% ncol; icons$y <- -(icons$i %/% ncol)
lab <- setNames(paste0(names(counts), " (", counts, ")"), names(counts))
ggsave("output/figures/fig_icon_array_profile.png",
  ggplot(icons, aes(x, y, fill = group, colour = group)) +
    geom_point(shape = 21, size = 3.6, stroke = 0.7) +
    scale_fill_manual(values = setNames(c(jama_cols[["navy"]], jama_cols[["steel"]], "white", jama_cols[["light"]]), names(counts)), labels = lab, name = NULL) +
    scale_colour_manual(values = setNames(c(jama_cols[["navy"]], jama_cols[["steel"]], jama_cols[["blue"]], "#B8C0C8"), names(counts)), labels = lab, name = NULL) +
    coord_equal(clip = "off") + guides(fill = guide_legend(nrow = 2, byrow = TRUE), colour = guide_legend(nrow = 2, byrow = TRUE)) +
    theme_void(base_size = 11) + theme(legend.position = "top", legend.justification = "left", legend.text = element_text(colour = jama_cols[["ink"]]), plot.margin = margin(8, 8, 8, 8)),
  width = 8.5, height = 3.6, dpi = 300)
message("icon array, per ", per, " admissions: ", paste(names(counts), counts, sep = " = ", collapse = "; "), "; extra false positives ", round(n[["extra_flagged"]] - n[["added"]]))

# ---- 3. How the bar moves --------------------------------------------------------------------------------------------------------------------------
raf0 <- 0.3
sw <- read.csv("output/calibration_check_aim4_sweep.csv") |> filter(comparison == "profile", threshold == 0.05, abs(raf - raf0) < 1e-9) |>
  mutate(band = cut(gain_oracle / mcid, c(-Inf, 0.25, 0.5, 1, Inf), labels = c("under a quarter of the MCID", "a quarter to half", "half to the MCID", "reaches the MCID"), right = FALSE))
st <- read.csv("output/aim4_main_stress_frontier_profile.csv") |> filter(threshold == 0.05, abs(raf - raf0) < 1e-9, sigma == 0.10, scenario %in% c("tag0.9", "tag0.7")) |>
  transmute(line = recode(scenario, tag0.9 = "tag of the causal variant, r = 0.9", tag0.7 = "tag of the causal variant, r = 0.7"), odds_ratio, n_loci = median_reaches)
pa <- read.csv("output/aim4_pathogen_frontier_profile.csv") |> filter(threshold == 0.05, abs(raf - raf0) < 1e-9, case_fraction == 0.8) |>
  transmute(line = "80% of cases due to the pathogen", odds_ratio, n_loci = min_loci_diluted)
moves <- bind_rows(st, pa) |> filter(!is.na(n_loci)) |> arrange(line, odds_ratio)
ends <- moves |> group_by(line) |> slice_min(odds_ratio, n = 1) |> ungroup()
exact_end <- sw |> filter(band == "reaches the MCID") |> group_by(odds_ratio) |> summarise(n_loci = min(n_loci), .groups = "drop") |> slice_min(odds_ratio, n = 1)
ggsave("output/figures/fig_frontier_moves_profile.png",
  ggplot(sw, aes(n_loci, odds_ratio)) +
    geom_tile(aes(fill = band), colour = "white", linewidth = 0.2) +
    geom_line(data = moves, aes(group = line), colour = "white", linewidth = 0.9) +
    geom_point(data = moves, shape = 21, size = 1.9, stroke = 0.7, fill = jama_cols[["navy"]], colour = "white") +
    geom_text(data = ends, aes(x = 20.9, label = line), hjust = 0, size = 3.1, colour = jama_cols[["ink"]]) +
    geom_segment(data = ends, aes(x = n_loci + 0.35, xend = 20.75, yend = odds_ratio), colour = jama_cols[["grey"]], linewidth = 0.3) +
    annotate("text", x = 20.9, y = exact_end$odds_ratio, label = "causal variants, every case responds\n(edge of the darkest band)", hjust = 0, vjust = 0.75, size = 3.1, lineheight = 0.95, colour = jama_cols[["ink"]]) +
    scale_fill_manual(values = unname(jama_cols[c("light", "pale", "steel", "navy")]), name = NULL, drop = FALSE) +
    scale_x_continuous(breaks = c(1, 5, 10, 15, 20), expand = c(0, 0)) + scale_y_continuous(breaks = seq(1.5, 3, 0.5), expand = c(0, 0)) +
    coord_cartesian(clip = "off") +
    labs(x = "Number of independent risk loci", y = "Odds ratio per risk allele") +
    theme_jama(grid = "none") + theme(plot.margin = margin(8, 190, 8, 8)),
  width = 8.5, height = 4.8, dpi = 300)
print(as.data.frame(moves))
