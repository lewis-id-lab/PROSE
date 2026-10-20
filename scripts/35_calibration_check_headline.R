# The compact comparison table under the three comparisons (addendum 24): one row per analysis and metric, columns current, marginal and
# profile. Run from the project root, after scripts 31, 32 and 33. It reads saved outputs only.
#
# Writes   output/calibration_check_headline.csv   (read by tbl-calibration in manuscript.qmd)

cmp <- c("current", "marginal", "profile")
row3 <- function(analysis, metric, d, value) data.frame(analysis = analysis, metric = metric, t(setNames(value[match(cmp, d$comparison)], cmp)), row.names = NULL)

lay <- read.csv("output/calibration_check_layer.csv")
panel_rows <- do.call(rbind, lapply(c("four terms", "five terms"), function(pn) {
  x <- lay[lay$panel == pn & lay$truth == "mean" & lay$score == "ceiling", ]
  a <- paste0(pn, " (ceiling, pooled mean)")
  at <- function(t) x[x$threshold == t, ]
  rbind(row3(a, "gain in NB at 5% (x1000)", at(0.05), 1000 * at(0.05)$gain),
        row3(a, "P(gain >= MCID) at 5%", at(0.05), at(0.05)$p_meets_mcid),
        row3(a, "gain in NB at 20% (x1000)", at(0.20), 1000 * at(0.20)$gain),
        row3(a, "P(gain >= MCID) at 20%", at(0.20), at(0.20)$p_meets_mcid))
}))

a4 <- read.csv("output/calibration_check_aim4.csv")
at4 <- function(t) a4[a4$threshold == t, ]
sweep_rows <- rbind(
  row3("Aim 4 sweep, exact weights", "share of 4,000 architectures reaching MCID at 5%", at4(0.05), at4(0.05)$share_reaching_mcid),
  row3("Aim 4 sweep, exact weights", "share of 4,000 architectures reaching MCID at 20%", at4(0.20), at4(0.20)$share_reaching_mcid),
  row3("Aim 4 sweep, exact weights", "verdicts changed vs current at 5%", at4(0.05), at4(0.05)$verdicts_changed))

fo <- read.csv("output/calibration_check_aim4_frontier.csv")
front_rows <- do.call(rbind, lapply(c(1.5, 1.8, 2.0), function(o) {
  x <- fo[fo$threshold == 0.05 & fo$raf == 0.30 & fo$odds_ratio == o, ]
  row3("Frontier, exact weights, 5%, freq 0.3", paste0("fewest loci at OR ", o, " (NA = none within 20)"), x, x$min_loci)
}))

rs <- read.csv("output/calibration_check_realised_summary.csv"); rf <- read.csv("output/calibration_check_realised_frontier.csv")
real_rows <- do.call(rbind, lapply(c(0.10, 0.40), function(s) {
  x <- rs[rs$sigma == s & rs$threshold == 0.05, ]
  y <- rf[rf$kind == "min_loci" & rf$sigma == s & rf$threshold == 0.05 & rf$raf == 0.3 & rf$odds_ratio == 1.8, ]
  rbind(row3(paste0("Realised score, weight error ", s), "share of architectures reaching MCID at 5% (median replicate)", x, x$share_reaching_mcid),
        row3(paste0("Realised score, weight error ", s), "fewest loci at OR 1.8, freq 0.3, 5% (NA = none within 20)", y, y$median_reaches))
}))

nr <- read.csv("output/calibration_check_nri.csv")
nri_row <- function(sc, m, label = m, scale = 1) { x <- nr[nr$scenario == sc & nr$metric == m, ]; row3(paste0("NRI, ", sc), label, x, scale * x$median) }
nri_rows <- rbind(
  do.call(rbind, lapply(c("four-term panel, ceiling", "four-term panel, realised"), function(sc) rbind(
    nri_row(sc, "nri_cat5"), nri_row(sc, "nri_cat3"), nri_row(sc, "nri_cont"), nri_row(sc, "dnb_05", "gain in NB at 5% (x1000)", 1000)))),
  nri_row("PTX3 h2 haplotype", "nri_cat5"), nri_row("PTX3 h2 haplotype", "nri_cont"))

pr <- read.csv("output/calibration_check_prs.csv"); pc <- pr[pr$part == "ceiling", ]
prs_rows <- do.call(rbind, lapply(c(0.05, 0.20), function(t) {
  x <- pc[pc$threshold == t, ]; row3("PRS benchmark, ceiling", paste0("OR per SD needed to reach MCID at ", t * 100, "%"), x, x$or_sd_needed)
}))

headline <- rbind(panel_rows, sweep_rows, front_rows, real_rows, nri_rows, prs_rows)
write.csv(headline, "output/calibration_check_headline.csv", row.names = FALSE)
h <- headline; h[cmp] <- lapply(h[cmp], function(v) ifelse(is.na(v), "NA", signif(v, 3)))
print(h, row.names = FALSE)
