# Decision-curve machinery on a discrete distribution of cells (here the 128 clinical
# profiles). Everything is an exact function of the cell probabilities: no patient-level
# sampling is involved. Decisions use the model's predicted risk; outcomes use the true risk,
# which equals the predicted risk for the baseline and will differ once genetic terms are
# added to the data-generating model.

# Net benefit and its components at each threshold.
#   w          cell probabilities (sum to 1)
#   risk_true  outcome probability in each cell
#   risk_pred  risk the decision rule uses; flag = risk_pred >= threshold
dca_exact <- function(w, risk_true, risk_pred = risk_true, thresholds = seq(0.01, 0.30, by = 0.001)) {
  stopifnot(length(w) == length(risk_true), length(w) == length(risk_pred),
            abs(sum(w) - 1) < 1e-6, all(thresholds > 0 & thresholds < 1))
  o <- order(risk_pred, decreasing = TRUE)
  w <- w[o]; rt <- risk_true[o]; rp <- risk_pred[o]
  cum <- list(tp = cumsum(w * rt), fp = cumsum(w * (1 - rt)), flagged = cumsum(w))
  k <- vapply(thresholds, function(t) sum(rp >= t), integer(1))
  pick <- function(x) ifelse(k == 0L, 0, x[pmax(k, 1L)])
  prev <- sum(w * rt)
  odds <- thresholds / (1 - thresholds)
  tp <- pick(cum$tp); fp <- pick(cum$fp)
  data.frame(
    threshold = thresholds,
    tp = tp, fp = fp, flagged = pick(cum$flagged),
    nb = tp - fp * odds,
    nb_all = prev - (1 - prev) * odds,
    nb_none = 0,
    snb = (tp - fp * odds) / prev,
    # MCID, reading (b): 20% more true cases at no increase in false positives
    mcid_b = 0.2 * tp,
    prevalence = prev
  )
}

# Compare a candidate with a baseline curve on the same threshold grid.
# The MCID is a property of the baseline, so it is taken from the baseline curve.
dca_compare <- function(base, new) {
  stopifnot(isTRUE(all.equal(base$threshold, new$threshold)))
  data.frame(
    threshold = base$threshold,
    delta_nb = new$nb - base$nb,
    mcid_b = base$mcid_b,
    meets_mcid = (new$nb - base$nb) >= base$mcid_b
  )
}

# Summary at chosen thresholds, per 1,000 admissions where useful
dca_summary <- function(curve, at = c(0.05, 0.10, 0.20)) {
  s <- curve[match(round(at, 6), round(curve$threshold, 6)), ]
  data.frame(
    threshold = s$threshold, nb = s$nb, nb_all = s$nb_all, snb = s$snb,
    tp_per_1000 = 1000 * s$tp, fp_per_1000 = 1000 * s$fp, flagged = s$flagged, mcid_b = s$mcid_b
  )
}

# Plots. `band` is an optional data frame with threshold, lo and hi for the model curve.
plot_dca <- function(curve, band = NULL, ylim = c(-0.02, 0.03), marks = c(0.05, 0.20)) {
  long <- rbind(
    data.frame(threshold = curve$threshold, strategy = "Baseline model", nb = curve$nb),
    data.frame(threshold = curve$threshold, strategy = "Treat all", nb = curve$nb_all),
    data.frame(threshold = curve$threshold, strategy = "Treat none", nb = curve$nb_none)
  )
  p <- ggplot2::ggplot(long, ggplot2::aes(threshold, nb))
  if (!is.null(band)) {
    p <- p + ggplot2::geom_ribbon(data = band, ggplot2::aes(threshold, ymin = lo, ymax = hi),
                                  inherit.aes = FALSE, alpha = 0.25)
  }
  p +
    ggplot2::geom_vline(xintercept = marks, linetype = "dashed", colour = "grey70") +
    ggplot2::geom_line(ggplot2::aes(colour = strategy)) +
    ggplot2::coord_cartesian(ylim = ylim) +
    ggplot2::scale_x_continuous(labels = scales::label_percent()) +
    ggplot2::labs(x = "Threshold probability", y = "Net benefit per admission", colour = NULL) +
    ggplot2::theme(legend.position = "bottom")
}

plot_dca_rates <- function(curve, marks = c(0.05, 0.20)) {
  long <- rbind(
    data.frame(threshold = curve$threshold, measure = "True positives per 1,000", value = 1000 * curve$tp),
    data.frame(threshold = curve$threshold, measure = "False positives per 1,000", value = 1000 * curve$fp)
  )
  ggplot2::ggplot(long, ggplot2::aes(threshold, value)) +
    ggplot2::geom_vline(xintercept = marks, linetype = "dashed", colour = "grey70") +
    ggplot2::geom_line() +
    ggplot2::facet_wrap(~measure, scales = "free_y") +
    ggplot2::scale_x_continuous(labels = scales::label_percent()) +
    ggplot2::labs(x = "Threshold probability", y = NULL)
}

plot_mcid <- function(curve, band = NULL, marks = c(0.05, 0.20)) {
  p <- ggplot2::ggplot(curve, ggplot2::aes(threshold, mcid_b))
  if (!is.null(band)) {
    p <- p + ggplot2::geom_ribbon(data = band, ggplot2::aes(threshold, ymin = lo, ymax = hi),
                                  inherit.aes = FALSE, alpha = 0.25)
  }
  p +
    ggplot2::geom_vline(xintercept = marks, linetype = "dashed", colour = "grey70") +
    ggplot2::geom_line() +
    ggplot2::scale_x_continuous(labels = scales::label_percent()) +
    ggplot2::labs(x = "Threshold probability",
                  y = "Net-benefit gain required (MCID, reading b)")
}
