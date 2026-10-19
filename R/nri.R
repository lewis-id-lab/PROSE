# Net reclassification improvement (NRI), computed exactly over weighted cells (addendum 22).
#
# NRI is a tertiary comparator only (project description, section 3). It counts movements across risk categories without weighing
# their consequences, so it is reported to show how the predecessor's headline figure (22.5%) arose, not as an endpoint.
#
# Requires R/baseline.R (expected_auc) and R/decision_curve.R (dca_exact).

# Exact NRI over cells with weights w, true event probability r_true, and old and new predicted risks.
#   cuts  category boundaries for the categorical NRI (a risk at or above a cut is in the higher category, as in dca_exact);
#         NULL gives the continuous, category-free NRI, in which any increase in predicted risk counts as a move up.
# Returns the event and non-event components and their sum (scale -2 to 2; the predecessor's 22.5% is 0.225 on this scale).
nri_exact <- function(w, r_true, r_old, r_new, cuts = NULL, tol = 1e-12) {
  stopifnot(length(w) == length(r_true), length(w) == length(r_old), length(w) == length(r_new))
  if (is.null(cuts)) {
    up <- r_new > r_old + tol
    down <- r_new < r_old - tol
  } else {
    c_old <- findInterval(r_old, sort(cuts))
    c_new <- findInterval(r_new, sort(cuts))
    up <- c_new > c_old
    down <- c_new < c_old
  }
  we <- w * r_true
  wn <- w * (1 - r_true)
  events <- (sum(we[up]) - sum(we[down])) / sum(we)
  nonevents <- (sum(wn[down]) - sum(wn[up])) / sum(wn)
  c(events = events, nonevents = nonevents, nri = events + nonevents)
}

# NRI (categorical at 5%, three categories at 5% and 20%, continuous), net-benefit gain and aROC change of a new score against the
# baseline, all under the same truth. Arguments are cell weights and logits.
reclassification_metrics <- function(w, lp_true, lp_base, lp_new, thresholds = c(0.05, 0.20)) {
  r_true <- stats::plogis(lp_true)
  r0 <- stats::plogis(lp_base)
  r1 <- stats::plogis(lp_new)
  d0 <- dca_exact(w, r_true, r0, thresholds = thresholds)
  d1 <- dca_exact(w, r_true, r1, thresholds = thresholds)
  n5 <- nri_exact(w, r_true, r0, r1, cuts = 0.05)
  n3 <- nri_exact(w, r_true, r0, r1, cuts = c(0.05, 0.20))
  nc <- nri_exact(w, r_true, r0, r1)
  auc <- function(lp) expected_auc(lp, w * r_true, w * (1 - r_true))
  c(nri_cat5 = n5[["nri"]], nri_cat5_events = n5[["events"]], nri_cat5_nonevents = n5[["nonevents"]],
    nri_cat3 = n3[["nri"]],
    nri_cont = nc[["nri"]], nri_cont_events = nc[["events"]], nri_cont_nonevents = nc[["nonevents"]],
    dnb_05 = d1$nb[1] - d0$nb[1], dnb_20 = d1$nb[2] - d0$nb[2],
    dauc = auc(lp_new) - auc(lp_base))
}

# The PTX3 trio as diplotypes of phased haplotypes (alleles in the order rs2305619, rs3816527, rs1840680; gnomAD NFE rescaled
# frequencies), with the h2 dosage (A at rs2305619, D11) and the dosage of each SNP's reported risk allele (rs2305619 A,
# rs3816527 C, rs1840680 A). The h2 and rs2305619 columns are equal by the definition of h2; both are kept so that a truth on h2
# and a score on the three SNPs can be written as separate effect vectors over the same grid.
ptx3_diplotype_grid <- function(path = "data/haplotype_frequencies_phased.csv") {
  h <- utils::read.csv(path)
  h <- h[h$block == "PTX3_trio" & h$freq_gnomAD_NFE_rescaled > 0, ]
  f <- h$freq_gnomAD_NFE_rescaled / sum(h$freq_gnomAD_NFE_rescaled)
  pr <- which(upper.tri(diag(nrow(h)), diag = TRUE), arr.ind = TRUE)
  i <- pr[, 1]
  j <- pr[, 2]
  cnt <- function(col, allele) as.numeric(h[[col]][i] == allele) + as.numeric(h[[col]][j] == allele)
  G <- cbind(h2 = cnt("allele_1", "A"), rs2305619 = cnt("allele_1", "A"),
             rs3816527 = cnt("allele_2", "C"), rs1840680 = cnt("allele_3", "A"))
  list(G = G, prob = ifelse(i == j, f[i]^2, 2 * f[i] * f[j]),
       diplotype = paste(h$haplotype[i], h$haplotype[j], sep = "/"))
}
