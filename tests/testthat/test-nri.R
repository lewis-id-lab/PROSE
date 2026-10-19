source(file.path(proj_root, "R", "decision_curve.R"))
source(file.path(proj_root, "R", "nri.R"))

test_that("NRI is zero when the new score equals the old one", {
  w <- c(0.3, 0.7); rt <- c(0.1, 0.02); r <- c(0.08, 0.03)
  expect_equal(unname(nri_exact(w, rt, r, r, cuts = 0.05)), c(0, 0, 0))
  expect_equal(unname(nri_exact(w, rt, r, r)), c(0, 0, 0))
})

test_that("NRI matches a hand calculation", {
  # cell 1 moves up across 5%, cell 2 moves down across 5%
  w <- c(0.5, 0.5); rt <- c(0.10, 0.02); r0 <- c(0.04, 0.06); r1 <- c(0.06, 0.03)
  ev <- (0.05 - 0.01) / 0.06
  ne <- (0.49 - 0.45) / 0.94
  expect_equal(unname(nri_exact(w, rt, r0, r1, cuts = 0.05)), c(ev, ne, ev + ne))
  expect_equal(unname(nri_exact(w, rt, r0, r1)), c(ev, ne, ev + ne))
  # a move that stays inside a category counts for the continuous NRI only
  expect_equal(unname(nri_exact(w, rt, c(0.06, 0.06), c(0.07, 0.06), cuts = 0.05)), c(0, 0, 0))
  expect_gt(nri_exact(w, rt, c(0.06, 0.06), c(0.07, 0.06))[["events"]], 0)
})

test_that("reclassification metrics are zero for an unchanged score", {
  w <- c(0.2, 0.3, 0.5); lp <- qlogis(c(0.15, 0.04, 0.01))
  m <- reclassification_metrics(w, lp, lp, lp)
  expect_equal(unname(m[c("nri_cat5", "nri_cat3", "nri_cont", "dnb_05", "dnb_20", "dauc")]), rep(0, 6))
})

test_that("the PTX3 diplotype grid reproduces the phased haplotype structure", {
  hap <- file.path(proj_root, "data", "haplotype_frequencies_phased.csv")
  skip_if_not(file.exists(hap))
  g <- ptx3_diplotype_grid(hap)
  expect_equal(sum(g$prob), 1)
  freq <- as.numeric(crossprod(g$G, g$prob)) / 2
  expect_equal(freq, c(0.482, 0.482, 0.4242, 0.4806), tolerance = 1e-3)
  expect_identical(g$G[, "h2"], g$G[, "rs2305619"])
  # linkage: rs3816527 C occurs only on h2 chromosomes, so its dosage never exceeds the h2 dosage
  expect_true(all(g$G[, "rs3816527"] <= g$G[, "h2"]))
})
