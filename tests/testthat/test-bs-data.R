source(file.path(proj_root, "R", "bs_data.R"))

# Synthetic data with the published structure, so the checking code is tested without patient data
make_fake_bs <- function(n = 4000, seed = 4719) {
  set.seed(seed)
  X <- sapply(PREDICTORS, function(p) rbinom(n, 1, pars$marginal[[p]]))
  lp <- pars$intercept + X %*% pars$coef[PREDICTORS]
  list(X = X, y = rbinom(n, 1, plogis(lp)), id = seq_len(n))
}

test_that("auc_rank agrees with a brute-force pairwise AUC, including ties", {
  set.seed(5503)
  score <- round(rnorm(60), 1); y <- rbinom(60, 1, 0.4)
  brute <- 0
  for (i in which(y == 1)) for (j in which(y == 0)) brute <- brute + (score[i] > score[j]) + 0.5 * (score[i] == score[j])
  expect_equal(auc_rank(score, y), brute / (sum(y == 1) * sum(y == 0)))
})

test_that("count and refit checks return the expected shapes on synthetic data", {
  bs <- make_fake_bs()
  counts <- check_bs_counts(bs, pars)
  expect_equal(nrow(counts), 7L)
  expect_equal(counts$cases_data + counts$non_cases_data, as.numeric(colSums(bs$X)))
  refit <- check_bs_refit(bs, pars)
  expect_equal(nrow(refit), 8L)
  expect_equal(refit$term[1], "(Intercept)")
  perf <- bs_apparent_performance(bs, pars)
  expect_gt(perf$apparent_auc, 0.7)
})

test_that("synthetic counts do not match the documented differences, but constructed ones do", {
  bs <- make_fake_bs()
  counts <- check_bs_counts(bs, pars)
  expect_false(bs_counts_as_expected(counts))

  exact <- counts
  exact$cases_data <- exact$cases_published; exact$non_cases_data <- exact$non_cases_published
  exact$diff_cases <- 0; exact$diff_non_cases <- 0
  expect_false(bs_counts_as_expected(exact))    # the known differences are required, not optional

  exact$diff_non_cases[exact$predictor %in% c("high_risk_chemo", "prolonged_neutropenia")] <- -1
  expect_true(bs_counts_as_expected(exact))

  exact$diff_cases[1] <- 1                      # one unexpected difference is a failure
  expect_false(bs_counts_as_expected(exact))
})

test_that("read_bs_data fails clearly when the file is absent", {
  expect_error(read_bs_data(file.path(tempdir(), "no-such-file.xlsx")), "gitignored")
})

bs_path <- file.path(proj_root, "bs.xlsx")

test_that("the derivation data reproduces the published counts, apart from the documented differences", {
  skip_if_not(file.exists(bs_path), "bs.xlsx is gitignored patient-level data and is not present")
  bs <- read_bs_data(bs_path)
  expect_equal(length(bs$y), pars$n_case + pars$n_non)
  expect_equal(sum(bs$y), pars$n_case)
  expect_false(anyDuplicated(bs$id) > 0)
  expect_true(bs_counts_as_expected(check_bs_counts(bs, pars)))
})

test_that("a refit on the derivation data recovers the published coefficients and aROC", {
  skip_if_not(file.exists(bs_path), "bs.xlsx is gitignored patient-level data and is not present")
  bs <- read_bs_data(bs_path)
  expect_lt(max(abs(check_bs_refit(bs, pars)$diff)), 0.01)
  perf <- bs_apparent_performance(bs, pars)
  expect_lt(abs(perf$apparent_auc - targets$aroc), 0.005)
  expect_lt(abs(perf$mean_predicted_risk - perf$observed_prevalence), 0.002)
})

test_that("empirical profile counts add up and index the right profiles", {
  bs <- make_fake_bs()
  cnt <- empirical_profile_counts(bs)
  expect_length(cnt, 128L)
  expect_equal(sum(cnt), nrow(bs$X))
  # A row with all predictors present lands in the last profile, none present in the first
  X <- rbind(rep(0, 7), rep(1, 7)); colnames(X) <- PREDICTORS
  cnt2 <- empirical_profile_counts(list(X = X))
  expect_equal(cnt2[c(1, 128)], c(1, 1))
  expect_equal(sum(cnt2), 2)
})

test_that("tetrachoric recovers the latent correlation of simulated binary variables", {
  set.seed(9036)
  for (rho in c(-0.3, 0, 0.5)) {
    z <- mvtnorm::rmvnorm(40000, sigma = matrix(c(1, rho, rho, 1), 2))
    x <- as.integer(z[, 1] > qnorm(0.7)); y <- as.integer(z[, 2] > qnorm(0.85))
    expect_equal(tetrachoric(x, y), rho, tolerance = 0.06, info = paste("rho =", rho))
  }
})

test_that("the real derivation data has 4,127 admissions across 75 observed profiles", {
  skip_if_not(file.exists(bs_path), "bs.xlsx is gitignored patient-level data and is not present")
  cnt <- empirical_profile_counts(read_bs_data(bs_path))
  expect_equal(sum(cnt), 4127)
  expect_equal(sum(cnt > 0), 75L)
})
