source(file.path(proj_root, "R", "pseudonymise.R"))

# Synthetic export with the same shape as the real one, including fake identifiers
make_fake_boscore <- function(n_patients = 40, seed = 7208) {
  set.seed(seed)
  pat <- sample(seq_len(n_patients), 150, replace = TRUE)
  adm <- as.integer(factor(paste(pat, sample(1:3, 150, replace = TRUE))))   # admissions nested in patients
  admit <- as.POSIXct("2010-01-01", tz = "UTC") + 86400 * sample(0:2000, 150, replace = TRUE)
  d <- data.frame(
    Reparto = sample(c("Ward A", "Ward B", "Ward C"), 150, replace = TRUE),
    IDRICO = 7000 + adm,
    IDIND = paste0("ID-", pat),
    Connome = paste("Surname", pat),
    Dnas = sprintf("%02d/%02d/%d", sample(1:28, 150, TRUE), sample(1:12, 150, TRUE), sample(1940:1995, 150, TRUE)),
    AGEAMM = sample(c(20:80, 91, 95), 150, replace = TRUE),
    DRICEM = admit,
    DDIMEM = admit + 86400 * 14,
    GGEM = 14,
    DTMO = admit + 86400 * 3,
    EMODES = sample(c("free text a", "free text b"), 150, replace = TRUE),
    DMOR = ifelse(runif(150) < 0.2, format(as.Date(admit) + 100, "%d/%m/%Y"), NA),
    DLFUP = format(as.Date(admit) + 200, "%d/%m/%Y"),
    DSTATUS = NA_character_,
    IDEPISOD = ifelse(runif(150) < 0.7, 500 + seq_len(150), NA),
    DEPISOD = admit + 86400 * 5,
    Expr1064 = as.POSIXct(NA, tz = "UTC"),
    TabMicetiIsolatiDISOL = as.POSIXct(NA, tz = "UTC"),
    TabBatteriIsolatiDISOL = admit + 86400 * 7,
    IDMICE = NA_real_, IDBACT = 1:150,
    neutropenia = rbinom(150, 1, 0.4),
    stringsAsFactors = FALSE
  )
  d
}

test_that("direct identifiers and calendar dates are removed", {
  d <- make_fake_boscore()
  ps <- pseudonymise_boscore(d, seed = 1)$data
  expect_true(all(unlist(check_pseudonymised(ps))))
  expect_false(any(c("Connome", "Dnas", "IDIND", "IDRICO", "DRICEM") %in% names(ps)))
  # No original identifier or name value survives anywhere in the extract
  everything <- unlist(lapply(ps, as.character))
  expect_false(any(grepl("Surname|ID-", everything)))
})

test_that("structure is preserved and every mapping is one to one", {
  d <- make_fake_boscore()
  ps <- pseudonymise_boscore(d, seed = 2)$data
  v <- verify_pseudonymisation(d, ps)
  expect_true(all(unlist(v)))
})

test_that("dates become offsets from admission, and the year is kept", {
  d <- make_fake_boscore()
  ps <- pseudonymise_boscore(d, seed = 3)$data
  expect_true(all(ps$days_from_admission_DTMO == 3))
  expect_true(all(ps$days_from_admission_DLFUP == 200))
  expect_true(all(ps$days_from_admission_DEPISOD == 5))
  expect_equal(ps$admission_year, as.integer(format(as.Date(d$DRICEM), "%Y")))
  expect_equal(sum(!is.na(ps$days_from_admission_DMOR)), sum(!is.na(d$DMOR)))
})

test_that("age is top-coded at 90 and clinical variables are untouched", {
  d <- make_fake_boscore()
  ps <- pseudonymise_boscore(d, seed = 4)$data
  expect_lte(max(ps$age_years), 90)
  expect_equal(ps$age_years, pmin(d$AGEAMM, 90))
  expect_equal(ps$neutropenia, d$neutropenia)
  expect_equal(ps$GGEM, d$GGEM)
})

test_that("codes are reproducible for a seed and differ between seeds", {
  d <- make_fake_boscore()
  a <- pseudonymise_boscore(d, seed = 5)$data$patient_code
  b <- pseudonymise_boscore(d, seed = 5)$data$patient_code
  c2 <- pseudonymise_boscore(d, seed = 6)$data$patient_code
  expect_identical(a, b)
  expect_false(identical(a, c2))
})

test_that("the checker catches a leaked date column, date string and identifier", {
  d <- make_fake_boscore()
  ps <- pseudonymise_boscore(d, seed = 7)$data
  leak1 <- ps; leak1$when <- as.Date("2010-05-05")
  expect_false(check_pseudonymised(leak1)$no_date_columns)
  leak2 <- ps; leak2$note <- "seen on 12/03/2011"
  expect_false(check_pseudonymised(leak2)$no_date_strings)
  for (s in c("12.03.2011", "12-03-2011", "2011-03-12", "on 3/4/11 seen")) {
    leak <- ps; leak$note <- s
    expect_false(check_pseudonymised(leak)$no_date_strings, info = s)
  }
  ok <- ps; ok$temp <- c("38.6-40°C", "38-38.5°C", ">40", "< 0,5", "0,5 - 0,7")[1:nrow(ps) %% 5 + 1]
  expect_true(check_pseudonymised(ok)$no_date_strings)    # category labels are not dates
  leak3 <- ps; leak3$IDIND <- "x"
  expect_false(check_pseudonymised(leak3)$no_identifier_columns)
  leak4 <- ps; leak4$age_years[1] <- 95
  expect_false(check_pseudonymised(leak4)$age_top_coded)
})

test_that("the manifest accounts for every original column", {
  d <- make_fake_boscore()
  m <- pseudonymise_boscore(d, seed = 8)$manifest
  expect_equal(m$column, names(d))
  expect_equal(m$treatment[m$column == "Connome"], "dropped: direct identifier")
  expect_equal(m$treatment[m$column == "neutropenia"], "kept unchanged")
})

extract_path <- file.path(proj_root, "data", "restricted", "boscore_pseudonymised.csv")

test_that("the real extract passes every check", {
  skip_if_not(file.exists(extract_path), "the pseudonymised extract is gitignored and not present")
  x <- read.csv(extract_path, stringsAsFactors = FALSE)
  expect_true(all(unlist(check_pseudonymised(x))))
})
