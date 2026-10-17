# Run from the project root: source("tests/run_tests.R")
testthat::test_dir("tests/testthat", reporter = "summary", stop_on_failure = FALSE)
