suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(gt)
})

source("R/baseline.R")
source("R/copula.R")

pars    <- load_baseline_params("data/baseline_model_parameters.csv")
targets <- load_validation_targets("data/baseline_validation_targets.csv")

theme_set(theme_minimal(base_size = 12))
