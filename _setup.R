suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(gt)
})

source("R/baseline.R")
source("R/copula.R")
# theme_jama() and the blue-grey palette (jama_cols, scale_colour_jama(), scale_fill_jama()); used by the manuscript figures
source("R/theme_jama.R")

pars    <- load_baseline_params("data/baseline_model_parameters.csv")
targets <- load_validation_targets("data/baseline_validation_targets.csv")

theme_set(theme_jama(base_size = 12))
