proj_root <- normalizePath(file.path("..", ".."))
source(file.path(proj_root, "R", "baseline.R"))
source(file.path(proj_root, "R", "copula.R"))

pars    <- load_baseline_params(file.path(proj_root, "data", "baseline_model_parameters.csv"))
targets <- load_validation_targets(file.path(proj_root, "data", "baseline_validation_targets.csv"))
profiles <- all_profiles()
