# Quarto starts R in this folder when it renders the submission; use the project's renv library (the project root is one level up)
local({
  owd <- setwd("..")
  on.exit(setwd(owd))
  source("renv/activate.R")
})
