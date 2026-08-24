#  Run_SizeCompetitionStrength.r   (how much colony SIZE decides contests)
#  is shifted in log-odds by beta * log(size_i / size_j); this sweeps the beta magnitude.
#  RUN:  Rscript Model/Extended_Sensitivity_Checks/Run_SizeCompetitionStrength.r
#  Output: Model/Results/SizeCompetitionStrength_<timestamp>_results.{rds,csv}

.root <- (function() {
  a <- commandArgs(FALSE)
  f <- sub("^--file=", "", grep("^--file=", a, value = TRUE))
  d <- normalizePath(if (length(f)) dirname(f[1]) else getwd(), mustWork = FALSE)
  while (basename(d) != "Model" &&
         !dir.exists(file.path(d, "Model")) &&
         dirname(d) != d) d <- dirname(d)
  if (basename(d) == "Model") dirname(d) else d
})()
setwd(.root)
cat("Project root:", getwd(), "\n")

SIM_CONFIG <- list(
  experiment    = paste0("SizeCompetitionStrength_", format(Sys.time(), "%Y%m%d_%H%M%S")),
  combinations  = c("RPS_sizecompImpact", "Linear_SizeImpactNormal",
                    "Linear_SizeImpactReverse", "Random_sizecompImpact"),
  n_species     = c(3, 7),
  sim_length    = 2500,
  replicates    = 30,
  reef          = 50,
  individuals   = 3,
  biases        = 0.9,
  intraspecific = 0.5,
  size_beta_max = c(0, 0.3, 0.6, 0.9),           # <-- the axis under test (off / weak / default / strong)
  disturbances  = c("off", "on"),
  dist_freq     = "often",
  dist_size     = "random",
  out_dir       = "Model/Results",
  checkpoint_every = 10,
  base_seed     = 1000,
  to_master     = FALSE
)

source("Model/Simulation_testing.r")
cat("\nSize-competition-strength run complete.\n")
