#  Run_MaturityAge.r   (life-history: age at first reproduction)
#  Sweeps maturity age: 10 (breeds young / fast life history), 30 (default), 60 (slow;
#  RUN:  Rscript Model/Extended_Sensitivity_Checks/Run_MaturityAge.r
#  Output: Model/Results/MaturityAge_<timestamp>_results.{rds,csv}

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
  experiment    = paste0("MaturityAge_", format(Sys.time(), "%Y%m%d_%H%M%S")),
  combinations  = c("RPS_reproductionEven", "Linear_reproductionEven",
                    "Neutral_reproductionEven", "Random_reproductionEven"),
  n_species     = c(3, 7),
  sim_length    = 2500,
  replicates    = 30,
  reef          = 50,
  individuals   = 3,
  biases        = 0.9,
  intraspecific = 0.5,
  maturity_age  = c(5, 30, 150),                 # <-- the axis under test, widened (very fast / default / very slow)
  disturbances  = c("off", "on"),
  dist_freq     = "often",
  dist_size     = "random",
  out_dir       = "Model/Results",
  checkpoint_every = 10,
  base_seed     = 1000,
  to_master     = FALSE
)

source("Model/Simulation_testing.r")
cat("\nMaturity-age (life-history) run complete.\n")
