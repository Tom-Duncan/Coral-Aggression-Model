#  Run_Robustness_Replicates.r
#  Re-uses the 4 base network models - NO new models needed - and varies ONLY the
#  Output: Model/Results/Robustness_Replicates_<timestamp>_results.{rds,csv}

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

# --- Robustness config: vary replicates, hold everything else at baseline ------------
SIM_CONFIG <- list(
  experiment    = paste0("Robustness_Replicates_",
                         format(Sys.time(), "%Y%m%d_%H%M%S")),
  combinations  = c("Classic_RPS", "Classic_Linear", "Classic_Neutral", "Classic_Random"),
  n_species     = c(3, 7),       # low/high richness endpoints (5 dropped)
  sim_length    = 1000,          # baseline
  replicates    = 100,           # <-- the axis under test (kept high; standard run is 30)
  reef          = 50,            # baseline
  individuals   = 3,             # baseline (middle founder count)
  biases        = 0.9,
  intraspecific = 0.5,
  growth        = 3,
  disturbances  = c("off", "on"),
  dist_freq     = "often",
  dist_size     = "random",
  out_dir       = "Model/Results",
  checkpoint_every = 10,
  base_seed     = 1000,
  to_master     = FALSE
)

# Simulation_testing.r sources the model, runs the grid from SIM_CONFIG, and saves.
source("Model/Simulation_testing.r")
cat("\nReplicate-adequacy robustness run complete.\n")
