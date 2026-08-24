#  Run_Sensitivity_DisturbanceRegime.r
#  PARAMETER SENSITIVITY: disturbance regime (frequency x intensity).
#  RUN:  Rscript Model/Sensitivity_Config_Checks/Run_Sensitivity_DisturbanceRegime.r
#  Output: Model/Results/Sensitivity_DisturbanceRegime_<timestamp>_results.{rds,csv}

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

# --- Sensitivity config: sweep disturbance frequency x intensity, baseline otherwise --
SIM_CONFIG <- list(
  experiment    = paste0("Sensitivity_DisturbanceRegime_",
                         format(Sys.time(), "%Y%m%d_%H%M%S")),
  combinations  = c("Classic_RPS", "Classic_Linear", "Classic_Neutral", "Classic_Random"),
  n_species     = c(3, 5, 7),
  sim_length    = 1000,                                   # standard
  replicates    = 30,                                     # standard
  reef          = 50,                                     # standard
  individuals   = 3,                                      # fixed (not an axis here)
  biases        = 0.9,                                    # baseline
  intraspecific = 0.5,
  growth        = 3,
  disturbances  = c("off", "on"),                         # off baseline + on regimes
  dist_freq     = c("rarely", "often", "very_often"),     # <-- frequency axis
  dist_size     = c("small", "medium", "large"),          # <-- intensity axis
  out_dir       = "Model/Results",
  checkpoint_every = 10,
  base_seed     = 1000,
  to_master     = FALSE
)

# Simulation_testing.r sources the model, runs the grid from SIM_CONFIG, and saves.
source("Model/Simulation_testing.r")
cat("\nDisturbance-regime sensitivity run complete.\n")
