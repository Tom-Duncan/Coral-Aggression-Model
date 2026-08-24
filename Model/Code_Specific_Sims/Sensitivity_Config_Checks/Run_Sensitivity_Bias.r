#  Run_Sensitivity_Bias.r
#  PARAMETER SENSITIVITY: competition bias.
#  RUN:  Rscript Model/Sensitivity_Config_Checks/Run_Sensitivity_Bias.r
#  Output: Model/Results/Sensitivity_Bias_<timestamp>_results.{rds,csv}

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

# --- Sensitivity config: sweep bias, hold everything else at baseline ----------------
SIM_CONFIG <- list(
  experiment    = paste0("Sensitivity_Bias_", format(Sys.time(), "%Y%m%d_%H%M%S")),
  combinations  = c("Classic_RPS", "Classic_Linear", "Classic_Neutral", "Classic_Random"),
  n_species     = c(3, 5, 7),
  sim_length    = 1000,                          # standard
  replicates    = 30,                            # standard
  reef          = 50,                            # standard
  individuals   = 3,                             # fixed (not the axis under test)
  biases        = c(0.5, 0.6, 0.7, 0.8, 0.9, 1.0),  # <-- the axis under test
  intraspecific = 0.5,
  growth        = 3,
  disturbances  = c("off", "on"),
  dist_freq     = "often",                       # baseline regime
  dist_size     = "random",                      # baseline regime
  out_dir       = "Model/Results",
  checkpoint_every = 10,
  base_seed     = 1000,
  to_master     = FALSE
)

# Simulation_testing.r sources the model, runs the grid from SIM_CONFIG, and saves.
source("Model/Simulation_testing.r")
cat("\nCompetition-bias sensitivity run complete.\n")
