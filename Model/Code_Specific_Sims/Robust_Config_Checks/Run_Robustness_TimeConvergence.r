# =============================================================================
#  Run_Robustness_TimeConvergence.r
# -----------------------------------------------------------------------------
#  ROBUSTNESS CHECK 1 of 3: has the model reached equilibrium by timestep 1000?
#  Re-runs a representative subset for 3000 timesteps (3x the standard length) so you
#  can confirm the metrics plateau well before 1000 (rather than still drifting).
#
#  Re-uses the 4 base network models - NO new models needed - and varies ONLY the
#  simulation length; every other setting is held at the standard baseline.
#
#  RUN (one command, no overrides needed):
#     Rscript Model/Run_Robustness_TimeConvergence.r
#  Output: Model/Results/Robustness_TimeConvergence_<timestamp>_results.{rds,csv}
# =============================================================================

# --- Locate the project root (the folder that CONTAINS Model) --------
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

# --- Robustness config: vary sim_length, hold everything else at baseline ------------
SIM_CONFIG <- list(
  experiment    = paste0("Robustness_TimeConvergence_",
                         format(Sys.time(), "%Y%m%d_%H%M%S")),
  combinations  = c("Classic_RPS", "Classic_Linear", "Classic_Neutral", "Classic_Random"),
  n_species     = c(3, 7),       # low/high richness endpoints (5 dropped)
  sim_length    = 3000,          # <-- the axis under test (standard run is 1000)
  replicates    = 15,            # standard 15 (only the replicate-adequacy check uses more)
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
cat("\nTime-convergence robustness run complete.\n")
