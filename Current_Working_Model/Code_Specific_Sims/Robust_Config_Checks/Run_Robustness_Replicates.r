# =============================================================================
#  Run_Robustness_Replicates.r
# -----------------------------------------------------------------------------
#  ROBUSTNESS CHECK 2 of 3: are 30 replicates enough for stable estimates?
#  Re-runs a representative subset with 100 replicates so you can compute the running
#  mean / CI over 1..100 reps and see where it stabilises. If it settles well before
#  30, the standard 30 is enough; if not, high-variance scenarios need more.
#
#  Re-uses the 4 base network models - NO new models needed - and varies ONLY the
#  replicate count; every other setting is held at the standard baseline.
#
#  RUN (one command, no overrides needed):
#     Rscript Current_Working_Model/Run_Robustness_Replicates.r
#  Output: Current_Working_Model/Results/Robustness_Replicates_<timestamp>_results.{rds,csv}
# =============================================================================

# --- Locate the project root (the folder that CONTAINS Current_Working_Model) --------
.root <- (function() {
  a <- commandArgs(FALSE)
  f <- sub("^--file=", "", grep("^--file=", a, value = TRUE))
  d <- normalizePath(if (length(f)) dirname(f[1]) else getwd(), mustWork = FALSE)
  while (basename(d) != "Current_Working_Model" &&
         !dir.exists(file.path(d, "Current_Working_Model")) &&
         dirname(d) != d) d <- dirname(d)
  if (basename(d) == "Current_Working_Model") dirname(d) else d
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
  out_dir       = "Current_Working_Model/Results",
  checkpoint_every = 10,
  base_seed     = 1000,
  to_master     = FALSE
)

# Simulation_testing.r sources the model, runs the grid from SIM_CONFIG, and saves.
source("Current_Working_Model/Simulation_testing.r")
cat("\nReplicate-adequacy robustness run complete.\n")
