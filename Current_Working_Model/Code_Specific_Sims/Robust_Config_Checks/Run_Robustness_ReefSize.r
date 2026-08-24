# =============================================================================
#  Run_Robustness_ReefSize.r
# -----------------------------------------------------------------------------
#  ROBUSTNESS CHECK 3 of 3: do the results depend on reef (domain) size?
#  Re-runs a representative subset across reef sizes 30, 50, 75 and 100 so you can
#  check whether coexistence / exclusion outcomes are a real property or a
#  finite-size artifact of the 50x50 grid.
#
#  Re-uses the 4 base network models - NO new models needed - and varies ONLY the
#  reef size; every other setting is held at the standard baseline.
#
#  NOTE: larger reefs are slower per simulation (reef 100 has 4x the cells of reef 50),
#  so this run is dominated by the big-reef scenarios.
#
#  RUN (one command, no overrides needed):
#     Rscript Current_Working_Model/Run_Robustness_ReefSize.r
#  Output: Current_Working_Model/Results/Robustness_ReefSize_<timestamp>_results.{rds,csv}
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

# --- Robustness config: vary reef size, hold everything else at baseline -------------
SIM_CONFIG <- list(
  experiment    = paste0("Robustness_ReefSize_",
                         format(Sys.time(), "%Y%m%d_%H%M%S")),
  combinations  = c("Classic_RPS", "Classic_Linear", "Classic_Neutral", "Classic_Random"),
  n_species     = c(3, 7),            # low/high richness endpoints (5 dropped)
  sim_length    = 1000,               # baseline
  replicates    = 15,
  reef          = c(30, 50, 100),     # <-- the axis under test (standard run is 50; 75 dropped)
  individuals   = 3,                  # baseline (middle founder count)
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
cat("\nReef-size robustness run complete.\n")
