# =============================================================================
#  Run_MaturityAge.r   (life-history: age at first reproduction)
# -----------------------------------------------------------------------------
#  Tests how the AGE AT FIRST REPRODUCTION shapes the community. A colony must survive
#  MATURITY_AGE timesteps before it can spawn any recruit. This is a life-history
#  bottleneck that interacts with disturbance: if die-offs kill colonies before they
#  mature, reproduction is effectively suppressed, so a species must reach reproductive
#  age FAST to persist under disturbance.
#
#  Sweeps maturity age: 10 (breeds young / fast life history), 30 (default), 60 (slow;
#  must survive long before contributing recruits). Maturity only bites when reproduction
#  is ON, so it is tested on the equal-reproduction models of the 4 networks (a clean
#  "recruitment present" set that isolates the maturity effect).
#
#  Prediction: high maturity age should hurt most under disturbance (the recruitment
#  bottleneck), and the OFF arm shows the effect without the disturbance interaction.
#
#  Run length 2500.
#
#  RUN:  Rscript Current_Working_Model/Extended_Sensitivity_Checks/Run_MaturityAge.r
#  Output: Current_Working_Model/Results/MaturityAge_<timestamp>_results.{rds,csv}
# =============================================================================

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
  out_dir       = "Current_Working_Model/Results",
  checkpoint_every = 10,
  base_seed     = 1000,
  to_master     = FALSE
)

source("Current_Working_Model/Simulation_testing.r")
cat("\nMaturity-age (life-history) run complete.\n")
