# =============================================================================
#  Run_Sensitivity_Intraspecific.r
# -----------------------------------------------------------------------------
#  PARAMETER SENSITIVITY: intraspecific competition strength.
#  Sweeps the matrix DIAGONAL (the same-species overgrowth probability) from 0 (a
#  colony never overgrows another colony of its OWN species) to 1 (always), on the
#  4 base network models. Modern coexistence theory predicts coexistence is stabilised
#  when intraspecific competition is strong relative to interspecific - this tests
#  whether the model behaves that way.
#
#  NOTE: for Classic_Neutral the diagonal is the ONLY non-zero interaction, so this
#  sweep is especially meaningful there (it is the only competition present).
#
#  Standard settings otherwise: 30 replicates, 1000 timesteps, reef 50, species {3,5,7},
#  bias 0.9, disturbance off/on at the baseline regime. Founders fixed at 3.
#
#  RUN:  Rscript Current_Working_Model/Sensitivity_Config_Checks/Run_Sensitivity_Intraspecific.r
#  Output: Current_Working_Model/Results/Sensitivity_Intraspecific_<timestamp>_results.{rds,csv}
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

# --- Sensitivity config: sweep intraspecific (diagonal), baseline otherwise ----------
SIM_CONFIG <- list(
  experiment    = paste0("Sensitivity_Intraspecific_", format(Sys.time(), "%Y%m%d_%H%M%S")),
  combinations  = c("Classic_RPS", "Classic_Linear", "Classic_Neutral", "Classic_Random"),
  n_species     = c(3, 5, 7),
  sim_length    = 1000,                          # standard
  replicates    = 30,                            # standard
  reef          = 50,                            # standard
  individuals   = 3,                             # fixed (not the axis under test)
  biases        = 0.9,                           # baseline
  intraspecific = c(0, 0.25, 0.5, 0.75, 1.0),    # <-- the axis under test (diagonal)
  growth        = 3,
  disturbances  = c("off", "on"),
  dist_freq     = "often",                       # baseline regime
  dist_size     = "random",                      # baseline regime
  out_dir       = "Current_Working_Model/Results",
  checkpoint_every = 10,
  base_seed     = 1000,
  to_master     = FALSE
)

# Simulation_testing.r sources the model, runs the grid from SIM_CONFIG, and saves.
source("Current_Working_Model/Simulation_testing.r")
cat("\nIntraspecific-competition sensitivity run complete.\n")
