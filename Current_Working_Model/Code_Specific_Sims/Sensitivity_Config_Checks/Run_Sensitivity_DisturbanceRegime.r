# =============================================================================
#  Run_Sensitivity_DisturbanceRegime.r
# -----------------------------------------------------------------------------
#  PARAMETER SENSITIVITY: disturbance regime (frequency x intensity).
#  Crosses disturbance FREQUENCY {rarely, often, very_often} with INTENSITY (patch
#  size) {small, medium, large} - 9 regimes - plus a disturbance-OFF baseline, on the
#  4 base network models. Tests how the frequency and intensity of mass die-offs shape
#  coexistence/richness (e.g. the Intermediate Disturbance Hypothesis: does a moderate
#  regime maximise richness?).
#
#  The engine now sweeps AND records dist_freq/dist_size, so every row is labelled with
#  its regime (off rows carry NA for both, and are not duplicated across regimes).
#
#  Standard settings otherwise: 30 replicates, 1000 timesteps, reef 50, species {3,5,7},
#  bias 0.9. Founders fixed at 3 (single value) so the regime axes are not multiplied
#  out - widen `individuals` if you want the founder sweep too.
#
#  Frequency: rarely ~1 / often ~2 / very_often ~4 events per 100 steps.
#  Intensity: small 10% / medium 20% / large 30% of the reef per event.
#
#  RUN:  Rscript Current_Working_Model/Sensitivity_Config_Checks/Run_Sensitivity_DisturbanceRegime.r
#  Output: Current_Working_Model/Results/Sensitivity_DisturbanceRegime_<timestamp>_results.{rds,csv}
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
  out_dir       = "Current_Working_Model/Results",
  checkpoint_every = 10,
  base_seed     = 1000,
  to_master     = FALSE
)

# Simulation_testing.r sources the model, runs the grid from SIM_CONFIG, and saves.
source("Current_Working_Model/Simulation_testing.r")
cat("\nDisturbance-regime sensitivity run complete.\n")
