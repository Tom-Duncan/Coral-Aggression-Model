# =============================================================================
#  Run_Reproduction_DisturbanceRegime.r  (Reproduction experiment B)
# -----------------------------------------------------------------------------
#  Tests whether recruitment's benefit depends on the DISTURBANCE REGIME - i.e. does
#  reproduction matter more when die-offs are frequent/severe, because recruitment
#  enables post-disturbance recolonisation (a recovery / storage-like effect)?
#
#  Contrast: reproduction OFF (Classic_*) vs ON (equal recruitment, *_reproductionEven)
#  for each network, crossed with the FULL disturbance regime grid (frequency x
#  intensity) plus the off baseline.
#
#  Run length 2500 (recruitment + recovery dynamics are slow).
#
#  RUN:  Rscript Current_Working_Model/Reproduction_Config_Checks/Run_Reproduction_DisturbanceRegime.r
#  Output: Current_Working_Model/Results/Reproduction_DisturbanceRegime_<timestamp>_results.{rds,csv}
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
  experiment    = paste0("Reproduction_DisturbanceRegime_", format(Sys.time(), "%Y%m%d_%H%M%S")),
  combinations  = c("Classic_RPS",     "RPS_reproductionEven",
                    "Classic_Linear",  "Linear_reproductionEven",
                    "Classic_Neutral", "Neutral_reproductionEven",
                    "Classic_Random",  "Random_reproductionEven"),   # OFF vs ON per network
  n_species     = c(3, 7),
  sim_length    = 2500,
  replicates    = 30,
  reef          = 50,
  individuals   = 3,
  biases        = 0.9,
  intraspecific = 0.5,
  growth        = 3,
  disturbances  = c("off", "on"),
  dist_freq     = c("rarely", "very_often"),            # <-- frequency extremes (middle "often" dropped)
  dist_size     = c("small", "large"),                  # <-- intensity extremes (medium dropped)
  out_dir       = "Current_Working_Model/Results",
  checkpoint_every = 10,
  base_seed     = 1000,
  to_master     = FALSE
)

source("Current_Working_Model/Simulation_testing.r")
cat("\nReproduction x disturbance-regime run complete.\n")
