# =============================================================================
#  Run_Reproduction_Fecundity.r   (Reproduction experiment A: dose-response)
# -----------------------------------------------------------------------------
#  Tests how the MAGNITUDE of reproduction (base per-timestep spawn chance) affects
#  the outcomes, holding the reproduction STRUCTURE (even vs graded) fixed. Sweeps the
#  base fecundity across a WIDE range so the dose-response is obvious: from
#  recruitment-limited (0.02) up to recruitment-dominated (0.80).
#
#  Uses only reproduction-ON models (even + graded per network, plus the Linear
#  competition-fecundity trade-off pair). Reproduction OFF is the Classic_* baseline,
#  already in the main 22-model run, so it is not repeated here.
#
#  Run length 2500 (not 1000): recruitment adds slow dynamics, so the longer runs let
#  the metrics reach equilibrium (verify with the time-convergence diagnostic).
#
#  RUN:  Rscript Current_Working_Model/Reproduction_Config_Checks/Run_Reproduction_Fecundity.r
#  Output: Current_Working_Model/Results/Reproduction_Fecundity_<timestamp>_results.{rds,csv}
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
  experiment    = paste0("Reproduction_Fecundity_", format(Sys.time(), "%Y%m%d_%H%M%S")),
  combinations  = c("RPS_reproductionEven",     "RPS_reproduction",
                    "Linear_reproductionEven",  "Linear_reproductionNormal", "Linear_reproductionOpposite",
                    "Neutral_reproductionEven", "Neutral_reproduction",
                    "Random_reproductionEven",  "Random_reproduction"),
  n_species     = c(3, 7),
  sim_length    = 2500,                                # longer: recruitment dynamics are slow
  replicates    = 30,
  reef          = 50,
  individuals   = 3,
  biases        = 0.9,
  intraspecific = 0.5,
  repro_base    = c(0.02, 0.08, 0.20, 0.80),          # <-- the axis under test (4 levels, wide range; 0.08 = default)
  growth        = 3,
  disturbances  = c("off", "on"),
  dist_freq     = "often",
  dist_size     = "random",
  out_dir       = "Current_Working_Model/Results",
  checkpoint_every = 10,
  base_seed     = 1000,
  to_master     = FALSE
)

source("Current_Working_Model/Simulation_testing.r")
cat("\nReproduction fecundity dose-response run complete.\n")
