#  Run_Reproduction_Fecundity.r   (Reproduction experiment A: dose-response)
#  the outcomes, holding the reproduction STRUCTURE (even vs graded) fixed. Sweeps the
#  RUN:  Rscript Model/Reproduction_Config_Checks/Run_Reproduction_Fecundity.r
#  Output: Model/Results/Reproduction_Fecundity_<timestamp>_results.{rds,csv}

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
  out_dir       = "Model/Results",
  checkpoint_every = 10,
  base_seed     = 1000,
  to_master     = FALSE
)

source("Model/Simulation_testing.r")
cat("\nReproduction fecundity dose-response run complete.\n")
