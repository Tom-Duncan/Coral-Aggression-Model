#  Run_Reproduction_DisturbanceRegime.r  (Reproduction experiment B)
#  RUN:  Rscript Model/Reproduction_Config_Checks/Run_Reproduction_DisturbanceRegime.r
#  Output: Model/Results/Reproduction_DisturbanceRegime_<timestamp>_results.{rds,csv}

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
  out_dir       = "Model/Results",
  checkpoint_every = 10,
  base_seed     = 1000,
  to_master     = FALSE
)

source("Model/Simulation_testing.r")
cat("\nReproduction x disturbance-regime run complete.\n")
