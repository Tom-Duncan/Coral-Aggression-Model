#  Run_BackgroundRecruitment.r   (metacommunity / rescue-effect experiment)
#  Sweeps the per-species external-larva chance: 0 (CLOSED community, the control),
#  RUN:  Rscript Model/Extended_Sensitivity_Checks/Run_BackgroundRecruitment.r
#  Output: Model/Results/BackgroundRecruitment_<timestamp>_results.{rds,csv}

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
  experiment    = paste0("BackgroundRecruitment_", format(Sys.time(), "%Y%m%d_%H%M%S")),
  combinations  = c("Classic_RPS", "Classic_Linear", "Classic_Neutral", "Classic_Random"),
  n_species     = c(3, 7),
  sim_length    = 2500,
  replicates    = 30,
  reef          = 50,
  individuals   = 3,
  biases        = 0.9,
  intraspecific = 0.5,
  background_chance   = c(0, 0.5, 1.0),          # <-- 0 = closed (control) / moderate / max external supply
  background_interval = 15,                       # external larvae every 15 steps (more frequent -> stronger, clearer signal)
  disturbances  = c("off", "on"),
  dist_freq     = "often",
  dist_size     = "random",
  out_dir       = "Model/Results",
  checkpoint_every = 10,
  base_seed     = 1000,
  to_master     = FALSE
)

source("Model/Simulation_testing.r")
cat("\nBackground-recruitment (rescue-effect) run complete.\n")
