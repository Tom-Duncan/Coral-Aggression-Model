#  Run_Growth_Rate.r   (Growth-rate experiment: uniform growth x disturbance)
#  here we sweep the extremes + midpoint: SLOW (1), MEDIUM (3, the default), FAST (5).
#  RUN:  Rscript Model/Growth_Config_Checks/Run_Growth_Rate.r
#  Output: Model/Results/Growth_Rate_<timestamp>_results.{rds,csv}

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
  experiment    = paste0("Growth_Rate_", format(Sys.time(), "%Y%m%d_%H%M%S")),
  combinations  = c("Classic_RPS", "Classic_Linear", "Classic_Neutral", "Classic_Random"),
  n_species     = c(3, 7),
  sim_length    = 2500,
  replicates    = 30,
  reef          = 50,
  individuals   = 3,
  biases        = 0.9,
  intraspecific = 0.5,
  repro_base    = 0.08,
  growth        = c(1, 3, 5),                    # <-- the axis under test (slow / medium / fast)
  disturbances  = c("off", "on"),                # off = time-rescaling control; on = the real test
  dist_freq     = "often",
  dist_size     = "random",
  out_dir       = "Model/Results",
  checkpoint_every = 10,
  base_seed     = 1000,
  to_master     = FALSE
)

source("Model/Simulation_testing.r")
cat("\nGrowth-rate sensitivity run complete.\n")
