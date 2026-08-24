# =============================================================================
#  Run_InitialPlacement.r   (structural robustness: founder spatial pattern)
# -----------------------------------------------------------------------------
#  Tests whether the results depend on the INITIAL SPATIAL ARRANGEMENT of the founding
#  colonies - the three canonical spatial point patterns:
#    * random    - uniform-random positions (the standard default)
#    * clustered - all founders packed into ONE sub-region (a single settlement patch;
#                  aggregated). Founders suffer CORRELATED mortality: one disturbance
#                  patch can wipe out several at once.
#    * spread    - founders placed as evenly as possible across the reef (over-dispersed;
#                  farthest-point sampling). Disturbance-buffered.
#
#  Robustness question: is an initial-condition difference transient (washes out to the
#  same end state) or persistent (a spatial priority effect -> different end state)?
#  Analyse the ENDPOINT and the TRAJECTORY together. The biggest placement effect is
#  expected in the disturbance-ON arm.
#
#  Run length 1000 (this validates the standard-length main results).
#
#  RUN:  Rscript Model/Structural_Config_Checks/Run_InitialPlacement.r
#  Output: Model/Results/InitialPlacement_<timestamp>_results.{rds,csv}
# =============================================================================

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
  experiment    = paste0("InitialPlacement_", format(Sys.time(), "%Y%m%d_%H%M%S")),
  combinations  = c("Classic_RPS", "Classic_Linear", "Classic_Neutral", "Classic_Random"),
  n_species     = c(3, 7),
  sim_length    = 1000,                          # validates the standard-length main results
  replicates    = 30,
  reef          = 50,
  individuals   = 3,
  biases        = 0.9,
  intraspecific = 0.5,
  placement     = c("random", "clustered", "spread"),   # <-- the axis under test
  disturbances  = c("off", "on"),
  dist_freq     = "often",
  dist_size     = "random",
  out_dir       = "Model/Results",
  checkpoint_every = 10,
  base_seed     = 1000,
  to_master     = FALSE
)

source("Model/Simulation_testing.r")
cat("\nInitial-placement (structural robustness) run complete.\n")
