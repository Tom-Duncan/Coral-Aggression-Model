# Run from a terminal (path resolved automatically):
#   Rscript Model/Run_Models.r --replicates=100 --run_name=batch1
# Each invocation writes its own standalone results file,
#   <out_dir>/<run_name>_<timestamp>_results.{rds,csv}
# It never appends to a shared master. Set unique_runs = FALSE to reuse a fixed file.


# Locate the project root (the directory containing "Model") - no hard-coded paths.
# Works under Rscript (script location) and when sourced interactively (walks up).
findProjectRoot <- function() {
  a <- commandArgs(FALSE)
  f <- sub("^--file=", "", grep("^--file=", a, value = TRUE))
  start <- if (length(f)) dirname(normalizePath(f[1], mustWork = FALSE)) else getwd()
  d <- normalizePath(start, mustWork = FALSE)
  repeat {
    if (basename(d) == "Model") return(dirname(d))
    if (dir.exists(file.path(d, "Model"))) return(d)
    p <- dirname(d)
    if (p == d) break
    d <- p
  }
  getwd()
}
setwd(findProjectRoot())
cat("Project root:", getwd(), "\n")


# ---- Edit settings here -----------------------------------------------------
CONFIG <- list(

  # all 22 model variants across the 4 network types (RPS needs odd n)
  models = c(
    # RPS (intransitive cycle, odd n)
    "Classic_RPS", "RPS_sizecompImpact", "RPS_growthImpact",
    "RPS_reproduction", "RPS_reproductionEven",
    # Linear (transitive hierarchy)
    "Classic_Linear", "Linear_GrowthImpact",
    "Linear_SizeImpactNormal", "Linear_SizeImpactReverse",
    "Linear_reproductionNormal", "Linear_reproductionOpposite",
    "Linear_reproductionEven",
    # Neutral (no overgrowth)
    "Classic_Neutral", "Neutral_reproductionEven", "Neutral_reproduction",
    "Neutral_sizecompImpact", "Neutral_growthImpact",
    # Random (arbitrary network)
    "Classic_Random", "Random_reproductionEven", "Random_reproduction",
    "Random_sizecompImpact", "Random_growthImpact"
  ),

  species     = c(3, 5, 7), # species counts (odd for RPS); each is run in turn
  timesteps   = 1000,
  replicates  = 30,         # runs per scenario

  reef        = 50,         # reef x reef cells
  founders    = c(2, 3, 5),  # starting colonies per species (swept)

  bias          = 0.9,      # stronger species' win chance in a decided pair (0.5-1)
  intraspecific = 0.5,      # same-species overgrowth (matrix diagonal)
  growth        = 3,        # base growth 1-5 (growth-rule models override)

  disturbance = c("off", "on"),   # "off", "on", or both
  dist_freq   = "often",    # per 100 steps: rarely(1) often(2) very_often(4)
  dist_size   = "random",   # small(10%) medium(20%) large(30%) random

  habitat           = "aquarium", # "aquarium"/"none" (closed) or "reef" (open); sweep with c(...)
  background_chance = 0.10,        # reef: external recruitment chance per species
  supply_mode       = "pulse",     # reef: "pulse" or "continuous"

  run_name    = "AllModels_sweep",
  out_dir     = "Model/Results",
  save_every  = 10,         # record metrics every N steps
  seed        = 1000,
  unique_runs = TRUE        # timestamp each run so batches accumulate
)



# Apply --key=value overrides, coercing to the existing entry's type (comma = vector).
.applyOverride <- function(cfg, key, val) {
  old <- cfg[[key]]
  cfg[[key]] <-
    if (is.null(old))         strsplit(val, ",")[[1]]
    else if (is.logical(old)) as.logical(val)
    else if (is.numeric(old)) as.numeric(strsplit(val, ",")[[1]])
    else                      strsplit(val, ",")[[1]]
  cfg
}
for (a in commandArgs(trailingOnly = TRUE)) {
  m <- regmatches(a, regexec("^--([^=]+)=(.*)$", a))[[1]]
  if (length(m) != 3) { warning("Ignoring unrecognised argument: ", a); next }
  key <- m[2]; val <- m[3]
  if (!key %in% names(CONFIG)) { warning("Ignoring unknown --", key); next }
  CONFIG <- .applyOverride(CONFIG, key, val)
  cat("Override:", key, "=", paste(CONFIG[[key]], collapse = ","), "\n")
}


# Unique label so each run writes its own file (unique_runs = FALSE reuses one).
experiment_label <- if (isTRUE(CONFIG$unique_runs)) {
  sprintf("%s_%s", CONFIG$run_name, format(Sys.time(), "%Y%m%d_%H%M%S"))
} else {
  CONFIG$run_name
}
cat("Experiment label:", experiment_label, "\n\n")


# Map CONFIG onto the run engine's settings.
SIM_CONFIG <- list(
  experiment       = experiment_label,
  combinations     = CONFIG$models,
  n_species        = CONFIG$species,
  sim_length       = CONFIG$timesteps,
  replicates       = CONFIG$replicates,
  reef             = CONFIG$reef,
  individuals      = CONFIG$founders,
  biases           = CONFIG$bias,
  intraspecific    = CONFIG$intraspecific,
  growth           = CONFIG$growth,
  disturbances     = CONFIG$disturbance,
  dist_freq        = CONFIG$dist_freq,
  dist_size        = CONFIG$dist_size,
  habitat          = CONFIG$habitat,
  background_chance = CONFIG$background_chance,
  supply_mode      = CONFIG$supply_mode,
  out_dir          = CONFIG$out_dir,
  checkpoint_every = CONFIG$save_every,
  base_seed        = CONFIG$seed,
  to_master        = FALSE          # each run stands alone
)

# Sources the engine + builders + run engine, runs SIM_CONFIG, saves the table.
source("Model/Simulation_testing.r")

cat("\nDone. Results saved to",
    file.path(CONFIG$out_dir, paste0(experiment_label, "_results.rds")), "\n")
