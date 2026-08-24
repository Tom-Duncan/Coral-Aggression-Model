
#  RUN FROM A TERMINAL (project root or anywhere - the path is resolved):
#     "C:\Program Files\R\R-4.5.1\bin\Rscript.exe" ^
#         "Current_Working_Model/Run_Models.r" --replicates=100 --run_name=batch1

#  Each invocation writes its OWN standalone results file,
#     <out_dir>/<run_name>_<timestamp>_results.{rds,csv}
#  It does NOT append to or overwrite any shared master table - every run stands
#  alone. The timestamp keeps repeated runs from overwriting each other; set
#  unique_runs = FALSE (or --unique_runs=FALSE) to reuse a fixed <run_name> file.
# =============================================================================


# --- Locate the project root robustly, so NO absolute path is hard-coded. ----
#  The root is the directory that CONTAINS "Current_Working_Model". Works when
#  launched with Rscript (uses the script's own location) and when sourced in an
#  interactive session (walks up from the working directory).
findProjectRoot <- function() {
  a <- commandArgs(FALSE)
  f <- sub("^--file=", "", grep("^--file=", a, value = TRUE))
  start <- if (length(f)) dirname(normalizePath(f[1], mustWork = FALSE)) else getwd()
  d <- normalizePath(start, mustWork = FALSE)
  repeat {
    if (basename(d) == "Current_Working_Model") return(dirname(d))
    if (dir.exists(file.path(d, "Current_Working_Model"))) return(d)
    p <- dirname(d)
    if (p == d) break            # reached the filesystem root
    d <- p
  }
  getwd()                        # fallback: assume we are already at the root
}
setwd(findProjectRoot())
cat("Project root:", getwd(), "\n")


#DO EDITS HERE
CONFIG <- list(

  # wHICH MODELS TO USE - all 22 combinations across the 4 network types.
  #   RPS models need ODD species counts (3,5,7 below are all odd, so all 22 run).
  models = c(
    # RPS (intransitive cycle, odd n)
    "Classic_RPS", "RPS_sizecompImpact", "RPS_growthImpact",
    "RPS_reproduction", "RPS_reproductionEven",
    # Linear (transitive hierarchy)
    "Classic_Linear", "Linear_GrowthImpact",
    "Linear_SizeImpactNormal", "Linear_SizeImpactReverse",
    "Linear_reproductionNormal", "Linear_reproductionOpposite",
    "Linear_reproductionEven",
    # Neutral (no overgrowth null)
    "Classic_Neutral", "Neutral_reproductionEven", "Neutral_reproduction",
    "Neutral_sizecompImpact", "Neutral_growthImpact",
    # Random (set arbitrary network)
    "Classic_Random", "Random_reproductionEven", "Random_reproduction",
    "Random_sizecompImpact", "Random_growthImpact"
  ),

  # ---- RUN SIZE ------------------------------------------------------------
  species     = c(3, 5, 7), # number of species (ODD for RPS models). A vector such
                            #   as c(3, 5, 7) runs each species count in turn.
  timesteps   = 1000,        # length of each simulation
  replicates  = 30,         # repeat runs (different random seed) per scenario

  # ---- REEF & STARTING COLONIES -------------------------------------------
  reef        = 50,         # reef is reef x reef cells
  founders    = c(2, 3, 5),  # starting colonies per species (swept)

  # ---- COMPETITION / GROWTH ------------------------------------------------
  bias          = 0.9,      # win chance of the stronger species in a decided pair (0.5-1)
  intraspecific = 0.5,      # same-species overgrowth probability (matrix diagonal)
  growth        = 3,        # base growth trait 1-5 (models with their own growth rule override this)

  # ---- DISTURBANCE ---------------------------------------------------------
  disturbance = c("off", "on"),   # "off", "on", or c("off","on") to run and compare both
  dist_freq   = "often",    # events per 100 steps: "rarely"(1)  "often"(2)  "very_often"(4)
  dist_size   = "random",   # patch size: "small"(10%)  "medium"(20%)  "large"(30%)  "random"

  # ---- HABITAT (open reef vs closed aquarium) ------------------------------
  habitat           = "aquarium", # "aquarium"/"none" (closed patch) or "reef" (open: external larval supply); c("aquarium","reef") sweeps both
  background_chance = 0.10,        # reef only: external recruitment chance per species (equal for all species)
  supply_mode       = "pulse",     # reef only: "pulse" (periodic event) or "continuous" (stochastic trickle, matched supply)

  # ---- OUTPUT --------------------------------------------------------------
  run_name    = "AllModels_sweep",               # labels the output files
  out_dir     = "Current_Working_Model/Results", # where results are saved
  save_every  = 10,                              # record metrics every N timesteps
  seed        = 1000,                            # base RNG seed (reproducible)
  unique_runs = TRUE                             # timestamp each run so batches accumulate
)



# --- Apply any --key=value command-line overrides onto CONFIG ---------------
#  Coerces each value to the type of the existing CONFIG entry; a comma splits a
#  value into a vector (e.g. --models=A,B  or  --species=3,5,7).
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


# --- Unique experiment label so each run writes its own separate results file
#  (rather than overwriting a previous one). A fixed label (unique_runs = FALSE)
#  reuses the same <run_name>_results file, letting a re-run overwrite itself.
experiment_label <- if (isTRUE(CONFIG$unique_runs)) {
  sprintf("%s_%s", CONFIG$run_name, format(Sys.time(), "%Y%m%d_%H%M%S"))
} else {
  CONFIG$run_name
}
cat("Experiment label:", experiment_label, "\n\n")


# --- Map the friendly CONFIG onto the run engine's settings and run ---------
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
  to_master        = FALSE          # each run stands alone; no shared master table
)

# Sources the model, the trait-combination builders and the run engine, then runs
# SIM_CONFIG and saves the results table. Afterwards `sim_results` holds the table.
source("Current_Working_Model/Simulation_testing.r")

cat("\nDone. Results saved to",
    file.path(CONFIG$out_dir, paste0(experiment_label, "_results.rds")), "\n")
