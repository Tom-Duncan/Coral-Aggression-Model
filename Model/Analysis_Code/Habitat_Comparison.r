# =============================================================================
#  HABITAT COMPARISON (Aim 2) - natural-reef vs aquarium configurations and the
#  drivers that run them (see the FINAL RUN sections below). The first section is
#  a coverage check that skips scenarios already present in the run index.
# -----------------------------------------------------------------------------
#  Given a PLAN (the scenarios you intend to run), this checks the combined,
#  de-duplicated run index for runs we ALREADY have that are equivalent, so those
#  cells can be skipped (or only topped up to the target replicate count).
#
#  The habitat/supply dimension needs care, because the existing data predates the
#  habitat feature (no `habitat`/`supply_mode` columns). But background recruitment
#  IS the natural-reef mechanism, in its original PULSE form, so we map:
#     background_chance = 0 or NA  -> "aquarium"        (closed, no external supply)
#     background_chance > 0        -> "reef_pulse_<chance>"  (existing reef runs)
#  A planned CONTINUOUS-supply scenario has no equivalent in the old data, so it is
#  always reported as needing to be run. Older files that did not record newer
#  parameters (maturity_age, size_beta_max, background_chance) are normalised to the
#  model DEFAULTS below, so a baseline plan matches the baseline runs and NOT the
#  parameter-sweep runs that varied those settings.
# =============================================================================

# data.table is only needed by the coverage tool below; the simulation run path
# does not use it. Load if available, else skip coverage so the run server (no
# data.table installed) can still source this file.
.have_dt <- requireNamespace("data.table", quietly = TRUE)
if (.have_dt) suppressMessages(library(data.table))

# --- Portable paths: resolve the engine dir (holds Sim_func.r) and the results
#     root (holds Simulation_Results) by searching from the script location or
#     the working dir, so this runs unchanged on the laptop (Windows) and the
#     server (Linux). Original absolute laptop paths kept as the fallback. --------
.scriptDir <- function() {
  a <- commandArgs(FALSE); f <- sub("^--file=", "", grep("^--file=", a, value = TRUE))
  if (length(f)) dirname(normalizePath(f[1], mustWork = FALSE)) else getwd()
}
.ancestorWith <- function(child, starts = unique(c(.scriptDir(), getwd()))) {
  for (start in starts) {
    d <- normalizePath(start, mustWork = FALSE)
    repeat { if (dir.exists(file.path(d, child))) return(d)
             p <- dirname(d); if (p == d) break; d <- p }
  }
  NA_character_
}
.cwm <- .ancestorWith("Model")
MODEL_DIR    <- if (!is.na(.cwm)) file.path(.cwm, "Model") else
  "C:/Users/dell/Documents/Research Project 2/Modelling/Model"
REAL_INT_MAT <- file.path(MODEL_DIR, "Analysis_Code", "Empirical_Matrices.r")
MIA_MATRIX   <- file.path(MODEL_DIR, "Interaction_Matrices", "Mia's_Matrix.r")
.resr <- .ancestorWith("Simulation_Results")
RESULTS_DIR  <- if (!is.na(.resr)) file.path(.resr, "Simulation_Results") else
  "C:/Users/dell/Documents/Research Project 2/Simulation_Results"

index_path <- file.path(RESULTS_DIR, "0_Combined_Master", "run_index.rds")
out_csv    <- file.path(RESULTS_DIR, "0_Combined_Master", "final_run_coverage.csv")
run_index  <- if (.have_dt && file.exists(index_path)) as.data.table(readRDS(index_path)) else NULL

# Model defaults (from the engine) used to normalise NA / unrecorded parameters.a plan that names that regime instead of reading as "NA".
DEFAULTS <- list(maturity_age = 30, repro_base = 0.08, size_beta = 0.5,
                 reef = 50, sim_length = 1000, dist_freq = "often", dist_size = "random")

# --- Map a habitat / supply setting to one comparable key --------------------
# Existing data has neither `habitat` nor `supply_mode`; its supply is fully
# described by background_chance (and was always PULSE). This makes plan and index
# comparable on the same axis.
supplyKey <- function(habitat, supply_mode, background_chance) {
  bc   <- ifelse(is.na(background_chance), 0, as.numeric(background_chance))
  aq   <- (!is.na(habitat) & habitat == "aquarium") | bc <= 0
  mode <- ifelse(is.na(supply_mode), "pulse", as.character(supply_mode))   # legacy data = pulse
  ifelse(aq, "aquarium", paste0("reef_", mode, "_", bc))
}

# --- Reduce a set of scenarios (plan OR index) to a canonical signature ------
# Adds the columns needed if missing, normalises the NA-prone / context-only
# parameters to their defaults, and builds one `sig` string per row. Two runs with
# the same `sig` are the same experiment for coverage purposes.
prepScenarios <- function(dt) {
  dt <- copy(as.data.table(dt))
  need <- c("combination", "n_species", "individuals", "reef_x", "sim_length", "bias",
            "intraspecific", "size_mode", "size_beta_max", "maturity_age",
            "reproduction_on", "repro_base", "disturbance_on", "dist_freq", "dist_size",
            "habitat", "supply_mode", "background_chance")
  for (cc in need) if (!cc %in% names(dt)) dt[[cc]] <- NA

  dt[, reef_k  := fifelse(is.na(reef_x), DEFAULTS$reef, as.numeric(reef_x))]
  dt[, len_k   := fifelse(is.na(sim_length), DEFAULTS$sim_length, as.numeric(sim_length))]
  dt[, mat_k   := fifelse(is.na(maturity_age), DEFAULTS$maturity_age, as.numeric(maturity_age))]
  # size effect: 0 when size_mode is off, else the recorded beta or the default
  dt[, sb_k    := fifelse(size_mode %in% "none", 0,
                          fifelse(is.na(size_beta_max), DEFAULTS$size_beta, as.numeric(size_beta_max)))]
  # reproduction: "off", else its base rate (context-only when off, so collapsed)
  dt[, repro_k := fifelse(reproduction_on %in% TRUE,
                          as.character(fifelse(is.na(repro_base), DEFAULTS$repro_base, as.numeric(repro_base))),
                          "off")]
  # disturbance: "off", else its regime (freq + size, NA -> the baseline default)
  dt[, dist_k  := fifelse(disturbance_on %in% TRUE,
                          paste(fifelse(is.na(dist_freq), DEFAULTS$dist_freq, as.character(dist_freq)),
                                fifelse(is.na(dist_size), DEFAULTS$dist_size, as.character(dist_size))),
                          "off")]
  dt[, supply_k := supplyKey(habitat, supply_mode, background_chance)]

  dt[, sig := paste(combination, n_species, individuals, reef_k, len_k, bias,
                    intraspecific, size_mode, sb_k, mat_k, repro_k, dist_k, supply_k,
                    sep = "|")]
  dt[]
}

# =============================================================================
#  checkCoverage(plan) - the main function
#    plan      : a data.frame, one row per scenario you intend to run. Specify the
#                factors you vary; unset parameters default to the model baseline.
#    min_reps  : target replicate count per scenario (data below this = TOP-UP).
#    Returns the plan annotated with n_existing / to_run / status / matching sources.
# =============================================================================
checkCoverage <- function(plan, index = run_index, min_reps = 30, exclude_dups = TRUE) {
  p <- prepScenarios(plan)
  x <- prepScenarios(index)
  if (exclude_dups && "is_duplicate" %in% names(x)) x <- x[is_duplicate %in% FALSE]

  agg <- x[, .(n_existing = .N,
               sources = paste(sort(unique(basename(source_file))), collapse = "; ")),
           by = sig]
  out <- merge(p, agg, by = "sig", all.x = TRUE, sort = FALSE)
  out[is.na(n_existing), `:=`(n_existing = 0L, sources = "")]
  out[, to_run  := pmax(0L, as.integer(min_reps) - n_existing)]
  out[, status  := fifelse(n_existing >= min_reps, "SKIP - already covered",
                    fifelse(n_existing > 0,         "TOP-UP - partial data",
                                                    "RUN - no prior data"))]
  out[]
}

# Compact, readable view of a coverage result.
coverageReport <- function(cov, show = c("combination", "n_species", "disturbance_on",
                                         "habitat", "supply_mode", "background_chance",
                                         "n_existing", "to_run", "status", "sources")) {
  show <- intersect(show, names(cov))
  as.data.frame(cov[, ..show])
}


# =============================================================================
#  TEMPLATE PLAN - the natural-vs-aquarium final design. EDIT to match your runs.
#  Only the varied factors + the fixed baseline are set here; anything omitted
#  falls back to the model default during matching.
# =============================================================================
if (!is.null(run_index)) {   # coverage check only runs when a master run_index exists
plan <- setDT(expand.grid(
  combination      = c("Classic_RPS", "Classic_Random", "Classic_Linear", "Classic_Neutral"),
  n_species        = c(3, 7),
  disturbance_on   = c(FALSE, TRUE),
  habitat          = c("aquarium", "reef"),
  supply_mode      = c("pulse", "continuous"),   # only meaningful for a reef
  background_chance = 0.5,                        # external supply level for reef runs
  stringsAsFactors = FALSE))

# Fixed baseline parameters for every planned scenario (edit to your final setup)
plan[, `:=`(individuals = 5, reef_x = 50, reef_y = 50, sim_length = 1000,
            bias = 0.9, intraspecific = 0.5, size_mode = "none", size_beta_max = NA,
            maturity_age = 30, reproduction_on = FALSE, repro_base = 0.08,
            dist_freq = "often", dist_size = "random")]

# Aquarium is a closed patch: no supply mode / chance, so collapse those duplicates
plan[habitat == "aquarium", `:=`(supply_mode = NA_character_, background_chance = 0)]
plan <- unique(plan)

# --- Run the check -----------------------------------------------------------
coverage <- checkCoverage(plan, min_reps = 30)

cat("Planned scenarios:", nrow(coverage), "\n")
cat("  already covered (>=30 reps):", sum(coverage$status == "SKIP - already covered"), "\n")
cat("  partial (top-up):           ", sum(coverage$status == "TOP-UP - partial data"), "\n")
cat("  to run (no data):           ", sum(coverage$status == "RUN - no prior data"), "\n\n")
print(coverageReport(coverage), row.names = FALSE)

# Save the coverage plan so the run script (and you) can see exactly what to run.
fwrite(coverageReport(coverage), out_csv)
cat("\nSaved coverage plan to:", out_csv, "\n")
} else {
  cat("[coverage] no run_index found at\n  ", index_path,
      "\n  -> skipping the coverage check (normal on the run server).\n")
}


# =============================================================================
#  MANUAL SCENARIO BUILDER - hand-specify the interaction matrix + every trait
# -----------------------------------------------------------------------------
#  Everything for a final run in one editable place. You type the values directly;
#  buildManualTraits() turns them into a `species_traits` object the engine accepts
#  (it mirrors the engine's attachInteraction / attachGrowthSize / setReproduction),
#  previewTraits() lets you eyeball it, and runManualScenario() runs it end to end.
#
#  Engine facts to respect:
#    * INTERACTION MATRIX: M[i, j] = probability species i (ROW) overgrows species j
#      (COLUMN) when they contest a cell; the DIAGONAL is intraspecific overgrowth.
#      Entries are probabilities in [0, 1].
#    * growth        : integer 1..GROWTH_MAX (= 5); the base growth rate trait.
#    * size_beta     : per-species size effect on COMPETITION (log-ratio mode). 0 =
#                      size-independent; + = larger colonies win more. Needs size_mode
#                      = "logratio" to have any effect.
#    * growth_gamma  : per-species size effect on GROWTH. 0 = size-independent; + =
#                      grows faster when large. Needs growth_size_on = TRUE.
#    * repro_chance  : per-species per-timestep spawn chance once mature. 0 = repro off.
#    * MATURITY AGE  : GLOBAL (one value for ALL species) - the engine has no
#                      per-species maturity, so it is set once in run_cfg below.
# =============================================================================

# MODEL_DIR and REAL_INT_MAT are resolved portably at the top of this file.
GROWTH_MAX_REF <- 5      # matches Sim_func.r GROWTH_MAX

# ---- 1. SPECIES + INTERACTION MATRIX (edit these directly) -------------------
species <- c("Sp1", "Sp2", "Sp3")

# Row beats Column. Diagonal = intraspecific. This example is a clean RPS cycle
# (each species reliably overgrows the next, loses to the previous); edit freely.
M <- rbind(
  Sp1 = c(0.50, 0.85, 0.15),
  Sp2 = c(0.15, 0.50, 0.85),
  Sp3 = c(0.85, 0.15, 0.50))

# ---- 2. PER-SPECIES TRAITS (one row per species; edit the columns) -----------
traits <- data.frame(
  species      = species,
  growth       = c(3,   3,   3),     # 1..5
  size_beta    = c(0,   0,   0),     # size effect on competition (needs size_mode = "logratio")
  growth_gamma = c(0,   0,   0),     # size effect on growth      (needs growth_size_on = TRUE)
  repro_chance = c(0,   0,   0),     # per-step spawn chance when mature; 0 = reproduction off
  repro_ref    = c(5,   5,   5),     # reference size (% cover) reproduction scales around
  maturity_age = c(30,  30,  30),    # PER-SPECIES age of maturity (timesteps); engine now supports this
  stringsAsFactors = FALSE)

# ---- 3. RUN CONFIG (global settings + the run environment) -------------------
run_cfg <- list(
  size_mode        = "none",         # "none" or "logratio" (whether size_beta bites)
  growth_size_on   = FALSE,          # whether growth_gamma bites
  maturity_age     = 30,             # FALLBACK maturity (timesteps) for any species whose per-species value is unset
  reef             = 50,             # reef is reef x reef cells
  individuals      = 5,              # founding colonies per species
  sim_length       = 1000,
  checkpoint_every = 10,
  disturbance_on   = FALSE, dist_freq = "often", dist_size = "random",
  habitat          = "aquarium",     # "aquarium"/"none" (closed) or "reef" (open supply)
  background_chance = 0.10, supply_mode = "pulse")   # reef only

# --- Validate the hand-entered inputs and assemble species_traits -------------
# IMPORTANT: the engine (createCorals) labels colonies "Sp1..SpN" by ORDER and
# looks up traits + the interaction matrix by that label, so species_traits$species
# and the matrix dimnames MUST be "Sp1..SpN". We therefore relabel the input
# species/matrix to Sp1..SpN (preserving order) and store the real names as
# attr(st, "species_map") = c(Sp1 = "<real name 1>", ...) for interpretation.
buildManualTraits <- function(species, M, traits, cfg) {
  stopifnot(length(species) == nrow(M), nrow(M) == ncol(M), nrow(traits) == length(species))
  if (!all(traits$species == species)) stop("traits$species must match `species`, in order.")
  if (any(M < 0 | M > 1, na.rm = TRUE)) stop("interaction entries must be probabilities in [0, 1].")
  if (any(traits$growth < 1 | traits$growth > GROWTH_MAX_REF))
    stop("growth must be an integer 1..", GROWTH_MAX_REF, ".")
  if (!cfg$size_mode %in% c("none", "logratio")) stop("size_mode must be 'none' or 'logratio'.")

  labels <- paste0("Sp", seq_along(species))   # engine's internal colony labels
  dimnames(M) <- list(labels, labels)          # matrix keyed to those labels

  # per-species maturity: from the traits table if given, else the run_cfg fallback
  mat <- if (!is.null(traits$maturity_age)) as.integer(traits$maturity_age)
         else rep(as.integer(cfg$maturity_age), length(species))
  st <- data.frame(
    species           = labels,
    growth            = as.integer(traits$growth),
    maturity_age      = mat,                               # read by the engine's maturityAge()
    reproduction      = traits$repro_chance > 0,
    repro_chance      = pmin(1, pmax(0, traits$repro_chance)),
    repro_ref_percent = ifelse(traits$repro_chance > 0, traits$repro_ref, NA_real_),
    stringsAsFactors  = FALSE)
  # attributes exactly as the engine's attach* functions store them (Sp-labelled)
  attr(st, "interaction") <- list(
    enabled = TRUE, matrix = M, size_mode = cfg$size_mode,
    size_beta = setNames(as.numeric(traits$size_beta), labels))
  attr(st, "growth_size") <- list(
    enabled = isTRUE(cfg$growth_size_on),
    gamma   = setNames(as.numeric(traits$growth_gamma), labels))
  attr(st, "species_map") <- setNames(species, labels)     # Sp<i> -> real name
  st
}

# --- Human-readable preview (matrix + traits + the switches) ------------------
previewTraits <- function(st, cfg) {
  it <- attr(st, "interaction"); gs <- attr(st, "growth_size")
  cat("Interaction matrix (row overgrows column):\n"); print(round(it$matrix, 2))
  cat("\nPer-species traits:\n")
  print(data.frame(st[c("species", "growth", "maturity_age", "reproduction", "repro_chance")],
                   size_beta = it$size_beta[st$species],
                   growth_gamma = gs$gamma[st$species], row.names = NULL))
  cat(sprintf("\nsize_mode = %s | growth_size = %s | maturity_age = per-species (fallback %d)\n",
              it$size_mode, gs$enabled, as.integer(cfg$maturity_age)))
  cat(sprintf("reef = %dx%d | founders/sp = %d | length = %d | habitat = %s (%s)\n",
              cfg$reef, cfg$reef, cfg$individuals, cfg$sim_length, cfg$habitat, cfg$supply_mode))
  invisible(st)
}

# --- Source the model engine ONCE (the same files the sweep runner sources; NOT
#     Simulation_testing.r, which would auto-launch a sweep) ---------------------
.sourceEngine <- function() {
  mf <- c("Size_Impact_Functions.r", "Intialisation_Functions.r", "Disturbance_Functions.r",
          "Sim_func.r", "Reproduction_Functions.r", "Visual_Functions.r", "Saving_Functions.r")
  for (f in mf) source(file.path(MODEL_DIR, f))
}

# --- Apply run-config overrides to the engine's global constants ---------------
# Called after .sourceEngine(); anything absent from cfg keeps the engine default.
#   maturity_age        -> MATURITY_AGE  (fallback; per-species column overrides it)
#   repro_interval      -> REPRO_INTERVAL (spawning events every N steps; 1 = every step)
#   repro_cooldown      -> REPRO_COOLDOWN (per-colony wait after spawning)
#   dist_size_fractions -> DISTURBANCE_SIZE_FRACTION (named small/medium/large)
#   dist_size_weights   -> DISTURBANCE_SIZE_WEIGHTS ("weighted" size draw probs)
.applyEngineOverrides <- function(cfg) {
  MATURITY_AGE <<- as.integer(cfg$maturity_age)
  if (!is.null(cfg$repro_interval))      REPRO_INTERVAL <<- as.integer(cfg$repro_interval)
  if (!is.null(cfg$repro_cooldown))      REPRO_COOLDOWN <<- as.integer(cfg$repro_cooldown)
  if (!is.null(cfg$dist_size_fractions)) DISTURBANCE_SIZE_FRACTION <<- cfg$dist_size_fractions
  if (!is.null(cfg$dist_size_weights))   DISTURBANCE_SIZE_WEIGHTS  <<- cfg$dist_size_weights
}

# --- Set up the reef + founding colonies and run one replicate (engine sourced) --
# Returns the per-timestep states AND the disturbance config actually used (its
# realised schedule is needed for the output row).
.setupAndRun <- function(st, cfg, seed) {
  n <- nrow(st); indiv <- rep(cfg$individuals, length.out = n)
  # rectangular reefs: reef_x (width) / reef_y (height); square cfg$reef as fallback
  rx <- if (!is.null(cfg$reef_x)) cfg$reef_x else cfg$reef
  ry <- if (!is.null(cfg$reef_y)) cfg$reef_y else cfg$reef
  set.seed(seed)
  # placement: "random" (default) or "spread" (farthest-point sampling - founders
  # placed to maximise the space between colonies, e.g. a stocked aquarium)
  plc <- if (!is.null(cfg$placement)) cfg$placement else "random"
  coords  <- if (identical(plc, "spread")) spreadCoordinates(sum(indiv), rx, ry)
             else                          getCoordinates(sum(indiv), rx, ry, 2)
  cd      <- createCorals(rx, ry, n, indiv, coords, st)
  disturb <- buildDisturbanceConfig(cfg$disturbance_on, cfg$dist_freq, cfg$dist_size, cfg$sim_length)
  bint    <- if (!is.null(cfg$background_interval)) cfg$background_interval else BACKGROUND_REPRO_INTERVAL
  blr     <- if (!is.null(cfg$local_retention)) cfg$local_retention else 0
  habitat <- buildHabitatConfig(cfg$habitat, cfg$background_chance, cfg$supply_mode, bint, blr)
  invisible(capture.output(
    states <- runSimulation(cd$reef, cd$corals, cfg$sim_length, cd$colony_species, st, disturb, habitat)))
  list(states = states, disturbance = disturb)
}

# --- Run ONE manual scenario end-to-end --------------------------------------
# Returns the per-timestep states + the checkpoint metrics the sweeps save.
runManualScenario <- function(st, cfg, seed = 1) {
  .sourceEngine()
  .applyEngineOverrides(cfg)
  run <- .setupAndRun(st, cfg, seed)
  list(states = run$states,
       metrics = mainCheckpointMetrics(run$states, cfg$checkpoint_every), traits = st)
}

# --- Replicate backend --------------------------------------------------------
# Each replicate is fully independent and sets its OWN seed inside .setupAndRun,
# so forking the replicate loop is byte-identical to running it serially. Fork on
# Linux/macOS (the run server); plain serial lapply on Windows (the laptop) or
# where `parallel` is unavailable. Override cores with options(final_sim_cores=N)
# or the FINAL_SIM_CORES env var.
.repApply <- function(x, fn) {
  cores <- getOption("final_sim_cores",
                     max(0L, suppressWarnings(as.integer(Sys.getenv("FINAL_SIM_CORES", "0")))))
  can_fork <- .Platform$OS.type != "windows" && requireNamespace("parallel", quietly = TRUE)
  if (cores <= 0L && can_fork) cores <- max(1L, parallel::detectCores() - 2L)
  if (isTRUE(cores > 1L) && can_fork) parallel::mclapply(x, fn, mc.cores = cores)
  else                                lapply(x, fn)
}

# --- Run a manual scenario across REPLICATES and return the tidy sweep-format ---
# results table (one row per checkpoint per replicate), so a hand-built scenario
# feeds straight into the existing plotting / summary / retention pipeline.
#   label     : names the experiment / scenario column (make it unique per design)
#   n_reps    : replicate runs (independent seeds)
#   base_seed : seeds are base_seed + 1 .. + n_reps (reproducible)
#   save      : also write <label>_results.csv/.rds to out_dir
#   to_master : additionally append to a shared final-runs master
runManualBatch <- function(st, cfg, n_reps = 30, base_seed = 1000, label = "manual_scenario",
                           out_dir = file.path(RESULTS_DIR, "6_Final_Runs"),
                           save = TRUE, to_master = FALSE) {
  .sourceEngine()
  .applyEngineOverrides(cfg)
  cat("Running", n_reps, "replicates of '", label, "'...\n", sep = "")

  rows <- .repApply(seq_len(n_reps), function(rep) {
    seed <- base_seed + rep
    run  <- .setupAndRun(st, cfg, seed)
    cp   <- mainCheckpointMetrics(run$states, cfg$checkpoint_every)
    dsteps <- if (isTRUE(run$disturbance$enabled) && length(run$disturbance$schedule) > 0)
                paste(run$disturbance$schedule, collapse = ";") else ""
    params <- list(
      combination       = label,
      disturbance_on    = cfg$disturbance_on,
      dist_freq         = if (isTRUE(cfg$disturbance_on)) cfg$dist_freq else NA_character_,
      dist_size         = if (isTRUE(cfg$disturbance_on)) cfg$dist_size else NA_character_,
      disturbance_steps = dsteps,
      reproduction_on   = any(st$reproduction),
      n_species         = nrow(st),
      individuals       = cfg$individuals,
      reef_x            = if (!is.null(cfg$reef_x)) cfg$reef_x else cfg$reef,
      reef_y            = if (!is.null(cfg$reef_y)) cfg$reef_y else cfg$reef,
      placement         = if (!is.null(cfg$placement)) cfg$placement else "random",
      sim_length        = cfg$sim_length,
      growth            = paste(st$growth, collapse = ";"),         # per-species vector
      size_mode         = attr(st, "interaction")$size_mode,
      maturity_age      = paste(st$maturity_age, collapse = ";"),   # per-species vector
      habitat           = cfg$habitat,
      supply_mode       = if (identical(cfg$habitat, "reef")) cfg$supply_mode else NA_character_,
      background_chance = if (identical(cfg$habitat, "reef")) cfg$background_chance else 0,
      local_retention   = if (identical(cfg$habitat, "reef") && !is.null(cfg$local_retention))
                            cfg$local_retention else 0)
    assembleRunResults(run_id = sprintf("%s_r%02d", label, rep), experiment = label,
                       seed = seed, replicate = rep, params = params, checkpoints = cp,
                       scenario = label)
  })
  results <- rbindUnion(rows)
  if (isTRUE(save)) {
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    saveResults(results, name = paste0(label, "_results"), dir = out_dir)
    if (isTRUE(to_master))
      appendToMaster(results, master = file.path(out_dir, "final_runs_master_results.rds"))
  }
  invisible(results)
}

# --- Save the exact spec (matrix + traits + config) for the record ------------
saveScenarioSpec <- function(st, cfg, name, out_dir = dirname(out_csv)) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  it <- attr(st, "interaction"); gs <- attr(st, "growth_size")
  write.csv(it$matrix, file.path(out_dir, paste0(name, "_matrix.csv")))
  write.csv(data.frame(st[c("species","growth","maturity_age","reproduction","repro_chance","repro_ref_percent")],
                       size_beta = it$size_beta[st$species], growth_gamma = gs$gamma[st$species],
                       row.names = NULL),
            file.path(out_dir, paste0(name, "_traits.csv")), row.names = FALSE)
  writeLines(vapply(names(cfg), function(k) paste0(k, " = ", paste(cfg[[k]], collapse = ",")), ""),
             file.path(out_dir, paste0(name, "_config.txt")))
  cat("Saved spec:", name, "(_matrix.csv, _traits.csv, _config.txt) to", out_dir, "\n")
}

# --- Build + preview now; run/save are opt-in (un-comment to use) -------------
manual_traits <- buildManualTraits(species, M, traits, run_cfg)
previewTraits(manual_traits, run_cfg)
# res  <- runManualScenario(manual_traits, run_cfg, seed = 1)  # one replicate (states + metrics)
# print(tail(res$metrics))
# saveScenarioSpec(manual_traits, run_cfg, name = "final_RPS_baseline")
#
# Full replicate batch -> tidy sweep-format table (feeds the analysis pipeline):
# out <- runManualBatch(manual_traits, run_cfg, n_reps = 30, base_seed = 1000,
#                       label = "final_RPS_baseline", save = TRUE, to_master = TRUE)


# =============================================================================
#  FINAL RUN CONFIGURATIONS - natural reef vs aquarium
# -----------------------------------------------------------------------------
#  One timestep = 2 WEEKS (26 steps = 1 year). Shared biology in both models:
#   * Growth-size effect REMOVED (growth_size_on = FALSE, gamma = 0 all species).
#   * Competition-size effect ON and EQUAL across species (logratio, beta = +0.35
#     for every species - slightly reduced from 0.5, still a real advantage to
#     larger colonies in contested cells).
#   * Maturity at 78 steps (= 3 years); only mature colonies spawn.
#   * Reproduction = synchronised spawning EVENTS every 13 steps (twice a year),
#     equal rate across all species (every mature colony participates,
#     repro_chance = 1 per event, size-scaled by the engine as usual);
#     cooldown = 13 so the event interval alone paces spawning. On the reef this
#     internal spawning stacks ON TOP of the external background recruitment
#     pulse (same 13-step cycle), so as more colonies mature, total recruitment
#     rises above the constant external supply.
#   * Disturbance ON in both: "often" (2 events / 100 steps ~ 1 per ~2 years),
#     size "weighted" - like natural reefs, LARGE disturbances are RARER than
#     small ones (small 50% / medium 30% / large 20%), and all sizes are
#     slightly larger than the old defaults (15% / 25% / 35% of the reef).
#  The two models differ ONLY in openness, arena geometry and initial placement:
#   * AQUARIUM     - closed: no external larval supply; 250 x 160 rectangle;
#     founders placed by farthest-point sampling to MAXIMISE spacing (a
#     deliberately stocked tank).
#   * NATURAL REEF - open: external recruitment pulse every 13 steps (with the
#     spawning cycle), equal chance per species (10% per event), so locally
#     extinct species can return; 200 x 200 square; founder positions fully
#     RANDOM (natural settlement).
#   Both arenas = 40,000 cells, so cover, carrying capacity and absolute
#   disturbance sizes are identical.
# =============================================================================

FINAL_SHARED <- list(
  size_mode           = "logratio",                 # comp size effect ON
  growth_size_on      = FALSE,                      # growth size effect REMOVED
  maturity_age        = 78,                         # mature at 78 steps = 3 years
  repro_interval      = 13,                         # spawning event every ~6 months
  repro_cooldown      = 13,                         # one event cycle
  individuals         = 5,
  sim_length          = 750,                        # ~ 29 years at 2 wk/step
  checkpoint_every    = 10,
  disturbance_on      = TRUE,
  dist_freq           = "often",                    # 2 events / 100 steps
  dist_size           = "weighted",                 # small common, large rare
  dist_size_fractions = c(small = 0.15, medium = 0.25, large = 0.35),
  dist_size_weights   = c(small = 0.50, medium = 0.30, large = 0.20))

#TOTAL external larval supply on the natural reef is CONSTANT across species
#counts (metacommunity convention: larval flux is a property of the patch and its
#connectivity, not of regional richness; a richer pool partitions the same larval
#rain). Expected recruits per pulse = 0.7 (anchored to the tested n = 7
#parameterisation, 7 x 0.10); the grid sets per-species chance = SUPPLY / n.
FINAL_TOTAL_SUPPLY <- 0.7

#LOCAL RETENTION: on top of the fixed base, each species' arrival chance rises by
#FINAL_LOCAL_RETENTION x (its own MATURE cover as a fraction of the reef) - the
#self-seeding of locally produced larvae. Very slight by design: even a reef
#FULLY covered by mature colonies adds at most +0.05 recruits/pulse (~7% of the
#base). Absent species keep the full base chance, so regional rescue is never
#weakened. (Defined here as reef_cfg references it.)
FINAL_LOCAL_RETENTION <- 0.05

#Both arenas hold EXACTLY 40,000 cells (identical area, carrying capacity and
#absolute disturbance sizes); only the geometry differs. 250x160 is the unique
#integer rectangle of that exact area with both sides <= 300 (aspect 25:16).
#Boundaries remain periodic in both, so shape affects wrap distances only.
aquarium_cfg <- modifyList(FINAL_SHARED, list(
  habitat = "aquarium",                             # closed system
  reef_x = 250, reef_y = 160,                       # rectangular tank, 40,000 cells
  placement = "spread",                             # founders placed to MAXIMISE spacing
  background_chance = 0, supply_mode = NA_character_))

reef_cfg <- modifyList(FINAL_SHARED, list(
  habitat = "reef",                                 # open system
  reef_x = 200, reef_y = 200,                       # square reef, 40,000 cells
  placement = "random",                             # fully random settlement positions
  supply_mode         = "pulse",
  background_interval = 13,                         # supply pulse with each spawning cycle
  background_chance   = 0.10,                       # base chance; grid sets TOTAL/n
  local_retention     = FINAL_LOCAL_RETENTION))     # slight self-seeding boost

#Trait table for any species set. All species-asymmetry comes from the MATRIX
#unless growth is varied. comp-size beta equal (+0.35), growth-size gamma 0,
#reproduction equal (chance 1 per event), maturity equal (78 steps = 3 years).
#  growth: a single value (uniform, default 3 = "growth trait OFF"), or a named
#  per-species vector (keyed by species code, e.g. logan_growth) to turn the
#  morphology-based growth trait ON. Values are aligned to `species` by name.
finalTraits <- function(species, growth = 3) {
  g <- if (length(growth) == 1L) rep(as.integer(growth), length(species))
       else if (!is.null(names(growth))) as.integer(growth[species])
       else as.integer(rep(growth, length.out = length(species)))
  if (any(is.na(g))) stop("growth vector is missing an entry for some species (check names).")
  data.frame(species = species, growth = g, size_beta = 0.35, growth_gamma = 0,
             repro_chance = 1.0, repro_ref = 5, maturity_age = 78, stringsAsFactors = FALSE)
}

# =============================================================================
#  THE FINAL SIMULATION RUN - 750 timesteps, 40 replicates per scenario
# -----------------------------------------------------------------------------
#  Two scenarios (the two configurations above) on the same empirical community
#  (default M_dai7, the complete 7-species Dai matrix). Each writes its own
#  <label>_results.{csv,rds} to 6_Final_Runs and appends to the shared
#  final_runs_master, so the analysis pipeline picks both up together.
#  Run with:  runFinalSimulations()          (opt-in; nothing runs on source)
# =============================================================================

FINAL_REPLICATES <- 40
FINAL_BASE_SEED  <- 5000

runFinalSimulations <- function(M = NULL, community = "dai7", vary_growth = FALSE,
                                n_reps = FINAL_REPLICATES, base_seed = FINAL_BASE_SEED) {
  if (is.null(M)) {
    source(REAL_INT_MAT)
    M <- M_dai7
  }
  sp <- rownames(M)
  # growth trait toggle: ON = morphology-based per-species vector; OFF = uniform 3
  growth <- if (isTRUE(vary_growth)) communityGrowth(community) else 3L
  gtag   <- if (isTRUE(vary_growth)) "_growth" else ""
  st <- buildManualTraits(sp, M, finalTraits(sp, growth), FINAL_SHARED)
  saveScenarioSpec(st, FINAL_SHARED, name = paste0("final_", community, gtag, "_shared_spec"))

  t0 <- proc.time()[3]
  aq  <- runManualBatch(st, aquarium_cfg, n_reps = n_reps, base_seed = base_seed,
                        label = paste0("final_", community, gtag, "_aquarium"),
                        save = TRUE, to_master = TRUE)
  nat <- runManualBatch(st, reef_cfg, n_reps = n_reps, base_seed = base_seed + n_reps,
                        label = paste0("final_", community, gtag, "_reef"),
                        save = TRUE, to_master = TRUE)
  cat(sprintf("\nFINAL RUN COMPLETE: 2 scenarios x %d replicates x %d steps in %.1f min\n",
              n_reps, FINAL_SHARED$sim_length, (proc.time()[3] - t0) / 60))
  invisible(list(aquarium = aq, reef = nat))
}

# --- To launch the final runs for ONE community: ------------------------------
# runFinalSimulations()                                   # Dai-7, 2 x 40 x 750


# =============================================================================
#  RUN BOTH EMPIRICAL COMMUNITIES - Logan (Bermuda) and Heron (GBR)
# -----------------------------------------------------------------------------
#  Registry of the empirical communities available for the final experiment.
#  Each name maps to a matrix object defined in Empirical_Matrices.r (10 species each,
#  richness-matched). Extend with dai7 / heron / logan10 / logan17 as needed.
# =============================================================================
FINAL_COMMUNITIES <- c(
  logan10 = "M_logan10", logan7 = "M_logan7", logan3 = "M_logan3",
  heron   = "M_heron10", heron7 = "M_heron7", heron3 = "M_heron3",
  mia10   = "M_mia",     mia7   = "M_mia7",   mia3   = "M_mia3")

#Run one or more empirical communities through both scenarios (aquarium + reef).
#  communities : names from FINAL_COMMUNITIES (default: both Logan and Heron)
#  growth      : "off", "on", or c("off","on") - runs the growth toggle each way
#  Each (community x growth) gets its own 1000-seed block so no replicate seed is
#  reused; results + master are written per scenario, labelled by community.
runFinalCommunities <- function(communities = c("logan10", "heron"),
                                growth = "off",
                                n_reps = FINAL_REPLICATES, base_seed = FINAL_BASE_SEED) {
  source(REAL_INT_MAT)                                   # brings in the matrices + growth vectors
  if (file.exists(MIA_MATRIX)) source(MIA_MATRIX)  # Mia's matrices (M_mia/M_mia7/M_mia3) + mia_growth
  growth <- match.arg(growth, c("off", "on"), several.ok = TRUE)
  block <- 0L; out <- list(); t0 <- proc.time()[3]
  for (comm in communities) {
    Mname <- FINAL_COMMUNITIES[[comm]]
    if (is.null(Mname) || !exists(Mname)) stop("unknown community: ", comm)
    M <- get(Mname)
    for (g in growth) {
      vg  <- identical(g, "on")
      key <- paste0(comm, if (vg) "_growth" else "")
      cat(sprintf("\n=== %s | %d species | growth %s ===\n", comm, nrow(M), g))
      out[[key]] <- runFinalSimulations(M = M, community = comm, vary_growth = vg,
                                        n_reps = n_reps, base_seed = base_seed + block * 1000L)
      block <- block + 1L
    }
  }
  cat(sprintf("\nALL COMMUNITIES DONE: %d community-runs (x2 scenarios x %d reps x %d steps) in %.1f min\n",
              block, n_reps, FINAL_SHARED$sim_length, (proc.time()[3] - t0) / 60))
  invisible(out)
}

# --- To launch communities (any registered name): -----------------------------
# runFinalCommunities()                                     # Logan-10 + Heron (default)
# runFinalCommunities(c("logan3","logan7","heron3","heron7"))   # the 3- and 7-species subsets
# runFinalCommunities("heron3", growth = "on")              # one subset, growth ON
# names(FINAL_COMMUNITIES)  # all runnable: logan10/logan7/logan3, heron/heron7/heron3


# =============================================================================
#  FINAL INITIALISATION GRID - species number x founder density x habitat
# -----------------------------------------------------------------------------
#  Crossed initial conditions for the final experiment:
#    n species : 3, 7, 11      (all odd - required for a balanced RPS cycle)
#    founders  : 1, 3, 5 starting colonies per species
#    habitat   : aquarium (closed) vs natural reef (open)
#  -> 18 scenarios per interaction matrix x FINAL_REPLICATES replicates.
#
#  MATRIX 1: a simple RPS-style cyclic matrix. With n species around a circle,
#  each species beats the (n-1)/2 species AHEAD of it (prob = bias) and loses to
#  the (n-1)/2 behind it (1 - bias); diagonal = intraspecific. Perfectly
#  balanced intransitive competition at any odd n.
#  MATRIX 2: to be chosen (slot ready in `final_matrices`).
# =============================================================================

FINAL_SPECIES  <- c(3, 7, 11)
FINAL_FOUNDERS <- c(1, 3, 5)

#Cyclic RPS-style matrix for odd n (bias = the model's standard decided-pair 0.9)
rpsMatrix <- function(n, bias = 0.9, intraspecific = 0.5) {
  if (n %% 2 == 0) stop("RPS cycle needs an odd species count")
  sp <- paste0("Sp", seq_len(n))
  M  <- matrix(1 - bias, n, n, dimnames = list(sp, sp))
  for (i in seq_len(n))
    for (k in seq_len((n - 1) / 2))
      M[i, ((i - 1 + k) %% n) + 1] <- bias        # beats the next (n-1)/2 species
  diag(M) <- intraspecific
  M
}

#Strict LINEAR competitive hierarchy: Sp1 beats every lower-ranked species (prob
#= bias), Sp_n loses to all; a perfectly transitive pecking order. Sp1 is the
#top competitor and Sp_n the bottom. diagonal = intraspecific.
linearMatrix <- function(n, bias = 0.9, intraspecific = 0.5) {
  sp <- paste0("Sp", seq_len(n))
  M  <- matrix(intraspecific, n, n, dimnames = list(sp, sp))
  M[upper.tri(M)] <- bias        # row beats higher-index (lower-ranked) column
  M[lower.tri(M)] <- 1 - bias
  M
}

#A growth vector illustrating the competition-colonization trade-off: the best
#competitor (Sp1) grows slowest and the weakest (Sp_n) fastest, spread over the
#2-4 tier. Use with the growth toggle to show growth's impact against a hierarchy.
linearGrowthTradeoff <- function(n) {
  sp <- paste0("Sp", seq_len(n))
  setNames(as.integer(round(2 + 2 * (seq_len(n) - 1) / (n - 1))), sp)   # rank1->2 ... rankN->4
}

#Matrices for the final grid: each entry is a function(n) -> matrix (so it can
#build all species counts), or a fixed matrix used only at its own size.
final_matrices <- list(
  RPS    = rpsMatrix,
  Linear = linearMatrix
)

#Run the full grid: every matrix x species count x founder density x habitat,
#FINAL_REPLICATES replicates each, saved + appended to the final-runs master.
#Seed blocks of 1000 per scenario keep every replicate's seed unique.
runFinalGrid <- function(matrices = final_matrices, species = FINAL_SPECIES,
                         founders = FINAL_FOUNDERS, vary_growth = FALSE,
                         n_reps = FINAL_REPLICATES, base_seed = FINAL_BASE_SEED) {
  block <- 0L
  t0 <- proc.time()[3]
  gtag <- if (isTRUE(vary_growth)) "_growth" else ""
  for (mname in names(matrices)) {
    mk <- matrices[[mname]]
    for (nsp in species) {
      M <- if (is.function(mk)) mk(nsp)
           else if (nrow(mk) == nsp) mk
           else next                                  # fixed matrix: only its own size
      sp <- rownames(M)
      # growth toggle: ON = community morphology vector; OFF (or RPS) = uniform 3
      growth <- if (isTRUE(vary_growth)) communityGrowth(mname) else 3L
      st <- buildManualTraits(sp, M, finalTraits(sp, growth), FINAL_SHARED)
      for (ind in founders) for (cfg_nm in c("aquarium", "reef")) {
        cfg <- modifyList(get(paste0(cfg_nm, "_cfg")), list(individuals = ind))
        #constant TOTAL supply: divide the fixed larval rain among the n species
        if (cfg_nm == "reef")
          cfg$background_chance <- FINAL_TOTAL_SUPPLY / nsp
        label <- sprintf("final_%s%s_n%d_i%d_%s", mname, gtag, nsp, ind, cfg_nm)
        runManualBatch(st, cfg, n_reps = n_reps,
                       base_seed = base_seed + block * 1000L,
                       label = label, save = TRUE, to_master = TRUE)
        block <- block + 1L
        cat(sprintf("  [%d scenarios done, %.1f min elapsed]\n",
                    block, (proc.time()[3] - t0) / 60))
      }
    }
  }
  cat(sprintf("\nGRID COMPLETE: %d scenarios x %d replicates x %d steps\n",
              block, n_reps, FINAL_SHARED$sim_length))
}

# --- To launch the full initialisation grid: ----------------------------------
# runFinalGrid()          # RPS only: 18 scenarios x 40 reps x 750 steps
