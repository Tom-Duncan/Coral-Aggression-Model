# =============================================================================
#  Simulation testing.r  -  automatically run many simulations of the CURRENT
#                           spatial model and save all their results.
# -----------------------------------------------------------------------------
#  reefSetUp() is INTERACTIVE (it uses readline), so it cannot drive a sweep. This
#  script builds every piece programmatically - the interaction matrix, the species
#  traits, the colony placement - runs runSimulation() directly, collects the same
#  checkpoint metrics the main model saves, and writes ONE tidy results table (plus a
#  growing master) through Results_Table.r, exactly like the earlier RPS exploration.
#
#  HOW TO RUN
#    Usually driven by Run_Models.r (which pre-defines SIM_CONFIG). Directly:
#      "C:/Program Files/R/R-4.5.1/bin/Rscript.exe" "Current_Working_Model/Simulation_testing.r"
#    R session:
#      source("Current_Working_Model/Simulation_testing.r")
#
#  Edit the SIM_CONFIG block to choose which combinations to sweep. Every vector
#  field is crossed into a full grid, and each combination runs for `replicates`
#  stochastic replicates. Results land in
#      Current_Working_Model/Results/<experiment>_results.{rds,csv}
# =============================================================================

# Set the working directory to the project root (the folder that CONTAINS
# "Current_Working_Model") unless a caller (e.g. Run_Models.r) already did. No
# absolute path is hard-coded, so this runs on any machine / over a VPN. When
# launched with Rscript the script's own location is used; otherwise we walk up
# from the current working directory.
if (!file.exists("Current_Working_Model/Saving_Functions.r")) {
  .findRoot <- function() {
    a <- commandArgs(FALSE)
    f <- sub("^--file=", "", grep("^--file=", a, value = TRUE))
    start <- if (length(f)) dirname(normalizePath(f[1], mustWork = FALSE)) else getwd()
    d <- normalizePath(start, mustWork = FALSE)
    repeat {
      if (basename(d) == "Current_Working_Model") return(dirname(d))
      if (dir.exists(file.path(d, "Current_Working_Model"))) return(d)
      p <- dirname(d)
      if (p == d) break
      d <- p
    }
    getwd()
  }
  setwd(.findRoot())
}

# ---- Source the model (the consolidated files Master_Script sources) ----
.model_files <- c(
  "Size_Impact_Functions.r", "Intialisation_Functions.r", "Disturbance_Functions.r",
  "Sim_func.r", "Reproduction_Functions.r", "Visual_Functions.r", "Saving_Functions.r"
)
for (f in .model_files) source(file.path("Current_Working_Model", f))

# ---- Trait-combination builders (the competition-network matrices) ----
source("Current_Working_Model/Trait_Combinations.r")


# =============================================================================
#  CONFIG  -  edit this to control the sweep. (A caller can pre-define SIM_CONFIG
#  before sourcing to override it, e.g. for a quick test.)
# =============================================================================
if (!exists("SIM_CONFIG")) SIM_CONFIG <- list(
  experiment    = "sim_test",             # label stored on every row / output file
  combinations  = TRAIT_COMBINATIONS,     # named trait models (see Trait_Combinations.r)
  biases        = c(0.7, 0.9),            # win prob of the stronger species in a pair
  disturbances  = c("off", "on"),         # "off", or "on" (periodic mass die-off)
  n_species     = 3,                      # ODD for RPS-based models (3-15); may be a vector to sweep
  individuals   = 5,                      # founding colonies per species
  reef          = 30,                     # square reef side (auto-grown if too small)
  sim_length    = 150,                    # timesteps per run
  replicates    = 5,                      # stochastic replicates per combination
  growth        = 3,                      # growth trait for every species (1-GROWTH_MAX)
  intraspecific = 0.5,                    # diagonal (same-species) overgrowth prob
  dist_freq     = "often",                # disturbance frequency when "on"
  dist_size     = "random",               # disturbance size when "on"
  habitat          = "aquarium",          # "aquarium"/"none" (closed) or "reef" (open: external larval supply); may be a vector to sweep
  background_chance = 0.10,               # reef only: background recruitment chance per species (per pulse)
  supply_mode      = "pulse",             # reef only: "pulse" (periodic) or "continuous" (stochastic trickle, matched supply)
  base_seed     = 1000,                   # seeds are derived from this (reproducible)
  out_dir       = "Current_Working_Model/Results",
  checkpoint_every = 10                   # metrics recorded every N steps (+ start/end)
)


# =============================================================================
#  Resolve per-species growth from a trait-combination spec:
#    growth_random = TRUE -> a random growth trait per species (drawn under the run
#                            seed, so it varies across replicates but is reproducible)
#    spec$growth (vector) -> a fixed per-species growth (e.g. graded by rank)
#    otherwise            -> the uniform config growth for every species
# =============================================================================
resolveGrowth <- function(cfg, spec, n) {
  if (isTRUE(spec$growth_random)) {
    sample(seq_len(GROWTH_MAX), n, replace = TRUE)
  } else if (!is.null(spec$growth)) {
    rep(spec$growth, length.out = n)
  } else {
    rep(cfg$growth, n)
  }
}


# =============================================================================
#  Build species_traits directly (no readline) from a trait-combination spec and a
#  resolved per-species growth vector: reproduction off, the spec's interaction matrix
#  + size mode + per-species size_beta attached, growth size effect off. A NULL
#  size_beta defaults to no size effect.
# =============================================================================
buildTraits <- function(n, spec, growth_vec) {
  species <- paste0("Sp", seq_len(n))
  size_beta <- if (is.null(spec$size_beta)) defaultSizeBeta(species) else spec$size_beta
  st <- data.frame(species = species, growth = rep(growth_vec, length.out = n),
                   stringsAsFactors = FALSE)
  st <- setReproduction(st, spec, n)   # reproduction traits from the spec (or off)
  st <- attachInteraction(st, TRUE, spec$matrix, size_mode = spec$size_mode,
                          size_beta = size_beta)
  st <- attachGrowthSize(st, FALSE, defaultGrowthGamma(species))
  st
}


# =============================================================================
#  Set the reproduction columns from a spec. A NULL spec$repro_chance means the
#  feature is OFF (matches the model's "no reproduction" defaults); a numeric vector
#  switches it ON with a per-species base spawn chance and a shared reference size.
# =============================================================================
setReproduction <- function(st, spec, n) {
  if (is.null(spec$repro_chance)) {
    st$reproduction      <- rep(FALSE, n)
    st$repro_chance      <- rep(0, n)
    st$repro_ref_percent <- rep(NA_real_, n)
  } else {
    chance <- pmin(1, pmax(0, rep(spec$repro_chance, length.out = n)))
    ref    <- if (is.null(spec$repro_ref)) 5 else spec$repro_ref
    st$reproduction      <- chance > 0
    st$repro_chance      <- chance
    st$repro_ref_percent <- rep(ref, length.out = n)
  }
  st
}


# =============================================================================
#  Run ONE replicate and return a result shaped like reefSetUp()'s output (only the
#  fields the metrics need). Simulation console output is suppressed for the sweep.
# =============================================================================
runOneSim <- function(cfg, spec, seed) {
  n         <- cfg$n_species
  indiv_vec <- rep(cfg$individuals, n)
  total     <- sum(indiv_vec)

  #Grow the reef if it is too small to place all colonies (avoids the crowding error)
  grown  <- growReefToFit(cfg$reef, cfg$reef, total)
  reef_x <- grown[1]; reef_y <- grown[2]

  #Seed FIRST so any random growth (RPS_growthImpact) and the colony placement are
  #both reproducible from this run's seed.
  set.seed(seed)
  growth_vec <- resolveGrowth(cfg, spec, n)
  st         <- buildTraits(n, spec, growth_vec)

  #Founder spatial pattern: random (uniform), clustered (one aggregated patch) or spread
  #(over-dispersed). Defaults to random, which is byte-identical to the previous behaviour
  #(same getCoordinates call, same RNG draws), so existing runs are unchanged.
  coords <- switch(if (is.null(cfg$placement)) "random" else cfg$placement,
                   random    = getCoordinates(total, reef_x, reef_y, 2),
                   clustered = clusteredCoordinates(total, reef_x, reef_y),
                   spread    = spreadCoordinates(total, reef_x, reef_y),
                   getCoordinates(total, reef_x, reef_y, 2))
  cd     <- createCorals(reef_x, reef_y, n, indiv_vec, coords, st)

  dist_on <- !identical(cfg$dist_setting, "off")
  disturb <- buildDisturbanceConfig(dist_on, cfg$dist_freq, cfg$dist_size, cfg$sim_length)

  #Suppress the per-step "Step: t" progress printing during the sweep
  invisible(capture.output(
    states <- runSimulation(cd$reef, cd$corals, cfg$sim_length,
                            cd$colony_species, st, disturb, habitat)
  ))

  list(states = states, reef_x = reef_x, reef_y = reef_y,
       sim_length = cfg$sim_length, disturbance = disturb, species_traits = st)
}

# =============================================================================
#  Run the full grid and save. Returns the combined long results table.
# =============================================================================
runSimGrid <- function(cfg = SIM_CONFIG) {
  #These axes default to the model constants, so configs that don't set them behave
  #exactly as before (single value -> no extra grid rows, no change, seeds intact).
  if (is.null(cfg$growth))           cfg$growth            <- 3                     # uniform base growth trait
  if (is.null(cfg$repro_base))       cfg$repro_base       <- REPRO_BASE_CHANCE
  if (is.null(cfg$background_chance)) cfg$background_chance <- 0                    # external larvae OFF
  if (is.null(cfg$maturity_age))     cfg$maturity_age      <- MATURITY_AGE         # age at first reproduction
  if (is.null(cfg$size_beta_max))    cfg$size_beta_max     <- SIZECOMP_BETA_MAX    # size-competition strength
  if (is.null(cfg$placement))        cfg$placement         <- "random"            # founder spatial pattern

  grid <- expand.grid(combination = cfg$combinations,
                      n_species   = cfg$n_species,
                      individuals = cfg$individuals,   # founders per species (swept)
                      reef        = cfg$reef,          # square reef side (swept)
                      bias        = cfg$biases,
                      disturbance = cfg$disturbances,
                      dist_freq   = cfg$dist_freq,     # disturbance frequency (swept)
                      dist_size   = cfg$dist_size,     # disturbance patch size (swept)
                      intraspecific = cfg$intraspecific,  # same-species (diagonal) overgrowth prob (swept)
                      repro_base  = cfg$repro_base,    # base fecundity / spawn chance (swept)
                      stringsAsFactors = FALSE)

  #Disturbance frequency/size only matter when disturbance is ON. For OFF rows collapse
  #them to a single representative value so the (unused) regime combinations don't
  #duplicate the off scenario. With single-value dist_freq/dist_size (the usual case)
  #this is a no-op, so existing configs - and their per-row seeds - are unchanged.
  off_rows <- grid$disturbance == "off"
  if (any(off_rows)) {
    grid$dist_freq[off_rows] <- cfg$dist_freq[1]
    grid$dist_size[off_rows] <- cfg$dist_size[1]
    grid <- unique(grid)
    rownames(grid) <- NULL
  }

  library(parallel)
  ncore   <- max(1L, parallel::detectCores() - 1L)
  logfile <- paste0(cfg$experiment, ".log")
  cl      <- parallel::makeCluster(ncore, type = "FORK", outfile = logfile)
  on.exit(parallel::stopCluster(cl), add = TRUE)   # always shut the cluster down

  total_runs <- nrow(grid) * cfg$replicates
  cat("Running", nrow(grid), "scenarios x", cfg$replicates,
      "replicates =", total_runs, "simulations, parallelised on", ncore,
      "cores; per-worker logs in:", logfile, "\n\n")

  #Run every grid row; within each, its `replicates` runs execute in parallel.
  #parLapply RETURNS a list of results - each worker holds a PRIVATE copy of memory and
  #cannot write back to a shared object, so its return value is the only way results
  #come back. We collect those returns (this is exactly what the old code dropped).
  grid_rows <- lapply(seq_len(nrow(grid)), function(g) {
    combo <- grid$combination[g]
    nsp   <- grid$n_species[g]
    indiv <- grid$individuals[g]
    reef  <- grid$reef[g]
    bias  <- grid$bias[g]
    dset  <- grid$disturbance[g]
    dfreq <- grid$dist_freq[g]
    dsize <- grid$dist_size[g]
    intra <- grid$intraspecific[g]
    fec   <- grid$repro_base[g]
    #Build the trait-combination spec (interaction matrix + size mode + size impact).
    #intra sets the matrix DIAGONAL (same-species overgrowth), so sweeping it varies
    #intraspecific competition strength.
    spec  <- makeTraitCombination(combo, nsp, bias, intra)

    #Fecundity-magnitude sweep: proportionally rescale the reproduction profile so the
    #base spawn chance becomes `fec`, preserving the even/graded STRUCTURE while changing
    #MAGNITUDE. No-op for reproduction-off models (repro_chance NULL) and when fec equals
    #the default base, so existing runs are unchanged.
    if (!is.null(spec$repro_chance) && fec != REPRO_BASE_CHANCE) {
      spec$repro_chance <- pmin(1, pmax(0, spec$repro_chance * (fec / REPRO_BASE_CHANCE)))
    }

    #Size-competition strength sweep: proportionally rescale the per-species size_beta so
    #its magnitude becomes `sbmax`, preserving the graded/split STRUCTURE. No-op for models
    #with no size effect (size_beta all 0) and when sbmax equals the default, so existing
    #runs are unchanged. sbmax = 0 switches the size effect fully off.
    if (!is.null(spec$size_beta) && any(spec$size_beta != 0) && sbmax != SIZECOMP_BETA_MAX) {
      spec$size_beta <- spec$size_beta * (sbmax / SIZECOMP_BETA_MAX)
    }

    #Scenario tag: off rows read "nodist"; on rows carry their frequency+size regime so
    #the disturbance regime is visible in the run label and output filenames. bias,
    #intraspecific and base fecundity are encoded so those sweeps are distinguishable.
    regime_tag   <- if (dset == "off") "nodist" else sprintf("dist_%s_%s", dfreq, dsize)
    scenario_tag <- sprintf("%s_n%d_i%d_reef%d_b%02d_intra%02d_fec%03d_%s", combo, nsp, indiv, reef,
                            round(bias * 100), round(intra * 100), round(fec * 1000), regime_tag)

    parallel::parLapply(cl, seq_len(cfg$replicates), function(rep) {
      run_i <- (g - 1L) * cfg$replicates + rep       # deterministic index (no shared counter)
      seed  <- cfg$base_seed + g * 100L + rep
      cat(sprintf("  [%d/%d] %s rep %d (seed %d)\n",
                  run_i, total_runs, scenario_tag, rep, seed))

      #Inject this row's swept values so runOneSim (which reads them as scalars) uses them
      run_cfg <- modifyList(cfg, list(dist_setting = dset, n_species = nsp,
                                      individuals = indiv, reef = reef,
                                      dist_freq = dfreq, dist_size = dsize))

      rr <- runOneSim(run_cfg, spec, seed)

      #Checkpoint metrics (richness, diversity, per-species cover, ...) - the SAME
      #function the main model's save path uses.
      cp <- mainCheckpointMetrics(rr$states, cfg$checkpoint_every)

      #How growth was assigned for this combination (uniform / random / by rank)
      growth_mode <- if (isTRUE(spec$growth_random)) "random"
                     else if (!is.null(spec$growth)) "by_rank" else "uniform"

      #Per-run parameters recorded on every checkpoint row
      #Timesteps disturbance events landed on (";"-joined), for recovery analysis
      dist_steps <- if (isTRUE(rr$disturbance$enabled) &&
                        length(rr$disturbance$schedule) > 0) {
        paste(rr$disturbance$schedule, collapse = ";")
      } else ""

      params <- list(
        combination       = combo,
        bias              = bias,
        disturbance_on    = (dset != "off"),
        dist_freq         = if (dset == "off") NA_character_ else dfreq,
        dist_size         = if (dset == "off") NA_character_ else dsize,
        disturbance_steps = dist_steps,
        reproduction_on   = !is.null(spec$repro_chance),
        n_species         = nsp,
        individuals     = indiv,
        reef_x          = rr$reef_x,
        reef_y          = rr$reef_y,
        sim_length      = cfg$sim_length,
        growth          = grw,
        growth_mode     = growth_mode,
        intraspecific   = intra,
        repro_base      = fec,
        size_mode       = spec$size_mode
      )

      #RETURN the assembled row - do NOT assign into an outer all_rows, which would
      #only touch this worker's private copy and then be discarded.
      assembleRunResults(
        run_id      = sprintf("%s_%s_r%02d", cfg$experiment, scenario_tag, rep),
        experiment  = cfg$experiment,
        seed        = seed,
        replicate   = rep,
        params      = params,
        checkpoints = cp,
        scenario    = scenario_tag
      )
    })
  })

  #grid_rows is a list (one per grid row) of lists (one per replicate); flatten one
  #level to a single flat list of per-run data.frames.
  all_rows <- unlist(grid_rows, recursive = FALSE)

  results <- rbindUnion(all_rows)
  saveResults(results, name = paste0(cfg$experiment, "_results"), dir = cfg$out_dir)
  #Append to the shared master too, UNLESS the caller opts out (to_master = FALSE),
  #in which case this run's <experiment>_results file stands alone. Default stays on
  #so existing callers (e.g. Overnight_run.r) are unchanged.
  if (is.null(cfg$to_master) || isTRUE(cfg$to_master)) {
    appendToMaster(results, master = file.path(cfg$out_dir, "master_results.rds"))
  }
  invisible(results)
}


# =============================================================================
#  A compact end-of-run summary: mean FINAL species richness per scenario.
# =============================================================================
summariseSims <- function(results) {
  fin <- results[results$phase == "final", , drop = FALSE]
  agg <- aggregate(species_richness ~ combination + bias + disturbance_on,
                   data = fin, FUN = function(x) round(mean(x), 2))
  names(agg)[names(agg) == "species_richness"] <- "mean_final_richness"
  agg <- agg[order(agg$combination, agg$bias, agg$disturbance_on), ]
  cat("\n--- Mean final species richness by scenario ---\n")
  print(agg, row.names = FALSE)
  invisible(agg)
}


# =============================================================================
#  RUN
# =============================================================================
sim_results <- runSimGrid(SIM_CONFIG)
summariseSims(sim_results)
