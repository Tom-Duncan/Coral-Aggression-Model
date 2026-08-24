# Runs many simulations of the spatial model non-interactively and saves the
# results. reefSetUp() is interactive, so this builds everything programmatically
# (interaction matrix, traits, placement), runs runSimulation() across a grid of
# combinations, and writes one tidy results table.
#
# Usually driven by Run_Models.r (which pre-defines SIM_CONFIG). Directly:
#   Rscript Model/Simulation_testing.r   or   source(...) in a session.
# Edit SIM_CONFIG to choose the sweep; every vector field is crossed into a grid.

# Set working directory to the project root (folder containing "Model") unless a
# caller already did. Uses the script location under Rscript, else walks up.
if (!file.exists("Model/Saving_Functions.r")) {
  .findRoot <- function() {
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
  setwd(.findRoot())
}

# Source the model
.model_files <- c(
  "Size_Impact_Functions.r", "Intialisation_Functions.r", "Disturbance_Functions.r",
  "Sim_func.r", "Reproduction_Functions.r", "Visual_Functions.r", "Saving_Functions.r"
)
for (f in .model_files) source(file.path("Model", f))

source("Model/Trait_Combinations.r")   # competition-network matrix builders


# Sweep settings (a caller can pre-define SIM_CONFIG to override).
if (!exists("SIM_CONFIG")) SIM_CONFIG <- list(
  experiment    = "sim_test",             # label on every row / output file
  combinations  = TRAIT_COMBINATIONS,     # named trait models
  biases        = c(0.7, 0.9),            # stronger species' win prob
  disturbances  = c("off", "on"),
  n_species     = 3,                      # odd for RPS models; may be a vector
  individuals   = 5,                      # founders per species
  reef          = 30,                     # square reef side (auto-grown if too small)
  sim_length    = 150,
  replicates    = 5,
  growth        = 3,                      # growth trait for every species (1-GROWTH_MAX)
  intraspecific = 0.5,                    # diagonal overgrowth prob
  dist_freq     = "often",
  dist_size     = "random",
  habitat          = "aquarium",          # "aquarium"/"none" (closed) or "reef" (open); may be a vector
  background_chance = 0.10,               # reef only: recruitment chance per species
  supply_mode      = "pulse",             # reef only: "pulse" or "continuous"
  base_seed     = 1000,
  out_dir       = "Model/Results",
  checkpoint_every = 10
)


# Per-species growth: random per species / fixed vector (by rank) / uniform config.
resolveGrowth <- function(cfg, spec, n) {
  if (isTRUE(spec$growth_random)) {
    sample(seq_len(GROWTH_MAX), n, replace = TRUE)
  } else if (!is.null(spec$growth)) {
    rep(spec$growth, length.out = n)
  } else {
    rep(cfg$growth, n)
  }
}


# Build species_traits from a spec + growth vector (reproduction, matrix, size mode).
buildTraits <- function(n, spec, growth_vec) {
  species <- paste0("Sp", seq_len(n))
  size_beta <- if (is.null(spec$size_beta)) defaultSizeBeta(species) else spec$size_beta
  st <- data.frame(species = species, growth = rep(growth_vec, length.out = n),
                   stringsAsFactors = FALSE)
  st <- setReproduction(st, spec, n)
  st <- attachInteraction(st, TRUE, spec$matrix, size_mode = spec$size_mode,
                          size_beta = size_beta)
  st <- attachGrowthSize(st, FALSE, defaultGrowthGamma(species))
  st
}


# Reproduction columns from a spec. NULL repro_chance = off.
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


# Run one replicate; returns a reefSetUp()-shaped result (metrics fields only).
runOneSim <- function(cfg, spec, seed) {
  n         <- cfg$n_species
  indiv_vec <- rep(cfg$individuals, n)
  total     <- sum(indiv_vec)

  grown  <- growReefToFit(cfg$reef, cfg$reef, total)   # grow reef if too small
  reef_x <- grown[1]; reef_y <- grown[2]

  # seed first so random growth and placement are both reproducible
  set.seed(seed)
  growth_vec <- resolveGrowth(cfg, spec, n)
  st         <- buildTraits(n, spec, growth_vec)

  # founder pattern: random / clustered / spread (default random)
  coords <- switch(if (is.null(cfg$placement)) "random" else cfg$placement,
                   random    = getCoordinates(total, reef_x, reef_y, 2),
                   clustered = clusteredCoordinates(total, reef_x, reef_y),
                   spread    = spreadCoordinates(total, reef_x, reef_y),
                   getCoordinates(total, reef_x, reef_y, 2))
  cd     <- createCorals(reef_x, reef_y, n, indiv_vec, coords, st)

  dist_on <- !identical(cfg$dist_setting, "off")
  disturb <- buildDisturbanceConfig(dist_on, cfg$dist_freq, cfg$dist_size, cfg$sim_length)

  invisible(capture.output(   # suppress per-step progress
    states <- runSimulation(cd$reef, cd$corals, cfg$sim_length,
                            cd$colony_species, st, disturb, habitat)
  ))

  list(states = states, reef_x = reef_x, reef_y = reef_y,
       sim_length = cfg$sim_length, disturbance = disturb, species_traits = st)
}

# Run the full grid (replicates parallelised) and save. Returns the long table.
runSimGrid <- function(cfg = SIM_CONFIG) {
  # default the newer axes so configs that omit them behave as before
  if (is.null(cfg$growth))           cfg$growth            <- 3
  if (is.null(cfg$repro_base))       cfg$repro_base       <- REPRO_BASE_CHANCE
  if (is.null(cfg$background_chance)) cfg$background_chance <- 0
  if (is.null(cfg$maturity_age))     cfg$maturity_age      <- MATURITY_AGE
  if (is.null(cfg$size_beta_max))    cfg$size_beta_max     <- SIZECOMP_BETA_MAX
  if (is.null(cfg$placement))        cfg$placement         <- "random"

  grid <- expand.grid(combination = cfg$combinations,
                      n_species   = cfg$n_species,
                      individuals = cfg$individuals,
                      reef        = cfg$reef,
                      bias        = cfg$biases,
                      disturbance = cfg$disturbances,
                      dist_freq   = cfg$dist_freq,
                      dist_size   = cfg$dist_size,
                      intraspecific = cfg$intraspecific,
                      repro_base  = cfg$repro_base,
                      stringsAsFactors = FALSE)

  # frequency/size only matter when disturbance is on; collapse off rows so they
  # don't duplicate the off scenario
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
  on.exit(parallel::stopCluster(cl), add = TRUE)

  total_runs <- nrow(grid) * cfg$replicates
  cat("Running", nrow(grid), "scenarios x", cfg$replicates,
      "replicates =", total_runs, "simulations, parallelised on", ncore,
      "cores; per-worker logs in:", logfile, "\n\n")

  # each grid row's replicates run in parallel; parLapply returns their results
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
    # spec = interaction matrix + size mode; intra sets the diagonal
    spec  <- makeTraitCombination(combo, nsp, bias, intra)

    # fecundity sweep: rescale the reproduction profile to base chance `fec`
    if (!is.null(spec$repro_chance) && fec != REPRO_BASE_CHANCE) {
      spec$repro_chance <- pmin(1, pmax(0, spec$repro_chance * (fec / REPRO_BASE_CHANCE)))
    }

    # size-competition sweep: rescale size_beta to magnitude `sbmax` (0 = off)
    if (!is.null(spec$size_beta) && any(spec$size_beta != 0) && sbmax != SIZECOMP_BETA_MAX) {
      spec$size_beta <- spec$size_beta * (sbmax / SIZECOMP_BETA_MAX)
    }

    # scenario tag encodes the swept values, for labels and filenames
    regime_tag   <- if (dset == "off") "nodist" else sprintf("dist_%s_%s", dfreq, dsize)
    scenario_tag <- sprintf("%s_n%d_i%d_reef%d_b%02d_intra%02d_fec%03d_%s", combo, nsp, indiv, reef,
                            round(bias * 100), round(intra * 100), round(fec * 1000), regime_tag)

    parallel::parLapply(cl, seq_len(cfg$replicates), function(rep) {
      run_i <- (g - 1L) * cfg$replicates + rep
      seed  <- cfg$base_seed + g * 100L + rep
      cat(sprintf("  [%d/%d] %s rep %d (seed %d)\n",
                  run_i, total_runs, scenario_tag, rep, seed))

      run_cfg <- modifyList(cfg, list(dist_setting = dset, n_species = nsp,
                                      individuals = indiv, reef = reef,
                                      dist_freq = dfreq, dist_size = dsize))

      rr <- runOneSim(run_cfg, spec, seed)

      cp <- mainCheckpointMetrics(rr$states, cfg$checkpoint_every)   # same metrics as main save path

      growth_mode <- if (isTRUE(spec$growth_random)) "random"
                     else if (!is.null(spec$growth)) "by_rank" else "uniform"

      # disturbance timesteps (";"-joined) for recovery analysis
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

  # grid_rows: list (per grid row) of lists (per replicate) -> flatten one level
  all_rows <- unlist(grid_rows, recursive = FALSE)

  results <- rbindUnion(all_rows)
  saveResults(results, name = paste0(cfg$experiment, "_results"), dir = cfg$out_dir)
  # also append to the shared master unless the caller opts out
  if (is.null(cfg$to_master) || isTRUE(cfg$to_master)) {
    appendToMaster(results, master = file.path(cfg$out_dir, "master_results.rds"))
  }
  invisible(results)
}


# Mean final species richness per scenario.
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


sim_results <- runSimGrid(SIM_CONFIG)
summariseSims(sim_results)
