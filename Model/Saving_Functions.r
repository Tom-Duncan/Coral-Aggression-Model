# =============================================================================
#  Results_Table.r  -  a single tidy results table for large simulation sweeps
# =============================================================================
#  MAIN-MODEL copy. Identical storage layer to the one used by the RPS-style
#  exploration, but kept separate so the main model has its own copy and its own
#  output folder. (The exploration keeps its version under RPS_Style_Exploration/.)
#
#
#  SHAPE (tidy / long): one row per (run x checkpoint). Every row carries:
#     run_id      unique id for the simulation run (traceable back to its setup)
#     experiment  label for the model variation/sweep, so many different
#                 explorations can live in one growing master table
#     replicate   stochastic replicate number within a parameter combination
#     seed        RNG seed used for the run's dynamics (reproducibility)
#     <params...> every setup parameter (repeated on each of the run's rows)
#     phase       "initial" (first checkpoint) / "checkpoint" / "final" (last)
#     timestep    the checkpoint's timestep
#     <metrics...> the checkpoint outputs (richness, diversity, cover, ...)
#
#  The INITIAL and FINAL states are simply the phase == "initial" / "final" rows,
#  so they are stored alongside every intermediate checkpoint with no special case.
#
#  Different experiments may record different parameter/metric columns. They are
#  combined with a column-union row-bind (rbindUnion): any column a run does not
#  have is filled with NA. This lets one master table span every model variation.
#
# Examples:
#     RPS : 3 species, 5 indiv, small 30x30, disturbance on, repro off -> 3i5Dr_30x30
#     Main: 5 species, big 60x60, disturbance off, repro on            -> 5dR_60x60


scenarioCode <- function(n_species, disturbance_on, reproduction_on,
                         reef_x, reef_y, individuals = NA) {
  onoff <- function(flag, up) if (isTRUE(flag)) up else tolower(up)
  ind   <- if (!is.na(individuals)) paste0("i", individuals) else ""
  paste0(n_species, ind,
         onoff(disturbance_on, "D"),
         onoff(reproduction_on, "R"),
         "_", reef_x, "x", reef_y)
}


#  Row-bind a list of data.frames that may have DIFFERENT columns.
#  The result has the union of all columns; missing cells become NA.
#  This is what lets many experiments (each with its own parameters) share one
#  master table without a fixed, ever-growing schema.
rbindUnion <- function(dfs) {
  dfs <- Filter(function(d) !is.null(d) && nrow(d) > 0, dfs)
  if (length(dfs) == 0) return(NULL)

  all_cols <- unique(unlist(lapply(dfs, names)))
  dfs2 <- lapply(dfs, function(d) {
    for (m in setdiff(all_cols, names(d))) d[[m]] <- NA
    d[all_cols]                       # reorder to a common column order
  })
  out <- do.call(rbind, dfs2)
  rownames(out) <- NULL
  out
}


#  Assemble ONE run's rows for the master table.
#    run_id      : unique run id (character)
#    experiment  : experiment/sweep label (character)
#    seed        : RNG seed for the run's dynamics
#    replicate   : replicate number within the parameter combination
#    params      : named list or 1-row data.frame of setup parameters
#    checkpoints : data.frame, one row per checkpoint (must contain `timestep`)
#    scenario    : compact setup code (see scenarioCode); same across replicates
#  Returns a long data.frame: run/experiment/scenario/replicate/seed + params +
#  phase + the checkpoint rows, with params/meta recycled across every row.
assembleRunResults <- function(run_id, experiment, seed, replicate, params, checkpoints,
                               scenario = NA) {
  if (is.null(checkpoints) || nrow(checkpoints) == 0) return(NULL)
  n <- nrow(checkpoints)

  meta <- data.frame(
    run_id     = run_id,
    experiment = experiment,
    scenario   = scenario,
    replicate  = replicate,
    seed       = seed,
    stringsAsFactors = FALSE
  )
  param_df <- as.data.frame(params, stringsAsFactors = FALSE)

  #Mark the first / last checkpoint as the initial / final state
  ts    <- checkpoints$timestep
  phase <- ifelse(ts == min(ts), "initial",
           ifelse(ts == max(ts), "final", "checkpoint"))

  out <- cbind(
    meta[rep(1, n), , drop = FALSE],
    param_df[rep(1, n), , drop = FALSE],
    data.frame(phase = phase, stringsAsFactors = FALSE),
    checkpoints
  )
  rownames(out) <- NULL
  out
}


#  Save a results table to disk as BOTH an .rds (typed master, fast to reload)
#  and a .csv (portable, human-readable). Creates `dir` if needed.
#  Returns the .rds path invisibly.
saveResults <- function(results, name, dir = "Results") {
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  rds <- file.path(dir, paste0(name, ".rds"))
  csv <- file.path(dir, paste0(name, ".csv"))
  saveRDS(results, rds)
  write.csv(results, csv, row.names = FALSE)
  cat("Saved", nrow(results), "rows to", rds, "and", csv, "\n")
  invisible(rds)
}


#  Append a results table to the growing MASTER table (one file for all sweeps).
#  By default, rows from an experiment already in the master are REPLACED, so
#  re-running an experiment overwrites its old rows rather than duplicating them
#  (other experiments are left untouched). Keeps an .rds master + .csv mirror.
#  Returns the combined table invisibly.
appendToMaster <- function(results, master = file.path("Results", "master_results.rds"),
                           replace_existing = TRUE) {
  dir.create(dirname(master), showWarnings = FALSE, recursive = TRUE)

  if (file.exists(master)) {
    old <- readRDS(master)
    if (replace_existing && "experiment" %in% names(old)) {
      old <- old[!(old$experiment %in% unique(results$experiment)), , drop = FALSE]
    }
    combined <- rbindUnion(list(old, results))
  } else {
    combined <- results
  }

  saveRDS(combined, master)
  write.csv(combined, sub("\\.rds$", ".csv", master), row.names = FALSE)
  cat("Master table now holds", nrow(combined), "rows across",
      length(unique(combined$experiment)), "experiment(s):", master, "\n")
  invisible(combined)
}


#  Load the master table (or any saved .rds results table).
loadResults <- function(file = file.path("Results", "master_results.rds")) {
  if (!file.exists(file)) {
    stop("No results file at: ", file)
  }
  readRDS(file)
}

# =============================================================================
#  Save_Simulation.r  -  turn a main-model run into the flat checkpoint table
# =============================================================================
#  Records the simulation at checkpoints - the start (t = 0), then every
#  `interval` timesteps, up to and including the end - and stores each checkpoint
#  as a row in the tidy results table (see Results_Table.r). This is the main
#  model's equivalent of the RPS exploration's rpsCheckpointMetrics + save step.
#
#  Requires Results_Table.r to be sourced first (assembleRunResults / saveResults
#  / appendToMaster). Master_Script.r sources both.
# =============================================================================


# --------------------------------------------------------------------
#  Community metrics at each checkpoint of ONE run's `states` list.
#  Samples t = 0, then every `interval` steps, up to the final step (always
#  including start and end). At each checkpoint, per species (columns named
#  <species>_cover), plus richness, Shannon diversity, Pielou evenness, total
#  cover and the live colony count. All cover values are % of the reef.
#  Works for any number of species (derived from the initial state).
#  Returns a tidy data.frame: one row per checkpoint.
# --------------------------------------------------------------------
# --- Colony size classes (cells): small < 20, medium 20-60, large > 60 -------
CSIZE_SMALL_MAX  <- 20    # small  : size  < 20 cells
CSIZE_MEDIUM_MAX <- 60    # medium : 20 <= size < 60 ; large : size >= 60

# Structural complexity: of all adjacent pairs where BOTH cells are occupied, the
# fraction that belong to DIFFERENT colonies (edge/interface density of the occupied
# mosaic). 0 = one uniform patch, ->1 = highly interdigitated. NA if no occupied pairs.
.structuralComplexity <- function(reef) {
  E <- 0L; Ed <- 0L
  if (ncol(reef) > 1) {                       # horizontal neighbours
    a <- reef[, -ncol(reef)]; b <- reef[, -1]
    both <- !is.na(a) & !is.na(b)
    E <- E + sum(both); Ed <- Ed + sum(both & a != b)
  }
  if (nrow(reef) > 1) {                        # vertical neighbours
    a <- reef[-nrow(reef), ]; b <- reef[-1, ]
    both <- !is.na(a) & !is.na(b)
    E <- E + sum(both); Ed <- Ed + sum(both & a != b)
  }
  if (E == 0) NA_real_ else Ed / E
}

# The SPECIES occupying each cell as a character matrix ("" for empty), used for
# turnover (colony ids are volatile through splits, so we compare at species level).
.speciesGrid <- function(state) {
  ids <- vapply(state$corals, function(x) x$id, character(1))
  sp  <- vapply(state$corals, function(x) x$species, character(1))
  out <- setNames(sp, ids)[state$reef]
  out[is.na(out)] <- ""
  matrix(out, nrow = nrow(state$reef), ncol = ncol(state$reef))
}

# Spatial turnover between two states: of the cells occupied at EITHER checkpoint, the
# fraction whose occupying species changed (a cell emptied, filled, or swapped species).
.turnover <- function(prev_state, cur_state) {
  a <- .speciesGrid(prev_state); b <- .speciesGrid(cur_state)
  occ <- (a != "") | (b != "")
  if (!any(occ)) return(NA_real_)
  sum(a[occ] != b[occ]) / sum(occ)
}

# Colony deaths between two states: colonies alive (size > 0) at the previous
# checkpoint that hold 0 cells (or are gone) at the current one.
.deathsBetween <- function(prev_state, cur_state) {
  pid <- vapply(prev_state$corals, function(x) x$id, character(1))
  psz <- vapply(prev_state$corals, function(x) x$size, numeric(1))
  cid <- vapply(cur_state$corals, function(x) x$id, character(1))
  csz <- setNames(vapply(cur_state$corals, function(x) x$size, numeric(1)), cid)
  alive_prev <- pid[psz > 0]
  cur_of <- csz[alive_prev]                    # NA if the colony no longer exists
  sum(is.na(cur_of) | cur_of <= 0)
}

mainCheckpointMetrics <- function(states, interval = 10) {
  sim_length  <- length(states) - 1            # states[[1]] is timestep 0
  checkpoints <- seq(0, sim_length, by = interval)
  if (checkpoints[length(checkpoints)] != sim_length) {
    checkpoints <- c(checkpoints, sim_length)   # always include the final step
  }

  #Master species list from the initial state, so an extinct species still shows
  #as a 0-cover column at every checkpoint
  species_list <- sort(unique(vapply(states[[1]]$corals,
                                     function(x) x$species, character(1))))

  rows <- lapply(seq_along(checkpoints), function(k) {
    t        <- checkpoints[k]
    state_t  <- states[[t + 1]]
    corals_t <- state_t$corals
    reef_t   <- state_t$reef
    #Colony $size is a cell count; convert to % of the reef so cover metrics stay a
    #reef-independent percentage (and comparable to earlier % cover datasets).
    sizes <- vapply(corals_t, function(x) coralCoverPercent(x, reef_t), numeric(1))
    specs <- vapply(corals_t, function(x) x$species, character(1))
    alive <- sizes > 0
    sizes <- sizes[alive]; specs <- specs[alive]

    #Per-species reef cover (abundance, as % of reef) across all its colonies
    cover       <- vapply(species_list, function(sp) sum(sizes[specs == sp]), numeric(1))
    total_cover <- sum(cover)
    present     <- cover > 0
    richness    <- sum(present)

    #Shannon diversity H' over species relative cover, and Pielou's evenness
    if (total_cover > 0) {
      p       <- cover[present] / total_cover
      shannon <- -sum(p * log(p))
    } else {
      shannon <- 0
    }
    evenness <- if (richness > 1) shannon / log(richness) else NA_real_

    #Colony sizes in CELLS (raw) for the size structure metrics
    cell_sizes <- vapply(corals_t, function(x) x$size, numeric(1))
    live       <- cell_sizes[cell_sizes > 0]
    #Founding colonies still alive (founders have no parent_id; kept at size 0 if dead)
    founders_alive <- sum(vapply(corals_t,
      function(x) is.null(x$parent_id) && x$size > 0, logical(1)))

    #Between-checkpoint metrics need the previous checkpoint (NA at the first)
    if (k == 1) {
      turnover <- NA_real_; deaths <- NA_real_
    } else {
      prev_state <- states[[checkpoints[k - 1] + 1]]
      turnover   <- .turnover(prev_state, state_t)
      deaths     <- .deathsBetween(prev_state, state_t)
    }

    row <- data.frame(
      timestep              = t,
      species_richness      = richness,
      shannon               = round(shannon, 4),
      evenness              = round(evenness, 4),
      total_cover           = round(total_cover, 3),
      n_colonies            = length(sizes),
      structural_complexity = round(.structuralComplexity(reef_t), 4),
      turnover              = round(turnover, 4),
      deaths_since_last     = deaths,
      founders_alive        = founders_alive,
      n_small               = sum(live <  CSIZE_SMALL_MAX),
      n_medium              = sum(live >= CSIZE_SMALL_MAX & live < CSIZE_MEDIUM_MAX),
      n_large               = sum(live >= CSIZE_MEDIUM_MAX),
      mean_colony_size      = if (length(live)) round(mean(live), 2) else 0,
      median_colony_size    = if (length(live)) median(live) else 0,
      max_colony_size       = if (length(live)) max(live) else 0,
      stringsAsFactors      = FALSE
    )
    for (sp in species_list) {
      row[[paste0(sp, "_cover")]] <- round(cover[[sp]], 3)
    }
    row
  })

  do.call(rbind, rows)
}


# --------------------------------------------------------------------
#  Was disturbance / reproduction switched on for this replicate?
# --------------------------------------------------------------------
mainDisturbanceOn <- function(rep_result) {
  !is.null(rep_result$disturbance) && isTRUE(rep_result$disturbance$enabled)
}
mainReproductionOn <- function(rep_result) {
  st <- rep_result$species_traits
  !is.null(st) && "reproduction" %in% names(st) && any(st$reproduction, na.rm = TRUE)
}


# --------------------------------------------------------------------
#  Basic setup parameters for one replicate result, stored on every row.
# --------------------------------------------------------------------
mainRunParams <- function(rep_result, interval) {
  init_corals <- rep_result$states[[1]]$corals
  #Timesteps disturbance events occurred on (";"-joined), for recovery-rate analysis
  dist_steps <- if (mainDisturbanceOn(rep_result) &&
                    length(rep_result$disturbance$schedule) > 0) {
    paste(rep_result$disturbance$schedule, collapse = ";")
  } else ""
  data.frame(
    reef_x            = rep_result$reef_x,
    reef_y            = rep_result$reef_y,
    sim_length        = rep_result$sim_length,
    n_species         = length(unique(vapply(init_corals, function(x) x$species, character(1)))),
    n_individuals     = length(init_corals),
    disturbance_on    = mainDisturbanceOn(rep_result),
    disturbance_steps = dist_steps,
    reproduction_on   = mainReproductionOn(rep_result),
    checkpoint_every  = interval,
    stringsAsFactors  = FALSE
  )
}


# --------------------------------------------------------------------
#  Build the flat checkpoint table for a whole simulation (WITHOUT saving).
#  Handles a single-replicate sim (returned as one replicate result) and a
#  multi-replicate sim (a list with $replicates). One run_id per replicate.
#    sim        : the object returned by reefSetUp()
#    experiment : label for this run/variation (goes in the `experiment` column)
#    seed       : optional RNG seed to record (NA if you did not set one)
#  Returns the long data.frame (one row per replicate per checkpoint).
# --------------------------------------------------------------------
collectSimTable <- function(sim, experiment = "main", seed = NA, interval = 10) {
  reps <- if (!is.null(sim$replicates)) sim$replicates else list(sim)

  rows <- lapply(seq_along(reps), function(i) {
    rep_result <- reps[[i]]
    cp     <- mainCheckpointMetrics(rep_result$states, interval)
    params <- mainRunParams(rep_result, interval)
    rep_no <- if (!is.null(rep_result$replicate)) rep_result$replicate else i

    #Compact setup code (reef dimensions at the end; no individuals segment, as
    #the main model's species count is the leading number)
    scenario <- scenarioCode(
      n_species        = params$n_species,
      disturbance_on   = params$disturbance_on,
      reproduction_on  = params$reproduction_on,
      reef_x           = params$reef_x,
      reef_y           = params$reef_y
    )

    assembleRunResults(
      run_id     = sprintf("%s_%04d", experiment, i),
      experiment = experiment,
      seed       = seed,
      replicate  = rep_no,
      params     = params,
      checkpoints = cp,
      scenario   = scenario
    )
  })

  rbindUnion(rows)
}


# --------------------------------------------------------------------
#  Collect the flat checkpoint table for a simulation AND save it: an
#  <experiment>_results .rds/.csv, and (by default) append into the master.
#    dir : output folder (defaults to "Results" in the working directory =
#          the model root when Master_Script.r sets it)
#  Returns the flat table invisibly, so you can View() it too.
# --------------------------------------------------------------------
saveSimulation <- function(sim, experiment = "main", seed = NA, interval = 10,
                           dir = "Results", to_master = TRUE) {
  tbl <- collectSimTable(sim, experiment, seed, interval)
  saveResults(tbl, name = paste0(experiment, "_results"), dir = dir)
  if (to_master) {
    appendToMaster(tbl, master = file.path(dir, "master_results.rds"))
  }
  invisible(tbl)
}

