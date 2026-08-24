# One tidy (long) results table for simulation sweeps: one row per (run x
# checkpoint), carrying run_id / experiment / replicate / seed / setup params /
# phase (initial|checkpoint|final) / timestep / metrics. Different experiments may
# record different columns; rbindUnion combines them (missing cells -> NA), so one
# master table can span every model variation.


# Compact scenario code, e.g. "3i5Dr_30x30".
scenarioCode <- function(n_species, disturbance_on, reproduction_on,
                         reef_x, reef_y, individuals = NA) {
  onoff <- function(flag, up) if (isTRUE(flag)) up else tolower(up)
  ind   <- if (!is.na(individuals)) paste0("i", individuals) else ""
  paste0(n_species, ind,
         onoff(disturbance_on, "D"),
         onoff(reproduction_on, "R"),
         "_", reef_x, "x", reef_y)
}


# Row-bind data.frames with different columns; result = union of columns, missing -> NA.
rbindUnion <- function(dfs) {
  dfs <- Filter(function(d) !is.null(d) && nrow(d) > 0, dfs)
  if (length(dfs) == 0) return(NULL)

  all_cols <- unique(unlist(lapply(dfs, names)))
  dfs2 <- lapply(dfs, function(d) {
    for (m in setdiff(all_cols, names(d))) d[[m]] <- NA
    d[all_cols]
  })
  out <- do.call(rbind, dfs2)
  rownames(out) <- NULL
  out
}


# Assemble one run's rows (meta + params + phase + checkpoints), params recycled.
# checkpoints must contain `timestep`; scenario is the compact setup code.
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


# Save a table as .rds (fast reload) + .csv (portable). Returns the .rds path.
saveResults <- function(results, name, dir = "Results") {
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  rds <- file.path(dir, paste0(name, ".rds"))
  csv <- file.path(dir, paste0(name, ".csv"))
  saveRDS(results, rds)
  write.csv(results, csv, row.names = FALSE)
  cat("Saved", nrow(results), "rows to", rds, "and", csv, "\n")
  invisible(rds)
}


# Append to the growing master (one file for all sweeps). By default rows from an
# experiment already present are replaced, so re-running overwrites its old rows.
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


# Load the master (or any saved .rds results table).
loadResults <- function(file = file.path("Results", "master_results.rds")) {
  if (!file.exists(file)) {
    stop("No results file at: ", file)
  }
  readRDS(file)
}


# ---- Checkpoint metrics -----------------------------------------------------

# Colony size classes (cells): small < 20, medium 20-60, large > 60.
CSIZE_SMALL_MAX  <- 20
CSIZE_MEDIUM_MAX <- 60

# Structural complexity: of occupied adjacent pairs, the fraction in DIFFERENT
# colonies (interface density). 0 = one uniform patch, ->1 = interdigitated. NA if none.
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

# Species occupying each cell ("" for empty); used for turnover (ids are volatile).
.speciesGrid <- function(state) {
  ids <- vapply(state$corals, function(x) x$id, character(1))
  sp  <- vapply(state$corals, function(x) x$species, character(1))
  out <- setNames(sp, ids)[state$reef]
  out[is.na(out)] <- ""
  matrix(out, nrow = nrow(state$reef), ncol = ncol(state$reef))
}

# Spatial turnover: of cells occupied at either checkpoint, the fraction whose species changed.
.turnover <- function(prev_state, cur_state) {
  a <- .speciesGrid(prev_state); b <- .speciesGrid(cur_state)
  occ <- (a != "") | (b != "")
  if (!any(occ)) return(NA_real_)
  sum(a[occ] != b[occ]) / sum(occ)
}

# Colony deaths: colonies alive at the previous checkpoint holding 0 cells (or gone) now.
.deathsBetween <- function(prev_state, cur_state) {
  pid <- vapply(prev_state$corals, function(x) x$id, character(1))
  psz <- vapply(prev_state$corals, function(x) x$size, numeric(1))
  cid <- vapply(cur_state$corals, function(x) x$id, character(1))
  csz <- setNames(vapply(cur_state$corals, function(x) x$size, numeric(1)), cid)
  alive_prev <- pid[psz > 0]
  cur_of <- csz[alive_prev]                    # NA if the colony no longer exists
  sum(is.na(cur_of) | cur_of <= 0)
}

# Community metrics at each checkpoint (t=0, every `interval`, and the final step).
# Per-species cover (% reef), richness, Shannon, evenness, cover, colony structure.
mainCheckpointMetrics <- function(states, interval = 10) {
  sim_length  <- length(states) - 1            # states[[1]] is timestep 0
  checkpoints <- seq(0, sim_length, by = interval)
  if (checkpoints[length(checkpoints)] != sim_length) {
    checkpoints <- c(checkpoints, sim_length)
  }

  # master species list from the initial state (extinct species stay as 0-cover cols)
  species_list <- sort(unique(vapply(states[[1]]$corals,
                                     function(x) x$species, character(1))))

  rows <- lapply(seq_along(checkpoints), function(k) {
    t        <- checkpoints[k]
    state_t  <- states[[t + 1]]
    corals_t <- state_t$corals
    reef_t   <- state_t$reef
    # $size (cells) -> % of reef so cover is reef-independent
    sizes <- vapply(corals_t, function(x) coralCoverPercent(x, reef_t), numeric(1))
    specs <- vapply(corals_t, function(x) x$species, character(1))
    alive <- sizes > 0
    sizes <- sizes[alive]; specs <- specs[alive]

    cover       <- vapply(species_list, function(sp) sum(sizes[specs == sp]), numeric(1))
    total_cover <- sum(cover)
    present     <- cover > 0
    richness    <- sum(present)

    if (total_cover > 0) {
      p       <- cover[present] / total_cover
      shannon <- -sum(p * log(p))
    } else {
      shannon <- 0
    }
    evenness <- if (richness > 1) shannon / log(richness) else NA_real_

    cell_sizes <- vapply(corals_t, function(x) x$size, numeric(1))   # raw cells for size structure
    live       <- cell_sizes[cell_sizes > 0]
    founders_alive <- sum(vapply(corals_t,   # founders have no parent_id
      function(x) is.null(x$parent_id) && x$size > 0, logical(1)))

    if (k == 1) {   # between-checkpoint metrics need a previous checkpoint
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


# Was disturbance / reproduction on for this replicate?
mainDisturbanceOn <- function(rep_result) {
  !is.null(rep_result$disturbance) && isTRUE(rep_result$disturbance$enabled)
}
mainReproductionOn <- function(rep_result) {
  st <- rep_result$species_traits
  !is.null(st) && "reproduction" %in% names(st) && any(st$reproduction, na.rm = TRUE)
}


# Setup parameters for one replicate, stored on every row.
mainRunParams <- function(rep_result, interval) {
  init_corals <- rep_result$states[[1]]$corals
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


# Flat checkpoint table for a whole sim (without saving). Handles single- and
# multi-replicate ($replicates) sims; one run_id per replicate.
collectSimTable <- function(sim, experiment = "main", seed = NA, interval = 10) {
  reps <- if (!is.null(sim$replicates)) sim$replicates else list(sim)

  rows <- lapply(seq_along(reps), function(i) {
    rep_result <- reps[[i]]
    cp     <- mainCheckpointMetrics(rep_result$states, interval)
    params <- mainRunParams(rep_result, interval)
    rep_no <- if (!is.null(rep_result$replicate)) rep_result$replicate else i

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


# Collect the checkpoint table and save it (<experiment>_results + append to master).
saveSimulation <- function(sim, experiment = "main", seed = NA, interval = 10,
                           dir = "Results", to_master = TRUE) {
  tbl <- collectSimTable(sim, experiment, seed, interval)
  saveResults(tbl, name = paste0(experiment, "_results"), dir = dir)
  if (to_master) {
    appendToMaster(tbl, master = file.path(dir, "master_results.rds"))
  }
  invisible(tbl)
}
