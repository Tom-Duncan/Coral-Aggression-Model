# Disturbance events: a random-shaped patch of a random size, in a random place,
# where every coral cell is killed. Split into one small function per step:
#   1. tuning numbers            - size / frequency, all in one place
#   2. getDisturbanceCellTarget  - how many cells an event covers
#   3. pickDisturbanceCentre     - where it starts
#   4. generateDisturbanceCells  - the patch shape
#   5. applyDisturbance          - kill the cells
#   6. triggerDisturbance        - one full event (2-5)
#   7. buildDisturbanceSchedule  - which timesteps events fire on
#   8. maybeApplyDisturbance     - called each step, fires if due
#   9. buildDisturbanceConfig    - bundle the settings
#  10. getDisturbanceSetup       - ask / randomise during setup


# ---- Tuning numbers ---------------------------------------------------------

# Event size as a fraction of the whole reef (small 10%, medium 20%, large 30%).
DISTURBANCE_SIZE_FRACTION <- c(
  small  = 0.10,
  medium = 0.20,
  large  = 0.30
)

# Frequency as events per 100 timesteps.
DISTURBANCE_RATE_PER_100 <- c(
  rarely     = 1,
  often      = 2,
  very_often = 4
)

# Weights for the "weighted" size option: small common, large rare, as on real
# reefs (inverse size-frequency relation).
DISTURBANCE_SIZE_WEIGHTS <- c(
  small  = 0.50,
  medium = 0.30,
  large  = 0.20
)

# In fully-random setup, the chance disturbances are on at all.
DISTURBANCE_RANDOM_ENABLE_CHANCE <- 0.5

# Patch shapes; one chosen at random per event. Set to a single value to force it.
DISTURBANCE_SHAPES <- c("disc", "square", "rectangle")


# STEP 2: cells covered = a fraction of the reef's total (any reef shape).
getDisturbanceCellTarget <- function(size_category, total_cells) {
  fraction <- DISTURBANCE_SIZE_FRACTION[[size_category]]
  max(1, round(fraction * total_cells))   # never empty
}


# STEP 3: random centre cell as c(row, col); n_rows/n_cols allow rectangular reefs.
pickDisturbanceCentre <- function(n_rows, n_cols) {
  c(
    sample(seq_len(n_rows), 1),
    sample(seq_len(n_cols), 1)
  )
}


# STEP 4: build the patch shape (disc / square / rectangle) directly at the target
# area - faster than growing a blob. The shape WRAPS around reef edges (toroidal)
# rather than being clipped, so every event kills the same number of cells wherever
# it lands. Returns a 2-column [row, col] matrix.
generateDisturbanceCells <- function(centre, target_cells, n_rows, n_cols) {
  cr <- centre[1]
  cc <- centre[2]

  shape <- if (length(DISTURBANCE_SHAPES) == 1) {
    DISTURBANCE_SHAPES
  } else {
    sample(DISTURBANCE_SHAPES, 1)
  }

  if (shape == "disc") {
    # radius giving area pi*r^2 == target_cells; bounding box may overhang (wrapped below)
    r     <- sqrt(target_cells / pi)
    grid  <- expand.grid(row = floor(cr - r):ceiling(cr + r),
                         col = floor(cc - r):ceiling(cc + r))
    inside <- (grid$row - cr)^2 + (grid$col - cc)^2 <= r^2   # centre within the circle
    cells  <- as.matrix(grid[inside, c("row", "col")])

  } else {
    # square, or rectangle with a random aspect ratio, both of area ~ target_cells
    if (shape == "rectangle") {
      aspect <- runif(1, 0.4, 2.5)
    } else {
      aspect <- 1
    }
    half_w <- (sqrt(target_cells * aspect) - 1) / 2   # half-width  (columns)
    half_h <- (sqrt(target_cells / aspect) - 1) / 2   # half-height (rows)

    cells  <- as.matrix(expand.grid(row = round(cr - half_h):round(cr + half_h),
                                    col = round(cc - half_w):round(cc + half_w)))
  }

  # wrap cells past the edge around to the opposite side (toroidal); dedupe overlaps
  cells[, 1] <- ((cells[, 1] - 1L) %% n_rows) + 1L
  cells[, 2] <- ((cells[, 2] - 1L) %% n_cols) + 1L
  cells <- unique(cells)

  dimnames(cells) <- NULL   # plain matrix for reef[cells] <- NA
  cells
}


# STEP 5: kill every cell in the patch (set to NA).
applyDisturbance <- function(reef, cells) {
  if (nrow(cells) == 0) {
    return(reef)
  }
  reef[cells] <- NA
  reef
}


# STEP 6: run one event - pick size, place and shape, then kill those cells.
triggerDisturbance <- function(reef, size_category) {

  n_rows <- nrow(reef)   # height (y)
  n_cols <- ncol(reef)   # width  (x)

  # "random" = fresh uniform size this event; "weighted" = draw by SIZE_WEIGHTS
  if (identical(size_category, "random")) {
    size_category <- sample(c("small", "medium", "large"), 1)
  } else if (identical(size_category, "weighted")) {
    size_category <- sample(names(DISTURBANCE_SIZE_WEIGHTS), 1,
                            prob = DISTURBANCE_SIZE_WEIGHTS)
  }

  target <- getDisturbanceCellTarget(size_category, n_rows * n_cols)
  centre <- pickDisturbanceCentre(n_rows, n_cols)
  cells  <- generateDisturbanceCells(centre, target, n_rows, n_cols)

  applyDisturbance(reef, cells)
}


# STEP 7: pick the timesteps events fire on, once up front. Returns a sorted vector.
buildDisturbanceSchedule <- function(n_steps, frequency) {
  rate     <- DISTURBANCE_RATE_PER_100[[frequency]]
  n_events <- round(rate * n_steps / 100)

  if (n_events < 1) {
    return(integer(0))
  }

  n_events <- min(n_events, n_steps)   # no more events than timesteps
  sort(sample(seq_len(n_steps), n_events))
}


# STEP 8: fire any events scheduled for this timestep (no-op if off / none due).
maybeApplyDisturbance <- function(reef, timestep_num, config) {

  if (is.null(config) || !isTRUE(config$enabled)) {
    return(reef)
  }

  n_events <- sum(config$schedule == timestep_num)
  if (n_events == 0) {
    return(reef)
  }

  for (e in seq_len(n_events)) {
    reef <- triggerDisturbance(reef, config$size)
  }

  reef
}


# Human-readable settings summary, e.g. "Large and rarely", or "Off".
describeDisturbance <- function(config) {
  if (is.null(config) || !isTRUE(config$enabled)) {
    return("Off")
  }

  size_label <- switch(config$size,
    small  = "Small",
    medium = "Medium",
    large  = "Large",
    random = "Random",
    config$size
  )

  freq_label <- switch(config$frequency,
    rarely     = "rarely",
    often      = "often",
    very_often = "very often",
    config$frequency
  )

  paste(size_label, "and", freq_label)
}


# Is a disturbance due this timestep? (printing only)
isDisturbanceStep <- function(timestep_num, config) {
  if (is.null(config) || !isTRUE(config$enabled)) {
    return(FALSE)
  }
  any(config$schedule == timestep_num)
}


# STEP 9: bundle the settings; the schedule is fixed here for the whole run.
buildDisturbanceConfig <- function(enabled, frequency, size, n_steps) {
  schedule <- if (enabled) buildDisturbanceSchedule(n_steps, frequency) else integer(0)

  # dimensions aren't stored - events read them from the live reef when they fire
  list(
    enabled   = enabled,
    frequency = frequency,
    size      = size,
    schedule  = schedule
  )
}


# Manual-mode prompts for frequency and size.
getDisturbanceFrequencyChoice <- function() {
  choice <- getMenuChoice(
    "Disturbance frequency:",
    c("Rarely", "Often", "Very often")
  )
  c("rarely", "often", "very_often")[choice]
}

getDisturbanceSizeChoice <- function() {
  choice <- getMenuChoice(
    "Disturbance size:",
    c("Small (10%)", "Medium (20%)", "Large (30%)")
  )
  c("small", "medium", "large")[choice]
}


# STEP 10: get the settings during setup, following the overall setup mode
# (3 = random, 2 = ask on/off only, 1 = full manual). The schedule is built later
# by buildDisturbanceConfig, once the run length is known.
getDisturbanceSetup <- function(setup_mode) {

  # fully random
  if (setup_mode == 3) {
    enabled   <- runif(1) < DISTURBANCE_RANDOM_ENABLE_CHANCE
    frequency <- sample(c("rarely", "often", "very_often"), 1)
    size      <- "random"

  # semi-manual: ask only whether they happen
  } else if (setup_mode == 2) {
    enabled   <- getYesNo("Include disturbance events?")
    frequency <- sample(c("rarely", "often", "very_often"), 1)
    size      <- "random"

  # manual: whether, how often, how big
  } else {
    enabled  <- getYesNo("Include disturbance events?")

    if (enabled) {
      frequency <- getDisturbanceFrequencyChoice()
      size      <- getDisturbanceSizeChoice()
    } else {
      frequency <- "rarely"
      size      <- "small"
    }
  }

  list(enabled = enabled, frequency = frequency, size = size)
}
