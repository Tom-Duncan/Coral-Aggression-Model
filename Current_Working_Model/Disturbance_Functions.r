#Functions related to generating disturbance events
#
#A disturbance event is a sudden mass die-off: a random shaped patch, of a
#random size, in a random place on the reef, where every coral cell is killed
#
#Everything here is deliberately split into one small function per step so each
#piece can be tuned or replaced on its own:
#   1. Tuning numbers          - how big / how often, all in one place
#   2. getDisturbanceCellTarget- how many cells a disturbance covers
#   3. pickDisturbanceCentre   - where it starts
#   4. generateDisturbanceCells- the random shape it takes
#   5. applyDisturbance        - actually kill the cells
#   6. triggerDisturbance      - run one full event (steps 2-5 together)
#   7. buildDisturbanceSchedule- which timesteps events happen on
#   8. maybeApplyDisturbance   - called each timestep, fires events if due
#   9. buildDisturbanceConfig  - bundle the settings into one list
#  10. getDisturbanceSetup     - ask the user / randomise during setup


# ====================================================================
#  TUNING NUMBERS  - change these to alter how disturbances behave.
#  Kept together and clearly named so they are easy to find and edit.
# ====================================================================

#Disturbance size as a FRACTION OF THE WHOLE REEF (number of cells covered).
#Small covers 10% of the reef, medium 20%, large 30%.
DISTURBANCE_SIZE_FRACTION <- c(
  small  = 0.10,
  medium = 0.20,
  large  = 0.30
)

#How often disturbances happen, expressed as events PER 100 TIMESTEPS.
#Rarely ~1 per 100 steps, often ~2 per 100, very often ~4 per 100.
DISTURBANCE_RATE_PER_100 <- c(
  rarely     = 1,
  often      = 2,
  very_often = 4
)

#Size-frequency weights for the "weighted" size option: like natural reefs, SMALL
#disturbances are common and LARGE ones rare (inverse size-frequency relation).
#Each event draws its size category with these probabilities.
DISTURBANCE_SIZE_WEIGHTS <- c(
  small  = 0.50,
  medium = 0.30,
  large  = 0.20
)

#In fully-random setup, the chance disturbances are switched on at all.
DISTURBANCE_RANDOM_ENABLE_CHANCE <- 0.5

#Shapes a disturbance patch may take; one is chosen at random per event. Set to a
#single value (e.g. c("disc")) to force one shape. Options: "disc", "square",
#"rectangle". Shapes are built directly at the target size (fast), rather than
#grown cell-by-cell.
DISTURBANCE_SHAPES <- c("disc", "square", "rectangle")


# --------------------------------------------------------------------
#  STEP 2: How many cells should this disturbance cover?
#  A fraction of the whole reef, so it takes the reef's total cell count
#  (works for square or rectangular reefs alike).
# --------------------------------------------------------------------
getDisturbanceCellTarget <- function(size_category, total_cells) {
  fraction <- DISTURBANCE_SIZE_FRACTION[[size_category]]
  #Always at least one cell so a disturbance is never empty
  max(1, round(fraction * total_cells))
}


# --------------------------------------------------------------------
#  STEP 3: Pick a random starting cell (the centre of the disturbance).
#  Takes the reef's number of rows (n_rows = height) and columns
#  (n_cols = width) so it works on rectangular reefs.
#  Returns a c(row, col) coordinate somewhere on the reef.
# --------------------------------------------------------------------
pickDisturbanceCentre <- function(n_rows, n_cols) {
  c(
    sample(seq_len(n_rows), 1),
    sample(seq_len(n_cols), 1)
  )
}


# --------------------------------------------------------------------
#  STEP 4: Build the SHAPE of the disturbance patch.
#  A pre-defined geometric shape of the target area is placed at the centre and its
#  cells are selected in one vectorised pass - far faster than growing a blob
#  cell-by-cell. A shape is chosen at random per event from DISTURBANCE_SHAPES, so
#  patches still vary in form:
#     disc      - a filled circle of radius sqrt(area / pi)
#     square    - a filled square of side sqrt(area)
#     rectangle - a filled rectangle of the target area with a random aspect ratio
#  The shape WRAPS around the reef edges (toroidal), rather than being clipped, so
#  that every event kills the same number of cells wherever its centre lands and no
#  part of the reef is disturbed more often than another. A patch whose centre is
#  near an edge simply continues onto the opposite edge instead of being cut short.
#  Returns a 2-column matrix of [row, col] cells to be killed.
# --------------------------------------------------------------------
generateDisturbanceCells <- function(centre, target_cells, n_rows, n_cols) {
  cr <- centre[1]
  cc <- centre[2]

  shape <- if (length(DISTURBANCE_SHAPES) == 1) {
    DISTURBANCE_SHAPES
  } else {
    sample(DISTURBANCE_SHAPES, 1)
  }

  if (shape == "disc") {
    #Radius giving area pi*r^2 == target_cells. The bounding box is taken at full
    #size (it may extend past the reef); the overhang is wrapped, not clipped.
    r     <- sqrt(target_cells / pi)
    grid  <- expand.grid(row = floor(cr - r):ceiling(cr + r),
                         col = floor(cc - r):ceiling(cc + r))
    #Keep cells whose centre lies within the circle
    inside <- (grid$row - cr)^2 + (grid$col - cc)^2 <= r^2
    cells  <- as.matrix(grid[inside, c("row", "col")])

  } else {
    #Square, or rectangle with a random aspect ratio, both of area ~ target_cells
    if (shape == "rectangle") {
      aspect <- runif(1, 0.4, 2.5)
    } else {
      aspect <- 1
    }
    half_w <- (sqrt(target_cells * aspect) - 1) / 2   #half-width  (columns)
    half_h <- (sqrt(target_cells / aspect) - 1) / 2   #half-height (rows)

    cells  <- as.matrix(expand.grid(row = round(cr - half_h):round(cr + half_h),
                                    col = round(cc - half_w):round(cc + half_w)))
  }

  #Wrap any cells that fell outside the reef around to the opposite edge (toroidal),
  #so the patch keeps its full size no matter where the centre is. Deduplicate in
  #case a very large patch overlaps itself after wrapping.
  cells[, 1] <- ((cells[, 1] - 1L) %% n_rows) + 1L
  cells[, 2] <- ((cells[, 2] - 1L) %% n_cols) + 1L
  cells <- unique(cells)

  #Return a plain 2-column [row, col] matrix for reef[cells] <- NA indexing
  dimnames(cells) <- NULL
  cells
}


# --------------------------------------------------------------------
#  STEP 5: Apply the disturbance - kill every cell in the patch.
#  Killing a cell just means setting it back to NA (empty).
#  Returns the updated reef.
# --------------------------------------------------------------------
applyDisturbance <- function(reef, cells) {
  if (nrow(cells) == 0) {
    return(reef)
  }
  #cells is a 2-column [row, col] matrix, which indexes the reef directly
  reef[cells] <- NA
  reef
}


# --------------------------------------------------------------------
#  STEP 6: Run one complete disturbance event.
#  Picks a size (random per event if requested), a place and a shape,
#  then kills those cells. Returns the updated reef.
# --------------------------------------------------------------------
triggerDisturbance <- function(reef, size_category) {

  #Read the reef's real dimensions so square and rectangular reefs both work
  n_rows <- nrow(reef)   #height (y)
  n_cols <- ncol(reef)   #width  (x)

  #"random" means choose a fresh size for this individual event (uniform);
  #"weighted" draws by DISTURBANCE_SIZE_WEIGHTS (small common, large rare -
  #the inverse size-frequency structure of natural disturbance regimes)
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


# --------------------------------------------------------------------
#  STEP 7: Decide which timesteps disturbances happen on.
#  Turns the chosen frequency into a number of events spread across the
#  whole run, then picks that many random timesteps. Done once up front
#  so the simulation loop just has to check the schedule each step.
#  Returns a sorted integer vector of timestep numbers.
# --------------------------------------------------------------------
buildDisturbanceSchedule <- function(n_steps, frequency) {
  rate     <- DISTURBANCE_RATE_PER_100[[frequency]]
  n_events <- round(rate * n_steps / 100)

  if (n_events < 1) {
    return(integer(0))
  }

  #Cannot have more events than there are timesteps
  n_events <- min(n_events, n_steps)
  sort(sample(seq_len(n_steps), n_events))
}


# --------------------------------------------------------------------
#  STEP 8: Called every timestep - fire a disturbance if one is due.
#  Does nothing if disturbances are switched off or none are scheduled
#  for this timestep. Returns the (possibly disturbed) reef.
# --------------------------------------------------------------------
maybeApplyDisturbance <- function(reef, timestep_num, config) {

  #No config, or disturbances turned off, so nothing happens
  if (is.null(config) || !isTRUE(config$enabled)) {
    return(reef)
  }

  #How many events land on this exact timestep
  n_events <- sum(config$schedule == timestep_num)
  if (n_events == 0) {
    return(reef)
  }

  for (e in seq_len(n_events)) {
    reef <- triggerDisturbance(reef, config$size)
  }

  reef
}


# --------------------------------------------------------------------
#  Short human-readable summary of the disturbance settings, e.g.
#  "Large and rarely", or "Off" when disturbances are disabled.
#  Used when printing the setup summary.
# --------------------------------------------------------------------
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


# --------------------------------------------------------------------
#  Quick check used only for printing: is a disturbance due this timestep?
#  Returns TRUE if disturbances are on and one is scheduled for this step.
# --------------------------------------------------------------------
isDisturbanceStep <- function(timestep_num, config) {
  if (is.null(config) || !isTRUE(config$enabled)) {
    return(FALSE)
  }
  any(config$schedule == timestep_num)
}


# --------------------------------------------------------------------
#  STEP 9: Bundle all the disturbance settings into one list.
#  The schedule is built here so it is fixed for the whole run.
# --------------------------------------------------------------------
buildDisturbanceConfig <- function(enabled, frequency, size, n_steps) {
  schedule <- if (enabled) buildDisturbanceSchedule(n_steps, frequency) else integer(0)

  #Reef dimensions are not stored: disturbances read them from the live reef when
  #they fire, so square and rectangular reefs are handled automatically.
  list(
    enabled   = enabled,
    frequency = frequency,
    size      = size,
    schedule  = schedule
  )
}


# --------------------------------------------------------------------
#  Small helpers to ask the user for the frequency and size in manual mode.
# --------------------------------------------------------------------
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


# --------------------------------------------------------------------
#  STEP 10: Ask for the disturbance settings during setup.
#  Behaviour follows the overall setup mode, matching the rest of the model:
#    Random (3)      - randomly toggled on/off, random frequency and size
#    Semi-manual (2) - ask only whether disturbances happen, rest randomised
#    Manual (1)      - ask whether, and if so the exact frequency and size
#  Returns just the settings (enabled/frequency/size); the per-run schedule is
#  built later by buildDisturbanceConfig once the simulation length is known.
# --------------------------------------------------------------------
getDisturbanceSetup <- function(setup_mode) {

  #FULLY RANDOM
  if (setup_mode == 3) {
    enabled   <- runif(1) < DISTURBANCE_RANDOM_ENABLE_CHANCE
    frequency <- sample(c("rarely", "often", "very_often"), 1)
    size      <- "random"

  #SEMI-MANUAL: only asks if they should happen, then randomises the details
  } else if (setup_mode == 2) {
    enabled   <- getYesNo("Include disturbance events?")
    frequency <- sample(c("rarely", "often", "very_often"), 1)
    size      <- "random"

  #MANUAL: direct control over whether, how often and how big
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
