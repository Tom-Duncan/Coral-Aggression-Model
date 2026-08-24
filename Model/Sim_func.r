
#Growth trait range. Growth spans 1-5, so GROWTH_MAX = 5 is the top of the growth
#trait range (growth is the only per-species trait; interactions live in the
#matrix). It is used to clamp the
#(effective) growth value before it is turned into a per-cell expansion probability
#by GROWTH_RATE_TABLE, and it caps growth wherever the trait is generated.
GROWTH_MAX <- 5


# ====================================================================
#  EXPRESSED GROWTH RATES  (trait value -> per-cell expansion probability)
# --------------------------------------------------------------------
#  Growth happens PER FRONTIER CELL: every cell on a colony's edge independently
#  rolls this probability. A big colony has a long perimeter = many rolls, so it
#  still adds many cells per step while a one-cell colony crawls - that natural
#  "big grows fast, small grows slow" behaviour is intentional and unchanged.
#
#  GROWTH_RATE_TABLE maps each growth trait value (1..5) to the per-cell chance,
#  applied to EVERY colony throughout the whole simulation. Growth rates are higher
#  overall than before, rising in steps of 15 percentage points from 15% (growth 1)
#  to 75% (growth 5). Edit the numbers to retune (keeping them increasing means a
#  higher growth trait still spreads faster).
#
#    growth trait:   1     2     3     4     5
#    per-cell chance:.15   .30   .45   .60   .75
GROWTH_RATE_TABLE <- c(0.15, 0.30, 0.45, 0.60, 0.75)


# --------------------------------------------------------------------
#  Per-cell expansion probability for a given (effective) growth trait value, read
#  from GROWTH_RATE_TABLE. The value is clamped to [1, GROWTH_MAX]; a fractional
#  value (e.g. from a size modifier) is linearly interpolated between table entries.
# --------------------------------------------------------------------
growthChance <- function(growth_rate) {
  g  <- max(1, min(GROWTH_MAX, growth_rate))
  lo <- floor(g); hi <- ceiling(g)
  if (lo == hi) {
    return(GROWTH_RATE_TABLE[lo])
  }
  frac <- g - lo
  GROWTH_RATE_TABLE[lo] * (1 - frac) + GROWTH_RATE_TABLE[hi] * frac
}


# ====================================================================
#  GROWTH SIZE EFFECT  (toggleable; log-ratio, the growth analogue of the competition
#  size effect)
# --------------------------------------------------------------------
#  For a size-affected species, a colony's per-cell growth chance is shifted in
#  LOG-ODDS space by the log of its OWN size (cell count), scaled by a per-species
#  sensitivity gamma:
#      p = logistic( logit(growthChance(trait)) + gamma * ln(size) )
#  A one-cell recruit (ln 1 = 0) grows at the base chance; the chance then shifts
#  gently as the colony gets bigger. gamma > 0 = grows relatively faster when large,
#  gamma < 0 = faster when small, gamma = 0 = unaffected. Deliberately gentle
#  (GROWTH_SIZE_GAMMA small) so the shift is moderate, not explosive. OFF unless
#  switched on during setup.
# --------------------------------------------------------------------
GROWTH_SIZE_GAMMA <- 0.1    #default per-species growth size sensitivity magnitude
GROWTH_SIZE_RANDOM_ENABLE_CHANCE <- 0.5

#Default growth gammas: one per species, no effect (all 0).
defaultGrowthGamma <- function(species_names) {
  setNames(rep(0, length(species_names)), species_names)
}

#Random growth gammas: each species gets a random sign x GROWTH_SIZE_GAMMA (some grow
#faster when large, some when small). Used when the effect is switched on at random.
randomGrowthGamma <- function(species_names) {
  setNames(sample(c(-1, 1), length(species_names), replace = TRUE) * GROWTH_SIZE_GAMMA,
           species_names)
}

#Attach / read the growth-size config on species_traits (an attribute, so it rides
#along with the traits through the simulation like the interaction config).
#Stored as list(enabled = logical, gamma = named per-species sensitivity vector).
attachGrowthSize <- function(species_traits, enabled, gamma = NULL) {
  if (is.null(gamma)) gamma <- defaultGrowthGamma(species_traits$species)
  attr(species_traits, "growth_size") <- list(
    enabled = isTRUE(enabled),
    gamma   = gamma
  )
  species_traits
}

getGrowthSize <- function(species_traits) {
  gs <- attr(species_traits, "growth_size", exact = TRUE)
  if (is.null(gs)) return(list(enabled = FALSE, gamma = NULL))
  if (is.null(gs$gamma)) gs$gamma <- defaultGrowthGamma(species_traits$species)
  gs
}

#Per-cell growth chance for a specific colony WITH the size effect applied: the base
#chance from its growth trait, shifted in log-odds space by gamma * ln(cell count).
#When the effect is off (or gamma = 0) this is exactly growthChance(coral$growth).
growthChanceSized <- function(coral, growth_size) {
  base <- growthChance(coral$growth)
  if (!isTRUE(growth_size$enabled) || is.null(growth_size$gamma)) return(base)
  g <- growth_size$gamma[coral$species]
  if (is.na(g) || g == 0) return(base)
  z <- coral$size
  if (is.null(z) || is.na(z) || z <= 0) return(base)
  #Hard 0/1 base chances stay hard (logit undefined at the bounds)
  if (base <= 0) return(0)
  if (base >= 1) return(1)
  x <- log(base / (1 - base)) + as.numeric(g) * log(z)
  1 / (1 + exp(-x))
}

#Ask during setup whether the growth size effect is on. Fully random switches it on
#or off at random; manual / semi-manual ask (1 = On, 2 = Off).
getGrowthSizeSetup <- function(setup_mode) {
  if (setup_mode == 3) return(runif(1) < GROWTH_SIZE_RANDOM_ENABLE_CHANCE)
  getMenuChoice("Growth size effect (growth rate shifts with colony size)?",
                c("On", "Off")) == 1
}


#Function to decide whether one coral may take a cell from another coral
#Single source of truth for every takeover rule, so the growth and resolveClaims
#functions always agree on what is allowed. The pairwise overgrowth matrix is the
#ONLY interaction mechanism:
#  * interspecific pairs use the off-diagonal probability,
#  * same-species (different colony) pairs use the DIAGONAL (intraspecific) value,
#  * a 0 entry means "never" - this is how neutrality is expressed.
#A takeover is POSSIBLE whenever that probability is above zero; the actual
#stochastic outcome is rolled once in resolveClaims. See Interaction_Matrix.r.
#Takes the attacking coral, the occupying coral, the reef and species_traits.
#Returns TRUE if the attacker is allowed to take the occupant's cell.
#`interaction` may be passed pre-resolved (getInteraction() is deterministic) so it
#is not re-fetched for every contested cell in the hot loop; NULL resolves it here.
canTakeCell <- function(attacker, occupant, reef, species_traits, interaction = NULL) {

  #A colony never attacks its own cells, only those of other colonies
  if(attacker$id == occupant$id) {
    return(FALSE)
  }

  #Isolated fragments (walled-in islands with no colony behind them) are eroded by
  #any adjacent colony, one cell at a time, until they disappear.
  if(isTRUE(occupant$isolated)) {
    return(TRUE)
  }

  #The interaction matrix governs every other case (inter- and intraspecific).
  if (is.null(interaction)) interaction <- getInteraction(species_traits)
  return(overgrowthProb(attacker$species, occupant$species, interaction) > 0)
}


# --------------------------------------------------------------------
#  Shared hot-loop helpers.
#  NEIGHBOUR_OFFSETS: the 4-connected N/S/W/E offsets, as a module-level constant so
#  it is not rebuilt on every growth() call (the row order is unchanged, so the
#  random direction draw sample(seq_len(nrow(.)), 1) is identical).
#  .coralsById: an id -> coral lookup built ONCE per timestep, replacing the repeated
#  O(n) `corals[sapply(corals, \(x) x$id == target)]` scans with an O(1) name lookup.
#  Ids are unique, so `map[[id]]` returns the same object the old first-match scan did.
# --------------------------------------------------------------------
NEIGHBOUR_OFFSETS <- rbind(
  c(-1, 0),
  c(1, 0),
  c(0, -1),
  c(0, 1)
)

.coralsById <- function(corals) {
  setNames(corals, vapply(corals, function(x) x$id, character(1)))
}


#Function to create list of potential cells for each coral to grow into
#corals_by_id (an id -> coral map), occupied_lin (this colony's cells as column-major
#linear indices) and interaction may be supplied by the caller so they are computed
#ONCE per timestep instead of once per colony; each defaults to being resolved here so
#the function still works when called on its own. None of this changes the order or
#count of random draws, so results are identical.
growth <- function(coral, reef, colony_species, corals, species_traits,
                   corals_by_id = NULL, occupied_lin = NULL, interaction = NULL,
                   reef_vec = NULL) {
  colony_id <- coral$id
  if (is.null(corals_by_id)) corals_by_id <- .coralsById(corals)
  if (is.null(interaction))  interaction  <- getInteraction(species_traits)
  if (is.null(reef_vec))     reef_vec     <- as.vector(reef)
  nr <- nrow(reef); nc <- ncol(reef)

  #Per-cell spread chance for THIS colony: its growth trait (1-GROWTH_MAX) turned into
  #a probability, then shifted by the optional growth size effect for the colony's
  #current size category. It is the same for every frontier cell of this colony (size
  #is fixed within a step), so it is computed once here.
  growth_size     <- getGrowthSize(species_traits)
  per_cell_chance <- growthChanceSized(coral, growth_size)

#The cells already occupied by this colony, as (row, col). Prefer the precomputed
#linear indices (from a single split() over the whole reef) and convert them - this is
#the same set, in the same column-major order, that which(reef == id, arr.ind = TRUE)
#returns, so the frontier order (and thus the RNG draws) is unchanged.
  if (is.null(occupied_lin)) {
    occupied_lin <- which(reef == colony_id)
  }
  n_occ <- length(occupied_lin)
  if (n_occ == 0L) {
    return(matrix(numeric(0), ncol = 2, nrow = 0))
  }
  occ_r <- ((occupied_lin - 1L) %% nr) + 1L
  occ_c <- ((occupied_lin - 1L) %/% nr) + 1L

#List of surrounding N,S,W,E offsets (module-level constant)
  neighbours <- NEIGHBOUR_OFFSETS

#Finding frontier cells: a cell is on the frontier if any in-bounds neighbour is empty
#or belongs to a DIFFERENT colony. Computed by vectorised, direction-at-a-time
#comparisons over ALL occupied cells at once (no per-cell R loop), then OR-ed across
#the four directions. This yields exactly the same mask as checking each cell's
#neighbours in turn - out-of-bounds neighbours are excluded, an empty (NA) or
#different-id neighbour marks the cell - and it preserves occupied order, so the
#frontier order (and thus every later random draw) is unchanged.
  is_frontier <- logical(n_occ)
  for (j in seq_len(nrow(neighbours))) {
    rr <- occ_r + neighbours[j, 1]
    cc <- occ_c + neighbours[j, 2]
    inb <- rr >= 1L & rr <= nr & cc >= 1L & cc <= nc
    if (!any(inb)) {
      next
    }
    nb_lin <- (cc[inb] - 1L) * nr + rr[inb]      #column-major linear index
    vals   <- reef_vec[nb_lin]
    #An empty (NA) or different-colony neighbour puts this occupied cell on the frontier
    makes_frontier <- is.na(vals) | vals != colony_id
    is_frontier[which(inb)[makes_frontier]] <- TRUE
  }

  front_r <- occ_r[is_frontier]
  front_c <- occ_c[is_frontier]
  n_front <- length(front_r)

  # ---------- GROW FROM FRONTIER ----------
  #Second part of growth function
  #Takes frontier cells, tries to grow from it into ONE neighbouring cell.
  #Candidate cells are written into a preallocated matrix (at most one per frontier
  #cell) and trimmed at the end, instead of grown by rbind.
  potential <- matrix(numeric(0), ncol = 2, nrow = n_front)
  np <- 0L

  #Loops every frontier cell found
  for(i in seq_len(n_front)) {
    #Gets coordinates of the cells
    r <- front_r[i]
    c <- front_c[i]

    #Draws random number 0-1; if greater than the per-cell growth chance (computed
    #once above from the growth trait and any size effect), skip this cell this step.
    if(runif(1) > per_cell_chance) {
      next
    }

    #Chooses on random neighbour direction from 4 possible moves
    j <- sample(
      seq_len(nrow(neighbours)),
      1
    )

    #Gets its coordinates
    r_new <- r + neighbours[j,1]
    c_new <- c + neighbours[j,2]

#If empty...
   if(
  r_new >= 1 &&
  r_new <= nr &&
  c_new >= 1 &&
  c_new <= nc
) {

  target <- reef[r_new,c_new]

  if(is.na(target)) {

    np <- np + 1L
    potential[np, ] <- c(r_new, c_new)

  } else {

    #Identify the coral currently occupying the target cell (O(1) id lookup)
     target_coral <- corals_by_id[[target]]

    #Use the shared rule (the interaction matrix) to decide if this coral may take
    #the cell: possible whenever the overgrowth probability is above zero.
    if(canTakeCell(coral, target_coral, reef, species_traits, interaction)){

      np <- np + 1L
      potential[np, ] <- c(r_new, c_new)
    }
  }
}
  }
#Returns list of unique candidate cells for potential growth
#species could expand into
  return(
    unique(potential[seq_len(np), , drop = FALSE])
  )
}


#Function that collects claims from growth function
collectClaims <- function(corals, reef, colony_species, species_traits) {

  #Resolve the interaction config and the id -> coral map ONCE for the whole step,
  #and index every colony's cells in a single pass over the reef (split() on the
  #column-major linear indices), rather than re-scanning the reef inside each growth()
  #call. These are passed straight through to growth().
  interaction  <- getInteraction(species_traits)
  corals_by_id <- .coralsById(corals)
  reef_vec     <- as.vector(reef)
  cells_by_id  <- split(seq_along(reef_vec), reef_vec)   # linear idx per colony id

  #Collect each coral's candidate cells into a list, combined in one rbind at the end
  #(no O(n^2) row-by-row growth of the claims table). Order is preserved, so the
  #duplicate removal below sees the same rows in the same order as before.
  claim_list <- vector("list", length(corals))
  for(i in seq_along(corals)) {
    id        <- corals[[i]]$id
    new_cells <- growth(corals[[i]], reef, colony_species, corals, species_traits,
                        corals_by_id, cells_by_id[[id]], interaction, reef_vec)

#Builds data.frame, one row per candidate cell containing:
#coral species and the cell's current occupant
    if(nrow(new_cells) > 0) {
      claim_list[[i]] <- data.frame(
        row = new_cells[,1],
        col = new_cells[,2],
        colony_id = corals[[i]]$id,
        species = corals[[i]]$species,
        occupant = reef[cbind(
          new_cells[,1],
          new_cells[,2])],
        stringsAsFactors = FALSE
      )
    }
  }

  claims <- do.call(rbind, claim_list)
  #All corals declined to grow: return the empty claims table with the right columns
  if (is.null(claims)) {
    return(data.frame(
      row = integer(0), col = integer(0), colony_id = character(0),
      species = character(0), occupant = character(0), stringsAsFactors = FALSE))
  }

  #Removes potential duplicated after processing all corals when identical
  claims <- claims[
    !duplicated(
      paste(
        claims$row,
        claims$col,
        claims$colony_id
      )
    ),
    ,
    drop = FALSE
  ]
  #Returns the claims data.frame
  #Contains cells neighbouring coral want, the cells current occupant
  #Removes duplicates and returns for use in resolveClaims
  return(claims)
}


#How strongly a claimant's SIZE weights the tie-break among successful challengers.
#The size weight is size^TIEBREAK_SIZE_EXPONENT, so only the size RATIO matters and
#the exponent controls how much: 0 = size ignored (tie-break by overgrowth prob only),
#1 = size in full (linear). A small value keeps the tie-break NEAR RANDOM with only a
#slight nudge from size - e.g. at 0.15 a 10x-larger claimant wins ~58% of an otherwise
#even tie (vs ~91% at 1.0), and a 2x-larger one ~53%.
TIEBREAK_SIZE_EXPONENT <- 0.15

#Function to resolve ownership of contested cells
#Candidate growth claims evaluated and reef map updated
#Takes reef, claims to cells and coral
resolveClaims <- function(claims, reef, corals, species_traits) {

  #Checks if any claims, returns if no claims by any coral
  if (nrow(claims) == 0) {
    return(reef)
  }

  #The pairwise overgrowth matrix decides every contested cell: each claimant rolls
  #its per-pair probability once here. See Interaction_Matrix.r.
  interaction <- getInteraction(species_traits)

  #Id -> coral lookup built once, replacing per-claimant O(n) scans of the corals list.
  corals_by_id <- .coralsById(corals)

  #Splits claims into groups of coral that want the same cell
  #So contested cells resolved one at a time
  claim_groups <- split(claims,
    paste(claims$row, claims$col))

  #Loop over each group of claims targeting same cell
  for (group in claim_groups) {
    #Target coordinates
    row <- group$row[1]
    col <- group$col[1]
    #Checks if target cell empty or held by other coral
    occupant <- reef[row, col]

    #An empty cell has no incumbent to overgrow, so settle it at random among the
    #claimants.
    if (is.na(occupant)) {
      winner <- group[sample(nrow(group), 1), , drop = FALSE]
      reef[row, col] <- winner$colony_id}

      #If cell is alreayd occupied
      else {

      #Search list for coral whose ID matches occupant
      #So its rules can be checked against each claimant
      occupant_coral <- corals_by_id[[occupant]]

      #Resolve the claimant coral objects once, then find which are even ALLOWED to
      #take this cell (a nonzero overgrowth probability, or an isolated occupant),
      #so resolution matches what growth permitted.
      claimant_corals <- lapply(seq_len(nrow(group)), function(k)
        corals_by_id[[group$colony_id[k]]])
      allowed <- vapply(claimant_corals, function(cl)
        canTakeCell(cl, occupant_coral, reef, species_traits, interaction), logical(1))

      #If no attacker is allowed to win, the occupant keeps the cell
      if (!any(allowed)) {
        next}

      #Roll each allowed claimant's overgrowth probability ONCE, KEEPING the
      #probability so it can weight the tie-break below. Isolated-fragment takeovers
      #are deterministic (p = 1); every other contest (interspecific or intraspecific =
      #diagonal) uses the matrix probability, size-adjusted per interaction$size_mode
      #(none / logratio).
      idx   <- which(allowed)
      probs <- vapply(idx, function(k) {
        if (isTRUE(occupant_coral$isolated)) 1
        else overgrowthProbSized(claimant_corals[[k]], occupant_coral, interaction)
      }, numeric(1))
      won   <- runif(length(idx)) < probs

      winners_idx <- idx[won]
      #No claimant won its roll: the occupant keeps the cell this step
      if (length(winners_idx) == 0) {
        next}

      #Winner among the claimants that overgrew, chosen with probability WEIGHTED by
      #overgrowth probability x size^TIEBREAK_SIZE_EXPONENT: the stronger challenger is
      #favoured, with only a SLIGHT extra nudge from being larger (the exponent keeps
      #the size effect near-random - see TIEBREAK_SIZE_EXPONENT).
      if (length(winners_idx) == 1) {
        pick <- winners_idx
      } else {
        w_prob   <- probs[won]
        w_size   <- vapply(winners_idx, function(k) claimant_corals[[k]]$size, numeric(1))
        weights  <- w_prob * w_size ^ TIEBREAK_SIZE_EXPONENT
        pick <- if (all(is.finite(weights)) && sum(weights) > 0) {
          sample(winners_idx, 1, prob = weights)
        } else {
          sample(winners_idx, 1)          #defensive: uniform if weights degenerate
        }
      }
      winner <- group[pick, , drop = FALSE]

      #Replace occupant cell with winner
      reef[row, col] <- winner$colony_id
    }
  }
  #Return changed reef with updated cells
  reef
}


#Function to update each coral's size after conflict and growth.
#Size is stored as the NUMBER OF CELLS a colony occupies (a raw cell count). Percent
#cover, where it is needed (size categories, reproduction scaling, reported metrics),
#is derived from this via coralCoverPercent() using the reef dimensions.
#Takes coral list and reef matrix as inputs
updateSize <- function(corals, reef) {

  #Tally every colony's cell count in ONE pass over the reef (empty cells are NA and
  #excluded by table), then read each colony's size from the tally - instead of
  #re-scanning the whole reef once per colony. An id absent from the tally has 0 cells.
  counts <- table(reef)

  #loops through each coral in coral list
  for(i in seq_along(corals)){
    n <- counts[corals[[i]]$id]
    corals[[i]]$size <- if (is.na(n)) 0L else as.integer(n)
  }
  #Returns modified corals list with updated size (as a cell count)
return(corals)
}


#Function for one complete simulation step, in order:
#  growth -> competition -> disturbance -> consolidation -> reproduction
#Disturbance is applied AFTER growth/combat but BEFORE consolidation so that a
#single consolidation reconciles the fragmentation caused by both combat and the
#disturbance - there is no need to consolidate twice. Reproduction happens last,
#after the disturbance, so a recruit cannot be killed by the same step's event and
#colonies reproduce at their post-disturbance size.
#Takes reef, corals, colony_species, running counters for split fragments and
#reproduction recruits (so every colony born this step gets a unique id), and the
#optional disturbance config plus the current timestep number.
timestep <- function(reef, corals, colony_species, species_traits,
                     fragment_counter, recruit_counter,
                     disturbance_config = NULL, timestep_num = NA, habitat = NULL) {

#Age every colony currently alive by one timestep. Done first so that this step's
#maturity checks (reproduction) see the updated age; colonies born LATER this step
#(split fragments inherit their parent's age, recruits start at 0) are appended
#afterwards and so are not aged until the next timestep. Age drives maturity.
  for (i in seq_along(corals)) {
    corals[[i]]$age <- (if (is.null(corals[[i]]$age)) 0L else corals[[i]]$age) + 1L
  }

#Asks each coral where it wants to grow, returning data.frame
#of all potential cells each coral claims, and current occupant of each cell

  claims <- collectClaims(corals, reef, colony_species, species_traits)


#Processes all claims for contested cells, updates reef with new occupant
  reef <- resolveClaims(claims, reef, corals, species_traits)

#Apply any disturbance event scheduled for this timestep (mass die-off), after the
#growth/combat above but before consolidation. Does nothing if disturbances are off
#or none are due. Defined in Disturbance.r
  reef <- maybeApplyDisturbance(reef, timestep_num, disturbance_config)

#Relabel any colony cut into disconnected pieces (by combat OR disturbance) so each
#contiguous blob is its own colony, refresh sizes, absorb islands and drop fully
#eroded fragments - all in one pass. Defined in Colony_Tracking.r
  consolidated <- consolidateColonies(reef, corals, colony_species, fragment_counter)

#Mature colonies may spawn new single-cell recruits on empty cells (Reproduction.r).
#Does nothing if reproduction is off or no species can reproduce. timestep_num
#drives the synchronised-spawning gate when REPRO_INTERVAL > 1.
  offspring <- reproduceColonies(
    consolidated$reef,
    consolidated$corals,
    consolidated$colony_species,
    species_traits,
    recruit_counter,
    timestep_num
  )

#On a natural reef, larvae also arrive from OUTSIDE the patch: every N steps each
#species in the pool may drop a background recruit anywhere on empty reef, so a
#locally extinct species can return. Does nothing for aquarium/none. (Reproduction.r)
  background <- backgroundReproduce(
    offspring$reef,
    offspring$corals,
    offspring$colony_species,
    species_traits,
    habitat,
    offspring$recruit_counter,
    timestep_num
  )

#Returns updated reef, coral list, colony_species and both id counters
  return(list(reef             = background$reef,
              corals           = background$corals,
              colony_species   = background$colony_species,
              fragment_counter = consolidated$fragment_counter,
              recruit_counter  = background$recruit_counter))
}

#Function to execute entire model and storing every state for n timesteps
#disturbance_config is optional: pass NULL (default) for no disturbances, or a
#config from getDisturbanceSetup() in Disturbance.r to enable mass die-off events
runSimulation <- function(reef, corals, n_steps, colony_species, species_traits,
                          disturbance_config = NULL, habitat = NULL) {

#Empty list to store simulation state at each timestep
  states <- vector("list", n_steps + 1)

#stores starting positions at position 1 in list
  states[[1]] <- list(
    reef = reef,
    corals = corals
  )

#Running counters that hand out unique ids to colonies born mid-run: split-off
#fragments (Colony_Tracking.r) and reproduction recruits (Reproduction.r).
#Shared across every timestep so ids never collide
  fragment_counter <- 0L
  recruit_counter  <- 0L

#Loops once each timestep from 1 to n timesteps
  for(t in 1:n_steps) {

#Prints current step number in console to see progress and when slows
#Also flags the step if a disturbance event (mass die-off) happens on it
  cat("Step:", t,
      if (isDisturbanceStep(t, disturbance_config)) "(Disturbance event)" else "",
      "\n")

#Run one complete step (growth, combat, disturbance, consolidation, reproduction),
#returning the updated state. Disturbance is handled inside timestep now, so only a
#single consolidation is done per step. colony_species and both counters carry
#forward the colonies born here
  state <- timestep(
    reef,
    corals,
    colony_species,
    species_traits,
    fragment_counter,
    recruit_counter,
    disturbance_config,
    t,
    habitat
  )
  reef <- state$reef
  corals <- state$corals
  colony_species <- state$colony_species
  fragment_counter <- state$fragment_counter
  recruit_counter <- state$recruit_counter

#Saves new state into states list at next position
  states[[t + 1]] <- list(
    reef = reef,
    corals = corals
  )
}
#Returns states after loop finishes
  return(states)
}

#Functions for tracking individual colonies as they split apart
#
#A "colony" is one contiguous (4-connected) blob of cells sharing a colony id.
#During growth and combat a colony can be cut into disconnected pieces. If those
#pieces kept the same id, the colony would count every scattered cell toward one
#size - inflating its size category and its size-adjusted overgrowth probability -
#even though it is really several smaller colonies.
#
#To prevent that, after every timestep we relabel the reef so that each id is a
#single contiguous blob again: the largest piece keeps the original id, and every
#other piece becomes a NEW colony (a fresh id) that inherits its parent's species
#and traits but is sized by its own blob. Because each broken-off fragment is now
#a small, weak colony, ordinary combat can slowly erode it (no special erosion
#rule is added - fragments simply compete at their true, small size).
#
#Everything is split into one small function per step, matching the rest of the
#model:
#   1. Tuning numbers            - the "isolated fragment" flag threshold
#   2. labelConnectedComponents  - label every contiguous same-id blob
#   3. splitColonies             - spin disconnected pieces off as new colonies
#   4. flagIsolatedFragments     - mark tiny enclosed fragments (informational)
#   5. consolidateColonies       - the full per-timestep pass (split + resize)
#   6. buildColonySpeciesFromStates - rebuild the id -> species map after a run


# ====================================================================
#  TUNING NUMBERS
# ====================================================================

#A broken-off fragment smaller than this (in cells) that is completely enclosed
#by other colonies is flagged as "isolated" on the coral object ($isolated).
#This is informational only - it does NOT change combat - so plots or exports can
#highlight little embedded islands. Raise/lower to taste, or ignore the flag.
ISOLATED_FRAGMENT_MAX <- 5

#Any colony piece this size or smaller (in cells) that is embedded inside a larger
#colony is absorbed into that colony during consolidation, so isolated cells never
#persist inside another colony. This is the decisive island cleanup - it applies to
#ALL small embedded pieces (fragments or whole colonies), unlike the flag above.
#Set to 0 to switch the cleanup off entirely.
ISLAND_ABSORB_MAX <- 5


# --------------------------------------------------------------------
#  STEP 2: Label every 4-connected blob of same-id cells.
#  Returns an integer matrix the same size as the reef: each contiguous run of
#  identical colony ids gets a unique label, empty (NA) cells get 0.
#  Uses linear indices and a preallocated stack so it stays O(number of cells).
# --------------------------------------------------------------------
labelConnectedComponents <- function(reef) {
  nr <- nrow(reef)
  nc <- ncol(reef)
  n  <- nr * nc

  reef_vec   <- as.vector(reef)     #column-major, matches linear indexing
  labels     <- integer(n)          #0 = unlabelled / empty
  next_label <- 0L
  stack      <- integer(n)          #preallocated flood-fill stack

  for (start in seq_len(n)) {
    #Skip empty cells and cells already assigned to a component
    if (is.na(reef_vec[start]) || labels[start] != 0L) {
      next
    }

    next_label <- next_label + 1L
    id <- reef_vec[start]
    sp <- 1L
    stack[sp]      <- start
    labels[start]  <- next_label

    #Flood fill outward from this seed cell across matching neighbours
    while (sp > 0L) {
      cur <- stack[sp]
      sp  <- sp - 1L

      #Recover row/col of this linear index (column-major layout)
      r <- ((cur - 1L) %% nr) + 1L
      c <- ((cur - 1L) %/% nr) + 1L

      #Up / down / left / right neighbours as linear-index offsets
      if (r > 1L)  { nb <- cur - 1L;  if (labels[nb] == 0L && !is.na(reef_vec[nb]) && reef_vec[nb] == id) { labels[nb] <- next_label; sp <- sp + 1L; stack[sp] <- nb } }
      if (r < nr)  { nb <- cur + 1L;  if (labels[nb] == 0L && !is.na(reef_vec[nb]) && reef_vec[nb] == id) { labels[nb] <- next_label; sp <- sp + 1L; stack[sp] <- nb } }
      if (c > 1L)  { nb <- cur - nr;  if (labels[nb] == 0L && !is.na(reef_vec[nb]) && reef_vec[nb] == id) { labels[nb] <- next_label; sp <- sp + 1L; stack[sp] <- nb } }
      if (c < nc)  { nb <- cur + nr;  if (labels[nb] == 0L && !is.na(reef_vec[nb]) && reef_vec[nb] == id) { labels[nb] <- next_label; sp <- sp + 1L; stack[sp] <- nb } }
    }
  }

  matrix(labels, nrow = nr, ncol = nc)
}


# --------------------------------------------------------------------
#  STEP 3: Split any colony that is in more than one piece.
#  For every id whose cells span several components, the largest component keeps
#  the id; each other component is turned into a new colony that inherits the
#  parent's species and traits (colour included) but is sized by its own blob.
#  fragment_counter guarantees globally unique new ids across the whole run.
#  Returns the updated reef, corals, colony_species and fragment_counter.
# --------------------------------------------------------------------
splitColonies <- function(reef, corals, colony_species, fragment_counter,
                          labels = NULL) {
  nr       <- nrow(reef)
  #Reuse a caller-supplied label matrix when the reef has not changed since it was
  #computed (labelConnectedComponents is a pure function of the reef, so a reused
  #matrix is identical to recomputing); otherwise label here.
  if (is.null(labels)) labels <- labelConnectedComponents(reef)
  reef_vec <- as.vector(reef)
  lab_vec  <- as.vector(labels)

  #Current ids in the corals list, kept in step as we append new fragments
  ids <- vapply(corals, function(x) x$id, character(1))

  present_ids <- unique(reef_vec[!is.na(reef_vec)])

  for (id in present_ids) {
    cell_idx <- which(reef_vec == id)   #linear indices of this colony's cells
    labs     <- lab_vec[cell_idx]
    ulabs    <- unique(labs)

    #Already a single contiguous blob - nothing to split
    if (length(ulabs) <= 1L) {
      next
    }

    #The largest piece keeps the original id
    counts     <- tabulate(match(labs, ulabs))
    keep_label <- ulabs[which.max(counts)]

    parent_pos <- which(ids == id)[1]
    parent     <- corals[[parent_pos]]

    #Every other piece becomes its own colony
    for (lab in ulabs[ulabs != keep_label]) {
      fragment_counter <- fragment_counter + 1L
      new_id   <- paste0(id, ".f", fragment_counter)
      frag_idx <- cell_idx[labs == lab]

      #Relabel this fragment's cells on the reef
      reef_vec[frag_idx] <- new_id

      #Fragment inherits the parent's identity/traits, but is its own colony
      frag <- parent
      frag$id        <- new_id
      frag$parent_id <- id
      #Keep the very first ancestor so lineage can be traced back to setup
      frag$origin_id <- if (!is.null(parent$origin_id)) parent$origin_id else id
      #Size as a cell count (updateSize recomputes this too, but keep it right)
      frag$size      <- length(frag_idx)

      #Centroid of the fragment, stored as coords c(x = col, y = row)
      rows <- ((frag_idx - 1L) %% nr) + 1L
      cols <- ((frag_idx - 1L) %/% nr) + 1L
      frag$coords <- c(round(mean(cols)), round(mean(rows)))

      corals[[length(corals) + 1L]] <- frag
      colony_species[new_id]        <- parent$species
      ids <- c(ids, new_id)
    }
  }

  list(
    reef             = matrix(reef_vec, nrow = nr),
    corals           = corals,
    colony_species   = colony_species,
    fragment_counter = fragment_counter
  )
}


# --------------------------------------------------------------------
#  STEP 4: Flag tiny fragments that are completely walled in by other colonies.
#  Sets $isolated = TRUE on any fragment (a colony with a parent) that is smaller
#  than ISOLATED_FRAGMENT_MAX and whose every edge touches another colony (never
#  an empty cell and never the reef border). Purely informational - combat is
#  untouched - so downstream plots/exports can pick out embedded islands.
# --------------------------------------------------------------------
flagIsolatedFragments <- function(reef, corals) {
  nr       <- nrow(reef)
  nc       <- ncol(reef)
  reef_vec <- as.vector(reef)

  for (k in seq_along(corals)) {
    co <- corals[[k]]
    corals[[k]]$isolated <- FALSE

    #Only fragments (not founder colonies) can be isolated islands
    if (is.null(co$parent_id)) {
      next
    }

    #Island cleanup is a grid-mechanics rule measured in CELLS. co$size is a cell
    #count now, but we recount from the reef here so the test never uses a stale size.
    #Only small, live pieces qualify.
    idx     <- which(reef_vec == co$id)
    n_cells <- length(idx)
    if (n_cells == 0L || n_cells >= ISOLATED_FRAGMENT_MAX) {
      next
    }

    enclosed <- TRUE

    for (cur in idx) {
      r <- ((cur - 1L) %% nr) + 1L
      c <- ((cur - 1L) %/% nr) + 1L

      #A cell on the reef border is not fully enclosed
      if (r == 1L || r == nr || c == 1L || c == nc) {
        enclosed <- FALSE
        break
      }

      nbrs <- c(cur - 1L, cur + 1L, cur - nr, cur + nr)
      #An empty neighbour means the fragment is open, not walled in
      if (any(is.na(reef_vec[nbrs]))) {
        enclosed <- FALSE
        break
      }
    }

    corals[[k]]$isolated <- enclosed
  }

  corals
}


# --------------------------------------------------------------------
#  Absorb small islands embedded inside a larger colony.
#  For every contiguous piece of ISLAND_ABSORB_MAX cells or fewer, the cells just
#  outside it are inspected. If a single OTHER colony borders the piece at least as
#  much as empty water does, and that colony is at least as big, the piece is an
#  island trapped inside it: all its cells are reassigned to that colony. This is
#  what finally clears stranded single cells regardless of whether they are a
#  tracked fragment, a whittled-down original colony, or a strong defender that
#  ordinary combat could never dislodge. A piece sitting mostly in open water is
#  left alone (it is a real small colony, not an embedded island).
#  Returns the updated reef matrix.
# --------------------------------------------------------------------
absorbIslands <- function(reef, corals, colony_species) {
  if (ISLAND_ABSORB_MAX < 1) {
    return(list(reef = reef, labels = NULL, changed = FALSE))
  }

  nr       <- nrow(reef)
  nc       <- ncol(reef)
  labels   <- labelConnectedComponents(reef)
  reef_vec <- as.vector(reef)
  lab_vec  <- as.vector(labels)
  changed  <- FALSE

  for (lab in setdiff(unique(lab_vec), 0L)) {
    comp_idx <- which(lab_vec == lab)
    if (length(comp_idx) > ISLAND_ABSORB_MAX) {
      next
    }

    #Tally the cells bordering this piece: empties vs each neighbouring colony
    empty_count <- 0L
    nbr_ids     <- character(0)
    for (cur in comp_idx) {
      r <- ((cur - 1L) %% nr) + 1L
      c <- ((cur - 1L) %/% nr) + 1L
      nbrs <- c(if (r > 1L) cur - 1L, if (r < nr) cur + 1L,
                if (c > 1L) cur - nr, if (c < nc) cur + nr)
      for (nb in nbrs) {
        if (lab_vec[nb] == lab) {
          next                       #still inside the same piece
        }
        if (is.na(reef_vec[nb])) {
          empty_count <- empty_count + 1L
        } else {
          nbr_ids <- c(nbr_ids, reef_vec[nb])
        }
      }
    }

    #No colony neighbours: floating in open water, not an embedded island
    if (length(nbr_ids) == 0L) {
      next
    }

    #The colony that surrounds it most
    counts   <- sort(table(nbr_ids), decreasing = TRUE)
    dominant <- names(counts)[1]
    dom_bord <- as.integer(counts[1])

    #Must be genuinely embedded: the dominant colony walls it in at least as much
    #as open water does
    if (dom_bord < empty_count) {
      next
    }

    #Must be swallowed by something at least as large (a "larger colony")
    dom_size <- sum(reef_vec == dominant, na.rm = TRUE)
    if (dom_size < length(comp_idx)) {
      next
    }

    #Absorb: hand the island's cells to the surrounding colony
    reef_vec[comp_idx] <- dominant
    changed <- TRUE
  }

  #Return the labels too: if nothing was absorbed the reef is unchanged, so these
  #labels are still valid and splitColonies can reuse them instead of re-labelling.
  list(reef    = matrix(reef_vec, nrow = nr, ncol = nc),
       labels  = if (changed) NULL else labels,
       changed = changed)
}


# --------------------------------------------------------------------
#  STEP 5: One full consolidation pass, run after growth/combat each timestep.
#  First absorbs small islands embedded in larger colonies, then splits any
#  remaining disconnected colonies, refreshes every colony's size from the reef,
#  drops fragments that have been fully eroded (size 0) while keeping the original
#  setup colonies so their extinction is still recorded, then flags tiny islands.
#  Returns the updated reef, corals, colony_species and fragment_counter.
# --------------------------------------------------------------------
consolidateColonies <- function(reef, corals, colony_species, fragment_counter) {

  #Clear embedded islands first so no isolated cells survive inside a colony.
  #absorbIslands returns the connected-component labels it computed; when it absorbed
  #nothing the reef is unchanged, so splitColonies reuses those labels rather than
  #labelling the reef a second time this step.
  absorbed <- absorbIslands(reef, corals, colony_species)
  reef     <- absorbed$reef

  res <- splitColonies(reef, corals, colony_species, fragment_counter,
                       labels = absorbed$labels)

  #Recount sizes now that every id is a single contiguous blob
  res$corals <- updateSize(res$corals, res$reef)

  #Remove dead fragments (0 cells) but always keep the original colonies, which
  #have no parent, so a species going extinct still shows up in the time series
  keep <- vapply(res$corals, function(x) x$size > 0 || is.null(x$parent_id), logical(1))
  res$corals <- res$corals[keep]

  #Mark any tiny enclosed islands (informational only)
  res$corals <- flagIsolatedFragments(res$reef, res$corals)

  res
}


# --------------------------------------------------------------------
#  STEP 6: Rebuild a complete id -> species lookup after a whole run.
#  Because colonies are created mid-run, the initial colony_species map misses
#  the fragment ids. This gathers every colony that ever existed across all saved
#  states so species-level plots can map every reef cell back to its species.
# --------------------------------------------------------------------
buildColonySpeciesFromStates <- function(states, base = NULL) {
  all_corals <- unlist(lapply(states, `[[`, "corals"), recursive = FALSE)

  species_map <- setNames(
    vapply(all_corals, function(x) x$species, character(1)),
    vapply(all_corals, function(x) x$id,      character(1))
  )
  species_map <- species_map[!duplicated(names(species_map))]

  #Fold in any base ids that (unexpectedly) never appeared in a saved state
  if (!is.null(base)) {
    extra <- base[!(names(base) %in% names(species_map))]
    species_map <- c(species_map, extra)
  }

  species_map
}
