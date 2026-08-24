
# Growth trait range (1-5); GROWTH_MAX caps the effective growth value before it is
# turned into a per-cell expansion probability.
GROWTH_MAX <- 5


# Expressed growth rates: growth trait -> per-cell expansion probability. Growth is
# PER FRONTIER CELL - every edge cell independently rolls this, so a long perimeter
# grows fast and a one-cell colony crawls. Rises 15% (growth 1) to 75% (growth 5).
#    growth trait:   1     2     3     4     5
#    per-cell chance:.15   .30   .45   .60   .75
GROWTH_RATE_TABLE <- c(0.15, 0.30, 0.45, 0.60, 0.75)


# Per-cell expansion probability from GROWTH_RATE_TABLE; clamped to [1, GROWTH_MAX],
# fractional values linearly interpolated.
growthChance <- function(growth_rate) {
  g  <- max(1, min(GROWTH_MAX, growth_rate))
  lo <- floor(g); hi <- ceiling(g)
  if (lo == hi) {
    return(GROWTH_RATE_TABLE[lo])
  }
  frac <- g - lo
  GROWTH_RATE_TABLE[lo] * (1 - frac) + GROWTH_RATE_TABLE[hi] * frac
}


# Growth size effect (toggleable; log-ratio). A colony's per-cell growth chance is
# shifted in log-odds by gamma * ln(size):
#   p = logistic( logit(growthChance(trait)) + gamma * ln(size) )
# gamma > 0 = faster when large, < 0 = faster when small, 0 = unaffected. Off by default.
GROWTH_SIZE_GAMMA <- 0.1
GROWTH_SIZE_RANDOM_ENABLE_CHANCE <- 0.5

# Default growth gammas: one per species, no effect.
defaultGrowthGamma <- function(species_names) {
  setNames(rep(0, length(species_names)), species_names)
}

# Random growth gammas: each species a random sign x GROWTH_SIZE_GAMMA.
randomGrowthGamma <- function(species_names) {
  setNames(sample(c(-1, 1), length(species_names), replace = TRUE) * GROWTH_SIZE_GAMMA,
           species_names)
}

# Attach/read the growth-size config on species_traits (as an attribute).
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

# Per-cell growth chance for a colony with the size effect applied. When off (or
# gamma = 0) this equals growthChance(coral$growth).
growthChanceSized <- function(coral, growth_size) {
  base <- growthChance(coral$growth)
  if (!isTRUE(growth_size$enabled) || is.null(growth_size$gamma)) return(base)
  g <- growth_size$gamma[coral$species]
  if (is.na(g) || g == 0) return(base)
  z <- coral$size
  if (is.null(z) || is.na(z) || z <= 0) return(base)
  # hard 0/1 base chances stay hard (logit undefined at the bounds)
  if (base <= 0) return(0)
  if (base >= 1) return(1)
  x <- log(base / (1 - base)) + as.numeric(g) * log(z)
  1 / (1 + exp(-x))
}

# Ask during setup whether the growth size effect is on.
getGrowthSizeSetup <- function(setup_mode) {
  if (setup_mode == 3) return(runif(1) < GROWTH_SIZE_RANDOM_ENABLE_CHANCE)
  getMenuChoice("Growth size effect (growth rate shifts with colony size)?",
                c("On", "Off")) == 1
}


# May one coral take a cell from another? Single source of truth for takeover rules.
# The overgrowth matrix is the only mechanism: interspecific = off-diagonal, same-
# species = diagonal, 0 = never (neutrality). Takeover POSSIBLE if that prob > 0; the
# stochastic outcome is rolled in resolveClaims. `interaction` may be pre-resolved.
canTakeCell <- function(attacker, occupant, reef, species_traits, interaction = NULL) {

  # a colony never attacks its own cells
  if(attacker$id == occupant$id) {
    return(FALSE)
  }

  # isolated fragments (walled-in islands) are eroded by any adjacent colony
  if(isTRUE(occupant$isolated)) {
    return(TRUE)
  }

  if (is.null(interaction)) interaction <- getInteraction(species_traits)
  return(overgrowthProb(attacker$species, occupant$species, interaction) > 0)
}


# Shared hot-loop helpers. NEIGHBOUR_OFFSETS: 4-connected N/S/W/E offsets as a
# constant (not rebuilt per growth() call). .coralsById: an O(1) id -> coral lookup
# built once per timestep.
NEIGHBOUR_OFFSETS <- rbind(
  c(-1, 0),
  c(1, 0),
  c(0, -1),
  c(0, 1)
)

.coralsById <- function(corals) {
  setNames(corals, vapply(corals, function(x) x$id, character(1)))
}


# Candidate cells one colony may grow into. corals_by_id / occupied_lin / interaction /
# reef_vec may be passed in (computed once per timestep) or resolved here; either way
# the order/count of random draws is identical.
growth <- function(coral, reef, colony_species, corals, species_traits,
                   corals_by_id = NULL, occupied_lin = NULL, interaction = NULL,
                   reef_vec = NULL) {
  colony_id <- coral$id
  if (is.null(corals_by_id)) corals_by_id <- .coralsById(corals)
  if (is.null(interaction))  interaction  <- getInteraction(species_traits)
  if (is.null(reef_vec))     reef_vec     <- as.vector(reef)
  nr <- nrow(reef); nc <- ncol(reef)

  # this colony's per-cell spread chance (same for every frontier cell this step)
  growth_size     <- getGrowthSize(species_traits)
  per_cell_chance <- growthChanceSized(coral, growth_size)

  # this colony's occupied cells (column-major linear indices -> row/col)
  if (is.null(occupied_lin)) {
    occupied_lin <- which(reef == colony_id)
  }
  n_occ <- length(occupied_lin)
  if (n_occ == 0L) {
    return(matrix(numeric(0), ncol = 2, nrow = 0))
  }
  occ_r <- ((occupied_lin - 1L) %% nr) + 1L
  occ_c <- ((occupied_lin - 1L) %/% nr) + 1L

  neighbours <- NEIGHBOUR_OFFSETS

  # frontier cells: an occupied cell with any in-bounds empty/foreign neighbour.
  # Vectorised direction-at-a-time over all occupied cells, OR-ed across directions.
  is_frontier <- logical(n_occ)
  for (j in seq_len(nrow(neighbours))) {
    rr <- occ_r + neighbours[j, 1]
    cc <- occ_c + neighbours[j, 2]
    inb <- rr >= 1L & rr <= nr & cc >= 1L & cc <= nc
    if (!any(inb)) {
      next
    }
    nb_lin <- (cc[inb] - 1L) * nr + rr[inb]
    vals   <- reef_vec[nb_lin]
    makes_frontier <- is.na(vals) | vals != colony_id
    is_frontier[which(inb)[makes_frontier]] <- TRUE
  }

  front_r <- occ_r[is_frontier]
  front_c <- occ_c[is_frontier]
  n_front <- length(front_r)

  # grow from each frontier cell into one random neighbour (candidates into a
  # preallocated matrix, trimmed at the end)
  potential <- matrix(numeric(0), ncol = 2, nrow = n_front)
  np <- 0L

  for(i in seq_len(n_front)) {
    r <- front_r[i]
    c <- front_c[i]

    # roll the per-cell growth chance
    if(runif(1) > per_cell_chance) {
      next
    }

    # random neighbour direction
    j <- sample(
      seq_len(nrow(neighbours)),
      1
    )

    r_new <- r + neighbours[j,1]
    c_new <- c + neighbours[j,2]

   if(
  r_new >= 1 &&
  r_new <= nr &&
  c_new >= 1 &&
  c_new <= nc
) {

  target <- reef[r_new,c_new]

  if(is.na(target)) {   # empty

    np <- np + 1L
    potential[np, ] <- c(r_new, c_new)

  } else {

     target_coral <- corals_by_id[[target]]   # O(1) id lookup

    # occupied: allowed only if the overgrowth prob is above zero
    if(canTakeCell(coral, target_coral, reef, species_traits, interaction)){

      np <- np + 1L
      potential[np, ] <- c(r_new, c_new)
    }
  }
}
  }
  return(
    unique(potential[seq_len(np), , drop = FALSE])
  )
}


# Collect every colony's growth claims into one table.
collectClaims <- function(corals, reef, colony_species, species_traits) {

  # resolve interaction + id map + per-colony cells ONCE for the step, passed to growth()
  interaction  <- getInteraction(species_traits)
  corals_by_id <- .coralsById(corals)
  reef_vec     <- as.vector(reef)
  cells_by_id  <- split(seq_along(reef_vec), reef_vec)   # linear idx per colony id

  # one list entry per colony, combined in a single rbind at the end
  claim_list <- vector("list", length(corals))
  for(i in seq_along(corals)) {
    id        <- corals[[i]]$id
    new_cells <- growth(corals[[i]], reef, colony_species, corals, species_traits,
                        corals_by_id, cells_by_id[[id]], interaction, reef_vec)

    # one row per candidate cell: claimant species + the cell's current occupant
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
  # no growth at all: return an empty table with the right columns
  if (is.null(claims)) {
    return(data.frame(
      row = integer(0), col = integer(0), colony_id = character(0),
      species = character(0), occupant = character(0), stringsAsFactors = FALSE))
  }

  # drop duplicate (cell, claimant) rows
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
  return(claims)
}


# How strongly size weights the tie-break among successful challengers: weight is
# size^exponent, so only the ratio matters. 0 = size ignored, 1 = linear. Small keeps
# it near-random (0.15: a 10x-larger claimant wins ~58% of an even tie).
TIEBREAK_SIZE_EXPONENT <- 0.15

# Resolve contested cells and update the reef.
resolveClaims <- function(claims, reef, corals, species_traits) {

  if (nrow(claims) == 0) {
    return(reef)
  }

  interaction <- getInteraction(species_traits)   # decides every contested cell

  corals_by_id <- .coralsById(corals)

  # group claims by target cell, resolve one cell at a time
  claim_groups <- split(claims,
    paste(claims$row, claims$col))

  for (group in claim_groups) {
    row <- group$row[1]
    col <- group$col[1]
    occupant <- reef[row, col]

    # empty cell: settle at random among claimants
    if (is.na(occupant)) {
      winner <- group[sample(nrow(group), 1), , drop = FALSE]
      reef[row, col] <- winner$colony_id}

      else {

      occupant_coral <- corals_by_id[[occupant]]

      # which claimants are ALLOWED to take this cell (matches what growth permitted)
      claimant_corals <- lapply(seq_len(nrow(group)), function(k)
        corals_by_id[[group$colony_id[k]]])
      allowed <- vapply(claimant_corals, function(cl)
        canTakeCell(cl, occupant_coral, reef, species_traits, interaction), logical(1))

      if (!any(allowed)) {   # occupant keeps the cell
        next}

      # roll each allowed claimant's overgrowth prob once (isolated occupants: p = 1;
      # otherwise the matrix prob, size-adjusted per interaction$size_mode)
      idx   <- which(allowed)
      probs <- vapply(idx, function(k) {
        if (isTRUE(occupant_coral$isolated)) 1
        else overgrowthProbSized(claimant_corals[[k]], occupant_coral, interaction)
      }, numeric(1))
      won   <- runif(length(idx)) < probs

      winners_idx <- idx[won]
      if (length(winners_idx) == 0) {   # nobody won their roll
        next}

      # winner weighted by overgrowth prob x size^exponent (slight size nudge)
      if (length(winners_idx) == 1) {
        pick <- winners_idx
      } else {
        w_prob   <- probs[won]
        w_size   <- vapply(winners_idx, function(k) claimant_corals[[k]]$size, numeric(1))
        weights  <- w_prob * w_size ^ TIEBREAK_SIZE_EXPONENT
        pick <- if (all(is.finite(weights)) && sum(weights) > 0) {
          sample(winners_idx, 1, prob = weights)
        } else {
          sample(winners_idx, 1)          # defensive: uniform if weights degenerate
        }
      }
      winner <- group[pick, , drop = FALSE]

      reef[row, col] <- winner$colony_id
    }
  }
  reef
}


# Update each colony's size (cell count) in one pass over the reef.
updateSize <- function(corals, reef) {

  counts <- table(reef)   # cell count per id (empty cells excluded)

  for(i in seq_along(corals)){
    n <- counts[corals[[i]]$id]
    corals[[i]]$size <- if (is.na(n)) 0L else as.integer(n)
  }
return(corals)
}


# One full simulation step: growth -> competition -> disturbance -> consolidation ->
# reproduction -> background recruitment. Disturbance is applied before consolidation
# so one consolidation reconciles both combat and disturbance; reproduction is last so
# a recruit cannot be killed by the same step's event.
timestep <- function(reef, corals, colony_species, species_traits,
                     fragment_counter, recruit_counter,
                     disturbance_config = NULL, timestep_num = NA, habitat = NULL) {

  # age every living colony by one (first, so this step's maturity checks see it;
  # colonies born later this step are appended afterwards and aged next step)
  for (i in seq_along(corals)) {
    corals[[i]]$age <- (if (is.null(corals[[i]]$age)) 0L else corals[[i]]$age) + 1L
  }

  claims <- collectClaims(corals, reef, colony_species, species_traits)

  reef <- resolveClaims(claims, reef, corals, species_traits)

  # disturbance (after combat, before consolidation)
  reef <- maybeApplyDisturbance(reef, timestep_num, disturbance_config)

  # relabel split colonies, refresh sizes, absorb islands, drop eroded fragments
  consolidated <- consolidateColonies(reef, corals, colony_species, fragment_counter)

  # mature colonies spawn recruits
  offspring <- reproduceColonies(
    consolidated$reef,
    consolidated$corals,
    consolidated$colony_species,
    species_traits,
    recruit_counter,
    timestep_num
  )

  # natural reef: external larvae arrive (no-op for aquarium/none)
  background <- backgroundReproduce(
    offspring$reef,
    offspring$corals,
    offspring$colony_species,
    species_traits,
    habitat,
    offspring$recruit_counter,
    timestep_num
  )

  return(list(reef             = background$reef,
              corals           = background$corals,
              colony_species   = background$colony_species,
              fragment_counter = consolidated$fragment_counter,
              recruit_counter  = background$recruit_counter))
}

# Run the whole model for n_steps, storing every state. disturbance_config: NULL for
# none, or a config from getDisturbanceSetup().
runSimulation <- function(reef, corals, n_steps, colony_species, species_traits,
                          disturbance_config = NULL, habitat = NULL) {

  states <- vector("list", n_steps + 1)

  states[[1]] <- list(   # starting state (timestep 0)
    reef = reef,
    corals = corals
  )

  # counters handing out unique ids to colonies born mid-run (fragments + recruits)
  fragment_counter <- 0L
  recruit_counter  <- 0L

  for(t in 1:n_steps) {

  cat("Step:", t,
      if (isDisturbanceStep(t, disturbance_config)) "(Disturbance event)" else "",
      "\n")

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

  states[[t + 1]] <- list(
    reef = reef,
    corals = corals
  )
}
  return(states)
}

# Colony split tracking. A colony is one contiguous (4-connected) blob. Combat and
# disturbance can cut it into pieces; if they kept one id, size would be inflated. So
# after each step the reef is relabelled: the largest piece keeps the id, each other
# piece becomes a new colony (inheriting species/traits, sized by its own blob) that
# ordinary combat can then erode.
# Steps: labelConnectedComponents -> splitColonies -> flagIsolatedFragments ->
#        consolidateColonies -> buildColonySpeciesFromStates.


# ---- Tuning numbers ----

# A walled-in fragment this small (cells) is flagged $isolated (informational only).
ISOLATED_FRAGMENT_MAX <- 5

# A piece this small (cells) embedded in a larger colony is absorbed into it during
# consolidation (the decisive island cleanup). 0 = off.
ISLAND_ABSORB_MAX <- 5


# Label every 4-connected same-id blob. Returns an integer matrix (empty cells = 0);
# flood-fill on linear indices with a preallocated stack, O(cells).
labelConnectedComponents <- function(reef) {
  nr <- nrow(reef)
  nc <- ncol(reef)
  n  <- nr * nc

  reef_vec   <- as.vector(reef)     # column-major
  labels     <- integer(n)          # 0 = unlabelled / empty
  next_label <- 0L
  stack      <- integer(n)

  for (start in seq_len(n)) {
    if (is.na(reef_vec[start]) || labels[start] != 0L) {
      next
    }

    next_label <- next_label + 1L
    id <- reef_vec[start]
    sp <- 1L
    stack[sp]      <- start
    labels[start]  <- next_label

    while (sp > 0L) {   # flood fill from this seed
      cur <- stack[sp]
      sp  <- sp - 1L

      r <- ((cur - 1L) %% nr) + 1L
      c <- ((cur - 1L) %/% nr) + 1L

      # up / down / left / right as linear-index offsets
      if (r > 1L)  { nb <- cur - 1L;  if (labels[nb] == 0L && !is.na(reef_vec[nb]) && reef_vec[nb] == id) { labels[nb] <- next_label; sp <- sp + 1L; stack[sp] <- nb } }
      if (r < nr)  { nb <- cur + 1L;  if (labels[nb] == 0L && !is.na(reef_vec[nb]) && reef_vec[nb] == id) { labels[nb] <- next_label; sp <- sp + 1L; stack[sp] <- nb } }
      if (c > 1L)  { nb <- cur - nr;  if (labels[nb] == 0L && !is.na(reef_vec[nb]) && reef_vec[nb] == id) { labels[nb] <- next_label; sp <- sp + 1L; stack[sp] <- nb } }
      if (c < nc)  { nb <- cur + nr;  if (labels[nb] == 0L && !is.na(reef_vec[nb]) && reef_vec[nb] == id) { labels[nb] <- next_label; sp <- sp + 1L; stack[sp] <- nb } }
    }
  }

  matrix(labels, nrow = nr, ncol = nc)
}


# Split any colony spanning several components: the largest piece keeps the id, each
# other becomes a new colony (inheriting species/traits). fragment_counter gives
# unique ids. A caller-supplied label matrix is reused when the reef is unchanged.
splitColonies <- function(reef, corals, colony_species, fragment_counter,
                          labels = NULL) {
  nr       <- nrow(reef)
  if (is.null(labels)) labels <- labelConnectedComponents(reef)
  reef_vec <- as.vector(reef)
  lab_vec  <- as.vector(labels)

  ids <- vapply(corals, function(x) x$id, character(1))

  present_ids <- unique(reef_vec[!is.na(reef_vec)])

  for (id in present_ids) {
    cell_idx <- which(reef_vec == id)
    labs     <- lab_vec[cell_idx]
    ulabs    <- unique(labs)

    if (length(ulabs) <= 1L) {   # already one blob
      next
    }

    counts     <- tabulate(match(labs, ulabs))
    keep_label <- ulabs[which.max(counts)]   # largest piece keeps the id

    parent_pos <- which(ids == id)[1]
    parent     <- corals[[parent_pos]]

    for (lab in ulabs[ulabs != keep_label]) {   # each other piece -> new colony
      fragment_counter <- fragment_counter + 1L
      new_id   <- paste0(id, ".f", fragment_counter)
      frag_idx <- cell_idx[labs == lab]

      reef_vec[frag_idx] <- new_id

      frag <- parent
      frag$id        <- new_id
      frag$parent_id <- id
      frag$origin_id <- if (!is.null(parent$origin_id)) parent$origin_id else id   # trace to setup
      frag$size      <- length(frag_idx)

      # centroid, stored as coords c(x = col, y = row)
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


# Flag tiny fragments walled in by other colonies ($isolated = TRUE). Informational
# only - combat is untouched.
flagIsolatedFragments <- function(reef, corals) {
  nr       <- nrow(reef)
  nc       <- ncol(reef)
  reef_vec <- as.vector(reef)

  for (k in seq_along(corals)) {
    co <- corals[[k]]
    corals[[k]]$isolated <- FALSE

    if (is.null(co$parent_id)) {   # only fragments can be isolated islands
      next
    }

    # recount from the reef so the test never uses a stale size; only small live pieces
    idx     <- which(reef_vec == co$id)
    n_cells <- length(idx)
    if (n_cells == 0L || n_cells >= ISOLATED_FRAGMENT_MAX) {
      next
    }

    enclosed <- TRUE

    for (cur in idx) {
      r <- ((cur - 1L) %% nr) + 1L
      c <- ((cur - 1L) %/% nr) + 1L

      if (r == 1L || r == nr || c == 1L || c == nc) {   # on the border = not enclosed
        enclosed <- FALSE
        break
      }

      nbrs <- c(cur - 1L, cur + 1L, cur - nr, cur + nr)
      if (any(is.na(reef_vec[nbrs]))) {   # an empty neighbour = open, not walled in
        enclosed <- FALSE
        break
      }
    }

    corals[[k]]$isolated <- enclosed
  }

  corals
}


# Absorb small islands embedded in a larger colony. For each piece <= ISLAND_ABSORB_MAX
# cells, if one other colony borders it at least as much as empty water does and is at
# least as big, the piece is reassigned to that colony. Pieces mostly in open water are
# left (real small colonies). Returns the updated reef.
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

    # tally cells bordering this piece: empties vs each neighbouring colony
    empty_count <- 0L
    nbr_ids     <- character(0)
    for (cur in comp_idx) {
      r <- ((cur - 1L) %% nr) + 1L
      c <- ((cur - 1L) %/% nr) + 1L
      nbrs <- c(if (r > 1L) cur - 1L, if (r < nr) cur + 1L,
                if (c > 1L) cur - nr, if (c < nc) cur + nr)
      for (nb in nbrs) {
        if (lab_vec[nb] == lab) {
          next
        }
        if (is.na(reef_vec[nb])) {
          empty_count <- empty_count + 1L
        } else {
          nbr_ids <- c(nbr_ids, reef_vec[nb])
        }
      }
    }

    if (length(nbr_ids) == 0L) {   # floating in open water
      next
    }

    counts   <- sort(table(nbr_ids), decreasing = TRUE)
    dominant <- names(counts)[1]
    dom_bord <- as.integer(counts[1])

    if (dom_bord < empty_count) {   # not genuinely embedded
      next
    }

    dom_size <- sum(reef_vec == dominant, na.rm = TRUE)
    if (dom_size < length(comp_idx)) {   # surrounder must be at least as large
      next
    }

    reef_vec[comp_idx] <- dominant   # absorb
    changed <- TRUE
  }

  # return the labels too: if nothing was absorbed they are still valid for reuse
  list(reef    = matrix(reef_vec, nrow = nr, ncol = nc),
       labels  = if (changed) NULL else labels,
       changed = changed)
}


# One consolidation pass per step: absorb islands, split disconnected colonies,
# refresh sizes, drop fully-eroded fragments (keeping founder colonies so extinction
# is recorded), then flag tiny islands.
consolidateColonies <- function(reef, corals, colony_species, fragment_counter) {

  # absorb embedded islands first; reuse its labels in splitColonies if unchanged
  absorbed <- absorbIslands(reef, corals, colony_species)
  reef     <- absorbed$reef

  res <- splitColonies(reef, corals, colony_species, fragment_counter,
                       labels = absorbed$labels)

  res$corals <- updateSize(res$corals, res$reef)

  # drop dead fragments (0 cells) but keep founder colonies (no parent)
  keep <- vapply(res$corals, function(x) x$size > 0 || is.null(x$parent_id), logical(1))
  res$corals <- res$corals[keep]

  res$corals <- flagIsolatedFragments(res$reef, res$corals)

  res
}


# Rebuild a complete id -> species lookup after a run (the initial map misses mid-run
# fragment ids), gathering every colony that ever existed across all saved states.
buildColonySpeciesFromStates <- function(states, base = NULL) {
  all_corals <- unlist(lapply(states, `[[`, "corals"), recursive = FALSE)

  species_map <- setNames(
    vapply(all_corals, function(x) x$species, character(1)),
    vapply(all_corals, function(x) x$id,      character(1))
  )
  species_map <- species_map[!duplicated(names(species_map))]

  # fold in any base ids that never appeared in a saved state
  if (!is.null(base)) {
    extra <- base[!(names(base) %in% names(species_map))]
    species_map <- c(species_map, extra)
  }

  species_map
}
