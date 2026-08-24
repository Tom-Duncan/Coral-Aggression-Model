# Coral reproduction (larval recruitment). When on, a mature colony can each
# timestep spawn a single-cell recruit of the same species onto a random empty cell.
# Maturity is age-based (MATURITY_AGE, or a per-species maturity_age column); once
# mature, spawn chance scales with size:
#   effective chance = min(1, repro_chance * cover% / repro_ref_percent)
# After spawning a colony waits REPRO_COOLDOWN steps. Off = no colony reproduces.
# Steps: getReproductionSetup -> assignReproductionToSpecies ->
#        colonyReproductionChance -> spawnRecruit -> reproduceColonies.


# ---- Tuning numbers ---------------------------------------------------------

# Steps a colony must be alive before it is mature (per-species maturity_age overrides).
MATURITY_AGE <- 30

# Steps a colony waits after spawning before it can spawn again.
REPRO_COOLDOWN <- 30

# Spawning event interval: colonies spawn only on multiples of this (1 = every step).
REPRO_INTERVAL <- 1

# Most transparent a recruit shade can be (alpha; recruits get a random shade).
REPRO_SHADE_ALPHA_MIN <- 0.45

# Random-setup: chance the feature is on, and chance a given species can reproduce.
REPRO_RANDOM_ENABLE_CHANCE <- 0.5
REPRO_SPECIES_CHANCE <- 0.6

# Random-setup ranges: base spawn chance, and reference size (% reef) at base chance.
REPRO_CHANCE_MIN <- 0.02
REPRO_CHANCE_MAX <- 0.15
REPRO_REF_PERCENT_MIN <- 1
REPRO_REF_PERCENT_MAX <- 15


# ---- Background recruitment (natural-reef habitat only) ---------------------
# Larvae also arrive from outside the patch: each species has an equal, independent
# chance to drop a recruit on an empty cell, so a locally extinct species can return.
# Aquarium/"none" gets no external supply. supply_mode sets the temporal pattern:
#   "pulse"      - every BACKGROUND_REPRO_INTERVAL steps, chance = background_chance
#   "continuous" - every step, chance = background_chance / interval (matched supply)

BACKGROUND_REPRO_INTERVAL <- 30
DEFAULT_SUPPLY_MODE <- "pulse"
BACKGROUND_REPRO_RANDOM_MIN <- 0.05   # random-setup range
BACKGROUND_REPRO_RANDOM_MAX <- 0.30


# Ask during setup whether reproduction is on (random on setup_mode 3, else prompt).
getReproductionSetup <- function(setup_mode) {
  if (setup_mode == 3) {
    return(runif(1) < REPRO_RANDOM_ENABLE_CHANCE)
  }
  getMenuChoice(
    "Include coral reproduction (colonies spawn new recruits)?",
    c("On", "Off")
  ) == 1
}


# Add reproduction / repro_chance / repro_ref_percent columns to species_traits
# (off = defaults; "manual" = prompt per species; "random" = random traits).
assignReproductionToSpecies <- function(species_traits, mode, enabled) {
  n <- nrow(species_traits)

  reproduction      <- logical(n)
  repro_chance      <- numeric(n)
  repro_ref_percent <- rep(NA_real_, n)

  if (!isTRUE(enabled)) {   # off: fill "no reproduction" defaults
    species_traits$reproduction      <- rep(FALSE, n)
    species_traits$repro_chance      <- rep(0, n)
    species_traits$repro_ref_percent <- rep(NA_real_, n)
    return(species_traits)
  }

  for (i in seq_len(n)) {
    sp <- species_traits$species[i]

    if (mode == "manual") {
      reproduction[i] <- getYesNo(paste0(sp, " able to reproduce?"))
      if (reproduction[i]) {
        cat("  ", sp, " reproduction: a colony's spawn chance scales with its size.\n", sep = "")
        cat("    You give a reference size and the spawn chance a colony that size has;\n")
        cat("    bigger colonies then get proportionally more chance, smaller ones less.\n")

        ref <- getValidInteger(   # reference size (% reef)
          paste0("    Reference size for ", sp, " (% of reef a colony covers, 1-100): "),
          min_val = 1, max_val = 100
        )
        pct <- getValidInteger(   # spawn chance at that size
          paste0("    Spawn chance per timestep for a colony at ", ref, "% cover (%, 0-100): "),
          min_val = 0, max_val = 100
        )

        repro_ref_percent[i] <- ref
        repro_chance[i]      <- pct / 100
      }

    } else {
      reproduction[i] <- runif(1) < REPRO_SPECIES_CHANCE
      if (reproduction[i]) {
        repro_chance[i]      <- runif(1, REPRO_CHANCE_MIN, REPRO_CHANCE_MAX)
        repro_ref_percent[i] <- sample(REPRO_REF_PERCENT_MIN:REPRO_REF_PERCENT_MAX, 1)
      }
    }
  }

  species_traits$reproduction      <- reproduction
  species_traits$repro_chance      <- repro_chance
  species_traits$repro_ref_percent <- repro_ref_percent
  species_traits
}


# The maturity age for a colony: per-species maturity_age column if present, else global.
maturityAge <- function(coral, species_traits = NULL) {
  if (!is.null(species_traits) && "maturity_age" %in% names(species_traits)) {
    i <- match(coral$species, species_traits$species)
    if (!is.na(i) && !is.na(species_traits$maturity_age[i]))
      return(species_traits$maturity_age[i])
  }
  MATURITY_AGE
}

# Is a colony old enough to reproduce? (age-based; missing age counts as 0)
isMature <- function(coral, species_traits = NULL) {
  age <- if (is.null(coral$age)) 0L else coral$age
  age >= maturityAge(coral, species_traits)
}


# A mature colony's spawn chance this step: base chance scaled by cover / reference,
# capped at 1. Returns 0 if the species cannot reproduce. (Maturity gated elsewhere.)
colonyReproductionChance <- function(coral, species_traits, reef) {
  sp_row <- which(species_traits$species == coral$species)
  if (length(sp_row) == 0 || !isTRUE(species_traits$reproduction[sp_row])) {
    return(0)
  }

  ref_percent <- species_traits$repro_ref_percent[sp_row]
  if (is.na(ref_percent) || ref_percent <= 0) {
    return(0)
  }

  cover_percent <- coralCoverPercent(coral, reef)

  base <- species_traits$repro_chance[sp_row]
  min(1, base * cover_percent / ref_percent)
}


# A recruit's shade: the parent's solid species colour + a fresh random transparency.
recruitShade <- function(parent_colour) {
  cc   <- col2rgb(parent_colour)   # solid RGB, ignoring alpha
  base <- rgb(cc[1, 1], cc[2, 1], cc[3, 1], maxColorValue = 255)
  adjustcolor(base, alpha.f = runif(1, REPRO_SHADE_ALPHA_MIN, 1))
}


# Place one recruit for a parent colony: a new single-cell colony inheriting the
# parent's species/traits but with its own id, parent_id, and recruit flag.
spawnRecruit <- function(reef, corals, colony_species, parent, cell_index, recruit_id) {
  nr <- nrow(reef)
  row <- ((cell_index - 1L) %% nr) + 1L
  col <- ((cell_index - 1L) %/% nr) + 1L

  reef[cell_index] <- recruit_id

  recruit <- parent
  recruit$id            <- recruit_id
  recruit$parent_id     <- parent$id
  recruit$origin_id     <- if (!is.null(parent$origin_id)) parent$origin_id else parent$id
  recruit$recruit       <- TRUE   # a reproduction recruit, not a split
  recruit$isolated      <- FALSE
  recruit$size          <- 1L
  recruit$age           <- 0L     # newborn: must mature before spawning
  recruit$repro_cooldown <- 0
  recruit$colour        <- recruitShade(parent$colour)
  recruit$coords        <- c(col, row)  # c(x, y)

  corals[[length(corals) + 1L]] <- recruit
  colony_species[recruit_id]    <- parent$species

  list(reef = reef, corals = corals, colony_species = colony_species)
}


# One reproduction pass per timestep: count down cooldowns; each mature, off-cooldown
# colony rolls its size-scaled chance and, on success, settles a recruit on a random
# empty cell and restarts its cooldown. Colonies are iterated in random order so none
# has priority for scarce empty cells; recruits made this step are not iterated.
reproduceColonies <- function(reef, corals, colony_species, species_traits, recruit_counter,
                              timestep_num = NA) {

  # synchronised-spawning gate (REPRO_INTERVAL > 1 = only on event steps)
  if (REPRO_INTERVAL > 1 &&
      (is.na(timestep_num) || timestep_num < REPRO_INTERVAL ||
       timestep_num %% REPRO_INTERVAL != 0)) {
    return(list(reef = reef, corals = corals,
                colony_species = colony_species, recruit_counter = recruit_counter))
  }

  if (!("reproduction" %in% names(species_traits)) || !any(species_traits$reproduction)) {
    return(list(reef = reef, corals = corals,
                colony_species = colony_species, recruit_counter = recruit_counter))
  }

  empties <- NULL   # found lazily the first time one is needed

  for (i in sample(seq_along(corals))) {
    co <- corals[[i]]

    # count down cooldown; skip if still waiting
    cd <- if (is.null(co$repro_cooldown)) 0 else co$repro_cooldown
    if (cd > 0) {
      corals[[i]]$repro_cooldown <- cd - 1L
      next
    }
    corals[[i]]$repro_cooldown <- 0

    if (co$size <= 0) {   # dead
      next
    }

    if (!isMature(co, species_traits)) {
      next
    }

    if (runif(1) >= colonyReproductionChance(co, species_traits, reef)) {
      next
    }

    if (is.null(empties)) {
      empties <- which(is.na(as.vector(reef)))
    }
    if (length(empties) == 0) {
      next   # reef full
    }

    pick <- sample(length(empties), 1)
    cell <- empties[pick]
    empties <- empties[-pick]

    recruit_counter <- recruit_counter + 1L
    recruit_id <- paste0(co$species, "_r", recruit_counter)

    placed <- spawnRecruit(reef, corals, colony_species, co, cell, recruit_id)
    reef           <- placed$reef
    corals         <- placed$corals
    colony_species <- placed$colony_species

    corals[[i]]$repro_cooldown <- REPRO_COOLDOWN
  }

  list(reef = reef, corals = corals,
       colony_species = colony_species, recruit_counter = recruit_counter)
}


# Ask which habitat this run is (controls background recruitment). Aquarium/none =
# no external supply; reef = supply on (also asks the per-species chance + pattern).
getHabitatSetup <- function(setup_mode) {

  if (setup_mode == 3) {   # fully random
    type <- sample(c("aquarium", "reef", "none"), 1)
    if (type == "reef") {
      return(reefHabitat(runif(1, BACKGROUND_REPRO_RANDOM_MIN, BACKGROUND_REPRO_RANDOM_MAX),
                         sample(c("pulse", "continuous"), 1)))
    }
    return(noBackground(type))
  }

  choice <- getMenuChoice(
    "Habitat type (natural reef adds background recruitment from outside)?",
    c("Aquarium", "Natural reef", "None")
  )
  type <- c("aquarium", "reef", "none")[choice]

  if (type == "reef") {
    cat("  Natural reef: each species may spawn a background recruit anywhere on the\n",
        "  reef (larvae arriving from outside), so a locally extinct species can return.\n", sep = "")
    pct <- getValidInteger(
      paste0("  Background recruitment chance per species (per pulse of ",
             BACKGROUND_REPRO_INTERVAL, " steps, %, 0-100): "),
      min_val = 0, max_val = 100
    )
    smode_choice <- getMenuChoice(
      "Larval supply pattern (same expected supply, different timing)?",
      c(paste0("Pulse (a periodic event every ", BACKGROUND_REPRO_INTERVAL, " steps)"),
        "Continuous (a stochastic trickle every step)")
    )
    smode <- c("pulse", "continuous")[smode_choice]
    return(reefHabitat(pct / 100, smode))
  }

  noBackground(type)
}


# Habitat config constructors (open reef = external supply; aquarium/none = closed).
noBackground <- function(type) {
  list(type = type, background_enabled = FALSE, background_chance = NA_real_,
       background_interval = BACKGROUND_REPRO_INTERVAL, supply_mode = NA_character_)
}
reefHabitat <- function(chance, supply_mode = DEFAULT_SUPPLY_MODE,
                        interval = BACKGROUND_REPRO_INTERVAL,
                        local_retention = 0) {
  list(type = "reef", background_enabled = TRUE, background_chance = chance,
       background_interval = interval, supply_mode = supply_mode,
       local_retention = local_retention)
}

# Non-interactive habitat builder for batch sweeps. local_retention adds self-seeding:
# a species' arrival chance rises by local_retention x (its mature cover fraction).
buildHabitatConfig <- function(habitat = "aquarium", chance = 0.1,
                               supply_mode = DEFAULT_SUPPLY_MODE,
                               interval = BACKGROUND_REPRO_INTERVAL,
                               local_retention = 0) {
  if (identical(habitat, "reef")) reefHabitat(chance, supply_mode, interval, local_retention)
  else                            noBackground(habitat)
}


# Each species' solid base colour (matches createCorals), so a background recruit can
# take its species colour even with no living colony to copy.
speciesColourMap <- function(species_traits) {
  setNames(hcl.colors(nrow(species_traits), "Dark 3"), species_traits$species)
}


# One-line habitat description for the setup summary.
describeHabitat <- function(habitat) {
  if (is.null(habitat)) return("None")
  if (isTRUE(habitat$background_enabled)) {
    mode <- if (is.null(habitat$supply_mode) || is.na(habitat$supply_mode)) "pulse"
            else habitat$supply_mode
    pattern <- if (identical(mode, "continuous"))
                 paste0("continuous trickle, ", round(100 * habitat$background_chance),
                        "% per species per ", habitat$background_interval, " steps")
               else
                 paste0("pulse every ", habitat$background_interval, " steps, ",
                        round(100 * habitat$background_chance), "% per species")
    return(paste0("Natural reef (background recruitment: ", pattern, ")"))
  }
  switch(habitat$type,
         aquarium = "Aquarium (no background recruitment)",
         none     = "None (no background recruitment)",
         "None")
}


# Place one background recruit: no parent (arrives from outside), so built from the
# species' trait row + base colour. $background = TRUE marks its larval-influx origin.
spawnBackgroundRecruit <- function(reef, corals, colony_species, trait_row,
                                   species_colour, cell_index, recruit_id) {
  nr  <- nrow(reef)
  row <- ((cell_index - 1L) %% nr) + 1L
  col <- ((cell_index - 1L) %/% nr) + 1L

  reef[cell_index] <- recruit_id

  recruit <- list(
    id       = recruit_id,
    species  = trait_row$species,
    growth   = trait_row$growth,
    colour   = recruitShade(species_colour),
    coords   = c(col, row),
    size     = 1L,
    age      = 0L,       # newborn: must mature before spawning
    recruit    = TRUE,
    background = TRUE,   # arrived from outside the reef
    isolated   = FALSE,
    parent_id  = NA_character_,
    origin_id  = recruit_id,
    repro_cooldown = 0
  )

  corals[[length(corals) + 1L]] <- recruit
  colony_species[recruit_id]    <- trait_row$species

  list(reef = reef, corals = corals, colony_species = colony_species)
}


# One background-recruitment pass per timestep (natural reef only). supply_mode sets
# whether it fires on interval pulses (chance) or every step (chance/interval). Every
# species rolls the same chance in random order, drawing from the whole pool.
backgroundReproduce <- function(reef, corals, colony_species, species_traits,
                                habitat, recruit_counter, timestep_num) {

  unchanged <- list(reef = reef, corals = corals,
                    colony_species = colony_species, recruit_counter = recruit_counter)

  if (is.null(habitat) || !isTRUE(habitat$background_enabled)) {
    return(unchanged)
  }

  chance   <- habitat$background_chance
  interval <- habitat$background_interval
  if (is.na(chance) || chance <= 0) {
    return(unchanged)
  }

  # per-species chance this step from the supply mode (default pulse for older configs)
  mode <- if (is.null(habitat$supply_mode)) "pulse" else habitat$supply_mode
  if (identical(mode, "continuous")) {
    step_chance <- chance / interval   # every step, matched expected supply
  } else {
    # pulse: only on the interval, after the first (t = 30, 60, ...)
    if (is.na(timestep_num) || timestep_num < interval || timestep_num %% interval != 0) {
      return(unchanged)
    }
    step_chance <- chance
  }

  colours <- speciesColourMap(species_traits)
  empties <- which(is.na(as.vector(reef)))

  # local retention (optional): raise each species' chance by its mature cover fraction
  # x local_retention (self-seeding). Absent species keep the full base chance.
  n_sp <- nrow(species_traits)
  sp_chance <- rep(step_chance, n_sp)
  lr <- habitat$local_retention
  if (!is.null(lr) && !is.na(lr) && lr > 0 && length(corals) > 0) {
    total_cells <- length(reef)
    mature_cells <- setNames(numeric(n_sp), species_traits$species)
    for (cc in corals) {
      if (is.null(cc$size) || cc$size <= 0) next
      if ((if (is.null(cc$age)) 0L else cc$age) >= maturityAge(cc, species_traits))
        mature_cells[cc$species] <- mature_cells[cc$species] + cc$size
    }
    boost <- lr * (mature_cells / total_cells)          # <= lr even at full cover
    if (identical(mode, "continuous")) boost <- boost / interval
    sp_chance <- pmin(1, sp_chance + as.numeric(boost[species_traits$species]))
  }

  # every species rolls independently, in random order (no priority for scarce cells)
  for (i in sample(seq_len(n_sp))) {
    if (length(empties) == 0) break
    if (runif(1) >= sp_chance[i]) next

    pick <- sample(length(empties), 1)
    cell <- empties[pick]
    empties <- empties[-pick]

    sp <- species_traits$species[i]
    recruit_counter <- recruit_counter + 1L
    recruit_id <- paste0(sp, "_b", recruit_counter)

    placed <- spawnBackgroundRecruit(reef, corals, colony_species,
                                     species_traits[i, ], colours[[sp]], cell, recruit_id)
    reef           <- placed$reef
    corals         <- placed$corals
    colony_species <- placed$colony_species
  }

  list(reef = reef, corals = corals,
       colony_species = colony_species, recruit_counter = recruit_counter)
}
