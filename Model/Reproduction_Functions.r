#Functions related to coral reproduction (larval recruitment)
#
#When reproduction is switched ON, a mature colony can, each timestep, spawn a
#new single-cell colony (a "recruit") of the SAME species onto a random empty
#cell anywhere on the reef. A recruit is a brand-new, separately tracked colony
#that then grows and competes like any other.
#
#A colony must first be MATURE to reproduce at all: maturity is age-based, reached
#once a colony has been alive for MATURITY_AGE timesteps (age is tracked per colony
#on $age, incremented each timestep in timestep()). Size no longer determines
#maturity. Once mature, a colony's spawning CHANCE still scales continuously with its
#size, so bigger mature colonies reproduce more and tiny mature ones almost never do.
#Whether/how a colony reproduces depends on per-species traits:
#   reproduction     - can this species reproduce at all?
#   repro_ref_percent- the reference size, as a PERCENT OF THE REEF, at which the
#                      base chance applies. It is a scaling reference, NOT a cutoff:
#                      a colony this big reproduces at repro_chance, smaller ones
#                      proportionally less, larger ones proportionally more. A low
#                      value means the species reproduces readily even when small; a
#                      high value means it only reproduces much once large.
#   repro_chance     - base per-timestep probability at the reference size.
#Percent cover is measured with coralCoverPercent() (Size_Impact.r), the SAME helper
#the small/medium/large size categories use, so everything agrees and stays cheap:
#   effective chance = min(1, repro_chance * cover% / repro_ref_percent)
#
#After a colony reproduces it must WAIT REPRO_COOLDOWN timesteps before it can roll
#again, so it cannot spawn every single step. The wait is tracked per colony on
#$repro_cooldown (timesteps remaining), counted down each step.
#
#When OFF, nothing changes - no colony reproduces.
#
#Everything is split into one small function per step, matching the rest of the
#model:
#   1. Tuning numbers            - rates and ranges, all in one place
#   2. getReproductionSetup      - ask during setup whether the feature is on
#   3. assignReproductionToSpecies - give each species its reproduction traits
#   4. colonyReproductionChance  - a colony's size-scaled spawning probability
#   5. spawnRecruit              - place one recruit on an empty cell
#   6. reproduceColonies         - the full per-timestep pass


# ====================================================================
#  TUNING NUMBERS  - change these to alter how reproduction behaves.
# ====================================================================

#Default timesteps a colony must be alive before it is MATURE and can reproduce at
#all. Maturity is age-based (not size-based): a colony younger than this cannot
#spawn, regardless of how large it is. This is the GLOBAL fallback; a per-species
#`maturity_age` column on species_traits overrides it per species (see maturityAge).
MATURITY_AGE <- 30

#Timesteps a colony must wait after reproducing before it can reproduce again.
#Change this ONE number to make reproduction more or less frequent per colony.
REPRO_COOLDOWN <- 30

#Reproduction EVENT interval: colony spawning only occurs on timesteps that are
#multiples of this (synchronised spawning events, like periodic mass spawning).
#1 = every timestep (the original behaviour; per-colony chance + cooldown govern).
REPRO_INTERVAL <- 1

#Recruits are drawn a random shade of the species colour (via transparency) so they
#are visibly distinct from their parent, exactly like multiple colonies of a species.
#This is the most transparent a recruit shade can be (alpha, 0 = clear, 1 = solid).
REPRO_SHADE_ALPHA_MIN <- 0.45

#In fully-random setup, the chance the whole feature is switched on.
REPRO_RANDOM_ENABLE_CHANCE <- 0.5

#When traits are randomised, the chance a given species can reproduce at all.
REPRO_SPECIES_CHANCE <- 0.6

#When randomised, the base per-timestep spawning probability is drawn from here.
REPRO_CHANCE_MIN <- 0.02
REPRO_CHANCE_MAX <- 0.15

#When randomised, the reference size - the PERCENT of the reef at which a colony
#reproduces at its base chance - is drawn from this range. It is a scaling point,
#not a cutoff. For reference, the size categories treat <5% as small and >15% as
#large, so a low value spans "reproduces readily while small" up to "only much once
#large".
REPRO_REF_PERCENT_MIN <- 1
REPRO_REF_PERCENT_MAX <- 15


# ====================================================================
#  BACKGROUND RECRUITMENT  (natural-reef habitat only)
# --------------------------------------------------------------------
#  Separate from colony reproduction above. On a NATURAL REEF, larvae also arrive
#  from OUTSIDE the modelled patch. Every species in the pool has an equal,
#  independent chance to drop one fresh recruit on a random empty cell - drawn from
#  the species pool that was created for the run, so a species that has gone locally
#  extinct can return. This runs ALONGSIDE colony reproduction (which still follows
#  the on/off reproduction toggle); an AQUARIUM or "none" habitat gets no external
#  supply. The TEMPORAL PATTERN of that supply is togglable (habitat$supply_mode):
#    "pulse"      - a periodic pulse: every BACKGROUND_REPRO_INTERVAL steps (after
#                   the first interval) each species has probability
#                   background_chance to settle one recruit. (The original mode.)
#    "continuous" - a stochastic trickle: EVERY step each species has probability
#                   background_chance / BACKGROUND_REPRO_INTERVAL. This holds the
#                   EXPECTED supply per species identical to the pulse mode, so the
#                   toggle isolates the temporal pattern of larval delivery, not its
#                   magnitude. In both modes the supply is species-symmetric (equal
#                   per-species chance), so no species is favoured by immigration.
# ====================================================================

#Pulse period (in timesteps): the pulse fires every interval; the continuous mode
#uses it to set the matched per-step rate (chance / interval).
BACKGROUND_REPRO_INTERVAL <- 30

#Default larval-supply temporal pattern for a natural reef ("pulse" or "continuous").
DEFAULT_SUPPLY_MODE <- "pulse"

#In fully-random setup only: the background chance per species is drawn from here.
BACKGROUND_REPRO_RANDOM_MIN <- 0.05
BACKGROUND_REPRO_RANDOM_MAX <- 0.30


# --------------------------------------------------------------------
#  STEP 2: Ask during setup whether reproduction is included.
#  Follows the overall setup mode, like disturbances:
#    Fully random (3) - switched on or off at random
#    Manual (1) / Semi-manual (2) - ask the user (1 = On, 2 = Off)
#  Returns TRUE if reproduction should be enabled.
# --------------------------------------------------------------------
getReproductionSetup <- function(setup_mode) {
  if (setup_mode == 3) {
    return(runif(1) < REPRO_RANDOM_ENABLE_CHANCE)
  }
  getMenuChoice(
    "Include coral reproduction (colonies spawn new recruits)?",
    c("On", "Off")
  ) == 1
}


# --------------------------------------------------------------------
#  STEP 3: Give every species its reproduction traits.
#  Adds three columns to species_traits: reproduction (logical),
#  repro_chance (base per-timestep probability) and repro_ref_percent (the percent
#  of reef cover at which that base chance applies - a scaling reference, not a
#  cutoff; uses the same cover measure as the size categories).
#    enabled = FALSE            - feature off, no species reproduces
#    enabled = TRUE, "manual"   - ask per species (can it? chance? reference %?)
#    enabled = TRUE, "random"   - assign each species random reproduction traits
#  Returns the updated species_traits.
# --------------------------------------------------------------------
assignReproductionToSpecies <- function(species_traits, mode, enabled) {
  n <- nrow(species_traits)

  reproduction      <- logical(n)
  repro_chance      <- numeric(n)
  repro_ref_percent <- rep(NA_real_, n)

  #Feature off: fill in "no reproduction" defaults so downstream code is uniform
  if (!isTRUE(enabled)) {
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
        #Explain the two-part model before asking, since the two numbers are easy
        #to confuse: one is a SIZE, the other is a PROBABILITY at that size.
        cat("  ", sp, " reproduction: a colony's spawn chance scales with its size.\n", sep = "")
        cat("    You give a reference size and the spawn chance a colony that size has;\n")
        cat("    bigger colonies then get proportionally more chance, smaller ones less.\n")

        #(a) the SIZE, as a percent of the reef a colony covers (like small/large)
        ref <- getValidInteger(
          paste0("    Reference size for ", sp, " (% of reef a colony covers, 1-100): "),
          min_val = 1, max_val = 100
        )
        #(b) the PROBABILITY per timestep when a colony is at that reference size
        pct <- getValidInteger(
          paste0("    Spawn chance per timestep for a colony at ", ref, "% cover (%, 0-100): "),
          min_val = 0, max_val = 100
        )

        repro_ref_percent[i] <- ref
        repro_chance[i]      <- pct / 100
      }

    } else {
      #Random assignment
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


# --------------------------------------------------------------------
#  The maturity age that applies to a given colony. PER-SPECIES if species_traits
#  carries a `maturity_age` column (a value per species, so species can mature at
#  different ages); otherwise the global MATURITY_AGE. With no such column - as in
#  every run made before this feature - it returns MATURITY_AGE unchanged, so old
#  behaviour is preserved exactly.
# --------------------------------------------------------------------
maturityAge <- function(coral, species_traits = NULL) {
  if (!is.null(species_traits) && "maturity_age" %in% names(species_traits)) {
    i <- match(coral$species, species_traits$species)
    if (!is.na(i) && !is.na(species_traits$maturity_age[i]))
      return(species_traits$maturity_age[i])
  }
  MATURITY_AGE
}

# --------------------------------------------------------------------
#  Is a colony old enough to reproduce? Maturity is purely age-based: a colony is
#  mature once it has been alive for at least its species' maturity age (see
#  maturityAge()). Age is tracked on $age (incremented each timestep in timestep());
#  a missing age counts as 0.
# --------------------------------------------------------------------
isMature <- function(coral, species_traits = NULL) {
  age <- if (is.null(coral$age)) 0L else coral$age
  age >= maturityAge(coral, species_traits)
}


# --------------------------------------------------------------------
#  STEP 4: A colony's actual chance of spawning this timestep.
#  Assumes the colony is already MATURE (the maturity gate is applied separately in
#  reproduceColonies). The chance then scales CONTINUOUSLY with size: a colony at its
#  species' reference cover reproduces at the base chance, a bigger one more, a
#  smaller one proportionally less (so tiny colonies almost never spawn). Percent
#  cover comes from the shared coralCoverPercent() helper. Returns 0 only if the
#  species cannot reproduce at all. Capped at 1.
# --------------------------------------------------------------------
colonyReproductionChance <- function(coral, species_traits, reef) {
  sp_row <- which(species_traits$species == coral$species)
  if (length(sp_row) == 0 || !isTRUE(species_traits$reproduction[sp_row])) {
    return(0)
  }

  ref_percent <- species_traits$repro_ref_percent[sp_row]
  if (is.na(ref_percent) || ref_percent <= 0) {
    return(0)
  }

  #Percent of reef covered, from the same helper the size categories use
  cover_percent <- coralCoverPercent(coral, reef)

  base <- species_traits$repro_chance[sp_row]
  min(1, base * cover_percent / ref_percent)
}


# --------------------------------------------------------------------
#  Give a recruit its own shade of the species colour.
#  Strips any transparency from the parent's colour to recover the solid species
#  colour, then applies a fresh random transparency - the same way createCorals
#  shades multiple colonies of one species. So a recruit is the same hue but a
#  slightly different shade, never identical to its parent and never fading further
#  each generation (it always starts from the solid colour).
# --------------------------------------------------------------------
recruitShade <- function(parent_colour) {
  cc   <- col2rgb(parent_colour)   #ignores any alpha, returns the solid RGB
  base <- rgb(cc[1, 1], cc[2, 1], cc[3, 1], maxColorValue = 255)
  adjustcolor(base, alpha.f = runif(1, REPRO_SHADE_ALPHA_MIN, 1))
}


# --------------------------------------------------------------------
#  STEP 5: Place one recruit for a parent colony on a given empty cell.
#  Builds a new single-cell colony that inherits the parent's species and traits
#  but is its own tracked individual (new id, parent_id for lineage, recruit flag
#  so it is not mistaken for a split fragment). Updates the reef and colony_species.
#  Returns a list with the updated reef, corals, colony_species and new_id.
# --------------------------------------------------------------------
spawnRecruit <- function(reef, corals, colony_species, parent, cell_index, recruit_id) {
  nr <- nrow(reef)
  row <- ((cell_index - 1L) %% nr) + 1L
  col <- ((cell_index - 1L) %/% nr) + 1L

  reef[cell_index] <- recruit_id

  recruit <- parent
  recruit$id            <- recruit_id
  recruit$parent_id     <- parent$id
  recruit$origin_id     <- if (!is.null(parent$origin_id)) parent$origin_id else parent$id
  recruit$recruit       <- TRUE   #marks this as a reproduction recruit, not a split
  recruit$isolated      <- FALSE
  recruit$size          <- 1L    #one cell (size is a cell count)
  recruit$age           <- 0L     #a newborn: reset age so it must mature before spawning
                                  #(recruit <- parent copied the parent's age above)
  recruit$repro_cooldown <- 0     #a newborn is only gated by maturity, not a cooldown
  recruit$colour        <- recruitShade(parent$colour)  #own shade, not identical to parent
  recruit$coords        <- c(col, row)  #stored as c(x, y)

  corals[[length(corals) + 1L]] <- recruit
  colony_species[recruit_id]    <- parent$species

  list(reef = reef, corals = corals, colony_species = colony_species)
}


# --------------------------------------------------------------------
#  STEP 6: One full reproduction pass, run once per timestep.
#  For every existing colony: count its cooldown down by one; if it is still
#  cooling down from a previous spawn it cannot reproduce this step. Otherwise, if
#  it is mature enough, it rolls its size-scaled chance and on success drops a
#  recruit on a random empty cell (larvae settle anywhere on the reef) and starts
#  a fresh REPRO_COOLDOWN wait. Reef space is limited, so if empty cells run out no
#  further recruits settle that step. recruit_counter hands out unique ids.
#  Returns the updated reef, corals, colony_species and recruit_counter.
# --------------------------------------------------------------------
reproduceColonies <- function(reef, corals, colony_species, species_traits, recruit_counter,
                              timestep_num = NA) {

  #Synchronised-spawning gate: with REPRO_INTERVAL > 1, spawning only happens on
  #event timesteps (t = interval, 2*interval, ...). Interval 1 = every step (original).
  if (REPRO_INTERVAL > 1 &&
      (is.na(timestep_num) || timestep_num < REPRO_INTERVAL ||
       timestep_num %% REPRO_INTERVAL != 0)) {
    return(list(reef = reef, corals = corals,
                colony_species = colony_species, recruit_counter = recruit_counter))
  }

  #Nothing to do if the feature is off or no species can reproduce
  if (!("reproduction" %in% names(species_traits)) || !any(species_traits$reproduction)) {
    return(list(reef = reef, corals = corals,
                colony_species = colony_species, recruit_counter = recruit_counter))
  }

  #Empty cells available for settlement, found lazily the first time one is needed
  #(larvae land anywhere on empty reef)
  empties <- NULL

  #Iterate the ORIGINAL colonies in RANDOM order, so when empty cells are scarce no
  #colony (or species) has a fixed priority for settling recruits. We shuffle the
  #INDICES (not the list), so cooldowns are still updated in place; the index set is
  #fixed at the start, so recruits appended below are not iterated this step (a newborn
  #cannot reproduce the same step it is created).
  for (i in sample(seq_along(corals))) {
    co <- corals[[i]]

    #Count down this colony's cooldown; if still waiting, it cannot reproduce yet
    cd <- if (is.null(co$repro_cooldown)) 0 else co$repro_cooldown
    if (cd > 0) {
      corals[[i]]$repro_cooldown <- cd - 1L
      next
    }
    corals[[i]]$repro_cooldown <- 0

    #Skip dead colonies (size is % cover, so a live colony has size > 0)
    if (co$size <= 0) {
      next
    }

    #Age-based maturity gate: a colony must have been alive at least its species'
    #maturity age before it can reproduce, regardless of its size.
    if (!isMature(co, species_traits)) {
      next
    }

    #Mature: roll the size-scaled chance
    if (runif(1) >= colonyReproductionChance(co, species_traits, reef)) {
      next
    }

    #Find empty cells the first time we actually need one
    if (is.null(empties)) {
      empties <- which(is.na(as.vector(reef)))
    }
    if (length(empties) == 0) {
      next   #reef is full, no room to settle a recruit (still count down others)
    }

    #Pick and consume a random empty cell
    pick <- sample(length(empties), 1)
    cell <- empties[pick]
    empties <- empties[-pick]

    recruit_counter <- recruit_counter + 1L
    recruit_id <- paste0(co$species, "_r", recruit_counter)

    placed <- spawnRecruit(reef, corals, colony_species, co, cell, recruit_id)
    reef           <- placed$reef
    corals         <- placed$corals
    colony_species <- placed$colony_species

    #This colony must now wait before it can reproduce again
    corals[[i]]$repro_cooldown <- REPRO_COOLDOWN
  }

  list(reef = reef, corals = corals,
       colony_species = colony_species, recruit_counter = recruit_counter)
}


# --------------------------------------------------------------------
#  HABITAT SETUP: ask which habitat this run is, which controls background
#  recruitment (see above). Asked during setup, like reproduction/disturbance.
#    Aquarium     - closed system; NO background recruitment (only living colonies
#                   reproduce, via the reproduction toggle).
#    Natural reef - background recruitment ON; also asks the per-species chance.
#    None         - no habitat effect; no background recruitment.
#  Fully-random setup (3) picks a habitat at random (and a random chance if reef).
#  Returns a habitat config list used by backgroundReproduce().
# --------------------------------------------------------------------
getHabitatSetup <- function(setup_mode) {

  #Fully random: choose a habitat at random; a reef gets a random background chance
  #and a random supply pattern (pulse / continuous)
  if (setup_mode == 3) {
    type <- sample(c("aquarium", "reef", "none"), 1)
    if (type == "reef") {
      return(reefHabitat(runif(1, BACKGROUND_REPRO_RANDOM_MIN, BACKGROUND_REPRO_RANDOM_MAX),
                         sample(c("pulse", "continuous"), 1)))
    }
    return(noBackground(type))
  }

  #Manual / semi-manual: ask which habitat
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
    #Temporal pattern of that supply: a periodic pulse or a continuous stochastic trickle
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


# --------------------------------------------------------------------
#  Habitat config constructors, shared by the interactive setup (getHabitatSetup)
#  and the non-interactive batch builder (buildHabitatConfig). A reef is an OPEN
#  patch (external larval supply); an aquarium/"none" is a CLOSED patch (no supply).
# --------------------------------------------------------------------
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

# --------------------------------------------------------------------
#  Non-interactive habitat builder for batch sweeps (parallels getHabitatSetup).
#    habitat     - "reef" (open: external larval supply) or "aquarium"/"none"
#                  (closed: no supply).
#    chance      - background recruitment chance per species (reef only).
#    supply_mode - "pulse" or "continuous" (reef only; ignored otherwise).
#    interval    - timesteps between supply pulses (reef only).
#    local_retention - LOCAL SELF-SEEDING: a species' arrival chance is raised by
#                  local_retention x (its MATURE cover as a fraction of the reef),
#                  on top of the fixed base chance. Locally produced larvae partly
#                  settle back on the natal reef, so supply is very slightly
#                  boosted by resident, reproducing species. 0 = off (default).
# --------------------------------------------------------------------
buildHabitatConfig <- function(habitat = "aquarium", chance = 0.1,
                               supply_mode = DEFAULT_SUPPLY_MODE,
                               interval = BACKGROUND_REPRO_INTERVAL,
                               local_retention = 0) {
  if (identical(habitat, "reef")) reefHabitat(chance, supply_mode, interval, local_retention)
  else                            noBackground(habitat)
}


# --------------------------------------------------------------------
#  The solid base colour of each species, matching how createCorals assigns them
#  (hcl.colors "Dark 3" in species order). Lets a background recruit take its
#  species' colour even when that species has no living colony left to copy.
# --------------------------------------------------------------------
speciesColourMap <- function(species_traits) {
  setNames(hcl.colors(nrow(species_traits), "Dark 3"), species_traits$species)
}


# --------------------------------------------------------------------
#  One-line description of a habitat config, for the setup summary.
# --------------------------------------------------------------------
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


# --------------------------------------------------------------------
#  Place one BACKGROUND recruit of a given species on an empty cell. Unlike a
#  colony recruit it has no parent (it arrives from outside), so it is built from
#  the species' trait row and base colour. Fields match a normal coral plus the
#  recruit flags; $background = TRUE marks its larval-influx origin.
# --------------------------------------------------------------------
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
    size     = 1L,       #one cell (size is a cell count)
    age      = 0L,       #newborn: must mature (MATURITY_AGE) before it can spawn
    recruit    = TRUE,   #a recruit (not a split fragment)
    background = TRUE,   #arrived from outside the reef, not from a parent colony
    isolated   = FALSE,
    parent_id  = NA_character_,
    origin_id  = recruit_id,
    repro_cooldown = 0
  )

  corals[[length(corals) + 1L]] <- recruit
  colony_species[recruit_id]    <- trait_row$species

  list(reef = reef, corals = corals, colony_species = colony_species)
}


# --------------------------------------------------------------------
#  One background-recruitment pass, run once per timestep. Does nothing unless the
#  run is a natural reef with background recruitment enabled. The temporal pattern of
#  supply follows habitat$supply_mode (see the section header):
#    "pulse"      - fires only on interval events (t = interval, 2*interval, ...),
#                   each species with probability background_chance.
#    "continuous" - fires EVERY step, each species with probability
#                   background_chance / interval (matched expected supply).
#  In either case every species in the pool gets the SAME independent chance, tried
#  in random order, drawing from the whole species pool so a locally extinct species
#  can reappear. Reef space is limited: if empty cells run out no more settle this
#  step. Returns the updated reef, corals, colony_species and recruit_counter.
# --------------------------------------------------------------------
backgroundReproduce <- function(reef, corals, colony_species, species_traits,
                                habitat, recruit_counter, timestep_num) {

  unchanged <- list(reef = reef, corals = corals,
                    colony_species = colony_species, recruit_counter = recruit_counter)

  #Off unless this is a natural-reef run with background recruitment enabled
  if (is.null(habitat) || !isTRUE(habitat$background_enabled)) {
    return(unchanged)
  }

  chance   <- habitat$background_chance
  interval <- habitat$background_interval
  if (is.na(chance) || chance <= 0) {
    return(unchanged)
  }

  #Resolve the per-species settlement probability for THIS step from the supply mode.
  #Continuous trickles the same expected supply across every step; pulse concentrates
  #it into periodic events. Default to "pulse" when the field is absent (older configs).
  mode <- if (is.null(habitat$supply_mode)) "pulse" else habitat$supply_mode
  if (identical(mode, "continuous")) {
    #Every step; rate matched to the pulse so expected supply per species is unchanged
    step_chance <- chance / interval
  } else {
    #Pulse: only on the interval, after the first interval (t = 30, 60, 90, ...)
    if (is.na(timestep_num) || timestep_num < interval || timestep_num %% interval != 0) {
      return(unchanged)
    }
    step_chance <- chance
  }

  colours <- speciesColourMap(species_traits)
  empties <- which(is.na(as.vector(reef)))

  #LOCAL RETENTION (optional): each species' arrival chance is raised slightly in
  #proportion to its own MATURE cover on the reef (self-seeding - locally produced
  #larvae partly settle back on the natal patch). Locally absent species keep the
  #full base chance, so rescue from the regional pool is never weakened.
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
    boost <- lr * (mature_cells / total_cells)          #<= lr even at full cover
    if (identical(mode, "continuous")) boost <- boost / interval
    sp_chance <- pmin(1, sp_chance + as.numeric(boost[species_traits$species]))
  }

  #Every species rolls independently, tried in RANDOM order so none has a fixed
  #priority for scarce empty cells.
  for (i in sample(seq_len(n_sp))) {
    if (length(empties) == 0) break        #reef full, no room to settle
    if (runif(1) >= sp_chance[i]) next     #this species: no recruit this step

    #Pick and consume a random empty cell
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
