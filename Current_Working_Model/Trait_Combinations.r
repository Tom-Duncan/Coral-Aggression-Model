# =============================================================================
#  Trait_Combinations.r  -  library of competition-network (interaction-matrix)
#                           builders for the model sweeps.
# -----------------------------------------------------------------------------
#  Every builder returns an N x N overgrowth-probability matrix M where
#     M[i, j] = P(species i overgrows species j)   (rows = attacker, cols = occupant)
#  with the DIAGONAL set to the intraspecific (same-species) probability. These plug
#  straight into the model via attachInteraction(). This file holds only builders -
#  no runs - and is sourced by "Simulation testing.r".
#
#  `bias` is the stronger species' win probability in a decided pair (0.5-1); set it
#  to 1 for a deterministic "always wins". `intra` is the same-species probability.
# =============================================================================

RPS_MIN_SPECIES <- 3    # generalized RPS: at least rock-paper-scissors
RPS_MAX_SPECIES <- 15   # upper bound for the sweeps


# --------------------------------------------------------------------
#  Generalized ROCK-PAPER-SCISSORS (intransitive cycle).
#  Requires an ODD number of species (3, 5, 7, ... up to RPS_MAX_SPECIES): each
#  species beats the next k = (n-1)/2 species around the cycle and loses to the
#  previous k, so every pair is decided and no species is best - the balanced
#  intransitive tournament that permits coexistence. An even n is rejected because
#  it would leave a diametrically-opposite pair with no clean winner.
#     e.g. n = 3  -> Sp1>Sp2, Sp2>Sp3, Sp3>Sp1              (classic RPS)
#          n = 5  -> each species beats the next 2, loses to the previous 2 (RPSLS)
# --------------------------------------------------------------------
rpsMatrix <- function(n, bias = 0.9, intra = 0.5) {
  if (n < RPS_MIN_SPECIES || n > RPS_MAX_SPECIES) {
    stop("RPS species count must be between ", RPS_MIN_SPECIES, " and ",
         RPS_MAX_SPECIES, "; got ", n)
  }
  if (n %% 2 == 0) {
    stop("RPS needs an ODD number of species (so the cycle is balanced); got ", n)
  }

  M <- matrix(NA_real_, n, n)
  k <- (n - 1) / 2                          # how many species each one beats
  for (i in seq_len(n)) {
    for (off in seq_len(k)) {
      win  <- ((i - 1 + off) %% n) + 1      # the off-th species clockwise: i beats it
      lose <- ((i - 1 - off) %% n) + 1      # the off-th species anticlockwise: i loses
      M[i, win]  <- bias
      M[i, lose] <- 1 - bias
    }
  }
  diag(M) <- intra
  species <- paste0("Sp", seq_len(n))
  dimnames(M) <- list(species, species)
  M
}


# --------------------------------------------------------------------
#  Linear HIERARCHY (transitive network).
#  A strict ranking Sp1 > Sp2 > ... > Spn: the top species overgrows every other, the
#  second overgrows all but the first, and so on. Predicts competitive exclusion down
#  to the single strongest species. Works for any n >= 2 (no odd/even restriction).
#     M[i, j] = bias        if i is ranked above j (i < j)
#             = 1 - bias    if i is ranked below j (i > j)
# --------------------------------------------------------------------
hierarchyMatrix <- function(n, bias = 0.9, intra = 0.5) {
  if (n < 2) stop("Hierarchy needs at least 2 species; got ", n)
  M <- matrix(NA_real_, n, n)
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      if (i < j)      M[i, j] <- bias
      else if (i > j) M[i, j] <- 1 - bias
    }
  }
  diag(M) <- intra
  species <- paste0("Sp", seq_len(n))
  dimnames(M) <- list(species, species)
  M
}


# --------------------------------------------------------------------
#  NEUTRAL network: no species ever overgrows another (all off-diagonal 0), so
#  colonies compete only for empty space. A useful null model for the sweeps.
# --------------------------------------------------------------------
neutralMatrix <- function(n, intra = 0.5) {
  M <- matrix(0, n, n)
  diag(M) <- intra
  species <- paste0("Sp", seq_len(n))
  dimnames(M) <- list(species, species)
  M
}


# --------------------------------------------------------------------
#  RANDOM pairwise network ("set random hierarchy").
#  Every off-diagonal entry is an independent random overgrowth probability in
#  [lo, hi], so the network is an arbitrary asymmetric tournament - neither a clean
#  intransitive cycle (RPS) nor a strict transitive rank (linear). Works for any n.
#
#  It is a SET network: drawn under a FIXED seed derived from n, so the SAME matrix is
#  produced on every run and is shared by all replicates of that species count (like a
#  frozen entry in Pairwise_Matrices.r). That keeps results reproducible and lets the
#  exact network be recovered from n. The private seed is applied and then the caller's
#  RNG stream is restored, so it never perturbs the rest of the simulation's draws.
#  `bias` is unused (entries are random, not a single win probability); `intra` sets the
#  diagonal (same-species) probability.
# --------------------------------------------------------------------
RANDOM_PROB_MIN  <- 0.10   # smallest random off-diagonal overgrowth probability
RANDOM_PROB_MAX  <- 0.90   # largest random off-diagonal overgrowth probability
RANDOM_SEED_BASE <- 7000L  # fixed base; the per-n seed is RANDOM_SEED_BASE + n

randomMatrix <- function(n, intra = 0.5, lo = RANDOM_PROB_MIN, hi = RANDOM_PROB_MAX,
                         seed = NULL) {
  if (n < 2) stop("Random network needs at least 2 species; got ", n)
  #Draw under a private fixed seed (identical every run for a given n) WITHOUT
  #disturbing the caller's RNG stream - saved and restored on exit.
  if (!is.null(seed)) {
    had <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    old <- if (had) get(".Random.seed", envir = .GlobalEnv) else NULL
    set.seed(seed)
    on.exit({
      if (had) assign(".Random.seed", old, envir = .GlobalEnv)
      else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
        rm(".Random.seed", envir = .GlobalEnv)
    }, add = TRUE)
  }
  M <- matrix(NA_real_, n, n)
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      if (i != j) M[i, j] <- round(runif(1, lo, hi), 2)
    }
  }
  diag(M) <- intra
  species <- paste0("Sp", seq_len(n))
  dimnames(M) <- list(species, species)
  M
}


# --------------------------------------------------------------------
#  Dispatcher: build a matrix by name. Central entry point used by the harness.
#    "rps"       - generalized rock-paper-scissors (odd n)
#    "hierarchy" - linear hierarchy / transitive ("linear" is an alias)
#    "neutral"   - no overgrowth
# --------------------------------------------------------------------
makeTraitMatrix <- function(type, n, bias = 0.9, intra = 0.5) {
  switch(type,
         rps       = rpsMatrix(n, bias, intra),
         hierarchy = hierarchyMatrix(n, bias, intra),
         linear    = hierarchyMatrix(n, bias, intra),
         neutral   = neutralMatrix(n, intra),
         random    = randomMatrix(n, intra, seed = RANDOM_SEED_BASE + n),
         stop("Unknown trait-matrix type: ", type))
}

#Names of the networks this library can build (handy for configuring a sweep).
TRAIT_MATRIX_TYPES <- c("rps", "hierarchy", "neutral", "random")


# =============================================================================
#  SIZE-COMPETITION BETA  (per-species sensitivity for the log-ratio size mode)
# =============================================================================
#  The size effect is a per-species log-ratio sensitivity beta (attachInteraction
#  size_beta): beta > 0 = the bigger colony wins more, beta < 0 = the smaller wins
#  more, 0 = size has no effect. These builders return a NAMED vector (species -> beta).
SIZECOMP_BETA_MIN <- 0.3   # smallest graded sensitivity
SIZECOMP_BETA_MAX <- 0.6   # largest graded sensitivity

#RPS split: exactly ONE species is neutral (beta 0); the rest split EQUALLY into a
#bigger-wins type (beta > 0) and a smaller-wins type (beta < 0), graded and mirrored
#so the set is balanced. Requires an ODD n (matches the RPS network). With n species
#there are k = (n-1)/2 of each type, taking sensitivities evenly spaced in
#[beta_min, beta_max] (a single level uses the midpoint).
#    e.g. n = 3 -> Sp1 beta 0; Sp2 +0.45; Sp3 -0.45
#         n = 5 -> Sp1 0; Sp2,Sp3 +0.3,+0.6; Sp4,Sp5 -0.3,-0.6
rpsSizeCompBeta <- function(n, beta_min = SIZECOMP_BETA_MIN, beta_max = SIZECOMP_BETA_MAX) {
  if (n %% 2 == 0) stop("RPS_sizecompImpact needs an ODD number of species; got ", n)
  species <- paste0("Sp", seq_len(n))
  beta <- setNames(rep(0, n), species)
  k <- (n - 1) / 2
  mags <- if (k == 1) mean(c(beta_min, beta_max)) else seq(beta_min, beta_max, length.out = k)
  for (j in seq_len(k)) {
    beta[1 + j]     <-  mags[j]   # bigger-wins type
    beta[1 + k + j] <- -mags[j]   # smaller-wins type
  }
  beta
}


# --------------------------------------------------------------------
#  LINEAR (rank-graded) size betas for a hierarchy. The effect scales with DOMINANCE:
#  the most dominant species (Sp1) is affected most, the least dominant (Spn) not at
#  all, with a linear gradient between.
#    reverse = FALSE ("normal")  - dominants gain from being BIGGER (beta > 0): they
#                                  must be large to express dominance.
#    reverse = TRUE              - dominants gain from being SMALLER (beta < 0): more
#                                  dominant as they get smaller.
#  w = (n - rank)/(n - 1) is the per-species weight (1 at the top, 0 at the bottom).
# --------------------------------------------------------------------
linearSizeBeta <- function(n, reverse = FALSE, beta_max = SIZECOMP_BETA_MAX) {
  if (n < 2) stop("Linear size beta needs at least 2 species; got ", n)
  species <- paste0("Sp", seq_len(n))
  w   <- (n - seq_len(n)) / (n - 1)      # 1 for Sp1 (top) down to 0 for Spn
  sgn <- if (reverse) -1 else 1
  setNames(sgn * beta_max * w, species)
}


#Top of the growth trait range (falls back to 5 if Sim_func.r is not yet sourced).
.growthMax <- function() if (exists("GROWTH_MAX")) GROWTH_MAX else 5


# =============================================================================
#  REPRODUCTION-RATE PROFILES  (per-species base spawn chance = "frequency")
# -----------------------------------------------------------------------------
#  The reproduction models vary each species' base per-timestep spawn chance
#  (repro_chance) around a common baseline. A fixed reference size (repro_ref) is
#  shared by all species so ONLY the frequency differs, not the size-scaling.
# =============================================================================
REPRO_BASE_CHANCE <- 0.08   # baseline per-timestep spawn chance (the "neutral" rate)
REPRO_SPREAD      <- 0.06   # max deviation from baseline for the linear gradients
REPRO_SPREAD_MIN  <- 0.02   # smallest graded step for the RPS split
REPRO_SPREAD_MAX  <- 0.06   # largest graded step for the RPS split
REPRO_REF_DEFAULT <- 5      # reference size (% reef cover) at which the base chance applies

#Linear (rank-graded) reproduction chances for a hierarchy.
#  reverse = FALSE ("normal")  - MORE dominant -> higher chance (Sp1 fastest breeder)
#  reverse = TRUE ("opposite") - LESS dominant -> higher chance (Spn fastest breeder)
#The gradient is symmetric around the baseline (top = base + spread, bottom = base -
#spread), clamped at 0.
linearReproChance <- function(n, reverse = FALSE, base = REPRO_BASE_CHANCE,
                              spread = REPRO_SPREAD) {
  if (n < 2) stop("Linear reproduction needs at least 2 species; got ", n)
  w   <- 1 - 2 * (seq_len(n) - 1) / (n - 1)     # +1 for Sp1 (top) down to -1 for Spn
  sgn <- if (reverse) -1 else 1
  pmax(0, base + sgn * spread * w)
}

#RPS reproduction chances: exactly ONE neutral species (Sp1 at the baseline) and the
#rest split into higher and lower breeders on a GRADED scale so no two species share a
#rate. k = (n-1)/2 species take base + graded steps (higher) and k take base - graded
#steps (lower). Requires an odd n (matching the RPS network). Clamped at 0.
rpsReproChance <- function(n, base = REPRO_BASE_CHANCE,
                           spread_min = REPRO_SPREAD_MIN, spread_max = REPRO_SPREAD_MAX) {
  if (n %% 2 == 0) stop("RPS_reproduction needs an ODD number of species; got ", n)
  chance <- rep(base, n)
  k <- (n - 1) / 2
  mags <- if (k == 1) mean(c(spread_min, spread_max)) else seq(spread_min, spread_max, length.out = k)
  for (j in seq_len(k)) {
    chance[1 + j]     <- base + mags[j]     # higher breeder
    chance[1 + k + j] <- base - mags[j]     # lower breeder
  }
  pmax(0, chance)
}


# =============================================================================
#  NAMED TRAIT COMBINATIONS (the "models" swept by the harness)
# -----------------------------------------------------------------------------
#  Each returns a spec: list(name, matrix, size_mode, size_beta, growth,
#  growth_random). The harness turns a spec into species_traits (matrix via
#  attachInteraction, size_mode + the per-species size_beta vector for logratio).
#    size_beta    = NULL  -> no size effect
#    growth       = NULL  -> uniform growth from the config; a numeric vector sets a
#                            fixed per-species growth (e.g. graded by rank)
#    growth_random= TRUE  -> the harness draws a random growth per species each run
#    repro_chance = NULL  -> reproduction OFF; a numeric vector switches it ON with a
#                            per-species base spawn chance (repro_ref = reference size)
# =============================================================================

#Classic_RPS - generalized rock-paper-scissors, no size effect.
Classic_RPS <- function(n, bias = 0.9, intra = 0.5) {
  list(name = "Classic_RPS", matrix = rpsMatrix(n, bias, intra),
       size_mode = "none", size_beta = NULL)
}

#Classic_Linear - linear hierarchy (transitive), no size effect.
Classic_Linear <- function(n, bias = 0.9, intra = 0.5) {
  list(name = "Classic_Linear", matrix = hierarchyMatrix(n, bias, intra),
       size_mode = "none", size_beta = NULL)
}

#RPS_sizecompImpact - the RPS network PLUS the additive size-competition effect: a
#balanced, graded split of size responses (see rpsSizeCompBeta).
RPS_sizecompImpact <- function(n, bias = 0.9, intra = 0.5,
                               beta_min = SIZECOMP_BETA_MIN, beta_max = SIZECOMP_BETA_MAX) {
  list(name = "RPS_sizecompImpact", matrix = rpsMatrix(n, bias, intra),
       size_mode = "logratio",
       size_beta = rpsSizeCompBeta(n, beta_min, beta_max))
}

#RPS_growthImpact - RPS network, but each species is given a RANDOM growth rate (drawn
#per run by the harness). No size effect.
RPS_growthImpact <- function(n, bias = 0.9, intra = 0.5) {
  list(name = "RPS_growthImpact", matrix = rpsMatrix(n, bias, intra),
       size_mode = "none", size_beta = NULL,
       growth = NULL, growth_random = TRUE)
}

#Linear_GrowthImpact - linear hierarchy with a competition/growth TRADE-OFF: the less
#dominant a species, the faster it grows. Sp1 (most dominant) is slowest (growth 1),
#Spn (least dominant) is fastest (growth GROWTH_MAX), graded linearly between.
Linear_GrowthImpact <- function(n, bias = 0.9, intra = 0.5) {
  growth <- seq(1, .growthMax(), length.out = n)   # Sp1 slow ... Spn fast
  list(name = "Linear_GrowthImpact", matrix = hierarchyMatrix(n, bias, intra),
       size_mode = "none", size_beta = NULL,
       growth = growth, growth_random = FALSE)
}

#Linear_SizeImpactNormal - linear hierarchy + additive size mode where the more
#dominant a species, the more it is HURT when small (and helped when large): a
#dominant colony must be big to express its dominance. Graded by rank.
Linear_SizeImpactNormal <- function(n, bias = 0.9, intra = 0.5, beta = SIZECOMP_BETA_MAX) {
  list(name = "Linear_SizeImpactNormal", matrix = hierarchyMatrix(n, bias, intra),
       size_mode = "logratio",
       size_beta = linearSizeBeta(n, reverse = FALSE, beta_max = beta))
}

#Linear_SizeImpactReverse - linear hierarchy + additive size mode with the OPPOSITE
#rule: the more dominant a species, the more it is HELPED when small (and hurt when
#large) - more dominant as it gets smaller. Graded by rank.
Linear_SizeImpactReverse <- function(n, bias = 0.9, intra = 0.5, beta = SIZECOMP_BETA_MAX) {
  list(name = "Linear_SizeImpactReverse", matrix = hierarchyMatrix(n, bias, intra),
       size_mode = "logratio",
       size_beta = linearSizeBeta(n, reverse = TRUE, beta_max = beta))
}

#Linear_reproductionNormal - linear hierarchy where the MORE dominant a species, the
#higher its reproduction frequency (Sp1 breeds fastest, graded down by rank).
Linear_reproductionNormal <- function(n, bias = 0.9, intra = 0.5) {
  list(name = "Linear_reproductionNormal", matrix = hierarchyMatrix(n, bias, intra),
       size_mode = "none", size_beta = NULL,
       repro_chance = linearReproChance(n, reverse = FALSE),
       repro_ref = REPRO_REF_DEFAULT)
}

#Linear_reproductionOpposite - linear hierarchy where the LESS dominant a species, the
#higher its reproduction frequency (Spn breeds fastest): a competition/fecundity
#trade-off.
Linear_reproductionOpposite <- function(n, bias = 0.9, intra = 0.5) {
  list(name = "Linear_reproductionOpposite", matrix = hierarchyMatrix(n, bias, intra),
       size_mode = "none", size_beta = NULL,
       repro_chance = linearReproChance(n, reverse = TRUE),
       repro_ref = REPRO_REF_DEFAULT)
}

#RPS_reproduction - RPS network where one species is NEUTRAL (baseline reproduction)
#and the rest are graded higher/lower breeders, no two the same (see rpsReproChance).
RPS_reproduction <- function(n, bias = 0.9, intra = 0.5) {
  list(name = "RPS_reproduction", matrix = rpsMatrix(n, bias, intra),
       size_mode = "none", size_beta = NULL,
       repro_chance = rpsReproChance(n),
       repro_ref = REPRO_REF_DEFAULT)
}

#RPS_reproductionEven - RPS network with reproduction ON but EQUAL for every species
#(all at the baseline rate). The CONTROL that separates "reproduction is present"
#from "reproduction is unequal": Classic_RPS (off) -> this (even) -> RPS_reproduction
#(unequal).
RPS_reproductionEven <- function(n, bias = 0.9, intra = 0.5) {
  list(name = "RPS_reproductionEven", matrix = rpsMatrix(n, bias, intra),
       size_mode = "none", size_beta = NULL,
       repro_chance = rep(REPRO_BASE_CHANCE, n),
       repro_ref = REPRO_REF_DEFAULT)
}

#Linear_reproductionEven - linear hierarchy with reproduction ON but EQUAL for every
#species (the even-reproduction control for the linear family).
Linear_reproductionEven <- function(n, bias = 0.9, intra = 0.5) {
  list(name = "Linear_reproductionEven", matrix = hierarchyMatrix(n, bias, intra),
       size_mode = "none", size_beta = NULL,
       repro_chance = rep(REPRO_BASE_CHANCE, n),
       repro_ref = REPRO_REF_DEFAULT)
}

#Classic_Neutral - NEUTRAL null model: no species overgrows another (all off-diagonal
#0), so colonies compete only for empty space. No size effect, no reproduction. The
#"coexistence by chance alone" baseline. Works for any n (bias is unused).
Classic_Neutral <- function(n, bias = 0.9, intra = 0.5) {
  list(name = "Classic_Neutral", matrix = neutralMatrix(n, intra),
       size_mode = "none", size_beta = NULL)
}

#Classic_Random - a SET random pairwise network (arbitrary asymmetric tournament),
#fixed per species count and reproducible, no size effect. The 4th network TYPE
#alongside RPS / linear / neutral - a middle ground between the structured cycle and
#rank and the empty neutral null. Works for any n (bias is unused; the matrix is
#recovered from the fixed per-n seed).
Classic_Random <- function(n, bias = 0.9, intra = 0.5) {
  list(name = "Classic_Random",
       matrix = randomMatrix(n, intra, seed = RANDOM_SEED_BASE + n),
       size_mode = "none", size_beta = NULL)
}

# --------------------------------------------------------------------
#  REPRODUCTION variants for the NEUTRAL and RANDOM networks - completing the
#  reproduction ladder (OFF -> EVEN -> UNEQUAL) for every network type, matching the
#  RPS/Linear families. The graded ("unequal") profile uses linearReproChance, which
#  grades fecundity by species index (Sp1 highest -> Spn lowest) and works for ANY n,
#  so these keep the neutral/random "any species count" property (no odd-n limit). The
#  reference size is shared, so ONLY the frequency differs between species.
# --------------------------------------------------------------------

#Neutral_reproductionEven - neutral null WITH reproduction on but EQUAL for all species.
#Isolates the effect of simply ADDING recruitment to the space-only null.
Neutral_reproductionEven <- function(n, bias = 0.9, intra = 0.5) {
  list(name = "Neutral_reproductionEven", matrix = neutralMatrix(n, intra),
       size_mode = "none", size_beta = NULL,
       repro_chance = rep(REPRO_BASE_CHANCE, n),
       repro_ref = REPRO_REF_DEFAULT)
}

#Neutral_reproduction - neutral null WITH UNEQUAL (graded) reproduction: no overgrowth
#at all, but species differ in fecundity. Isolates whether fecundity differences ALONE
#- with zero competition - can structure the community.
Neutral_reproduction <- function(n, bias = 0.9, intra = 0.5) {
  list(name = "Neutral_reproduction", matrix = neutralMatrix(n, intra),
       size_mode = "none", size_beta = NULL,
       repro_chance = linearReproChance(n, reverse = FALSE),
       repro_ref = REPRO_REF_DEFAULT)
}

#Random_reproductionEven - set random network WITH reproduction on but EQUAL for all.
Random_reproductionEven <- function(n, bias = 0.9, intra = 0.5) {
  list(name = "Random_reproductionEven",
       matrix = randomMatrix(n, intra, seed = RANDOM_SEED_BASE + n),
       size_mode = "none", size_beta = NULL,
       repro_chance = rep(REPRO_BASE_CHANCE, n),
       repro_ref = REPRO_REF_DEFAULT)
}

#Random_reproduction - set random network WITH UNEQUAL (graded) reproduction: fecundity
#differences layered on an arbitrary competition network.
Random_reproduction <- function(n, bias = 0.9, intra = 0.5) {
  list(name = "Random_reproduction",
       matrix = randomMatrix(n, intra, seed = RANDOM_SEED_BASE + n),
       size_mode = "none", size_beta = NULL,
       repro_chance = linearReproChance(n, reverse = FALSE),
       repro_ref = REPRO_REF_DEFAULT)
}

# --------------------------------------------------------------------
#  SIZE and GROWTH variants for the NEUTRAL and RANDOM networks - completing the
#  size x network and growth x network grids (RPS and Linear already have theirs).
#  Size uses linearSizeBeta (graded by species index, any n); growth uses a random
#  per-species growth rate drawn per run, matching RPS_growthImpact.
# --------------------------------------------------------------------

#Neutral_sizecompImpact - neutral network + the size-competition (logratio) mode.
#NOTE: with no interspecific overgrowth, the size effect only modulates SAME-species
#(intraspecific / diagonal) contests, so between-species dynamics stay close to the
#neutral null. Included for a complete size x network grid.
Neutral_sizecompImpact <- function(n, bias = 0.9, intra = 0.5, beta = SIZECOMP_BETA_MAX) {
  list(name = "Neutral_sizecompImpact", matrix = neutralMatrix(n, intra),
       size_mode = "logratio",
       size_beta = linearSizeBeta(n, reverse = FALSE, beta_max = beta))
}

#Neutral_growthImpact - neutral network + a RANDOM per-species growth rate (drawn per
#run). With no overgrowth, growth-rate differences ALONE decide who claims empty space
#- a clean test of whether growth variation structures a community without competition.
Neutral_growthImpact <- function(n, bias = 0.9, intra = 0.5) {
  list(name = "Neutral_growthImpact", matrix = neutralMatrix(n, intra),
       size_mode = "none", size_beta = NULL,
       growth = NULL, growth_random = TRUE)
}

#Random_sizecompImpact - set random network + the size-competition (logratio) mode:
#size-modulated overgrowth on an arbitrary competition network.
Random_sizecompImpact <- function(n, bias = 0.9, intra = 0.5, beta = SIZECOMP_BETA_MAX) {
  list(name = "Random_sizecompImpact",
       matrix = randomMatrix(n, intra, seed = RANDOM_SEED_BASE + n),
       size_mode = "logratio",
       size_beta = linearSizeBeta(n, reverse = FALSE, beta_max = beta))
}

#Random_growthImpact - set random network + a RANDOM per-species growth rate (drawn per
#run): growth variation layered on an arbitrary competition network.
Random_growthImpact <- function(n, bias = 0.9, intra = 0.5) {
  list(name = "Random_growthImpact",
       matrix = randomMatrix(n, intra, seed = RANDOM_SEED_BASE + n),
       size_mode = "none", size_beta = NULL,
       growth = NULL, growth_random = TRUE)
}


# --------------------------------------------------------------------
#  Dispatcher: build a named trait-combination spec. Central entry point used by the
#  harness so a sweep just lists combination names.
# --------------------------------------------------------------------
makeTraitCombination <- function(name, n, bias = 0.9, intra = 0.5) {
  switch(name,
         Classic_RPS                 = Classic_RPS(n, bias, intra),
         Classic_Linear              = Classic_Linear(n, bias, intra),
         RPS_sizecompImpact          = RPS_sizecompImpact(n, bias, intra),
         RPS_growthImpact            = RPS_growthImpact(n, bias, intra),
         Linear_GrowthImpact         = Linear_GrowthImpact(n, bias, intra),
         Linear_SizeImpactNormal     = Linear_SizeImpactNormal(n, bias, intra),
         Linear_SizeImpactReverse    = Linear_SizeImpactReverse(n, bias, intra),
         Linear_reproductionNormal   = Linear_reproductionNormal(n, bias, intra),
         Linear_reproductionOpposite = Linear_reproductionOpposite(n, bias, intra),
         RPS_reproduction            = RPS_reproduction(n, bias, intra),
         RPS_reproductionEven        = RPS_reproductionEven(n, bias, intra),
         Linear_reproductionEven     = Linear_reproductionEven(n, bias, intra),
         Classic_Neutral             = Classic_Neutral(n, bias, intra),
         Classic_Random              = Classic_Random(n, bias, intra),
         Neutral_reproductionEven    = Neutral_reproductionEven(n, bias, intra),
         Neutral_reproduction        = Neutral_reproduction(n, bias, intra),
         Random_reproductionEven     = Random_reproductionEven(n, bias, intra),
         Random_reproduction         = Random_reproduction(n, bias, intra),
         Neutral_sizecompImpact      = Neutral_sizecompImpact(n, bias, intra),
         Neutral_growthImpact        = Neutral_growthImpact(n, bias, intra),
         Random_sizecompImpact       = Random_sizecompImpact(n, bias, intra),
         Random_growthImpact         = Random_growthImpact(n, bias, intra),
         stop("Unknown trait combination: ", name))
}

#The named combinations this library provides.
TRAIT_COMBINATIONS <- c("Classic_RPS", "Classic_Linear", "RPS_sizecompImpact",
                        "RPS_growthImpact", "Linear_GrowthImpact",
                        "Linear_SizeImpactNormal", "Linear_SizeImpactReverse",
                        "Linear_reproductionNormal", "Linear_reproductionOpposite",
                        "RPS_reproduction",
                        "RPS_reproductionEven", "Linear_reproductionEven",
                        "Classic_Neutral", "Classic_Random",
                        "Neutral_reproductionEven", "Neutral_reproduction",
                        "Random_reproductionEven", "Random_reproduction",
                        "Neutral_sizecompImpact", "Neutral_growthImpact",
                        "Random_sizecompImpact", "Random_growthImpact")

#The models added to fill dataset gaps (run at the ORIGINAL config, appended to the
#main dataframe): the even-reproduction controls and the neutral null model.
CONTROL_COMBINATIONS <- c("RPS_reproductionEven", "Linear_reproductionEven",
                          "Classic_Neutral")
