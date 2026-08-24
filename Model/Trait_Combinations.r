# Builders for the competition-network interaction matrices, plus the named
# "trait combinations" (models) swept by the harness. Every builder returns an
# N x N overgrowth matrix M (M[i,j] = P(i overgrows j); diagonal = intraspecific).
# bias = the stronger species' win probability (0.5-1); intra = diagonal probability.

RPS_MIN_SPECIES <- 3
RPS_MAX_SPECIES <- 15


# Generalized rock-paper-scissors (intransitive cycle). Needs ODD n: each species
# beats the next (n-1)/2 and loses to the previous (n-1)/2, so no species is best.
rpsMatrix <- function(n, bias = 0.9, intra = 0.5) {
  if (n < RPS_MIN_SPECIES || n > RPS_MAX_SPECIES) {
    stop("RPS species count must be between ", RPS_MIN_SPECIES, " and ",
         RPS_MAX_SPECIES, "; got ", n)
  }
  if (n %% 2 == 0) {
    stop("RPS needs an ODD number of species (so the cycle is balanced); got ", n)
  }

  M <- matrix(NA_real_, n, n)
  k <- (n - 1) / 2
  for (i in seq_len(n)) {
    for (off in seq_len(k)) {
      win  <- ((i - 1 + off) %% n) + 1      # i beats this one
      lose <- ((i - 1 - off) %% n) + 1      # i loses to this one
      M[i, win]  <- bias
      M[i, lose] <- 1 - bias
    }
  }
  diag(M) <- intra
  species <- paste0("Sp", seq_len(n))
  dimnames(M) <- list(species, species)
  M
}


# Linear (transitive) hierarchy: Sp1 > Sp2 > ... > Spn. Any n >= 2.
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


# Neutral null: no overgrowth (off-diagonal 0), so colonies compete only for space.
neutralMatrix <- function(n, intra = 0.5) {
  M <- matrix(0, n, n)
  diag(M) <- intra
  species <- paste0("Sp", seq_len(n))
  dimnames(M) <- list(species, species)
  M
}


# Random pairwise network: each off-diagonal an independent draw in [lo, hi]. Drawn
# under a fixed per-n seed, so the same matrix recurs every run (RNG stream restored
# afterwards, so it never perturbs the simulation). bias is unused; intra = diagonal.
RANDOM_PROB_MIN  <- 0.10
RANDOM_PROB_MAX  <- 0.90
RANDOM_SEED_BASE <- 7000L   # per-n seed is RANDOM_SEED_BASE + n

randomMatrix <- function(n, intra = 0.5, lo = RANDOM_PROB_MIN, hi = RANDOM_PROB_MAX,
                         seed = NULL) {
  if (n < 2) stop("Random network needs at least 2 species; got ", n)
  # draw under a private fixed seed, saving/restoring the caller's RNG stream
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


# Build a matrix by type ("rps" / "hierarchy"(="linear") / "neutral" / "random").
makeTraitMatrix <- function(type, n, bias = 0.9, intra = 0.5) {
  switch(type,
         rps       = rpsMatrix(n, bias, intra),
         hierarchy = hierarchyMatrix(n, bias, intra),
         linear    = hierarchyMatrix(n, bias, intra),
         neutral   = neutralMatrix(n, intra),
         random    = randomMatrix(n, intra, seed = RANDOM_SEED_BASE + n),
         stop("Unknown trait-matrix type: ", type))
}

TRAIT_MATRIX_TYPES <- c("rps", "hierarchy", "neutral", "random")


# ---- Size-competition beta (per-species log-ratio sensitivity) --------------
# beta > 0 = bigger colony wins more, < 0 = smaller wins more, 0 = no effect.
SIZECOMP_BETA_MIN <- 0.3
SIZECOMP_BETA_MAX <- 0.6

# RPS split: one neutral species (beta 0), the rest split evenly into bigger-wins
# (+) and smaller-wins (-) types, graded and mirrored. Odd n.
rpsSizeCompBeta <- function(n, beta_min = SIZECOMP_BETA_MIN, beta_max = SIZECOMP_BETA_MAX) {
  if (n %% 2 == 0) stop("RPS_sizecompImpact needs an ODD number of species; got ", n)
  species <- paste0("Sp", seq_len(n))
  beta <- setNames(rep(0, n), species)
  k <- (n - 1) / 2
  mags <- if (k == 1) mean(c(beta_min, beta_max)) else seq(beta_min, beta_max, length.out = k)
  for (j in seq_len(k)) {
    beta[1 + j]     <-  mags[j]   # bigger-wins
    beta[1 + k + j] <- -mags[j]   # smaller-wins
  }
  beta
}


# Rank-graded size betas: effect scales with dominance (Sp1 most, Spn none).
# reverse = FALSE: dominants gain from being bigger; TRUE: from being smaller.
linearSizeBeta <- function(n, reverse = FALSE, beta_max = SIZECOMP_BETA_MAX) {
  if (n < 2) stop("Linear size beta needs at least 2 species; got ", n)
  species <- paste0("Sp", seq_len(n))
  w   <- (n - seq_len(n)) / (n - 1)      # 1 at top, 0 at bottom
  sgn <- if (reverse) -1 else 1
  setNames(sgn * beta_max * w, species)
}


# Top of the growth range (5 if Sim_func.r not yet sourced).
.growthMax <- function() if (exists("GROWTH_MAX")) GROWTH_MAX else 5


# ---- Reproduction-rate profiles (per-species base spawn chance) -------------
# A shared reference size means ONLY the frequency differs between species.
REPRO_BASE_CHANCE <- 0.08
REPRO_SPREAD      <- 0.06
REPRO_SPREAD_MIN  <- 0.02
REPRO_SPREAD_MAX  <- 0.06
REPRO_REF_DEFAULT <- 5      # reference size (% cover) at which the base chance applies

# Rank-graded reproduction: reverse = FALSE -> more dominant breeds faster; TRUE ->
# less dominant breeds faster. Symmetric around baseline, clamped at 0.
linearReproChance <- function(n, reverse = FALSE, base = REPRO_BASE_CHANCE,
                              spread = REPRO_SPREAD) {
  if (n < 2) stop("Linear reproduction needs at least 2 species; got ", n)
  w   <- 1 - 2 * (seq_len(n) - 1) / (n - 1)     # +1 at top to -1 at bottom
  sgn <- if (reverse) -1 else 1
  pmax(0, base + sgn * spread * w)
}

# RPS reproduction: one neutral species (baseline), rest graded higher/lower. Odd n.
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


# ---- Named trait combinations (the swept "models") --------------------------
# Each returns a spec: list(name, matrix, size_mode, size_beta, [growth,
# growth_random, repro_chance, repro_ref]). size_beta NULL = no size effect;
# growth NULL = uniform (config), vector = fixed per-species, growth_random = drawn
# per run; repro_chance NULL = reproduction off, vector = on.

# RPS, no size effect.
Classic_RPS <- function(n, bias = 0.9, intra = 0.5) {
  list(name = "Classic_RPS", matrix = rpsMatrix(n, bias, intra),
       size_mode = "none", size_beta = NULL)
}

# Linear hierarchy, no size effect.
Classic_Linear <- function(n, bias = 0.9, intra = 0.5) {
  list(name = "Classic_Linear", matrix = hierarchyMatrix(n, bias, intra),
       size_mode = "none", size_beta = NULL)
}

# RPS + graded size-competition split.
RPS_sizecompImpact <- function(n, bias = 0.9, intra = 0.5,
                               beta_min = SIZECOMP_BETA_MIN, beta_max = SIZECOMP_BETA_MAX) {
  list(name = "RPS_sizecompImpact", matrix = rpsMatrix(n, bias, intra),
       size_mode = "logratio",
       size_beta = rpsSizeCompBeta(n, beta_min, beta_max))
}

# RPS + random per-species growth.
RPS_growthImpact <- function(n, bias = 0.9, intra = 0.5) {
  list(name = "RPS_growthImpact", matrix = rpsMatrix(n, bias, intra),
       size_mode = "none", size_beta = NULL,
       growth = NULL, growth_random = TRUE)
}

# Linear + competition/growth trade-off: less dominant grows faster (Sp1 slow, Spn fast).
Linear_GrowthImpact <- function(n, bias = 0.9, intra = 0.5) {
  growth <- seq(1, .growthMax(), length.out = n)
  list(name = "Linear_GrowthImpact", matrix = hierarchyMatrix(n, bias, intra),
       size_mode = "none", size_beta = NULL,
       growth = growth, growth_random = FALSE)
}

# Linear + size mode where dominants must be big to express dominance (graded by rank).
Linear_SizeImpactNormal <- function(n, bias = 0.9, intra = 0.5, beta = SIZECOMP_BETA_MAX) {
  list(name = "Linear_SizeImpactNormal", matrix = hierarchyMatrix(n, bias, intra),
       size_mode = "logratio",
       size_beta = linearSizeBeta(n, reverse = FALSE, beta_max = beta))
}

# Linear + reversed size rule: dominants gain from being small.
Linear_SizeImpactReverse <- function(n, bias = 0.9, intra = 0.5, beta = SIZECOMP_BETA_MAX) {
  list(name = "Linear_SizeImpactReverse", matrix = hierarchyMatrix(n, bias, intra),
       size_mode = "logratio",
       size_beta = linearSizeBeta(n, reverse = TRUE, beta_max = beta))
}

# Linear where more dominant breeds faster (Sp1 fastest).
Linear_reproductionNormal <- function(n, bias = 0.9, intra = 0.5) {
  list(name = "Linear_reproductionNormal", matrix = hierarchyMatrix(n, bias, intra),
       size_mode = "none", size_beta = NULL,
       repro_chance = linearReproChance(n, reverse = FALSE),
       repro_ref = REPRO_REF_DEFAULT)
}

# Linear competition/fecundity trade-off: less dominant breeds faster (Spn fastest).
Linear_reproductionOpposite <- function(n, bias = 0.9, intra = 0.5) {
  list(name = "Linear_reproductionOpposite", matrix = hierarchyMatrix(n, bias, intra),
       size_mode = "none", size_beta = NULL,
       repro_chance = linearReproChance(n, reverse = TRUE),
       repro_ref = REPRO_REF_DEFAULT)
}

# RPS with graded (unequal) reproduction.
RPS_reproduction <- function(n, bias = 0.9, intra = 0.5) {
  list(name = "RPS_reproduction", matrix = rpsMatrix(n, bias, intra),
       size_mode = "none", size_beta = NULL,
       repro_chance = rpsReproChance(n),
       repro_ref = REPRO_REF_DEFAULT)
}

# RPS with equal reproduction - the control isolating "reproduction present" from "unequal".
RPS_reproductionEven <- function(n, bias = 0.9, intra = 0.5) {
  list(name = "RPS_reproductionEven", matrix = rpsMatrix(n, bias, intra),
       size_mode = "none", size_beta = NULL,
       repro_chance = rep(REPRO_BASE_CHANCE, n),
       repro_ref = REPRO_REF_DEFAULT)
}

# Linear with equal reproduction (the even-reproduction control).
Linear_reproductionEven <- function(n, bias = 0.9, intra = 0.5) {
  list(name = "Linear_reproductionEven", matrix = hierarchyMatrix(n, bias, intra),
       size_mode = "none", size_beta = NULL,
       repro_chance = rep(REPRO_BASE_CHANCE, n),
       repro_ref = REPRO_REF_DEFAULT)
}

# Neutral null: no overgrowth, no size effect, no reproduction. Any n.
Classic_Neutral <- function(n, bias = 0.9, intra = 0.5) {
  list(name = "Classic_Neutral", matrix = neutralMatrix(n, intra),
       size_mode = "none", size_beta = NULL)
}

# Set random network (arbitrary asymmetric tournament), reproducible per n. Any n.
Classic_Random <- function(n, bias = 0.9, intra = 0.5) {
  list(name = "Classic_Random",
       matrix = randomMatrix(n, intra, seed = RANDOM_SEED_BASE + n),
       size_mode = "none", size_beta = NULL)
}

# ---- Reproduction / size / growth variants for neutral and random networks --
# (complete the ladders for every network type; graded profiles use linearReproChance,
#  which works for any n, so no odd-n limit.)

# Neutral + equal reproduction.
Neutral_reproductionEven <- function(n, bias = 0.9, intra = 0.5) {
  list(name = "Neutral_reproductionEven", matrix = neutralMatrix(n, intra),
       size_mode = "none", size_beta = NULL,
       repro_chance = rep(REPRO_BASE_CHANCE, n),
       repro_ref = REPRO_REF_DEFAULT)
}

# Neutral + graded reproduction: fecundity differences with zero competition.
Neutral_reproduction <- function(n, bias = 0.9, intra = 0.5) {
  list(name = "Neutral_reproduction", matrix = neutralMatrix(n, intra),
       size_mode = "none", size_beta = NULL,
       repro_chance = linearReproChance(n, reverse = FALSE),
       repro_ref = REPRO_REF_DEFAULT)
}

# Random + equal reproduction.
Random_reproductionEven <- function(n, bias = 0.9, intra = 0.5) {
  list(name = "Random_reproductionEven",
       matrix = randomMatrix(n, intra, seed = RANDOM_SEED_BASE + n),
       size_mode = "none", size_beta = NULL,
       repro_chance = rep(REPRO_BASE_CHANCE, n),
       repro_ref = REPRO_REF_DEFAULT)
}

# Random + graded reproduction.
Random_reproduction <- function(n, bias = 0.9, intra = 0.5) {
  list(name = "Random_reproduction",
       matrix = randomMatrix(n, intra, seed = RANDOM_SEED_BASE + n),
       size_mode = "none", size_beta = NULL,
       repro_chance = linearReproChance(n, reverse = FALSE),
       repro_ref = REPRO_REF_DEFAULT)
}

# Neutral + size mode. NOTE: with no interspecific overgrowth this only affects
# same-species contests, so between-species dynamics stay near neutral.
Neutral_sizecompImpact <- function(n, bias = 0.9, intra = 0.5, beta = SIZECOMP_BETA_MAX) {
  list(name = "Neutral_sizecompImpact", matrix = neutralMatrix(n, intra),
       size_mode = "logratio",
       size_beta = linearSizeBeta(n, reverse = FALSE, beta_max = beta))
}

# Neutral + random growth: growth differences alone decide who claims empty space.
Neutral_growthImpact <- function(n, bias = 0.9, intra = 0.5) {
  list(name = "Neutral_growthImpact", matrix = neutralMatrix(n, intra),
       size_mode = "none", size_beta = NULL,
       growth = NULL, growth_random = TRUE)
}

# Random + size mode.
Random_sizecompImpact <- function(n, bias = 0.9, intra = 0.5, beta = SIZECOMP_BETA_MAX) {
  list(name = "Random_sizecompImpact",
       matrix = randomMatrix(n, intra, seed = RANDOM_SEED_BASE + n),
       size_mode = "logratio",
       size_beta = linearSizeBeta(n, reverse = FALSE, beta_max = beta))
}

# Random + random growth.
Random_growthImpact <- function(n, bias = 0.9, intra = 0.5) {
  list(name = "Random_growthImpact",
       matrix = randomMatrix(n, intra, seed = RANDOM_SEED_BASE + n),
       size_mode = "none", size_beta = NULL,
       growth = NULL, growth_random = TRUE)
}


# Build a named trait-combination spec.
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

# Control models run at the original config (even-reproduction controls + neutral null).
CONTROL_COMBINATIONS <- c("RPS_reproductionEven", "Linear_reproductionEven",
                          "Classic_Neutral")
