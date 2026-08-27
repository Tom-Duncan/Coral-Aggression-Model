# =============================================================================
#  RUN_ME.R  -  quick start / smoke test for the coral competition ABM.
#
#  Runs one small simulation end-to-end (a 3-species rock-paper-scissors community
#  on a 40x40 reef for 120 timesteps) and prints the resulting diversity metrics.
#  Takes a few seconds and needs BASE R ONLY - no extra packages.
#
#  Run it from a terminal:   Rscript RUN_ME.R
#  or in R / RStudio:        source("RUN_ME.R")
#
#  Dependencies:
#   * The demo below           : base R only.
#   * Plots / GIFs             : install.packages("animation")
#   * Batch runs & analysis    : install.packages(c("data.table", "ggplot2"))
#  Developed on R 4.6.0.
# =============================================================================

# --- Resolve the repo root (folder containing "Model/") so this runs anywhere ---
setwd((function() {
  a <- commandArgs(FALSE); f <- sub("^--file=", "", grep("^--file=", a, value = TRUE))
  d <- normalizePath(if (length(f)) dirname(f[1]) else getwd(), mustWork = FALSE)
  while (basename(d) != "Model" && !dir.exists(file.path(d, "Model")) && dirname(d) != d) d <- dirname(d)
  if (basename(d) == "Model") dirname(d) else d
})())
cat("Project root:", getwd(), "\n")

# --- Load the model engine ---------------------------------------------------
for (f in c("Size_Impact_Functions.r", "Intialisation_Functions.r",
            "Disturbance_Functions.r", "Sim_func.r", "Reproduction_Functions.r",
            "Visual_Functions.r", "Saving_Functions.r"))
  source(file.path("Model", f))
source("Model/Trait_Combinations.r")   # competition-network builders
cat("Model loaded.\n\n")

# --- Build one small scenario ------------------------------------------------
set.seed(1)                                   # reproducible
n_species <- 3
reef_x <- reef_y <- 40
founders  <- 3                                # founding colonies per species
n_steps   <- 120

# 3-species intransitive (rock-paper-scissors) interaction matrix
spec    <- makeTraitCombination("Classic_RPS", n_species, bias = 0.9, intra = 0.5)
species <- paste0("Sp", seq_len(n_species))

# species_traits: uniform growth, reproduction off, the RPS matrix attached
traits <- data.frame(species = species, growth = 3, stringsAsFactors = FALSE)
traits <- assignReproductionToSpecies(traits, "random", enabled = FALSE)
traits <- attachInteraction(traits, TRUE, spec$matrix,
                            size_mode = "none", size_beta = defaultSizeBeta(species))
traits <- attachGrowthSize(traits, FALSE, defaultGrowthGamma(species))

# place founders at random and create the colonies + reef
indiv  <- rep(founders, n_species)
coords <- getCoordinates(sum(indiv), reef_x, reef_y, placement_type = 2)
colony <- createCorals(reef_x, reef_y, n_species, indiv, coords, traits)

# no disturbance in this quick demo
disturb <- buildDisturbanceConfig(FALSE, "often", "random", n_steps)

# --- Run it (step output suppressed for a clean demo) ------------------------
cat("Running a", n_species, "-species RPS simulation for", n_steps, "timesteps...\n")
invisible(capture.output(
  states <- runSimulation(colony$reef, colony$corals, n_steps,
                          colony$colony_species, traits, disturb, habitat = NULL)
))

# --- Report ------------------------------------------------------------------
cp  <- mainCheckpointMetrics(states, interval = 20)
cat("\nDiversity through time (every 20 steps):\n")
print(cp[, c("timestep", "species_richness", "shannon", "evenness", "total_cover")],
      row.names = FALSE)

final <- cp[nrow(cp), ]
cat(sprintf("\nFinal state: %d of %d species coexisting (%.0f%% of the pool), %.0f%% reef cover.\n",
            final$species_richness, n_species,
            100 * final$species_richness / n_species, final$total_cover))
cat("\nSuccess - the model ran end to end. See README.md and Model/ for full runs.\n")
