# Run this script from the repository root (the folder containing "Model/").
# Source paths below are relative to that directory.

source("Model/Size_Impact_Functions.r")     # cell count -> % cover
source("Model/Intialisation_Functions.r")   # setup, colony creation, overgrowth matrix
source("Model/Disturbance_Functions.r")     # disturbance events
source("Model/Sim_func.r")                  # simulation core + colony split/merge
source("Model/Reproduction_Functions.r")    # reproduction + background recruitment
source("Model/Visual_Functions.r")          # plots, gifs, colony-size summaries
source("Model/Saving_Functions.r")          # checkpoint results table + saving

# -------- Example run --------

# Setup + simulation in one call
sim <- reefSetUp()
states <- sim$states

states[[1]]$corals                    # details of coral 1
states[[25]]$corals[[1]]$size         # coral 1's cover at timestep 25

plotReefAtTimestep(states, 250)       # reef at a timestep

plotInteractionMatrix(sim$species_traits)   # the figure
printInteractionMatrix(sim$species_traits)  # text version

make_coral_gif(states)                                        # growth gif
plotIndiAbundance(states, states[[1]]$corals)                 # per-colony abundance
plotSpeciesAbundance(states, states[[1]]$corals, sim$colony_species)  # per-species abundance
make_abundance_gif(
  states,
  states[[1]]$corals,
  file = "coral_abundance.gif",
  interval = 0.1
)

# Colony sizes over time + summary table
colony_sizes <- getColonySizeTimeSeries(states, sim$corals)
head(colony_sizes, 20)
size_summary <- summarizeColonySizes(colony_sizes)
print(size_summary)
