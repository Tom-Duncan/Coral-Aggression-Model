#Source files containing necessary functions

# Run this script from the repository root (the folder that contains "Model/").
# All source() paths below are relative to that directory.

#Size helper: colony cell count -> % reef cover (for reproduction + metrics)
source("Model/Size_Impact_Functions.r")

#Setup, colony creation, coordinate shapes, and the pairwise overgrowth matrix
source("Model/Intialisation_Functions.r")

#Disturbance events
source("Model/Disturbance_Functions.r")

#Simulation core (growth, competition, timestep) + colony split/merge tracking
source("Model/Sim_func.r")

#Coral reproduction (colonies spawn new recruits) + habitat / background recruitment
source("Model/Reproduction_Functions.r")

#Plots, gifs and colony-size summaries
source("Model/Visual_Functions.r")

#Per-run data collection, the flat checkpoint results table, and saving
source("Model/Saving_Functions.r")

#--------RUNNING CODE----------

#Runs intialisation as well as the simulation in one line
sim <- reefSetUp()

#Assigns values of the simulation to other variable for ease later
states <- sim$states

#Able to view the details of coral "1"
states[[1]]$corals

#Able to view the % of the reef covered by coral "1" at a given timestep
states[[25]]$corals[[1]]$size

#Visualise the reef at a given timestep 
plotReefAtTimestep(states, 250)

#Table of each species traits
plotInteractionMatrix(sim$species_traits)   # the figure
printInteractionMatrix(sim$species_traits)  # text version

#To get a gif of the coral growth
make_coral_gif(states)

#To get graph of changes in abundance for all individuals
plotIndiAbundance(states, states[[1]]$corals)

#To get pop change for species
plotSpeciesAbundance(states, states[[1]]$corals, sim$colony_species)

#To make a gif of the population changes in species abundance graph above
make_abundance_gif(
  states,
  states[[1]]$corals,
  file = "coral_abundance.gif",
  interval = 0.1
)


# Extract size data for all colonies across all timesteps
colony_sizes <- getColonySizeTimeSeries(states, sim$corals)

# View the data
head(colony_sizes, 20)

# Create summary table
size_summary <- summarizeColonySizes(colony_sizes)

# View summary table
print(size_summary)
