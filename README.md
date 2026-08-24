# Coral-Aggression-Model

A spatially explicit **agent-based model (ABM)** of coral growth, competition,
life history and disturbance, used to study how community diversity is maintained
across different competitive interaction networks and between natural-reef and
closed (aquarium) conditions.

Competitive outcomes between colonies are governed by a species **interaction
matrix**; each timestep applies growth, competitive resolution, disturbance,
reproduction and external recruitment on a 2-D lattice "reef". The model
underlies the analyses reported in the associated paper.

## Requirements

- **R** (developed and run on version 4.6.0)
- Base R only for the core engine; the plotting/analysis scripts use common CRAN
  packages (e.g. `data.table`, `ggplot2`) where noted at the top of each script.

## Repository structure

```
Current_Working_Model/
├── Sim_func.r               # core simulation engine (growth, competition, timestep loop)
├── *_Functions.r            # process modules (disturbance, reproduction, size effects, saving, ...)
├── Trait_Combinations.r     # builds the theoretical networks (RPS, Linear, Neutral, Random)
├── Pairwise_Matrices.r      # interaction-matrix helpers
├── Master_Script.r / Run_Models.r   # entry points for running scenarios
├── Interaction_Matrices/    # empirical / derived interaction matrices
└── Plotting_Code/           # analysis, summary and figure scripts
    ├── Real_Int_Mat.r        # empirical matrix builders (Logan, Connell/Heron, ...)
    └── Final_Sim_Prep.r      # natural-reef vs aquarium configuration + run drivers
```

## Running the model

The scripts resolve their own paths, so the model runs unchanged on Windows or
Linux. A minimal run sources the engine and executes a scenario via the entry
scripts (`Master_Script.r` / `Run_Models.r`); the per-experiment scripts in
`Current_Working_Model/Code_Specific_Sims/` reproduce the sensitivity, robustness
and configuration runs described in the paper.

## Data availability

Simulation **output** (result datasets) is archived separately (see the paper's
data-availability statement) rather than in this code repository. The empirical
interaction data underlying the Logan and Connell (Heron Island) networks derive
from published/archived sources cited in the paper.

**Note:** the Webden (London Zoo) interaction matrix is unpublished and is
therefore not included in this repository; it is available directly from its
author (M. Webden) on reasonable request.

## Licence

Released under the MIT Licence (see `LICENSE`).
