# Coral-Aggression-Model

A spatially explicit **agent-based model (ABM)** of coral growth, competition,
life history and disturbance, used to study how community diversity is maintained
across different competitive interaction networks and between natural-reef and
closed (aquarium) conditions.

Competitive outcomes between colonies are governed by a species **interaction
matrix**; each timestep applies growth, competitive resolution, disturbance,
reproduction and external recruitment on a 2-D lattice "reef". The model
underlies the analyses reported in the associated MSc project.

## Quick start

Run one small simulation end to end (a few seconds, base R only):

```
Rscript RUN_ME.R          # from a terminal
```
or in R/RStudio: `source("RUN_ME.R")`. It runs a 3-species rock-paper-scissors
community and prints the resulting diversity metrics — a smoke test that the model
works. `RUN_ME.R` resolves the repo location itself, so it can be run from anywhere.

## Requirements

- **R** (developed and run on version 4.6.0)
- **The core engine and `RUN_ME.R` need base R only** — no packages to install.
- Optional, only for extra features:
  - `animation` — reef/abundance GIFs (`Model/Visual_Functions.r`)
  - `data.table`, `ggplot2` — batch runs, dataset consolidation and analysis tables
  ```r
  install.packages(c("animation", "data.table", "ggplot2"))
  ```
- All scripts resolve the repository root themselves, so they run from any working
  directory when launched with `Rscript`.

## Repository structure

```
Model/
├── Sim_func.r               # core simulation engine (growth, competition, timestep loop)
├── *_Functions.r            # process modules (disturbance, reproduction, size effects, saving, ...)
├── Trait_Combinations.r     # builds the theoretical networks (RPS, Linear, Neutral, Random)
├── Pairwise_Matrices.r      # interaction-matrix helpers
├── Master_Script.r / Run_Models.r   # entry points for running scenarios
├── Simulation_testing.r     # sources the engine and runs a scenario
├── Code_Specific_Sims/      # per-experiment run scripts (sensitivity, robustness, ...)
├── Interaction_Matrices/    # empirical / derived interaction matrices
└── Analysis_Code/           # analysis and summary scripts
    ├── Empirical_Matrices.r  # builds the empirical networks (Logan, Connell/Heron, ...)
    └── Habitat_Comparison.r  # natural-reef vs aquarium configurations + run drivers (Aim 2)
```

## Running the model

The scripts resolve their own paths, so the model runs unchanged on Windows or
Linux. A minimal run sources the engine and executes a scenario via the entry
scripts (`Master_Script.r` / `Run_Models.r`); the per-experiment scripts in
`Model/Code_Specific_Sims/` reproduce the sensitivity and robustness runs
described in the project.

## Data availability

Simulation **output** (result datasets) is archived separately (see the project's
data-availability statement) rather than in this code repository. The empirical
interaction data underlying the Logan and Connell (Heron Island) networks derive
from published/archived sources cited in the project.

**Note:** the Webden (London Zoo) interaction matrix is unpublished and is
therefore not included in this repository; it is available directly from its
author (M. Webden) on reasonable request.
