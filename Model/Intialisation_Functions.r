#Functions related to intialisation process of the model

#Function checking inputs are valid from intialisation setup answers
#Prevents invalid numbers, allow reusable input
#prompt = text shown to user, min and max values able to be inputted by user
getValidInteger <- function(prompt, min_val = -Inf, max_val = Inf) {

   #Infinite loop displaying prompt
  repeat {value <- suppressWarnings(as.integer(readline(prompt)))
    #Checks if real integer within allowed range, so both valid and appropriate
    if(
      !is.na(value) &&
      value >= min_val &&
      value <= max_val) {
      return(value)}

    #If input fails, prints this error message
    cat("Invalid input. Please enter a number between",
      min_val,"and",max_val,"\n")
      }
}

#Helper to show a numbered menu and return the chosen option number.
#The menu is printed with cat first, then only a short "Enter choice:" prompt is
#given to readline. Passing long multi-line text straight to readline is what
#causes the console to show just a "$" with no visible typing, so we avoid it.
#Every setup choice goes through here so the numbering style stays consistent.
getMenuChoice <- function(title, options) {
  cat("\n", title, "\n", sep = "")
  for (i in seq_along(options)) {
    cat("  ", i, ") ", options[i], "\n", sep = "")
  }
  getValidInteger("Enter choice: ", min_val = 1, max_val = length(options))
}

#Helper to ask a yes/no question as a numbered choice (1 = Yes, 2 = No).
#Keeps every yes/no input numeric and consistent. Returns TRUE for yes.
getYesNo <- function(question) {
  getMenuChoice(question, c("Yes", "No")) == 1
}

#Helper to select which parameters should remain fixed across replicates
getKeepSameOptions <- function() {
  prompt_text <- paste0(
    "Choose which setup parameters should stay the same across all replicates:\n",
    " 1) Reef size\n",
    " 2) Number of species\n",
    " 3) Individuals per species\n",
    " 4) Simulation length\n",
    " 5) Species traits and identity\n",
    " 0) None\n"
  )

  repeat {
    cat(prompt_text)
    choice <- readline("Enter numbers separated by commas (e.g. 1,3,5): ")
    choice <- gsub("[^0-9,]", "", choice)
    selected <- as.integer(unlist(strsplit(choice, ",")))
    if (length(selected) == 0 || all(is.na(selected))) {
      cat("Invalid input. Please enter one or more digits separated by commas (e.g. 1,3,5).\n")
      next
    }
    selected <- unique(selected[!is.na(selected) & selected != 0])
    if (length(selected) == 0) {
      return(list(
        keep_reef_size = FALSE,
        keep_n_species = FALSE,
        keep_individuals = FALSE,
        keep_sim_length = FALSE,
        keep_traits = FALSE
      ))
    }
    if (any(selected < 0 | selected > 5)) {
      cat("Invalid option. Enter values from 0 to 5.\n")
      next
    }
    keep_traits <- 5 %in% selected
    keep_n_species <- 2 %in% selected || keep_traits
    keep_individuals <- 3 %in% selected || keep_traits
    return(list(
      keep_reef_size = 1 %in% selected,
      keep_n_species = keep_n_species,
      keep_individuals = keep_individuals,
      keep_sim_length = 4 %in% selected,
      keep_traits = keep_traits
    ))
  }
}

# Gather all setup choices from the user, in three sections: Environment, Community,
# Run. Fully-random setup samples everything. Returns the info list.
getSpeciesInfo <- function() {

   #Setup mode decides how every following question is asked
   #Manual = all own input on every aspect
   #Semi-manual = input on a range of values for each aspect, sampled from
   #Random = all values random from a sensible default range, no questions
  setup_mode <- getMenuChoice(
    "Setup mode:",
    c("Manual", "Semi-Manual", "Fully Random")
  )

  #Fully random: ask nothing about the reef or coral, just sample everything.
  #Only the cross-cutting choices (replicates, disturbance) are gathered.
  if (setup_mode == 3) {
    info <- getRandomSpeciesInfo()
    info$n_replicates <- getValidInteger("How many replicates / simulation runs? ", min_val = 1)
    info$randomize_per_replicate <- TRUE
    info$keep_same <- list(
      keep_reef_size = FALSE, keep_n_species = FALSE, keep_individuals = FALSE,
      keep_sim_length = FALSE, keep_traits = FALSE
    )
    info$disturbance     <- getDisturbanceSetup(setup_mode)
    info$reproduction    <- getReproductionSetup(setup_mode)
    info$habitat         <- getHabitatSetup(setup_mode)
    #Interaction matrix is always used; fill it randomly in fully-random setup
    info$interaction_mode <- "random"
    info$interaction_size_mode <- "none"
    #Growth size effect switched on/off at random
    info$growth_size_enabled <- getGrowthSizeSetup(setup_mode)
    return(info)
  }

  #Manual (1) and semi-manual (2): ask questions in three grouped sections.
  #Semi-manual asks for a min/max range for each value and samples from it.
  info <- list(setup_mode = setup_mode)

  # ================= ENVIRONMENT =================
  cat("\n=== Environment ===\n")

  #Reef dimensions. Width (x) and height (y) are asked separately so the reef
  #can be a rectangle (enter the same value for both to keep it square). Fully
  #random setup always makes a square - see getRandomSpeciesInfo().
  if (setup_mode == 1) {
    info$reef_x <- getValidInteger("Reef width (x)? ", min_val = 5)
    info$reef_y <- getValidInteger("Reef height (y)? ", min_val = 5)
  } else {
    info$reef_x_min <- getValidInteger("Minimum reef width (x): ", min_val = 5)
    info$reef_x_max <- getValidInteger("Maximum reef width (x): ", min_val = info$reef_x_min)
    info$reef_y_min <- getValidInteger("Minimum reef height (y): ", min_val = 5)
    info$reef_y_max <- getValidInteger("Maximum reef height (y): ", min_val = info$reef_y_min)
    info$reef_x <- sample(info$reef_x_min:info$reef_x_max, 1)
    info$reef_y <- sample(info$reef_y_min:info$reef_y_max, 1)
  }

  #Simulation length
  if (setup_mode == 1) {
    info$Sim_length <- getValidInteger("Simulation length (timesteps)? ", min_val = 1)
  } else {
    info$sim_min    <- getValidInteger("Minimum simulation length: ", min_val = 1)
    info$sim_max    <- getValidInteger("Maximum simulation length: ", min_val = info$sim_min)
    info$Sim_length <- sample(info$sim_min:info$sim_max, 1)
  }

  #Disturbance events (reef-level). Only the settings are gathered here; the
  #per-run schedule is built later once the simulation length is known.
  info$disturbance <- getDisturbanceSetup(setup_mode)

  #Habitat type. Controls background recruitment: a natural reef receives larvae
  #from outside the patch, an aquarium (or none) does not. See getHabitatSetup /
  #backgroundReproduce in Reproduction.r.
  info$habitat <- getHabitatSetup(setup_mode)

  # ================= COMMUNITY =================
  cat("\n=== Community ===\n")

  #Number of species
  if (setup_mode == 1) {
    info$n_species <- getValidInteger("Number of coral species? ", min_val = 1)
  } else {
    info$species_min <- getValidInteger("Minimum number of species: ", min_val = 1)
    info$species_max <- getValidInteger("Maximum number of species: ", min_val = info$species_min)
    info$n_species   <- sample(info$species_min:info$species_max, 1)
  }

  #Individuals per species
  if (setup_mode == 1) {
    individuals <- integer(info$n_species)
    for (i in seq_len(info$n_species)) {
      individuals[i] <- getValidInteger(paste("Individuals for species", i, ": "), min_val = 1)
    }
    info$individuals_per_species <- individuals
  } else {
    info$indiv_min <- getValidInteger("Minimum individuals per species: ", min_val = 1)
    info$indiv_max <- getValidInteger("Maximum individuals per species: ", min_val = info$indiv_min)
    info$individuals_per_species <- sample(info$indiv_min:info$indiv_max, info$n_species, replace = TRUE)
  }
  info$total_individuals <- sum(info$individuals_per_species)

  #Grow the reef if it is too small to place all the colonies. Works on both
  #dimensions (growing the shorter side first) so a rectangle stays a rectangle.
  grown <- growReefToFit(info$reef_x, info$reef_y, info$total_individuals)
  if (grown[1] != info$reef_x || grown[2] != info$reef_y) {
    if (setup_mode == 1) {
      cat("Reef enlarged to", grown[1], "x", grown[2], "to fit all",
          info$total_individuals, "colonies.\n")
    }
    info$reef_x <- grown[1]
    info$reef_y <- grown[2]
  }

  #Species growth values: enter each species by hand, or randomise them. Offered
  #in semi-manual too, so you can hand-pick the species you want and still have
  #the rest of the reef (size, length, individuals) randomised from ranges.
  info$trait_mode <- getMenuChoice("Species growth:", c("Manual", "Random"))

  #Growth size effect: does a colony's growth rate shift with its OWN size? A flat
  #per-species small/large adjustment in growth-trait points, mirroring the interaction
  #size effect (see growthChanceSized). Toggleable; adjustments are randomised when on.
  info$growth_size_enabled <- getGrowthSizeSetup(setup_mode)

  #Reproduction: whether colonies can spawn new recruits. Toggled on/off here; if
  #on, each species' reproduction chance and reference size are set later (asked in
  #manual mode, randomised otherwise) by assignReproductionToSpecies. (Maturity
  #itself is age-based, not size-based - see MATURITY_AGE in Reproduction.r.)
  info$reproduction <- getReproductionSetup(setup_mode)

  #Pairwise interaction matrix: the sole interaction mechanism. Every species pair
  #(and each species with itself) has an overgrowth PROBABILITY; a 0 means "never"
  #(neutrality). Choose how the matrix is filled and whether colony size shifts the
  #probabilities. See Interaction_Matrix.r.
  info$interaction_mode <-
    c("manual", "random")[getMenuChoice("Interaction matrix values:", c("Manual", "Random"))]
  info$interaction_size_mode <- getInteractionSizeSetup(setup_mode)

  # ================= RUN =================
  cat("\n=== Run ===\n")

  #Placement of colonies (manual chooses; semi-manual is random)
  if (setup_mode == 1) {
    info$placement_type <- getMenuChoice("Placement type:", c("Manual", "Random", "Shape"))
  } else {
    info$placement_type <- 2
  }

  #Number of replicates, and (semi-manual only) whether to re-randomise each one
  info$n_replicates <- getValidInteger("How many replicates / simulation runs? ", min_val = 1)
  if (setup_mode == 1) {
    info$randomize_per_replicate <- FALSE
    info$keep_same <- list(
      keep_reef_size = TRUE, keep_n_species = TRUE, keep_individuals = TRUE,
      keep_sim_length = TRUE, keep_traits = TRUE
    )
  } else {
    info$randomize_per_replicate <- getYesNo("Randomize setup parameters for each replicate?")
    if (info$randomize_per_replicate) {
      info$keep_same <- getKeepSameOptions()
      #If the species were entered by hand, keep them (and the species count) fixed
      #across replicates so they are not re-typed every run; only the rest of the
      #reef setup is re-randomised.
      if (info$trait_mode == 1) {
        info$keep_same$keep_traits    <- TRUE
        info$keep_same$keep_n_species <- TRUE
      }
    } else {
      info$keep_same <- list(
        keep_reef_size = FALSE, keep_n_species = FALSE, keep_individuals = FALSE,
        keep_sim_length = FALSE, keep_traits = FALSE
      )
    }
  }

  info
}




#List of range for random intialisation to sample from to ensure reasonable range
#Can be changed 
DEFAULTS <- list(
  reef_min = 5,
  reef_max = 100,
  species_min = 1,
  species_max = 15,
  indiv_min = 1,
  indiv_max = 5,
  sim_min = 50,
  sim_max = 250)


# Smallest reef size holding at least 2x the colonies (so it stays ~half full and
# random placement is fast). Never below absolute_min.
minReefSizeForIndividuals <- function(total_individuals, absolute_min = 5) {
  reef_size <- absolute_min
  while (ceiling(reef_size / 2)^2 < 2 * total_individuals) {
    reef_size <- reef_size + 1
  }
  reef_size
}


# Spaced-apart colony capacity of an X by Y reef (~one per 2x2 block).
reefCapacity <- function(reef_x, reef_y) {
  ceiling(reef_x / 2) * ceiling(reef_y / 2)
}


# Grow a reef until it holds 2x the colonies (shorter side first). Returns c(x, y).
growReefToFit <- function(reef_x, reef_y, total_individuals, absolute_min = 5) {
  reef_x <- max(reef_x, absolute_min)
  reef_y <- max(reef_y, absolute_min)
  while (reefCapacity(reef_x, reef_y) < 2 * total_individuals) {
    if (reef_x <= reef_y) {
      reef_x <- reef_x + 1
    } else {
      reef_y <- reef_y + 1
    }
  }
  c(reef_x, reef_y)
}


#Function if Random intialisation selected
#Samples from Default range to give necessary data to be used later
#No questions asked — all values drawn automatically from DEFAULTS
getRandomSpeciesInfo <- function() {
  reef_size <- sample(DEFAULTS$reef_min:DEFAULTS$reef_max,1)
  n_species <- sample(DEFAULTS$species_min:DEFAULTS$species_max,1)
  individuals_per_species <- sample(DEFAULTS$indiv_min:DEFAULTS$indiv_max, n_species, replace = TRUE)
  Sim_length <- sample(DEFAULTS$sim_min:DEFAULTS$sim_max,1)
  #Automatic random placement (placement: 1 = Manual, 2 = Random, 3 = Shape)
  placement_type <- 2
  total_individuals <- sum(individuals_per_species)
  #Grow the reef if it is too small for all the colonies, so random placement
  #never runs out of room (prevents the "reef too crowded" failure)
  reef_size <- max(reef_size, minReefSizeForIndividuals(total_individuals))
  #Automatic random growth (1 = Manual, 2 = Random)
  trait_mode <- 2
  list(
    #Fully random reefs are always square, so width and height are equal
    reef_x = reef_size,
    reef_y = reef_size,
    n_species = n_species,
    individuals_per_species = individuals_per_species,
    total_individuals = total_individuals,
    placement_type = placement_type,
    Sim_length = Sim_length,
    trait_mode = 2,
    setup_mode = 3
  )
}

#Function to print summary of all intialisation data
printSetupSummary <- function(info) {

   #Simply calls and prints previously stored data 
  cat("\n--- Simulation Setup Summary ---\n")
  if (info$setup_mode == 1) {
    cat("Setup mode: Manual\n")
  } else if (info$setup_mode == 2) {
    cat("Setup mode: Semi-Manual\n")
  } else if (info$setup_mode == 3) {
    cat("Setup mode: Fully Random\n")
  }

  #Grouped to match the setup question order: Environment, Community, Run.
  cat("\n[Environment]\n")
  cat("Reef size:", info$reef_x, "x", info$reef_y,
      if (info$reef_x == info$reef_y) "(square)" else "(rectangle)", "\n")
  cat("Simulation length:", info$Sim_length, "timesteps\n")
  cat("Disturbances:", describeDisturbance(info$disturbance), "\n")
  cat("Habitat:", describeHabitat(info$habitat), "\n")

  cat("\n[Community]\n")
  cat("Number of species:", info$n_species, "\n")
  cat("Individuals per species:", paste(info$individuals_per_species, collapse = ", "), "\n")
  cat("Total individuals:", info$total_individuals, "\n")
  if (info$trait_mode == 1) {
    cat("Growth mode: Manual\n")
  } else if (info$trait_mode == 2) {
    cat("Growth mode: Random\n")
  }
  cat("Growth size effect:", ifelse(isTRUE(info$growth_size_enabled), "On", "Off"), "\n")
  cat("Reproduction:", ifelse(isTRUE(info$reproduction), "On", "Off"), "\n")
  cat("Interaction matrix: ",
      ifelse(identical(info$interaction_mode, "manual"), "manual", "random"),
      " values, ",
      switch(if (is.null(info$interaction_size_mode)) "none" else info$interaction_size_mode,
             ratio = "ratio size effect", logodds = "log-odds size effect",
             additive = "additive size effect", "no size effect"),
      "\n", sep = "")

  cat("\n[Run]\n")
  if (info$placement_type == 1) {
    cat("Placement type: Manual\n")
  } else if (info$placement_type == 2) {
    cat("Placement type: Random\n")
  } else if (info$placement_type == 3) {
    cat("Placement type: Shape\n")
  }
  cat("Number of replicates:", info$n_replicates, "\n")
  cat("Randomize setup per replicate:", ifelse(info$randomize_per_replicate, "Yes", "No"), "\n")
  cat("--------------------------------\n\n")
}

#Function that assigns coral individuals
#coordinates in reef based on number, reef size and placement type
#Follows simple IF function for each type of placement
getCoordinates <- function(total_individuals,
                            reef_x,
                            reef_y,
                            placement_type){
    
        #If RANDOM placement (placement_type 2)
        #Coordinates passed as an empty list to be filled
        #Stores coordinates of all individuals
        if (placement_type == 2){
            coords <- list()
            attempts <- 0
            max_attempts <- 100000
                    
                    #Keep generating coordinates until enough for all individuals
                    while(length(coords) < total_individuals){

                        attempts <- attempts + 1

                        if(attempts > max_attempts){
                            stop(
                            "Unable to place all colonies. Reef too crowded."
                            )
                        }

                        #Generates a random "candidate" coordinate as c(x, y),
                        #with x drawn across the width and y across the height so
                        #a rectangular reef is filled correctly
                        candidate <- c(
                            sample(1:reef_x, 1),
                            sample(1:reef_y, 1)
                        )
                        #Assumes candidate is valid
                        valid <- TRUE

                        #Checks if any coordinates already exist
                        #If list empty, nothing to compare to so doesnt execute
                        if(length(coords) > 0){
                            
                            #Loops through existing coordinates
                            for(existing in coords){
                                
                                #Calculates vertical and horizontal distance of two coordinates
                                dx <- abs(candidate[1] - existing[1])
                                dy <- abs(candidate[2] - existing[2])

                                #Checks if too close/overlapping intial squares 
                                #Candidate rejected if within a 5x5 square centred on existing individual
                                if(dx <= 1 && dy <= 1){

                                    #If not spaced enough, invalid
                                    #Exit for loop, don't have to compare to others yet
                                    valid <- FALSE
                                    break
                                }
                            }
                        }
                        #Adds candidate to list if valid
                        if(valid){

                            coords[[length(coords)+1]] <- candidate
                        }
                    }

        #If MANUAL placement (placement_type 1)
        } else if (placement_type == 1){

            #Creates empty list to hold all individuals, initally empty
            coords <- vector("list", total_individuals)

            #Loops through all individuals
            #Asks for a certain corals x and y coordinates
            for(i in seq_len(total_individuals)){

                x <- as.integer(
                    readline(
                        paste("Coral",i,"x coordinates: ")))

                y <- as.integer(
                    readline(
                        paste("Coral",i,"y coordinates: ")))

                #Stores a coordinate pair for a given coral
                coords[[i]] <- c(x,y)}

        #If SHAPE placement
        } else if(placement_type == 3) {

            #Shape layouts are square by design, so fit them inside the largest
            #square the reef can hold (the shorter side of a rectangle)
            shape_size <- min(reef_x, reef_y)

            #Runs through possible number of total individuals
            #Call the shape coordinate functions in other document
            if(total_individuals == 2){

                coords <- linecoordinates(shape_size)

            } else if(total_individuals == 3){

                coords <- trianglecoordinates(shape_size)

            } else if(total_individuals == 4){

                coords <- squarecoordinates(shape_size)

            } else if(total_individuals >= 5) {

                coords <- polygoncoordinates(
                    shape_size,
                    total_individuals
                )


}
        }
return(coords)
}


# Clustered placement: all founders packed into one ~1/3-size box around a random
# centre, with the same minimum spacing as random placement. Returns c(x, y) pairs.
clusteredCoordinates <- function(total_individuals, reef_x, reef_y) {
  cw <- max(3L, floor(reef_x / 3)); ch <- max(3L, floor(reef_y / 3))
  cx <- sample(seq_len(reef_x), 1); cy <- sample(seq_len(reef_y), 1)
  x_lo <- max(1L, cx - cw %/% 2L); x_hi <- min(reef_x, x_lo + cw - 1L)
  y_lo <- max(1L, cy - ch %/% 2L); y_hi <- min(reef_y, y_lo + ch - 1L)

  coords <- list(); attempts <- 0L; max_attempts <- 100000L
  while (length(coords) < total_individuals) {
    attempts <- attempts + 1L
    if (attempts > max_attempts) {
      stop("Unable to place all colonies in the cluster region (too crowded).")
    }
    candidate <- c(sample(x_lo:x_hi, 1), sample(y_lo:y_hi, 1))
    valid <- TRUE
    for (existing in coords) {
      if (abs(candidate[1] - existing[1]) <= 1 && abs(candidate[2] - existing[2]) <= 1) {
        valid <- FALSE; break
      }
    }
    if (valid) coords[[length(coords) + 1L]] <- candidate
  }
  coords
}


# Spread placement: over-dispersed by farthest-point sampling (each founder is the
# best of a candidate batch, maximising min distance to those placed). Returns c(x, y).
spreadCoordinates <- function(total_individuals, reef_x, reef_y, n_candidates = 60L) {
  coords <- list(c(sample(seq_len(reef_x), 1), sample(seq_len(reef_y), 1)))

  while (length(coords) < total_individuals) {
    best <- NULL; best_d <- -1
    for (k in seq_len(n_candidates)) {
      cand <- c(sample(seq_len(reef_x), 1), sample(seq_len(reef_y), 1))
      dmin <- Inf; ok <- TRUE
      for (existing in coords) {
        dx <- cand[1] - existing[1]; dy <- cand[2] - existing[2]
        if (abs(dx) <= 1 && abs(dy) <= 1) { ok <- FALSE; break }
        d <- dx * dx + dy * dy
        if (d < dmin) dmin <- d
      }
      if (ok && dmin > best_d) { best_d <- dmin; best <- cand }
    }
    if (is.null(best)) {           #fallback: take any spacing-valid random cell
      repeat {
        cand <- c(sample(seq_len(reef_x), 1), sample(seq_len(reef_y), 1)); ok <- TRUE
        for (existing in coords) if (abs(cand[1]-existing[1]) <= 1 && abs(cand[2]-existing[2]) <= 1) { ok <- FALSE; break }
        if (ok) { best <- cand; break }
      }
    }
    coords[[length(coords) + 1L]] <- best
  }
  coords
}


#Function for creating individual corals
#Takes all previous data generated to do so
createCorals<- function(reef_x,
                        reef_y,
                        n_species,
                        individuals_per_species,
                        coords_sample,
                        species_traits){

    #Create empty reef matrix (may be rectangular). Rows are the y-axis (height)
    #and columns are the x-axis (width), so cells are indexed reef[y, x].
    reef <- matrix(NA,
    nrow = reef_y,
    ncol = reef_x)

    #Species named after letters in order (species 1 = A)
    species_names <- paste0("Sp", seq_len(n_species))

    #Define species colours
    species_colour <- setNames(
    hcl.colors(n_species, "Dark 3"),
    species_names)

    #Creates individual colours for multiple individuals of same species
    #Empty character vector will eventially contain one colour for every 
    #individual coral
    individual_colours <- character(0)

    #Loops through each species
    for(sp in species_names){

        #Number individuals per species
        n_ind <- individuals_per_species[
            match(sp, species_names)]
        
        #Generates transparency values
        alpha_values <- seq(
            0.4,1,
            length.out = n_ind
        )

        #Creates shades for each colour for each species based on transparency values
        individual_colours <- c(
            individual_colours,
            sapply(alpha_values,(function(a)
            adjustcolor(
                species_colour[sp],
                alpha.f = a
            ))))
    }

    #Creates an ID for each species
    #Species vector becomes an ordered list of all individual coral species
    #In the form of A,A,A,B,B,C etc
    species_vector <- rep(
        species_names,
        times = individuals_per_species
    )
        #Makes a numeric vector of 0s, equal length to number of species
        #Labels 0 with species names
        #Prepares counter for number of individuals per species
        species_counts <- setNames(
        rep(0, n_species),
        species_names
    )
        #Creates empty character vector sized to number of coral individuals
        #Stores unique ids in the form of "A1, A2, B1" etc
        id_vector <- character(length(species_vector))

        #Loop assigns unique ID to each individual coral
        for(i in seq_along(species_vector)){

            #Gets species letter for current coral
            sp <- species_vector[i]

            #Increments count for that species
            species_counts[sp] <- species_counts[sp] + 1

            #Build string of IDs 
            id_vector[i] <-paste0(sp, "_", species_counts[sp])
        }

        colony_species <- setNames(species_vector, id_vector)

    #Counts number of individuals expected, should match coordinate pairs
    #If not, stops, prevents mismatch of number of coral and placement data
    total_individuals <- length(species_vector)
    if(length(coords_sample) != total_individuals){
        stop("Number of coordinates must equal total number of individuals")}

    #Creates empty list assigned corals to store coral data
    corals <- vector("list", total_individuals)

    #Loop to process all coral, placing in reef and creating its object
    for(i in seq_len(total_individuals)){

        #Species of coral, and its coordinates
        sp <- species_vector[i]
        x <- coords_sample[[i]][1]
        y <- coords_sample[[i]][2]

        trait_row <- species_traits[
             species_traits$species == sp,
            ]

        #Boundary check, makes sure the colony sits within the reef (x = width,
        #y = height, which may differ on a rectangular reef)
        if(
            x < 1 ||
            x > reef_x  ||
            y < 1 ||
            y > reef_y
        ){
            stop(
                paste(
                    "Coral",
                    i,
                    "is too close to reef boundary"
                )
            )
        }

        #Overlap check to see if two coral not touching
        if(!is.na(reef[y, x])) stop("Overlap detected")

        #Places colony on reef, filling 3x3 space with that corals ID
        #Its physical placement within the matrix
        colony_id <- id_vector[i]
                reef[y, x] <- colony_id

            #Create coral object and assign characteristics: its unique ID, species,
            #growth trait, colour, initial coordinates and size. Interactions are
            #governed entirely by the species interaction matrix (Interaction_Matrix.r),
            #so no per-colony combat traits are stored.
            corals[[i]] <- list(
                id = colony_id,
                species = sp,

                growth = trait_row$growth,

                colour = individual_colours[i],
                coords = coords_sample[[i]],
                #Size is a cell count; a founding colony is a single cell
                size = 1L,
                #Age in timesteps alive; drives maturity (see MATURITY_AGE). Founders
                #start at 0 and are aged one step per timestep in timestep().
                age = 0L
            )
            }

    #Returns both corals creates and the reef matrix
    list(
        corals = corals,
        reef = reef,
        colony_species = colony_species
    )
}


# Build species_traits: a growth trait per species (manual or random), reproduction
# traits, and the interaction matrix (the sole interaction mechanism) attached.
createSpeciesTraits <- function(n_species, trait_mode,
                                setup_mode = 3, reproduction_enabled = FALSE,
                                interaction_mode = "random",
                                interaction_matrix = NULL,
                                interaction_size_mode = "none",
                                growth_size_enabled = FALSE){

    species_names <- paste0("Sp", seq_len(n_species))
    growth        <- integer(n_species)

    #Growth trait per species: entered by hand (1-GROWTH_MAX) or drawn at random.
    for (i in seq_len(n_species)) {
        if (trait_mode == 1) {
            growth[i] <- getValidInteger(
                paste0("Growth for ", species_names[i], " (1-", GROWTH_MAX, "): "),
                min_val = 1, max_val = GROWTH_MAX)
        } else {
            growth[i] <- sample(1:GROWTH_MAX, 1)
        }
    }

    species_traits <- data.frame(
        species = species_names,
        growth  = growth,
        stringsAsFactors = FALSE
    )

    #Give each species its reproduction traits (can it reproduce, how often, and
    #the reference size that scales its spawn chance). Maturity is age-based, set by
    #MATURITY_AGE. Only asks/assigns if reproduction is enabled; manual growth mode
    #asks per species, otherwise randomised. See Reproduction.r.
    repro_label <- if (trait_mode == 1) "manual" else "random"
    species_traits <- assignReproductionToSpecies(species_traits, repro_label, reproduction_enabled)

    #The pairwise interaction matrix is the sole interaction mechanism and is always
    #built. A supplied matrix (e.g. empirical overgrowth data) is passed straight
    #through. It rides along on species_traits as an attribute so it reaches
    #canTakeCell/resolveClaims without changing any function signatures.
    M <- buildInteractionMatrix(
        species_names,
        mode     = if (!is.null(interaction_matrix)) "supplied" else interaction_mode,
        supplied = interaction_matrix
    )
    #For the log-ratio size mode, give each species a size sensitivity beta (random
    #here; set precise values programmatically via attachInteraction).
    size_beta <- if (identical(interaction_size_mode, "logratio")) {
        randomSizeBeta(species_names)
    } else {
        defaultSizeBeta(species_names)
    }
    species_traits <- attachInteraction(species_traits, TRUE, M,
                                         size_mode = interaction_size_mode,
                                         size_beta = size_beta)

    #Optional growth size effect: a per-species log-ratio shift to the growth chance
    #based on the colony's own size (see attachGrowthSize / growthChanceSized).
    #Randomised when enabled, otherwise no effect.
    growth_gamma <- if (isTRUE(growth_size_enabled)) {
        randomGrowthGamma(species_names)
    } else {
        defaultGrowthGamma(species_names)
    }
    species_traits <- attachGrowthSize(species_traits, growth_size_enabled,
                                       growth_gamma)

    species_traits
}


# Run the entire reef setup, across one or more replicates.
reefSetUp <- function(){
    #Gather every setup choice from the user (reef parameters first, then coral)
    info <- getSpeciesInfo()

    #Cross-cutting choices decided during setup, reused for every replicate
    disturbance_settings <- info$disturbance
    reproduction_enabled <- isTRUE(info$reproduction)
    habitat_settings     <- info$habitat   #controls background recruitment (reef only)
    interaction_mode     <- if (is.null(info$interaction_mode)) "random" else info$interaction_mode
    interaction_size_mode <- if (is.null(info$interaction_size_mode)) "none" else info$interaction_size_mode
    growth_size_enabled  <- isTRUE(info$growth_size_enabled)

    #Build the base species traits. For manual mode this is where the per-species
    #growth and reproduction questions are asked, so it runs before the single setup
    #summary below.
    species_traits <- createSpeciesTraits(
        info$n_species,
        info$trait_mode,
        info$setup_mode,
        reproduction_enabled,
        interaction_mode,
        interaction_matrix = NULL,
        interaction_size_mode = interaction_size_mode,
        growth_size_enabled = growth_size_enabled
    )

    #Print the full setup summary ONCE, now that every question has been answered
    printSetupSummary(info)

    #Initialize list to store all replicates
    all_replicates <- vector("list", info$n_replicates)
    
    #Run simulation for each replicate
    for (rep in 1:info$n_replicates) {
        
        cat("\n--- Running Replicate", rep, "of", info$n_replicates, "---\n")
        
        # Start with the base setup info
        current_info <- info
        current_traits <- species_traits
        
        # Re-sample parameters for each replicate if requested.
        # Replicate 1 always uses the setup that was shown in the summary; only
        # replicates 2+ are re-randomised, so the printed summary matches the run.
        if (info$randomize_per_replicate && rep > 1) {
            if (info$setup_mode == 2) {

                if (info$keep_same$keep_reef_size) {
                    new_reef_x <- info$reef_x
                    new_reef_y <- info$reef_y
                } else {
                    new_reef_x <- sample(info$reef_x_min:info$reef_x_max, 1)
                    new_reef_y <- sample(info$reef_y_min:info$reef_y_max, 1)
                }

                new_n_species <- if (info$keep_same$keep_n_species) {
                    info$n_species
                } else {
                    sample(info$species_min:info$species_max, 1)
                }

                new_individuals_per_species <- if (info$keep_same$keep_individuals) {
                    if (length(info$individuals_per_species) == new_n_species) {
                        info$individuals_per_species
                    } else {
                        sample(info$indiv_min:info$indiv_max, new_n_species, replace = TRUE)
                    }
                } else {
                    sample(info$indiv_min:info$indiv_max, new_n_species, replace = TRUE)
                }

                new_Sim_length <- if (info$keep_same$keep_sim_length) {
                    info$Sim_length
                } else {
                    sample(info$sim_min:info$sim_max, 1)
                }

                #Grow the reef if it is too small for this replicate's colonies
                #(prevents the "reef too crowded" placement failure)
                grown <- growReefToFit(
                    new_reef_x, new_reef_y, sum(new_individuals_per_species)
                )
                new_reef_x <- grown[1]
                new_reef_y <- grown[2]

                current_info <- list(
                    reef_x = new_reef_x,
                    reef_y = new_reef_y,
                    n_species = new_n_species,
                    individuals_per_species = new_individuals_per_species,
                    total_individuals = sum(new_individuals_per_species),
                    placement_type = info$placement_type,
                    Sim_length = new_Sim_length,
                    trait_mode = info$trait_mode,
                    setup_mode = 2,
                    reef_x_min = info$reef_x_min,
                    reef_x_max = info$reef_x_max,
                    reef_y_min = info$reef_y_min,
                    reef_y_max = info$reef_y_max,
                    species_min = info$species_min,
                    species_max = info$species_max,
                    indiv_min = info$indiv_min,
                    indiv_max = info$indiv_max,
                    sim_min = info$sim_min,
                    sim_max = info$sim_max,
                    n_replicates = info$n_replicates,
                    randomize_per_replicate = info$randomize_per_replicate,
                    keep_same = info$keep_same
                )

                if (info$keep_same$keep_traits) {
                    current_traits <- species_traits
                } else {
                    current_traits <- createSpeciesTraits(
                        current_info$n_species,
                        current_info$trait_mode,
                        current_info$setup_mode,
                        reproduction_enabled,
                        interaction_mode,
                        interaction_matrix = NULL,
                        interaction_size_mode = interaction_size_mode,
                        growth_size_enabled = growth_size_enabled
                    )
                }

            #Fully random — re-sample everything silently from DEFAULTS
            } else if (info$setup_mode == 3) {
                current_info <- getRandomSpeciesInfo()
                current_info$n_replicates <- info$n_replicates
                current_info$randomize_per_replicate <- info$randomize_per_replicate
                current_traits <- createSpeciesTraits(
                    current_info$n_species,
                    current_info$trait_mode,
                    current_info$setup_mode,
                    reproduction_enabled,
                    interaction_mode,
                    interaction_matrix = NULL,
                    interaction_size_mode = interaction_size_mode,
                    growth_size_enabled = growth_size_enabled
                )
            }
        }

        #Generates placement coordinates for each coral (different each replicate)
        coords <- getCoordinates(
            total_individuals = current_info$total_individuals,
            reef_x = current_info$reef_x,
            reef_y = current_info$reef_y,
            placement_type = current_info$placement_type
        )

        #Generates coral objects and places them on reef using previous data
        coral_data <- createCorals(
            reef_x = current_info$reef_x,
            reef_y = current_info$reef_y,
            n_species = current_info$n_species,
            individuals_per_species = current_info$individuals_per_species,
            coords_sample = coords,
            species_traits = current_traits)

        #Build this replicate's disturbance schedule from the chosen settings,
        #using this replicate's own length (new random timings). Disturbance sizes
        #are taken from the live reef each time, so no reef dimensions are needed.
        current_disturbance <- buildDisturbanceConfig(
            disturbance_settings$enabled,
            disturbance_settings$frequency,
            disturbance_settings$size,
            current_info$Sim_length
        )

        #Creates states, which is the variable containing all the simulations
        states <- runSimulation(
            reef = coral_data$reef,
            corals = coral_data$corals,
            n_steps = current_info$Sim_length,
            colony_species = coral_data$colony_species,
            species_traits = current_traits,
            disturbance_config = current_disturbance,
            habitat = habitat_settings)

        #Colonies can split mid-run, so rebuild the id -> species map from every
        #colony that ever existed; species-level plots need the fragment ids too
        full_colony_species <- buildColonySpeciesFromStates(
            states, coral_data$colony_species
        )

        #Store this replicate's results
        all_replicates[[rep]] <- list(
            reef_x = current_info$reef_x,
            reef_y = current_info$reef_y,
            species_traits = current_traits,
            sim_length = current_info$Sim_length,
            corals = coral_data$corals,
            reef = coral_data$reef,
            states = states,
            colony_species = full_colony_species,
            disturbance = current_disturbance,
            habitat = habitat_settings,
            replicate = rep
        )
    }
    
    #returns and stores all important data
    #If only 1 replicate, return in same format as before for backward compatibility
    if (info$n_replicates == 1) {
        return(all_replicates[[1]])
    } else {
        return(list(
            n_replicates = info$n_replicates,
            replicates = all_replicates,
            species_traits = species_traits,
            corals = all_replicates[[1]]$corals,
            sim_length = info$Sim_length,
            reef_x = info$reef_x,
            reef_y = info$reef_y,
            disturbance = disturbance_settings,
            habitat = habitat_settings,
            randomize_per_replicate = info$randomize_per_replicate
        ))
    }
}


## Shape placement: equally-spaced coordinates, chosen by individual count.

# 2 individuals: a central line.
linecoordinates <- function(reef_size){

  side <- reef_size * 0.4

  centre_x <- reef_size / 2
  centre_y <- reef_size / 2

  list(
    c(
      round(centre_x - side/2),
      round(centre_y)
    ),
    c(
      round(centre_x + side/2),
      round(centre_y)
    )
  )
}

#Function if 3 individuals given
#Based on older model
#Basic equilateral triangle in centre of reef
trianglecoordinates <- function(reef_size){

  #Creates variables used in triangle, side of triangle = 40* of reef
    side <-reef_size*0.4
    centre_x <- reef_size/2
    centre_y <- reef_size/2

  #creates value for height of equilateral triangle
    height <- sqrt(3)/2*side
    

  #list of coords for three points in triangle
    list(
    c(
      round(centre_x - side/2),
      round(centre_y - height/3)),
    c(
      round(centre_x + side/2),
      round(centre_y - height/3)),
    c(
      round(centre_x),
      round(centre_y + 2*height/3))
  )}


#Function if 4 individuals given
#Generate square in reef centre
squarecoordinates <- function(reef_size){
  
  #Generates central point and ensures square appropriate distance from sides
  side <-reef_size*0.4
  centre_x <- reef_size/2
  centre_y <- reef_size/2

  half_side <- side /2

  #four corners of square
  list(
    c(
      round(centre_x - half_side),
      round(centre_y - half_side)
    ),
    c(
      round(centre_x + half_side),
      round(centre_y - half_side)
    ),
    c(
      round(centre_x + half_side),
      round(centre_y + half_side)
    ),
    c(
      round(centre_x - half_side),
      round(centre_y + half_side)
    )
  )
}


#Function if >4 individuals given
#Generates a circle from centre of reef
#Appropriately distanced from edges
#Then places individuals equally spaced around circle
#Based on number of total individuals
polygoncoordinates <- function(reef_size, n_vertices){

  side <- reef_size * 0.6

  centre_x <- reef_size / 2
  centre_y <- reef_size / 2

  radius <- side / 2

  angles <- seq(
    pi/2,
    pi/2 + 2*pi,
    length.out = n_vertices + 1
  )[-(n_vertices + 1)]

  lapply(angles, function(theta){

    c(
      round(centre_x + radius * cos(theta)),
      round(centre_y + radius * sin(theta))
    )

  })
}

# ---- Pairwise interaction matrix (the sole interaction mechanism) -----------
# M[i, j] = P(species i overgrows species j) on contact (rows = attacker, cols =
# occupant); diagonal = intraspecific. Any network is expressible; 0 = never
# (neutrality). canTakeCell allows a takeover when M[i,j] > 0; resolveClaims rolls it
# once per contested cell (size-adjusted). The matrix rides on species_traits as an
# attribute (attachInteraction). Empirical matrices plug in via the "supplied" mode.

# ---- Tuning numbers ----

#Range that RANDOM matrix entries are drawn from (a probability per ordered pair).
#Kept away from 0 and 1 so random networks are genuinely stochastic contests.
INTERACTION_PROB_MIN <- 0.10
INTERACTION_PROB_MAX <- 0.90

#Default INTRASPECIFIC overgrowth probability - the value placed on the DIAGONAL of
#the matrix, used when two colonies of the SAME species meet. Kept on the one matrix
#(rather than a separate rule) so everything lives in one table. 0.5 = an even
#contest; set to 0 for "same-species colonies never overgrow each other".
INTRASPECIFIC_PROB_DEFAULT <- 0.5

# Size influence on the overgrowth probability. "none" = matrix as-is; "logratio" =
# shift in log-odds by beta_i * ln(size_i / size_j) (attacker sensitivity beta; only
# the ratio matters, log-scaled, hard 0/1 entries stay). Default mode:
INTERACTION_SIZE_MODE_DEFAULT <- "none"

#Default per-species size-sensitivity magnitude (beta) for a size-affected species.
#Deliberately gentle: at beta = 0.5 a 3x-larger attacker's win prob rises from 0.50
#to ~0.63 (noticeable, not drastic). 0 = size has no effect on that species.
INTERACTION_SIZE_BETA <- 0.5


# Ask whether colony size influences the matrix probabilities (off in random setup).
# Returns "none" or "logratio".
getInteractionSizeSetup <- function(setup_mode) {
  if (setup_mode == 3) {
    return("none")
  }
  on <- getMenuChoice(
    "Size effect on overgrowth (log-ratio of the two colonies' sizes)?",
    c("On", "Off")
  ) == 1
  if (on) "logratio" else "none"
}


# Validate/normalise a supplied matrix against the run's species: reorder named rows/
# cols to match, clamp to [0, 1], fill any NA diagonal with the intraspecific default.
validateInteractionMatrix <- function(M, species_names) {
  if (is.null(M)) {
    stop("Interaction matrix mode is 'supplied' but no matrix was given.")
  }
  M <- as.matrix(M)
  n <- length(species_names)
  if (nrow(M) != n || ncol(M) != n) {
    stop(sprintf("Interaction matrix must be %d x %d to match the %d species.",
                 n, n, n))
  }

  #Unnamed: assume rows/cols are already in species order. Named: reorder to match.
  if (is.null(rownames(M))) rownames(M) <- species_names
  if (is.null(colnames(M))) colnames(M) <- species_names
  if (!all(rownames(M) == species_names) && all(species_names %in% rownames(M))) {
    M <- M[species_names, , drop = FALSE]
  }
  if (!all(colnames(M) == species_names) && all(species_names %in% colnames(M))) {
    M <- M[, species_names, drop = FALSE]
  }

  storage.mode(M) <- "double"
  M[] <- pmin(1, pmax(0, M))   #keep every entry a valid probability
  #Diagonal = intraspecific probability. Fill any blank diagonal with the default,
  #but keep an explicitly supplied same-species value.
  d <- diag(M)
  d[is.na(d)] <- INTRASPECIFIC_PROB_DEFAULT
  diag(M) <- d
  M
}


# Build the N x N overgrowth matrix: "random" (draw off-diagonals), "manual" (prompt
# each pair), or "supplied" (validate a passed-in matrix). Diagonal = intraspecific.
buildInteractionMatrix <- function(species_names, mode = "random", supplied = NULL) {
  if (mode == "supplied") {
    return(validateInteractionMatrix(supplied, species_names))
  }

  n <- length(species_names)
  M <- matrix(NA_real_, n, n, dimnames = list(species_names, species_names))

  if (mode == "manual") {
    cat("\nPairwise overgrowth matrix: enter P(row overgrows column) for each pair.\n",
        "  0 = never overgrows (a wall, like neutrality), 100 = always overgrows.\n",
        sep = "")
  }

  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      if (i == j) next   #diagonal set to the intraspecific probability below
      if (mode == "manual") {
        M[i, j] <- getValidInteger(
          paste0("  P(", species_names[i], " overgrows ",
                 species_names[j], ") (%, 0-100): "),
          min_val = 0, max_val = 100
        ) / 100
      } else {
        M[i, j] <- round(runif(1, INTERACTION_PROB_MIN, INTERACTION_PROB_MAX), 2)
      }
    }
  }
  #Diagonal = intraspecific probability (same for every species by default)
  diag(M) <- INTRASPECIFIC_PROB_DEFAULT
  M
}


# Default per-species size betas: no effect (all 0).
defaultSizeBeta <- function(species_names) {
  setNames(rep(0, length(species_names)), species_names)
}


# Random per-species size betas (log-ratio mode): random sign x INTERACTION_SIZE_BETA.
randomSizeBeta <- function(species_names) {
  setNames(sample(c(-1, 1), length(species_names), replace = TRUE) * INTERACTION_SIZE_BETA,
           species_names)
}


# Read an empirical overgrowth matrix from a CSV (species names in first col/header;
# cells as probabilities or percentages, auto-detected). Validated against species_names.
readInteractionMatrixCSV <- function(path, species_names) {
  raw <- read.csv(path, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
  M   <- as.matrix(raw)
  storage.mode(M) <- "double"
  #Auto-detect percentages (any entry > 1) and rescale to probabilities
  if (any(M > 1, na.rm = TRUE)) {
    M <- M / 100
  }
  validateInteractionMatrix(M, species_names)
}


# Attach the interaction config on species_traits (an attribute) so it travels through
# the combat functions. list(enabled, matrix, size_mode, size_beta).
attachInteraction <- function(species_traits, enabled, matrix = NULL,
                              size_mode = INTERACTION_SIZE_MODE_DEFAULT,
                              size_beta = NULL) {
  #Default the per-species size sensitivities to no effect
  if (is.null(size_beta) && !is.null(matrix)) {
    size_beta <- defaultSizeBeta(rownames(matrix))
  }
  attr(species_traits, "interaction") <- list(
    enabled   = isTRUE(enabled),
    matrix    = matrix,
    size_mode = size_mode,   #"none" or "logratio"
    size_beta = size_beta    #named per-species size sensitivity (0 = unaffected)
  )
  species_traits
}

getInteraction <- function(species_traits) {
  it <- attr(species_traits, "interaction", exact = TRUE)
  if (is.null(it)) {
    return(list(enabled = FALSE, matrix = NULL, size_mode = "none", size_beta = NULL))
  }
  #Backfill defaults for configs made before these fields existed
  if (is.null(it$size_mode)) it$size_mode <- "none"
  if (is.null(it$size_beta) && !is.null(it$matrix)) {
    it$size_beta <- defaultSizeBeta(rownames(it$matrix))
  }
  #Cache species-name -> row/col index maps once, so overgrowthProb can index the
  #matrix by integer position (M[i, j]) instead of doing two `%in%` name checks plus a
  #name-based lookup on every contested cell. An unknown name indexes to NA -> prob 0,
  #exactly as the old `%in%` guard returned. Purely a speed cache; values are unchanged.
  if (is.null(it$row_idx) && !is.null(it$matrix)) {
    it$row_idx <- setNames(seq_len(nrow(it$matrix)), rownames(it$matrix))
    it$col_idx <- setNames(seq_len(ncol(it$matrix)), colnames(it$matrix))
  }
  it
}

# Attacker species' size sensitivity (beta) from the interaction config; 0 if unset.
sizeBetaFor <- function(species, interaction) {
  sb <- interaction$size_beta
  if (is.null(sb)) return(0)
  b <- unname(sb[species])          #NA if this species has no entry
  if (is.na(b)) 0 else as.numeric(b)
}


# P(attacker overgrows occupant) on contact; 0 if off/missing/unknown/0-entry. A same-
# species pair reads the diagonal. (Fast path via cached integer indices.)
overgrowthProb <- function(attacker_species, occupant_species, interaction) {
  if (is.null(interaction) || !isTRUE(interaction$enabled) ||
      is.null(interaction$matrix)) {
    return(0)
  }
  M <- interaction$matrix
  #Fast path: index by cached integer position. row_idx/col_idx are named lookups, so
  #an unknown species name yields NA (-> prob 0), matching the old `%in%` guard.
  if (!is.null(interaction$row_idx) && !is.null(interaction$col_idx)) {
    i <- interaction$row_idx[attacker_species]
    j <- interaction$col_idx[occupant_species]
    if (is.na(i) || is.na(j)) {
      return(0)
    }
    p <- M[i, j]
    return(if (is.na(p)) 0 else p)
  }
  #Fallback (interaction built without the index cache): original name-based lookup.
  if (!(attacker_species %in% rownames(M)) ||
      !(occupant_species %in% colnames(M))) {
    return(0)
  }
  #Diagonal (attacker == occupant) is the intraspecific probability; off-diagonal is
  #the interspecific probability. Both are just a matrix lookup now.
  p <- M[attacker_species, occupant_species]
  if (is.na(p)) 0 else p
}


# Overgrowth probability for a specific attacker vs occupant colony, size-adjusted
# (the value rolled in resolveClaims). "logratio" shifts the base by beta * ln(ratio).
overgrowthProbSized <- function(attacker, occupant, interaction) {
  base <- overgrowthProb(attacker$species, occupant$species, interaction)
  #Hard 0/1 entries stay hard (neutrality / certainty are absolute; logit undefined)
  if (base <= 0) return(0)
  if (base >= 1) return(1)

  mode <- if (is.null(interaction$size_mode)) "none" else interaction$size_mode
  if (mode != "logratio") return(base)

  beta <- sizeBetaFor(attacker$species, interaction)
  if (beta == 0) return(base)

  z_i <- attacker$size; z_j <- occupant$size
  if (is.null(z_i) || is.na(z_i) || z_i <= 0 ||
      is.null(z_j) || is.na(z_j) || z_j <= 0) {
    return(base)
  }
  x <- log(base / (1 - base)) + beta * log(z_i / z_j)
  1 / (1 + exp(-x))
}


# One-line description of the interaction config, for the setup summary.
describeInteraction <- function(species_traits) {
  it <- getInteraction(species_traits)
  if (!isTRUE(it$enabled) || is.null(it$matrix)) {
    return("No interaction matrix set")
  }
  n <- if (is.null(it$matrix)) 0L else nrow(it$matrix)
  size_desc <- switch(it$size_mode,
    logratio = "log-ratio size effect",
    "no size effect")
  paste0("On (", n, "x", n, " pairwise overgrowth matrix, ", size_desc, ")")
}


# Print the interaction matrix as a readable table (rows overgrow columns; diagonal =
# intraspecific; 0 = never). Returns the matrix invisibly.
printInteractionMatrix <- function(species_traits, digits = 2) {
  it <- getInteraction(species_traits)
  if (!isTRUE(it$enabled) || is.null(it$matrix)) {
    cat("No interaction matrix is set for these species.\n")
    return(invisible(NULL))
  }

  M    <- it$matrix
  sp   <- rownames(M)
  n    <- length(sp)
  wsp  <- max(nchar(sp), 5)                 #width for the row/column labels
  wnum <- max(digits + 3, 6)                #width for each probability cell

  fmtnum <- function(x) {
    if (is.na(x)) formatC("-", width = wnum) else
      formatC(x, format = "f", digits = digits, width = wnum)
  }

  cat("Species interaction matrix   -   P(row overgrows column)\n")
  cat("  Rows = attacker (taking a cell); Columns = occupant (being grown over).\n")
  cat("  Diagonal = intraspecific (two colonies of the SAME species).  0 = never (neutral).\n")
  cat("  Size effect:", it$size_mode, "\n\n")

  #Header row of column-species labels, aligned over the number cells
  cat(strrep(" ", wsp + 3),
      paste(vapply(sp, function(s) formatC(s, width = wnum), character(1)),
            collapse = " "), "\n", sep = "")

  sb <- it$size_beta
  for (i in seq_len(n)) {
    cells <- paste(vapply(M[i, ], fmtnum, character(1)), collapse = " ")
    cat("  ", formatC(sp[i], width = wsp, flag = "-"), " ", cells, "\n", sep = "")

    #One small line of size sensitivity beneath the row (log-ratio mode only)
    if (identical(it$size_mode, "logratio") && !is.null(sb)) {
      b <- sb[sp[i]]
      if (!is.na(b) && b != 0) {
        pad <- strrep(" ", wsp + 5)
        cat(pad, "size sensitivity beta = ", sprintf("%+.2f", b),
            " (+ = bigger wins more)\n", sep = "")
      }
    }
  }
  invisible(M)
}


# plotmath label: species name with size beta as a superscript (omitted if 0).
speciesSizeLabel <- function(sp_name, size_beta) {
  b <- 0
  if (!is.null(size_beta) && sp_name %in% names(size_beta)) b <- size_beta[[sp_name]]
  sup <- if (b != 0) sprintf('^{"b=%+g"}', b) else ""
  parse(text = paste0('"', sp_name, '"', sup))[[1]]
}


# plotmath label: growth number with growth-size gamma as a superscript (omitted if 0/off).
growthSizeLabel <- function(growth_value, sp_name, growth_size) {
  g <- 0
  if (!is.null(growth_size) && isTRUE(growth_size$enabled) &&
      !is.null(growth_size$gamma) && sp_name %in% names(growth_size$gamma)) {
    g <- growth_size$gamma[[sp_name]]
  }
  sup <- if (g != 0) sprintf('^{"g=%+g"}', g) else ""
  parse(text = paste0('"', as.character(growth_value), '"', sup))[[1]]
}


# Plot the interaction matrix as a table figure: species-name boxes (filled with the
# species colour, row headers also showing size beta + growth rate) and white P(row
# overgrows column) cells. Pass `corals` for live colony colours.
plotInteractionMatrix <- function(species_traits, corals = NULL,
                                  title = "Species interaction matrix") {
  it <- getInteraction(species_traits)
  if (!isTRUE(it$enabled) || is.null(it$matrix)) {
    plot.new(); title("No interaction matrix set for these species")
    return(invisible(NULL))
  }

  M  <- it$matrix
  sp <- rownames(M)
  n  <- length(sp)
  si <- it$size_beta
  #Growth-size config (may be absent) for the growth number's superscript. Read
  #straight off the attribute so the plot does not depend on Sim_func being sourced.
  gsz <- attr(species_traits, "growth_size", exact = TRUE)

  #One colour per species: from the coral list if given (as plotSpeciesTraits does),
  #otherwise the base species palette (hcl "Dark 3") at the same alpha = 0.4 that
  #plotSpeciesTraits' first-individual colours use, so the two figures match.
  species_colours <- if (!is.null(corals)) {
    vapply(sp, function(s) {
      m <- which(vapply(corals, function(x) x$species, character(1)) == s)
      if (length(m) == 0) return("grey80")
      col <- corals[[m[1]]]$colour
      if (is.null(col) || is.na(col)) "grey80" else col
    }, character(1))
  } else {
    adjustcolor(hcl.colors(n, "Dark 3"), alpha.f = 0.4)
  }

  #Growth rate per species (blank if there is no growth column)
  growth_of <- function(s) {
    if (is.null(species_traits$growth)) return(NA)
    g <- species_traits$growth[species_traits$species == s]
    if (length(g) == 0) NA else g[1]
  }

  op <- par(mar = c(0.5, 0.5, 7, 0.5)); on.exit(par(op))
  plot(NA, xlim = c(0, n + 1), ylim = c(0, n + 1),
       axes = FALSE, xlab = "", ylab = "")

  #A species-name box: species colour fill and the species name (larger). When
  #full = TRUE (the LEFT row headers only) it also shows the size-impact
  #super/subscript on the name and the growth rate underneath (smaller). The top
  #column headers use full = FALSE, so they carry only the plain species name.
  drawSpeciesBox <- function(x_left, y_bottom, k, full = TRUE) {
    rect(x_left, y_bottom, x_left + 1, y_bottom + 1,
         col = species_colours[k], border = "black")
    if (full) {
      text(x_left + 0.5, y_bottom + 0.58, labels = speciesSizeLabel(sp[k], si),
           font = 2, cex = 1.5)
      g <- growth_of(sp[k])
      if (!is.na(g)) {
        text(x_left + 0.5, y_bottom + 0.20,
             labels = growthSizeLabel(g, sp[k], gsz), font = 2, cex = 0.85)
      }
    } else {
      text(x_left + 0.5, y_bottom + 0.5, labels = sp[k], font = 2, cex = 1.5)
    }
  }

  #Top-left corner: blank white cell
  rect(0, n, 1, n + 1, col = "white", border = "black")

  #Column headers (occupant species) across the top - name only
  for (j in seq_len(n)) drawSpeciesBox(j, n, j, full = FALSE)

  #Row headers (attacker species) down the left - full detail (size-impact
  #super/subscript on the name, growth rate beneath it), then the interaction cells
  for (i in seq_len(n)) {
    yb <- n - i
    drawSpeciesBox(0, yb, i, full = TRUE)
    for (j in seq_len(n)) {
      rect(j, yb, j + 1, yb + 1, col = "white", border = "black")
      p   <- M[i, j]
      txt <- if (is.na(p)) "-" else formatC(p, format = "f", digits = 2)
      text(j + 0.5, yb + 0.5, labels = txt, cex = 1.05)
    }
  }

  #Title, then the key (same wording as plotSpeciesTraits), the orientation line,
  #and the growth-rate note - each on its own line above the grid.
  mtext(title, side = 3, line = 5.4, font = 2, cex = 1.2)
  mtext(paste0("Superscript = effect when large   Subscript = effect when small",
               "   ( +/- = increase / decrease )"),
        side = 3, line = 4.0, cex = 0.7)
  mtext("rows overgrow columns   ·   diagonal = intraspecific   ·   0 = neutral",
        side = 3, line = 3.0, cex = 0.7)
  mtext(paste0("Number below species = growth rate",
               "  (its super/subscript = growth-trait shift when large / small)"),
        side = 3, line = 2.0, cex = 0.7)
  invisible(M)
}
