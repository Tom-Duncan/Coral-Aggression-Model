#File containing all functions related to visualisations of the model


#Function to create plot of reef matrix with visible cells for each species
#Converts species letters into numeric values then draws them as different colours on grid
plotReef <- function(reef, corals, timestep_num = NULL) {
  if (is.null(reef)) {
    stop("reef is NULL.")
  }
  if (!is.matrix(reef)) {
    stop(paste("reef must be a matrix, but is:", class(reef)[1]))
  }
  #Reef may be rectangular: columns are the x-axis (width), rows the y-axis (height)
  reef_x <- ncol(reef)
  reef_y <- nrow(reef)
  if (length(reef_x) != 1 || is.na(reef_x) || length(reef_y) != 1 || is.na(reef_y)) {
    stop("reef dimensions are not valid numbers.")
  }
  colony_ids <- sapply(corals, function(x) x$id)
  colony_cols <- sapply(corals, function(x) x$colour)
  names(colony_cols) <- colony_ids

  tick_x <- seq(0, reef_x, length.out = 6)
  tick_y <- seq(0, reef_y, length.out = 6)

  old_par <- par(
    xaxs = "i",
    yaxs = "i",
    bty = "n"
  )
  on.exit(par(old_par))

  plot(
    x = c(0, reef_x),
    y = c(0, reef_y),
    type = "n",
    xlim = c(0, reef_x),
    ylim = c(0, reef_y),
    asp = 1,
    axes = FALSE,
    ann = FALSE
  )

  # Reef background
  rect(
    0, 0, reef_x, reef_y,
    col = "white",
    border = NA
  )

  # Coral cells
  for (id in colony_ids) {
    cells <- which(reef == id, arr.ind = TRUE)

    if (nrow(cells) > 0) {
      for (k in seq_len(nrow(cells))) {
        row <- cells[k, 1]
        col <- cells[k, 2]

        rect(
          xleft = col - 1,
          xright = col,
          ybottom = row - 1,
          ytop = row,
          col = colony_cols[id],
          border = NA
        )
      }
    }
  }

  # Grid lines INCLUDING the full reef boundary (vertical every x, horizontal every y)
  for (i in 0:reef_x) {
    segments(i, 0, i, reef_y, col = "grey90", lwd = 1)
  }
  for (i in 0:reef_y) {
    segments(0, i, reef_x, i, col = "grey90", lwd = 1)
  }

  # Draw axes exactly on the reef boundary
  axis(1, at = tick_x, labels = round(tick_x), pos = 0)
  axis(2, at = tick_y, labels = round(tick_y), pos = 0, las = 1)

  # Overdraw left/bottom axes in black so they sit on top of pale grid
  segments(0, 0, reef_x, 0, col = "black", lwd = 1)
  segments(0, 0, 0, reef_y, col = "black", lwd = 1)

  title(
    main = "Coral Growth and Interaction Simulation",
    xlab = "Reef X coordinate",
    ylab = "Reef Y coordinate"
  )

  if (!is.null(timestep_num)) {
    mtext(
      paste("Timestep", timestep_num),
      side = 3,
      line = 0.5,
      cex = 0.9
    )
  }
}


#Convenience wrapper: plot the reef at a single given timestep.
#Enter ONE number (the timestep) and it finds the matching state and labels it.
#states[[1]] is the initial state (timestep 0), so timestep t lives in states[[t + 1]].
plotReefAtTimestep <- function(states, timestep_num) {

  state_index <- timestep_num + 1

  if (state_index < 1 || state_index > length(states)) {
    stop(paste0(
      "Timestep ", timestep_num, " is out of range. ",
      "Valid timesteps are 0 to ", length(states) - 1, "."
    ))
  }

  state <- states[[state_index]]

  plotReef(
    state$reef,
    state$corals,
    timestep_num = timestep_num
  )
}



#Function to create plot of population change over timesteps
#Tracks each colony individually. Because colonies can split into new colonies
#mid-run, the set of colonies is gathered from EVERY timestep (not a single fixed
#list), so fragments born partway through are still counted from their first step.
#The `corals` argument is accepted for backward compatibility but no longer needed.
plotIndiAbundance <- function(states, corals = NULL) {

  # Every colony that ever existed, each with its (first-seen) colour
  colony_cols <- character(0)
  for (s in states) {
    for (co in s$corals) {
      if (!(co$id %in% names(colony_cols))) {
        colony_cols[co$id] <- co$colour
      }
    }
  }

  colony_ids  <- names(colony_cols)
  colony_cols <- unname(colony_cols)

  n_steps <- length(states)

  # Create abundance table, one column per colony
  abundance <- data.frame(timestep = 0:(n_steps - 1))
  for (id in colony_ids) {
    abundance[[id]] <- numeric(n_steps)
  }

  # Each colony's reef cover (% of reef) at every timestep (0 before it is born)
  for (i in seq_along(states)) {
    reef <- states[[i]]$reef
    total_cells <- nrow(reef) * ncol(reef)
    for (id in colony_ids) {
      abundance[[id]][i] <- sum(reef == id, na.rm = TRUE) / total_cells * 100
    }
  }

  # Determine y-axis limit (guard against an all-zero / empty case)
  ymax <- max(as.matrix(abundance[, colony_ids, drop = FALSE]))
  if (!is.finite(ymax) || ymax == 0) ymax <- 1

  # Plot first colony
  plot(
    abundance$timestep,
    abundance[[colony_ids[1]]],
    type = "l",
    col = colony_cols[1],
    lwd = 2,
    ylim = c(0, ymax),
    xlab = "Timestep",
    ylab = "Reef cover (%)",
    main = "Coral abundance through time (individual colonies)"
  )

  # Add remaining colonies
  if (length(colony_ids) > 1) {
    for (j in 2:length(colony_ids)) {
      lines(
        abundance$timestep,
        abundance[[colony_ids[j]]],
        col = colony_cols[j],
        lwd = 2
      )
    }
  }

  # A full legend is only legible for a handful of colonies; with many
  # fragments it would swamp the plot, so summarise the count instead
  if (length(colony_ids) <= 12) {
    legend("topright", legend = colony_ids, col = colony_cols, lwd = 2, cex = 0.8)
  } else {
    legend(
      "topright",
      legend = paste0(length(colony_ids), " colonies (incl. split-off fragments)"),
      bty = "n", cex = 0.8
    )
  }
}

#Function to plot population by species instead of individuals
plotSpeciesAbundance <- function(states, corals, colony_species) {

  #Lookup table
  coral_df <- data.frame(
    id = sapply(corals, function(x) x$id),
    colour = sapply(corals, function(x) x$colour),
    species <- sapply(corals, function(x) x$species)
  )

  #Species-level lookup 
  species_lookup <- aggregate(
    colour ~ species,
    data = coral_df,
    FUN = function(x) x[1]
  )

  species_ids <- species_lookup$species
  species_cols <- species_lookup$colour
  names(species_cols) <- species_ids  

  n_steps <- length(states)

  abundance <- data.frame(
    timestep = 0:(n_steps - 1)
  )

  for (sp in species_ids) {
    abundance[[sp]] <- numeric(n_steps)
  }


  for (i in seq_along(states)) {

    # Extract reef state 
    reef <- states[[i]]$reef

    # Convert to species labels 
    species_matrix <- matrix(
    colony_species[reef],
    nrow = nrow(reef),
    ncol = ncol(reef)
)

    # Each species' reef cover (% of reef) at this timestep
    total_cells <- nrow(reef) * ncol(reef)
    for (sp in species_ids) {
      abundance[[sp]][i] <- sum(species_matrix == sp, na.rm = TRUE) / total_cells * 100
    }
  }


  #Plotting


  ymax <- max(abundance[, species_ids, drop = FALSE])

  # First species line
  plot(
    abundance$timestep,
    abundance[[species_ids[1]]],
    type = "l",
    col = species_cols[species_ids[1]],
    lwd = 2,
    ylim = c(0, ymax),
    xlab = "Timestep",
    ylab = "Reef cover (%)",
    main = "Coral abundance through time (by species)"
  )

  # Remaining species lines
  if (length(species_ids) > 1) {
    for (j in 2:length(species_ids)) {

      sp <- species_ids[j]

      lines(
        abundance$timestep,
        abundance[[sp]],
        col = species_cols[sp],
        lwd = 2
      )
    }
  }

  # Legend
  legend(
    "topright",
    legend = species_ids,
    col = species_cols[species_ids],
    lwd = 2
  )
}


#Function to make gif of abundance change graph
make_abundance_gif <- function(states,
                               corals,
                               file = "coral_abundance.gif",
                               interval = 0.1) {

  library(animation)

  # Ensure previous graphics devices are closed safely
  if (dev.cur() > 1) dev.off()

  on.exit({
    while (dev.cur() > 1) dev.off()
  }, add = TRUE)

  # Every colony that ever existed, each with its (first-seen) colour, since
  # colonies can split into new colonies partway through the run
  colony_cols <- character(0)
  for (s in states) {
    for (co in s$corals) {
      if (!(co$id %in% names(colony_cols))) {
        colony_cols[co$id] <- co$colour
      }
    }
  }
  species_ids  <- names(colony_cols)
  species_cols <- unname(colony_cols)

  n_steps <- length(states)

  # Build abundance table once
  abundance <- data.frame(
    timestep = 0:(n_steps - 1)
  )

  for (species in species_ids) {
    abundance[[species]] <- numeric(n_steps)
  }

  for (i in seq_along(states)) {
    reef <- states[[i]]$reef
    total_cells <- nrow(reef) * ncol(reef)

    for (species in species_ids) {
      abundance[[species]][i] <- sum(reef == species, na.rm = TRUE) / total_cells * 100
    }
  }

  ymax <- max(as.matrix(abundance[, species_ids, drop = FALSE]))
  if (!is.finite(ymax) || ymax == 0) ymax <- 1

  saveGIF({

    for (i in seq_len(n_steps)) {

      current_steps <- 1:i

      plot(
        abundance$timestep[current_steps],
        abundance[[species_ids[1]]][current_steps],
        type = "l",
        col = species_cols[1],
        lwd = 2,
        ylim = c(0, ymax),
        xlim = c(0, n_steps - 1),
        xlab = "Timestep",
        ylab = "Reef cover (%)",
        main = "Coral abundance through time"
      )

      if (length(species_ids) > 1) {
        for (j in 2:length(species_ids)) {
          lines(
            abundance$timestep[current_steps],
            abundance[[species_ids[j]]][current_steps],
            col = species_cols[j],
            lwd = 2
          )
        }
      }

      if (length(species_ids) <= 12) {
        legend("topright", legend = species_ids, col = species_cols, lwd = 2, cex = 0.8)
      } else {
        legend(
          "topright",
          legend = paste0(length(species_ids), " colonies (incl. split-off fragments)"),
          bty = "n", cex = 0.8
        )
      }

      title(sub = paste("Timestep", i - 1))
    }

  },
  movie.name = file,
  interval = interval,
  ani.dev = "png",
  clean = TRUE
  )
}



#--------COLONY SIZE TRACKING----------

#Function to extract colony size time series across all timesteps
getColonySizeTimeSeries <- function(states, corals) {
  
  # Create data frame to store size data
  size_data <- data.frame(
    timestep = integer(0),
    colony_id = character(0),
    species = character(0),
    size = integer(0),
    stringsAsFactors = FALSE
  )
  
  # Loop through each timestep
  for (t in seq_along(states)) {
    timestep_corals <- states[[t]]$corals
    
    # For each coral, record its size
    for (coral in timestep_corals) {
      size_data <- rbind(
        size_data,
        data.frame(
          timestep = t - 1,  # Adjust to 0-based indexing
          colony_id = coral$id,
          species = coral$species,
          size = coral$size,
          stringsAsFactors = FALSE
        )
      )
    }
  }
  
  return(size_data)
}

#Function to create a summary table of colony sizes
summarizeColonySizes <- function(size_data) {
  
  # Get unique colonies
  unique_colonies <- unique(size_data$colony_id)
  
  # Initialize summary table
  summary_table <- data.frame(
    colony_id = character(0),
    species = character(0),
    initial_size = integer(0),
    final_size = integer(0),
    max_size = integer(0),
    min_size = integer(0),
    mean_size = numeric(0),
    stringsAsFactors = FALSE
  )
  
  # Loop through each colony and calculate stats
  for (col_id in unique_colonies) {
    colony_data <- size_data[size_data$colony_id == col_id, ]
    
    summary_table <- rbind(
      summary_table,
      data.frame(
        colony_id = col_id,
        species = colony_data$species[1],
        initial_size = colony_data$size[1],
        final_size = colony_data$size[nrow(colony_data)],
        max_size = max(colony_data$size),
        min_size = min(colony_data$size),
        mean_size = round(mean(colony_data$size), 2),
        stringsAsFactors = FALSE
      )
    )
  }
  
  return(summary_table)
}


make_coral_gif <- function(states,
                              file = "coral_growth.gif",
            interval = 0.1) {

    library(animation)

    # Ensure previous graphics devices are closed safely
    if (dev.cur() > 1) dev.off()

    on.exit({
      while (dev.cur() > 1) dev.off()
    }, add = TRUE)

    saveGIF({

      for (i in seq_along(states)) {

        plotReef(
          states[[i]]$reef,
          states[[i]]$corals,
          timestep_num = i - 1
        )
      }

    },
    movie.name = file,
    interval = interval,
    ani.dev = "png",
    clean = TRUE
    )
}

