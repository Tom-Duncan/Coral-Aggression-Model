# Visualisation functions for the model.


# Draw the reef grid, colouring each cell by its colony's colour.
plotReef <- function(reef, corals, timestep_num = NULL) {
  if (is.null(reef)) {
    stop("reef is NULL.")
  }
  if (!is.matrix(reef)) {
    stop(paste("reef must be a matrix, but is:", class(reef)[1]))
  }
  reef_x <- ncol(reef)   # width
  reef_y <- nrow(reef)   # height
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

  rect(
    0, 0, reef_x, reef_y,
    col = "white",
    border = NA
  )

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

  # grid lines including the reef boundary
  for (i in 0:reef_x) {
    segments(i, 0, i, reef_y, col = "grey90", lwd = 1)
  }
  for (i in 0:reef_y) {
    segments(0, i, reef_x, i, col = "grey90", lwd = 1)
  }

  axis(1, at = tick_x, labels = round(tick_x), pos = 0)
  axis(2, at = tick_y, labels = round(tick_y), pos = 0, las = 1)

  # overdraw left/bottom axes in black, above the pale grid
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


# Plot the reef at a single timestep (state t is states[[t + 1]]).
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



# Per-colony cover through time. Colonies are gathered from every timestep (not a
# fixed list) so fragments born mid-run are counted. `corals` is kept for compatibility.
plotIndiAbundance <- function(states, corals = NULL) {

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

  abundance <- data.frame(timestep = 0:(n_steps - 1))
  for (id in colony_ids) {
    abundance[[id]] <- numeric(n_steps)
  }

  # each colony's cover (% reef) at every step (0 before it is born)
  for (i in seq_along(states)) {
    reef <- states[[i]]$reef
    total_cells <- nrow(reef) * ncol(reef)
    for (id in colony_ids) {
      abundance[[id]][i] <- sum(reef == id, na.rm = TRUE) / total_cells * 100
    }
  }

  ymax <- max(as.matrix(abundance[, colony_ids, drop = FALSE]))
  if (!is.finite(ymax) || ymax == 0) ymax <- 1

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

  # full legend only for a few colonies; otherwise just the count
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

# Per-species cover through time.
plotSpeciesAbundance <- function(states, corals, colony_species) {

  coral_df <- data.frame(
    id = sapply(corals, function(x) x$id),
    colour = sapply(corals, function(x) x$colour),
    species <- sapply(corals, function(x) x$species)
  )

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

    reef <- states[[i]]$reef

    species_matrix <- matrix(
    colony_species[reef],
    nrow = nrow(reef),
    ncol = ncol(reef)
)

    total_cells <- nrow(reef) * ncol(reef)
    for (sp in species_ids) {
      abundance[[sp]][i] <- sum(species_matrix == sp, na.rm = TRUE) / total_cells * 100
    }
  }


  ymax <- max(abundance[, species_ids, drop = FALSE])

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

  legend(
    "topright",
    legend = species_ids,
    col = species_cols[species_ids],
    lwd = 2
  )
}


# GIF of the abundance-through-time plot.
make_abundance_gif <- function(states,
                               corals,
                               file = "coral_abundance.gif",
                               interval = 0.1) {

  library(animation)

  if (dev.cur() > 1) dev.off()

  on.exit({
    while (dev.cur() > 1) dev.off()
  }, add = TRUE)

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



# -------- Colony size tracking --------

# Colony sizes at every timestep as a long data.frame.
getColonySizeTimeSeries <- function(states, corals) {

  size_data <- data.frame(
    timestep = integer(0),
    colony_id = character(0),
    species = character(0),
    size = integer(0),
    stringsAsFactors = FALSE
  )

  for (t in seq_along(states)) {
    timestep_corals <- states[[t]]$corals

    for (coral in timestep_corals) {
      size_data <- rbind(
        size_data,
        data.frame(
          timestep = t - 1,   # 0-based
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

# Per-colony size summary (initial/final/max/min/mean).
summarizeColonySizes <- function(size_data) {

  unique_colonies <- unique(size_data$colony_id)

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


# GIF of the reef growing over time.
make_coral_gif <- function(states,
                              file = "coral_growth.gif",
            interval = 0.1) {

    library(animation)

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
