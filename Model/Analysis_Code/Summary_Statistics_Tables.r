#All the summary-statistics table builders in one file (formerly
#Summary_tables.r + Calculated_Stats.r). Shared helpers (ID codebook, table
#writers, metric list) come from Analysis_Utils.r. Definitions first; every
#table-generating call is in the run section at the bottom, grouped by
#simulation, so source the top half and run only the blocks you need.
#The retention/loss analyses need run_index.rds from Consolidate_Datasets.r.

# Resolve the repo root (folder containing "Model/") so this works from any directory.
setwd((function() {
  a <- commandArgs(FALSE); f <- sub("^--file=", "", grep("^--file=", a, value = TRUE))
  d <- normalizePath(if (length(f)) dirname(f[1]) else getwd(), mustWork = FALSE)
  while (basename(d) != "Model" && !dir.exists(file.path(d, "Model")) && dirname(d) != d) d <- dirname(d)
  if (basename(d) == "Model") dirname(d) else d
})())

source("Model/Analysis_Code/Analysis_Utils.r")


## Scenario design tables ------------------------------------------------------
#one row per scenario describing its design/settings, written as CSV + Word

# --- Build the one-row-per-scenario design table ----------------------------
# Columns: ID, Scenario, Matrix, N0 Species, N0 Colonies, Disturbance,
#          Reproduction, Reef Size, Length, Replicates, Bias
#   ID       : short permanent number from the codebook; differs per simulation
#              run (unique, never reused)
#   Scenario : the scenario column from the master file (same config -> same
#              string across studies; covers all replicates)
buildScenarioTable <- function(results) {
  scenes <- unique(results$scenario)

  rows <- lapply(scenes, function(s) {
    d <- results[results$scenario == s, ]

    data.frame(
      run_id         = sub("_r[0-9]+$", "", d$run_id[1]),        # unique per run-scenario (keys the codebook)
      Scenario       = d$scenario[1],                            # scenario column, as-is
      Matrix         = sub("^Classic_", "", d$combination[1]),   # RPS / Linear / Neutral / Random
      `N0 Species`   = d$n_species[1],
      `N0 Colonies`  = d$individuals[1],          # STARTING colonies per species
      Disturbance    = ifelse(d$disturbance_on[1],  "On", "Off"),
      Reproduction   = ifelse(d$reproduction_on[1], "On", "Off"),
      `Reef Size`    = paste0(d$reef_x[1], "x", d$reef_y[1]),
      Length         = d$sim_length[1],
      Replicates     = length(unique(d$replicate)),
      Bias           = d$bias[1],
      check.names = FALSE, stringsAsFactors = FALSE
    )
  })
  tab <- do.call(rbind, rows)
  tab <- tab[order(tab$Matrix, tab$`N0 Species`, tab$Disturbance), ]

  # look up / assign the short numeric ID, then drop the long run_id column
  tab$ID     <- assignScenarioIDs(tab$run_id, tab$Scenario)
  tab$run_id <- NULL
  tab <- tab[, c("ID", setdiff(names(tab), "ID"))]       # ID as the first column
  rownames(tab) <- NULL
  tab
}

# writeWordTable() and writeTablePair() are defined in Analysis_Utils.r.

# --- One call: read a master CSV, write CSV + Word table --------------------
#   csv_path : a master results file
#   out_dir  : where to save (defaults to summary_out_dir)
makeSummaryTable <- function(csv_path, out_dir = summary_out_dir) {
  tab <- buildScenarioTable(read.csv(csv_path, stringsAsFactors = FALSE))
  cat("Writing", nrow(tab), "scenarios:\n")
  writeTablePair(tab, tools::file_path_sans_ext(basename(csv_path)), out_dir)
  invisible(tab)
}

# --- Combined: several master files -> ONE stacked scenario table -----------
#   csv_paths : vector of master results files
#   out_name  : base name for the combined output files
# Builds each file's scenario table (IDs come from the shared codebook, so they
# stay unique across files), stacks them into one table, and writes one CSV +
# Word table. Empty/unreadable files are skipped with a note.
makeCombinedSummary <- function(csv_paths, out_name, out_dir = summary_out_dir) {
  tabs <- lapply(csv_paths, function(p) {
    tryCatch({
      t <- buildScenarioTable(read.csv(p, stringsAsFactors = FALSE))
      cat("  ", nrow(t), "scenarios from", basename(p), "\n")
      t
    }, error = function(e) {
      cat("  SKIPPED (empty/unreadable):", basename(p), "-", conditionMessage(e), "\n")
      NULL
    })
  })
  tab <- do.call(rbind, tabs[!vapply(tabs, is.null, logical(1))])
  if (is.null(tab) || nrow(tab) == 0) { cat("No readable scenarios - nothing written.\n"); return(invisible(NULL)) }

  tab <- tab[order(tab$ID), ]
  rownames(tab) <- NULL
  cat("Combined", nrow(tab), "scenarios from", length(csv_paths), "files:\n")
  writeTablePair(tab, out_name, out_dir)
  invisible(tab)
}


## Statistics tables -----------------------------------------------------------
########Statistical summary tables
# One row per scenario (model run). For each community metric, the mean +/- SD
# accumulated over the LAST n_check checkpoints (timesteps) pooled across all of
# that scenario's replicates. Per-species abundance columns are excluded.
# Uses the same ID codebook as the design table, so the two tables join on ID.

# stat_metrics / stat_labels (the community-level metrics) are in Analysis_Utils.r.

buildStatsTable <- function(results, n_check = 10, metrics = stat_metrics, digits = 3) {
  late <- tail(sort(unique(results$timestep)), n_check)     # last n_check checkpoints
  d    <- results[results$timestep %in% late, ]

  rows <- lapply(unique(d$scenario), function(s) {
    v   <- d[d$scenario == s, ]
    row <- data.frame(run_id   = sub("_r[0-9]+$", "", v$run_id[1]),
                      Scenario = v$scenario[1],
                      check.names = FALSE, stringsAsFactors = FALSE)
    for (m in metrics) {                                    # mean +/- SD over reps x checkpoints
      mu  <- mean(v[[m]], na.rm = TRUE)
      sdv <- sd(v[[m]],   na.rm = TRUE)
      row[[stat_labels[[m]]]] <-
        if (!is.finite(mu))       "n/a"                     # undefined (e.g. evenness at richness 1)
        else if (!is.finite(sdv)) sprintf("%.*f", digits, mu)
        else                      sprintf("%.*f ± %.*f", digits, mu, digits, sdv)
    }
    row
  })
  tab <- do.call(rbind, rows)

  tab$ID     <- assignScenarioIDs(tab$run_id, tab$Scenario)   # same codebook as design table
  tab$run_id <- NULL
  tab <- tab[order(tab$ID), c("ID", "Scenario", unname(stat_labels[metrics]))]
  rownames(tab) <- NULL
  tab
}

# Read a master CSV, build the stats table, write CSV + Word (_statistics.*).
makeStatsTable <- function(csv_path, n_check = 10, out_dir = summary_out_dir) {
  results <- read.csv(csv_path, stringsAsFactors = FALSE)
  tab     <- buildStatsTable(results, n_check = n_check)
  cat("Writing statistics for", nrow(tab), "scenarios (last", n_check, "checkpoints):\n")
  writeTablePair(tab, tools::file_path_sans_ext(basename(csv_path)), out_dir,
                 suffix = "_statistics")
  invisible(tab)
}


## Calculated (derived) statistics ---------------------------------------------
# Derived / post-hoc summary statistics: measures CALCULATED (not logged directly)
# from a master results CSV. Each is a small function of one replicate's time series
# (add one via the `measures` list); the driver runs them per replicate and reports
# mean +/- SD per scenario. Covers Simpson, surviving fraction, coexistence, stability
# (inverse CV), dominance, mortality/exclusion, turnover, time-to events, colonisation
# and recovery rate. Colony-level measures need the saved grid, not the metrics table.
# --- small helpers ----------------------------------------------------------
sp_cover_cols <- function(d) grep("^Sp[0-9]+_cover$", names(d), value = TRUE)

# Simpson index D = 1 - sum(p^2) from a vector of per-species covers (NAs = absent)
simpson <- function(covers) {
  p <- covers[!is.na(covers)]; s <- sum(p)
  if (s <= 0) return(NA_real_)
  1 - sum((p / s)^2)
}

# --- the measures: each takes (r, ctx) and returns ONE number per replicate --
#   r   = one replicate's rows, ordered by timestep (the FULL time series)
#   ctx = list(rep = row indices of the last n_report checkpoints,
#              spcols, cover_thresh, dom_thresh, eq_tol)
# "state" measures average over ctx$rep (the end window); "dynamics"/"time-to"
# measures scan the whole series.
measures <- list(
  # --- end state: MEAN over the last n_report checkpoints -------------------
  "Richness"          = function(r, ctx) mean(r$species_richness[ctx$rep], na.rm = TRUE),
  "Simpson (D)"       = function(r, ctx) mean(vapply(ctx$rep, function(i)
                          simpson(unlist(r[i, ctx$spcols])), numeric(1)), na.rm = TRUE),
  "Surviving frac"    = function(r, ctx) mean(r$species_richness[ctx$rep], na.rm = TRUE) / r$n_species[1],
  "Coexistence"       = function(r, ctx) as.numeric(mean(r$species_richness[ctx$rep], na.rm = TRUE) >= r$n_species[1]),
  "Stability (cover)" = function(r, ctx) { x <- r$total_cover[ctx$rep]; s <- sd(x, na.rm = TRUE)
                          if (is.na(s) || s == 0) NA else mean(x, na.rm = TRUE) / s },
  "Stability (rich)"  = function(r, ctx) { x <- r$species_richness[ctx$rep]; s <- sd(x, na.rm = TRUE)
                          if (is.na(s) || s == 0) NA else mean(x, na.rm = TRUE) / s },
  "Small colonies"    = function(r, ctx) mean(r$n_small[ctx$rep], na.rm = TRUE),
  "Large colonies"    = function(r, ctx) mean(r$n_large[ctx$rep], na.rm = TRUE),

  # --- dynamics over the whole run ------------------------------------------
  # Dominance frequency: fraction of timesteps the single most-often dominant
  # species leads (1 = one species monopolises; low = dominance keeps shifting)
  "Dominance freq"    = function(r, ctx) { M <- as.matrix(r[ctx$spcols]); M[is.na(M)] <- -Inf
                          keep <- is.finite(apply(M, 1, max)); if (!any(keep)) return(NA_real_)
                          max(tabulate(max.col(M, "first")[keep])) / sum(keep) },
  # Dominance duration: longest unbroken run (checkpoints) one species stays top
  "Dominance dur"     = function(r, ctx) { M <- as.matrix(r[ctx$spcols]); M[is.na(M)] <- -Inf
                          keep <- is.finite(apply(M, 1, max)); if (!any(keep)) return(NA_real_)
                          max(rle(max.col(M, "first")[keep])$lengths) },
  "Turnover (mean)"   = function(r, ctx) mean(r$turnover, na.rm = TRUE),
  "Mortality %/step"  = function(r, ctx) 100 * mean(r$deaths_since_last / pmax(r$n_colonies, 1), na.rm = TRUE),
  # Colonisation rate: mean recruits per step. Uses a logged births/recruits
  # column if a later model provides one; else a PROXY = gross colony gains
  # (positive changes in n_colonies), which also counts splits, not just births.
  "Colonisation /step"= function(r, ctx) {
                          bc <- intersect(c("recruits_since_last", "births_since_last",
                                            "n_recruits", "recruits", "births"), names(r))
                          if (length(bc)) mean(r[[bc[1]]], na.rm = TRUE)
                          else            mean(pmax(diff(r$n_colonies), 0), na.rm = TRUE) },
  # Exclusion rate: species lost per 1000 timesteps (S0 - end richness)/T
  "Exclusion /1000t"  = function(r, ctx)
                          (r$n_species[1] - mean(r$species_richness[ctx$rep], na.rm = TRUE)) / max(r$timestep) * 1000,
  # Recovery rate (post-disturbance, PER SPECIES): mean regrowth slope (% cover
  # per timestep) of each species' cover over recovery_window steps after each
  # disturbance event (parsed from disturbance_steps). NA if no disturbance.
  "Recovery /step"    = function(r, ctx) {
                          ds <- r$disturbance_steps[1]
                          if (is.na(ds) || !nzchar(as.character(ds))) return(NA_real_)
                          ev <- suppressWarnings(as.numeric(strsplit(as.character(ds), ";")[[1]]))
                          ev <- ev[!is.na(ev)]; if (!length(ev)) return(NA_real_)
                          w  <- ctx$recovery_window
                          sl <- unlist(lapply(ev, function(td) vapply(ctx$spcols, function(sc) {
                            i0 <- which.min(abs(r$timestep - td))
                            i1 <- which.min(abs(r$timestep - (td + w)))
                            if (i1 <= i0 || is.na(r[[sc]][i0]) || is.na(r[[sc]][i1])) return(NA_real_)
                            (r[[sc]][i1] - r[[sc]][i0]) / (r$timestep[i1] - r$timestep[i0])
                          }, numeric(1))))
                          mean(sl, na.rm = TRUE) },

  # --- times to events (NA if never reached) --------------------------------
  "T full cover"      = function(r, ctx) { i <- which(r$total_cover >= ctx$cover_thresh)[1]
                          if (is.na(i)) NA else r$timestep[i] },
  "T 1st extinction"  = function(r, ctx) { i <- which(r$species_richness < r$n_species[1])[1]
                          if (is.na(i)) NA else r$timestep[i] },
  # complete dominance = one species holds >= dom_thresh of the occupied cover
  "T complete dom"    = function(r, ctx) { M <- as.matrix(r[ctx$spcols])
                          tot <- rowSums(M, na.rm = TRUE); mx <- suppressWarnings(apply(M, 1, max, na.rm = TRUE))
                          share <- ifelse(tot > 0, mx / tot, NA)
                          i <- which(share >= ctx$dom_thresh)[1]; if (is.na(i)) NA else r$timestep[i] },
  # persistence = last timestep at which ALL starting species are still present
  "Persistence t"     = function(r, ctx) { full <- which(r$species_richness >= r$n_species[1])
                          if (!length(full)) 0 else r$timestep[max(full)] },
  # equilibrium = first timestep after which total cover stays within eq_tol of
  # its end value (cover-based; a settling time)
  "T equilibrium"     = function(r, ctx) { final <- mean(r$total_cover[ctx$rep], na.rm = TRUE)
                          bad <- which(abs(r$total_cover - final) > ctx$eq_tol * max(abs(final), 1))
                          i <- if (length(bad)) max(bad) + 1 else 1
                          if (i > nrow(r)) NA else r$timestep[i] }
)

# --- driver: per-replicate values -> per-scenario mean +/- SD ----------------
#   n_report : the metrics are reported as the MEAN over the last n_report
#              checkpoints (default 10); time-to measures use the full series
#   cover_thresh : % cover counted as "full"; dom_thresh : share for dominance;
#   eq_tol : fractional tolerance for the equilibrium settling time;
#   recovery_window : steps after each disturbance over which recovery is measured
calcSummaryStats <- function(results, n_report = 10, cover_thresh = 90,
                             dom_thresh = 0.95, eq_tol = 0.05, recovery_window = 50,
                             digits = 3) {
  spcols <- sp_cover_cols(results)

  # one row of measures per replicate (keyed by run_id). split() once (fast) then
  # iterate the sub-frames, rather than re-subsetting the big table per group.
  per <- lapply(split(results, results$run_id), function(r) {
    r <- r[order(r$timestep), ]
    late <- tail(sort(unique(r$timestep)), n_report)      # last n_report checkpoints
    ctx  <- list(rep = which(r$timestep %in% late), spcols = spcols,
                 cover_thresh = cover_thresh, dom_thresh = dom_thresh, eq_tol = eq_tol,
                 recovery_window = recovery_window)
    vals <- vapply(measures, function(f) tryCatch(as.numeric(f(r, ctx)), error = function(e) NA_real_), numeric(1))
    data.frame(scenario = r$scenario[1], t(vals), check.names = FALSE)
  })
  vals <- do.call(rbind, per)

  # aggregate to per-scenario mean +/- SD (SD across replicates = Repeatability)
  fmt <- function(mu, sdv)
    ifelse(!is.finite(mu), "n/a",
    ifelse(!is.finite(sdv), sprintf("%.*f", digits, mu),
                            sprintf("%.*f ± %.*f", digits, mu, digits, sdv)))
  mnames <- names(measures)
  mus <- aggregate(vals[mnames], by = list(Scenario = vals$scenario), FUN = mean, na.rm = TRUE)
  sds <- aggregate(vals[mnames], by = list(Scenario = vals$scenario), FUN = sd,   na.rm = TRUE)
  out <- data.frame(Scenario = mus$Scenario, Replicates = as.integer(table(vals$scenario)[mus$Scenario]),
                    check.names = FALSE, stringsAsFactors = FALSE)
  for (m in mnames) out[[m]] <- fmt(mus[[m]], sds[[m]])
  rownames(out) <- NULL
  out
}

# --- one call: read one OR MORE master CSVs, COMBINE, compute, write CSV ----
#   csv_paths : a single path or a vector of masters to pool into one table
#   out_name  : output base name (defaults to the first file's name)
# Only the columns the measures need are read, so the big sweep stays manageable.
makeCalculatedStats <- function(csv_paths, out_name = NULL, out_dir = calc_out_dir, ...) {
  has_dt <- requireNamespace("data.table", quietly = TRUE)
  want   <- c("run_id", "scenario", "timestep", "n_species", "species_richness",
              "total_cover", "turnover", "n_colonies", "deaths_since_last",
              "disturbance_steps", "n_small", "n_large")
  read_needed <- function(p) {
    hdr  <- if (has_dt) names(data.table::fread(p, nrows = 0)) else names(read.csv(p, nrows = 1))
    keep <- intersect(c(want, grep("^Sp[0-9]+_cover$", hdr, value = TRUE),
                        grep("recruit|birth", hdr, value = TRUE, ignore.case = TRUE)), hdr)
    if (has_dt) as.data.frame(data.table::fread(p, select = keep))
    else        read.csv(p, stringsAsFactors = FALSE)[keep]
  }
  dl <- lapply(csv_paths, function(p) { cat("reading", basename(p), "...\n"); read_needed(p) })
  cols    <- Reduce(intersect, lapply(dl, names))            # bind on shared columns only
  results <- do.call(rbind, lapply(dl, function(d) d[cols]))
  cat("combined", nrow(results), "rows,", length(unique(results$run_id)), "replicates from",
      length(csv_paths), "file(s)\n")

  tab  <- calcSummaryStats(results, ...)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  base <- if (!is.null(out_name)) out_name else tools::file_path_sans_ext(basename(csv_paths[1]))
  out  <- file.path(out_dir, paste0(base, "_calculated.csv"))
  tryCatch({ write.csv(tab, out, row.names = FALSE); cat("wrote", nrow(tab), "scenarios ->", out, "\n") },
           error = function(e) cat("COULD NOT write (open in Excel?):", out, "\n"))
  invisible(tab)
}


## Sweep stats matrices (big sweep, one matrix per metric) ----------------------
#Sim 4 - Full sweep of parameters
# 396 scenarios is unreadable as a flat table, so instead we make ONE compact
# matrix PER METRIC: rows = base model, columns = manipulation type, each cell
# = mean +/- SD over the last n_check checkpoints (pooled across replicates,
# species, founders and sub-variants). One set of tables per disturbance state.

# decomposeCombo() is defined in Analysis_Utils.r.

# list of per-metric matrices (base_model rows x manip_type cols), cells "m +/- sd"
# `window`: apply the last-n_check-checkpoints filter (skip it if the data was
# already windowed per-file, as makeSweepStatsTable does when combining files).
buildSweepStatsMatrices <- function(results, n_check = 10, disturbance = FALSE,
                                    species = NULL, metrics = stat_metrics, digits = 3,
                                    window = TRUE) {
  d <- results
  if (window) {
    late <- tail(sort(unique(d$timestep)), n_check)
    d <- d[d$timestep %in% late, ]
  }
  d <- d[d$disturbance_on == disturbance, ]
  if (!is.null(species)) d <- d[d$n_species %in% species, ]
  dc <- decomposeCombo(d$combination)
  d$base_model <- dc$base_model; d$manip_type <- dc$manip_type

  models <- sort(unique(d$base_model))
  types  <- intersect(c("none", "size", "growth", "reproduction", "other"), unique(d$manip_type))

  out <- lapply(metrics, function(m) {
    cells <- outer(models, types, Vectorize(function(bm, ty) {
      v <- d[d$base_model == bm & d$manip_type == ty, m]
      if (length(v) == 0) return("")
      mu <- mean(v, na.rm = TRUE); sdv <- sd(v, na.rm = TRUE)
      if      (!is.finite(mu))  "n/a"
      else if (!is.finite(sdv)) sprintf("%.*f", digits, mu)
      else                      sprintf("%.*f ± %.*f", digits, mu, digits, sdv)
    }))
    df <- data.frame(`Base model` = models, cells, check.names = FALSE, stringsAsFactors = FALSE)
    names(df)[-1] <- types
    df
  })
  names(out) <- unname(stat_labels[metrics])
  out
}

# writeWordTables() is defined in Analysis_Utils.r.

# Read one OR MORE sweep masters, COMBINE their late-window rows into one pool,
# and write per-metric matrix tables for BOTH disturbance states: a multi-table
# Word .doc + a stacked CSV, per state. Each file is windowed to its own last
# n_check checkpoints before pooling, so files with different timestep spacing
# still each contribute their end state.
makeSweepStatsTable <- function(csv_paths, out_name, n_check = 10, species = NULL,
                                out_dir = summary_out_dir) {
  reader <- if (requireNamespace("data.table", quietly = TRUE))
              function(p) as.data.frame(data.table::fread(p))
            else function(p) read.csv(p, stringsAsFactors = FALSE)
  need <- c("timestep", "disturbance_on", "combination", "n_species", stat_metrics)

  # read each file, keep only the needed columns and its last n_check checkpoints
  results <- do.call(rbind, lapply(csv_paths, function(p) {
    cat("reading", basename(p), "...\n")
    d    <- reader(p)[need]
    late <- tail(sort(unique(d$timestep)), n_check)
    d[d$timestep %in% late, ]
  }))
  cat("combined", nrow(results), "late-window rows from", length(csv_paths), "file(s)\n")

  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  sp_tag <- if (is.null(species)) "" else paste0("_n", paste(species, collapse = "-"))
  safe <- function(path, fun)
    tryCatch({ fun(path); cat("  wrote:", path, "\n") },
             error = function(e) cat("  COULD NOT write:", path, "(open in Word/Excel?)\n"))

  for (dist in c(FALSE, TRUE)) {
    mats <- buildSweepStatsMatrices(results, disturbance = dist, species = species, window = FALSE)
    tag  <- paste0(if (dist) "_dist_on" else "_dist_off", sp_tag)
    stacked <- do.call(rbind, Map(function(m, nm) cbind(Metric = nm, m), mats, names(mats)))
    ttl <- sprintf("Combined sweep statistics (last %d checkpoints, disturbance %s%s)",
                   n_check, if (dist) "on" else "off",
                   if (is.null(species)) "" else paste0(", species ", paste(species, collapse = "/")))
    safe(file.path(out_dir, paste0(out_name, "_sweepstats", tag, ".csv")),
         function(p) write.csv(stacked, p, row.names = FALSE))
    safe(file.path(out_dir, paste0(out_name, "_sweepstats", tag, ".doc")),
         function(p) writeWordTables(mats, p, title = ttl))
  }
  invisible(NULL)
}


## Grouped stats summary --------------------------------------------------------
########Grouped stats summary for the big sweep (ONE compact table)
# Instead of one row per scenario (hundreds), group the sweep and give each
# metric as mean +/- SD per group. Default rows = Matrix x Manipulation x
# N0 Species x Disturbance (~90 rows), columns = the six community metrics.
# Stats pooled over the last n_check checkpoints x replicates within each group.
buildGroupedStats <- function(results, n_check = 10,
                              group_by = c("Matrix", "Manipulation", "N0 Species", "Disturbance"),
                              metrics = stat_metrics, digits = 3) {
  late <- tail(sort(unique(results$timestep)), n_check)
  d <- results[results$timestep %in% late, ]
  dc <- decomposeCombo(d$combination)
  d$Matrix        <- dc$base_model
  d$Manipulation  <- dc$manip_type
  d$`N0 Species`  <- d$n_species
  d$`N0 Colonies` <- d$individuals
  d$Disturbance   <- ifelse(d$disturbance_on,  "On", "Off")
  d$Reproduction  <- ifelse(d$reproduction_on, "On", "Off")

  by <- as.list(d[group_by])                              # aggregate keeps this row order for both
  mus <- aggregate(d[metrics], by, FUN = mean, na.rm = TRUE)
  sds <- aggregate(d[metrics], by, FUN = sd,   na.rm = TRUE)

  tab <- mus[group_by]
  for (m in metrics) {
    mu <- mus[[m]]; sdv <- sds[[m]]
    tab[[stat_labels[[m]]]] <-
      ifelse(!is.finite(mu),  "n/a",
      ifelse(!is.finite(sdv), sprintf("%.*f", digits, mu),
                              sprintf("%.*f ± %.*f", digits, mu, digits, sdv)))
  }
  rownames(tab) <- NULL
  tab
}

# Read one OR MORE masters, window each to its last n_check checkpoints, pool,
# and write ONE grouped stats table (CSV + Word, _statssummary.*).
makeSweepStatsSummary <- function(csv_paths, out_name, n_check = 10,
                                  group_by = c("Matrix", "Manipulation", "N0 Species", "Disturbance"),
                                  out_dir = summary_out_dir) {
  reader <- if (requireNamespace("data.table", quietly = TRUE))
              function(p, cols) as.data.frame(data.table::fread(p, select = cols))
            else function(p, cols) read.csv(p, stringsAsFactors = FALSE)[cols]
  need <- c("timestep", "combination", "n_species", "individuals",
            "disturbance_on", "reproduction_on", stat_metrics)
  results <- do.call(rbind, lapply(csv_paths, function(p) {
    cat("reading", basename(p), "...\n")
    d <- reader(p, need); late <- tail(sort(unique(d$timestep)), n_check)
    d[d$timestep %in% late, ]
  }))
  tab <- buildGroupedStats(results, n_check = n_check, group_by = group_by)
  cat("grouped stats:", nrow(tab), "rows\n")
  writeTablePair(tab, out_name, out_dir, suffix = "_statssummary")
  invisible(tab)
}


## Condition stats table ---------------------------------------------------------
########Condition stats table (breaks EVERY factor into its own column)
# For files that cross several conditions (e.g. reproduction x disturbance
# regime), one row per unique CONDITION with each metric as mean +/- SD. Factor
# columns: Matrix, Reproduction, N0 Species, Disturbance, Frequency, Size.
# Disturbance frequency/size show "-" for no-disturbance rows (and for files
# without those columns). n_check = last checkpoints pooled with the replicates.
makeConditionStatsTable <- function(csv_path, out_name = NULL, n_check = 10,
                                    out_dir = summary_out_dir, digits = 3) {
  has_dt <- requireNamespace("data.table", quietly = TRUE)
  want   <- c("timestep", "combination", "n_species", "disturbance_on",
              "reproduction_on", "dist_freq", "dist_size", stat_metrics)
  hdr  <- if (has_dt) names(data.table::fread(csv_path, nrows = 0)) else names(read.csv(csv_path, nrows = 1))
  keep <- intersect(want, hdr)
  d <- if (has_dt) as.data.frame(data.table::fread(csv_path, select = keep))
       else        read.csv(csv_path, stringsAsFactors = FALSE)[keep]

  d <- d[d$timestep %in% tail(sort(unique(d$timestep)), n_check), ]   # end window
  blank <- function(x) is.na(x) | !nzchar(as.character(x))
  d$Matrix        <- decomposeCombo(d$combination)$base_model
  d$Reproduction  <- ifelse(d$reproduction_on, "On", "Off")
  d$`N0 Species`  <- d$n_species
  d$Disturbance   <- ifelse(d$disturbance_on, "On", "Off")
  d$Frequency     <- if ("dist_freq" %in% names(d)) ifelse(blank(d$dist_freq), "-", as.character(d$dist_freq)) else "-"
  d$Size          <- if ("dist_size" %in% names(d)) ifelse(blank(d$dist_size), "-", as.character(d$dist_size)) else "-"

  gcols <- c("Matrix", "Reproduction", "N0 Species", "Disturbance", "Frequency", "Size")
  by    <- as.list(d[gcols])
  mus   <- aggregate(d[stat_metrics], by, FUN = mean, na.rm = TRUE)
  sds   <- aggregate(d[stat_metrics], by, FUN = sd,   na.rm = TRUE)
  tab <- mus[gcols]
  for (m in stat_metrics) {
    mu <- mus[[m]]; sdv <- sds[[m]]
    tab[[stat_labels[[m]]]] <-
      ifelse(!is.finite(mu),  "n/a",
      ifelse(!is.finite(sdv), sprintf("%.*f", digits, mu),
                              sprintf("%.*f ± %.*f", digits, mu, digits, sdv)))
  }
  tab <- tab[do.call(order, tab[gcols]), ]; rownames(tab) <- NULL
  cat("condition stats:", nrow(tab), "rows\n")
  writeTablePair(tab, if (!is.null(out_name)) out_name else tools::file_path_sans_ext(basename(csv_path)),
                 out_dir, suffix = "_conditions")
  invisible(tab)
}

# Diversity retention, replicate by replicate (from the de-duplicated index; RPS and
# Random only, as Linear/Neutral collapse deterministically). Because a scenario's
# replicates are often bimodal (full coexistence vs collapse), the mean is misleading,
# so we summarise the distribution: p_retain (with Wilson 95% CI), the two modes
# (p_full / p_collapse), median, bimodality, and flags (bimodal / mean-misleading /
# rescue = a rare survivor amid collapse). Retention is relative to the initial
# community via Hill numbers (div_ret = exp(Shannon)_end / n_species; retains if
# >= ret_thresh, default 0.5); n_species and disturbance are in the grouping.
# The retention/loss functions below need data.table; load it if available so sourcing
# this file never fails (they will error clearly if called without it).
if (requireNamespace("data.table", quietly = TRUE)) suppressMessages(library(data.table)) else
  message("Note: retention/loss table functions require the 'data.table' package.")

combined_index_path <- "Simulation_Results/0_Combined_Master/run_index.rds"

# Parameter columns that define one "combination setting" (a factorial cell).
RETENTION_GROUP_COLS <- c("combination", "n_species", "disturbance_on", "individuals",
                          "bias", "intraspecific", "repro_base", "background_chance",
                          "maturity_age", "size_beta_max", "placement", "size_mode")

# Wilson 95% CI for a proportion - honest at extreme p and small replicate counts,
# where the normal approximation would give impossible (<0 or >1) bounds.
wilsonCI <- function(x, n, z = 1.96) {
  if (n == 0) return(c(NA_real_, NA_real_))
  p <- x / n; d <- 1 + z^2 / n
  half <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / d
  c(max(0, (p + z^2 / (2 * n)) / d - half), min(1, (p + z^2 / (2 * n)) / d + half))
}

# Sarle's bimodality coefficient; BC > 5/9 (~0.555) indicates a likely bimodal /
# heavily two-moded distribution (uniform = 5/9). NA when too few replicates.
bimodCoef <- function(x) {
  x <- x[is.finite(x)]; n <- length(x)
  if (n < 4) return(NA_real_)
  s <- sd(x); if (s == 0) return(0)
  m <- mean(x)
  g <- mean((x - m)^3) / s^3                       # skewness
  k <- mean((x - m)^4) / s^4 - 3                   # excess kurtosis
  (g^2 + 1) / (k + 3 * (n - 1)^2 / ((n - 2) * (n - 3)))
}

# ---- Per-replicate retention metrics (RPS/Random only, de-duplicated) --------
# (every retention function calls this first, so attaching data.table here means the
#  whole analysis works even if you run the functions individually.)
retentionReplicates <- function(index, exclude = "Linear|Neutral", drop_dups = TRUE) {
  suppressMessages(library(data.table))
  d <- as.data.table(index)
  d <- d[!grepl(exclude, combination)]
  if (drop_dups && "is_duplicate" %in% names(d)) d <- d[is_duplicate == FALSE]
  d[, hill1    := exp(shannon)]                        # effective no. of species (q = 1)
  d[, div_ret  := pmin(1, hill1 / n_species)]          # relative Hill-1 retention (0-1)
  d[, rich_ret := species_richness / n_species]        # proportion of species persisting
  d[]
}

# ---- Per-combination distribution summary + flags ----------------------------
#   exclude : networks dropped before analysis (regex on `combination`).
#             "Linear|Neutral" = RPS + Random (this retention/rescue analysis).
retentionSummary <- function(index, ret_thresh = 0.5, rescue_thresh = 0.25,
                             min_rep = 10, group_cols = RETENTION_GROUP_COLS,
                             exclude = "Linear|Neutral") {
  d  <- retentionReplicates(index, exclude = exclude)
  gc <- intersect(group_cols, names(d))
  d[, `:=`(retained = div_ret >= ret_thresh,
           full     = species_richness == n_species,
           collapse = species_richness <= 1)]
  s <- d[, {
    ci <- wilsonCI(sum(retained), .N)
    list(n_rep = .N, n_retain = sum(retained),
         p_retain = mean(retained), ci_lo = ci[1], ci_hi = ci[2],
         p_full = mean(full), p_collapse = mean(collapse),
         div_mean = mean(div_ret), div_median = as.numeric(median(div_ret)),
         div_sd = sd(div_ret), div_min = min(div_ret), div_max = max(div_ret),
         bimod = bimodCoef(div_ret))
  }, by = gc]
  s <- s[n_rep >= min_rep]
  s[, flag_heterogeneous   := n_retain > 0 & n_retain < n_rep]
  # Genuine two-attractor bimodality: BOTH the coexistence and the collapse mode
  # are populated (>=15% each). This - not the raw bimodality coefficient, which
  # also fires on tight discrete distributions - is when the mean is an artefact.
  s[, flag_bimodal         := p_retain >= 0.15 & p_collapse >= 0.15]
  # The mean cannot be trusted if bimodal, if it diverges from the median, or if
  # replicates straddle the retention threshold with real spread.
  s[, flag_mean_misleading := flag_bimodal | abs(div_mean - div_median) >= 0.15 |
                              (div_min < ret_thresh & div_max >= ret_thresh & div_sd >= 0.20)]
  # Rare "rescue": at least one replicate retains diversity while MOST collapse -
  # the survivor run a scenario mean would erase.
  s[, flag_rescue          := n_retain >= 1 & p_retain <= rescue_thresh]
  setorder(s, -p_retain, -div_median)
  s[]
}

# ---- The individual "survivor" replicates inside flagged scenarios -----------
# For every scenario flagged rescue OR bimodal, list the actual replicates that
# retained diversity - the runs a scenario mean would hide.
rescueReplicates <- function(index, ret_thresh = 0.5, rescue_thresh = 0.25,
                             min_rep = 10, group_cols = RETENTION_GROUP_COLS) {
  d  <- retentionReplicates(index)
  gc <- intersect(group_cols, names(d))
  s  <- retentionSummary(index, ret_thresh, rescue_thresh, min_rep, group_cols)
  ctx <- c(gc, "n_rep", "p_retain", "flag_rescue", "flag_bimodal")
  flagged <- s[flag_rescue == TRUE | flag_bimodal == TRUE, ..ctx]
  # survivor replicates = the runs that retained diversity inside those scenarios,
  # each tagged with its scenario's retention probability and which flag it triggered
  keep <- d[div_ret >= ret_thresh][flagged, on = gc, nomatch = 0]
  keep[, outcome := ifelse(flag_rescue, "rescue", "bimodal")]
  cols <- c(gc, "n_rep", "p_retain", "outcome", "run_id", "replicate", "seed",
            "species_richness", "div_ret", "rich_ret", "source_file")
  keep <- keep[, intersect(cols, names(keep)), with = FALSE]
  setorderv(keep, c("p_retain", gc, "div_ret"), c(1L, rep(1L, length(gc)), -1L))
  keep[]
}

# ---- Word/CSV-friendly formatting -------------------------------------------
formatRetention <- function(s, top = NULL) {
  s <- copy(s)
  if (!is.null(top)) s <- head(s, top)
  yn <- function(x) ifelse(x, "Y", "")
  data.frame(
    combination = s$combination, n = s$n_species,
    disturbance = ifelse(s$disturbance_on, "on", "off"), founders = s$individuals,
    bias = s$bias, intra = s$intraspecific, repro = s$repro_base,
    backg = s$background_chance,
    reps = s$n_rep,
    P_retain = sprintf("%.2f [%.2f-%.2f]", s$p_retain, s$ci_lo, s$ci_hi),
    P_full = round(s$p_full, 2), P_collapse = round(s$p_collapse, 2),
    div_median = round(s$div_median, 2), div_mean = round(s$div_mean, 2),
    bimodal = yn(s$flag_bimodal), mean_misleading = yn(s$flag_mean_misleading),
    rescue = yn(s$flag_rescue), stringsAsFactors = FALSE)
}

# ---- One call: build + save the retention tables -----------------------------
makeRetentionTables <- function(index_path = combined_index_path,
                                out_dir = summary_out_dir, out_name = "RPS_Random_diversity_retention",
                                ret_thresh = 0.5, rescue_thresh = 0.25, min_rep = 10) {
  idx <- readRDS(index_path)
  s   <- retentionSummary(idx, ret_thresh, rescue_thresh, min_rep)
  resc <- rescueReplicates(idx, ret_thresh, rescue_thresh, min_rep)

  # 1. Full ranked leaderboard (by probability of retention, mode-robust)
  writeTablePair(formatRetention(s), out_name, out_dir, suffix = "_ranked")
  # 2. Just the scenarios where the mean is untrustworthy (bimodal / rescue)
  flagged <- s[flag_mean_misleading == TRUE | flag_rescue == TRUE]
  writeTablePair(formatRetention(flagged), out_name, out_dir, suffix = "_flagged")
  # 3. The individual survivor replicates inside those flagged scenarios
  writeTablePair(as.data.frame(resc), out_name, out_dir, suffix = "_survivor_replicates")

  cat("\nRetention analysis (RPS + Random,", nrow(s), "combinations, >=", min_rep, "reps):\n")
  cat("  consistently retaining (p_retain >= 0.8):", sum(s$p_retain >= 0.8), "\n")
  cat("  bimodal / mean-misleading:", sum(s$flag_mean_misleading), "\n")
  cat("  rescue (rare survivor amid collapse):", sum(s$flag_rescue), "\n")
  cat("\nTop 10 by probability of retention:\n")
  print(head(formatRetention(s)[, c("combination","n","disturbance","P_retain",
                                    "P_full","P_collapse","div_median","bimodal","rescue")], 10))
  invisible(list(summary = s, rescue = resc))
}

# Anomalous richness loss - the mirror of the rescue analysis (networks RPS, Random,
# Linear; Neutral dropped as it retains everything). Within a scenario that normally
# retains (p_retain >= loss_norm_thresh, default 0.75), flags the replicates that
# nonetheless lost diversity. metric = "richness" (species kept, default) or
# "diversity" (Hill-1). Severity = species dropped and fraction of the metric lost.

# Per-combination retention summary for the loss analysis (richness or diversity),
# with a Wilson CI and the "normally retains yet loses" flag. Self-contained so the
# loss analysis does not depend on the RPS/Random retention summary above.
lossSummary <- function(index, ret_thresh = 0.5, loss_norm_thresh = 0.75, min_rep = 10,
                        group_cols = RETENTION_GROUP_COLS, exclude = "Neutral",
                        metric = c("richness", "diversity")) {
  metric <- match.arg(metric)
  d  <- retentionReplicates(index, exclude = exclude)
  gc <- intersect(group_cols, names(d))
  d[, val := if (metric == "richness") rich_ret else div_ret]
  d[, retained := val >= ret_thresh]
  s <- d[, {
    ci <- wilsonCI(sum(retained), .N)
    list(n_rep = .N, n_retain = sum(retained), p_retain = mean(retained),
         ci_lo = ci[1], ci_hi = ci[2],
         p_collapse = mean(species_richness <= 1),
         val_median = as.numeric(median(val)), val_min = min(val), val_max = max(val))
  }, by = gc]
  s <- s[n_rep >= min_rep]
  # NORMALLY retains yet AT LEAST ONE replicate lost it -> the surprising cases
  s[, flag_anomalous_loss := (n_rep - n_retain) >= 1 & p_retain >= loss_norm_thresh]
  setorder(s, -p_retain, -val_median)
  s[]
}

# Word/CSV-friendly view of the scenarios that normally retain but occasionally lose.
formatLoss <- function(s, top = NULL) {
  s <- copy(s)
  if (!is.null(top)) s <- head(s, top)
  data.frame(
    combination = s$combination, n = s$n_species,
    disturbance = ifelse(s$disturbance_on, "on", "off"), founders = s$individuals,
    bias = s$bias, intra = s$intraspecific, repro = s$repro_base, backg = s$background_chance,
    reps = s$n_rep, n_lost = s$n_rep - s$n_retain,
    P_retain = sprintf("%.2f [%.2f-%.2f]", s$p_retain, s$ci_lo, s$ci_hi),
    P_collapse = round(s$p_collapse, 2), median = round(s$val_median, 2),
    stringsAsFactors = FALSE)
}

# The individual replicates that LOST diversity inside normally-retaining scenarios.
# Mirror of rescueReplicates(); includes Linear (exclude = "Neutral").
lossReplicates <- function(index, ret_thresh = 0.5, loss_norm_thresh = 0.75,
                           min_rep = 10, group_cols = RETENTION_GROUP_COLS,
                           exclude = "Neutral", metric = c("richness", "diversity")) {
  metric <- match.arg(metric)
  d  <- retentionReplicates(index, exclude = exclude)
  gc <- intersect(group_cols, names(d))
  d[, val := if (metric == "richness") rich_ret else div_ret]
  s  <- lossSummary(index, ret_thresh, loss_norm_thresh, min_rep, group_cols, exclude, metric)
  ctx <- c(gc, "n_rep", "p_retain", "flag_anomalous_loss")
  flagged <- s[flag_anomalous_loss == TRUE, ..ctx]
  lost <- d[val < ret_thresh][flagged, on = gc, nomatch = 0]
  lost[, lost_species := n_species - species_richness]      # species dropped vs the start
  lost[, frac_lost    := round(1 - val, 3)]                 # fraction of the metric lost
  cols <- c(gc, "n_rep", "p_retain", "run_id", "replicate", "seed",
            "species_richness", "lost_species", "rich_ret", "div_ret", "frac_lost",
            "source_file")
  lost <- lost[, intersect(cols, names(lost)), with = FALSE]
  # most-reliable scenarios first (highest p_retain = most surprising loss), worst loss first
  setorderv(lost, c("p_retain", gc, "frac_lost"), c(-1L, rep(1L, length(gc)), -1L))
  lost[]
}

# ---- One call: build + save the anomalous-loss tables ------------------------
makeLossTables <- function(index_path = combined_index_path, out_dir = summary_out_dir,
                           out_name = "RPS_Random_Linear_richness_loss",
                           ret_thresh = 0.5, loss_norm_thresh = 0.75, min_rep = 10,
                           metric = "richness") {
  idx  <- readRDS(index_path)
  s    <- lossSummary(idx, ret_thresh, loss_norm_thresh, min_rep, exclude = "Neutral", metric = metric)
  lost <- lossReplicates(idx, ret_thresh, loss_norm_thresh, min_rep, exclude = "Neutral", metric = metric)

  # 1. Scenarios that normally retain yet suffer occasional loss
  flagged <- s[flag_anomalous_loss == TRUE][order(-p_retain, -p_collapse)]
  writeTablePair(formatLoss(flagged), out_name, out_dir, suffix = "_loss_scenarios")
  # 2. The individual replicates that lost diversity inside those scenarios
  writeTablePair(as.data.frame(lost), out_name, out_dir, suffix = "_loss_replicates")

  cat("\nAnomalous-loss analysis (RPS + Random + Linear; metric =", metric, "):\n")
  cat("  scenarios normally retaining (p_retain >=", loss_norm_thresh, "):",
      sum(s$p_retain >= loss_norm_thresh), "\n")
  cat("  of those, WITH >=1 anomalous loss:", sum(s$flag_anomalous_loss), "\n")
  cat("  individual anomalous-loss replicates:", nrow(lost), "\n")
  cat("\nMost surprising losses (most-reliable scenario first, biggest drop):\n")
  print(head(lost[, c("combination", "n_species", "disturbance_on", "p_retain",
                      "species_richness", "lost_species", "frac_lost")], 10))
  invisible(list(scenarios = s, loss = lost))
}


## Run section ------------------------------------------------------------------
#everything below regenerates tables - run the blocks you need
#(design tables for sims 1 and 2 were already done and registered)

#Sim 2, robustness checks: length (timesteps), reef size, replicates
# NOTE: these regenerate tables from result files (not shipped with the code).
# Sourcing this file only DEFINES the functions above; to regenerate a table,
# set RUN_TABLES <- TRUE (with the result files present) or run a block manually.
RUN_TABLES <- FALSE
if (RUN_TABLES) {
makeStatsTable("Simulation_Results/2_Robustness_Checks/TimeConvergence(Timesteps)/Master_Robustness_TimeConvergence_20260727_110034_results.csv")
makeStatsTable("Simulation_Results/2_Robustness_Checks/ReefSize/Master_Robustness_ReefSize_20260727_132335_results.csv")
makeStatsTable("Simulation_Results/2_Robustness_Checks/NumReplicates/Master_Robustness_Replicates_20260727_115226_results.csv")
makeCalculatedStats("Simulation_Results/2_Robustness_Checks/TimeConvergence(Timesteps)/Master_Robustness_TimeConvergence_20260727_110034_results.csv")
makeCalculatedStats("Simulation_Results/2_Robustness_Checks/ReefSize/Master_Robustness_ReefSize_20260727_132335_results.csv")
makeCalculatedStats("Simulation_Results/2_Robustness_Checks/NumReplicates/Master_Robustness_Replicates_20260727_115226_results.csv")

#Sim 3, disturbance regime sensitivity
makeStatsTable("Simulation_Results/3_Parameter_Sensitivity/Run_20260727_Sensitivity_DisturbanceRegime/Master_Sensitivity_DisturbanceRegime_20260727_213926_results.csv")
makeCalculatedStats("Simulation_Results/3_Parameter_Sensitivity/Dist_Sens_Test/Master_Sensitivity_DisturbanceRegime_20260727_213926_results.csv")

#Sim 3, bias sensitivity
makeSummaryTable("Simulation_Results/3_Parameter_Sensitivity/Bias_Sensitivity_Checks/Master_Sensitivity_Bias_20260731_010903_results.csv")
makeStatsTable("Simulation_Results/3_Parameter_Sensitivity/Bias_Sensitivity_Checks/Master_Sensitivity_Bias_20260731_010903_results.csv")
makeCalculatedStats("Simulation_Results/3_Parameter_Sensitivity/Bias_Sensitivity_Checks/Master_Sensitivity_Bias_20260731_010903_results.csv")

#Sim 3, intraspecific competition sensitivity
makeSummaryTable("Simulation_Results/3_Parameter_Sensitivity/Intra_Sens_Test/Master_Sensitivity_Intraspecific_20260731_054833_results.csv")
makeStatsTable("Simulation_Results/3_Parameter_Sensitivity/Intra_Sens_Test/Master_Sensitivity_Intraspecific_20260731_054833_results.csv")
makeCalculatedStats("Simulation_Results/3_Parameter_Sensitivity/Intra_Sens_Test/Master_Sensitivity_Intraspecific_20260731_054833_results.csv")

#Sim 3, reproduction sensitivity
makeSummaryTable("Simulation_Results/3_Parameter_Sensitivity/Repro_Sens_Test/Reproduction_Fecundity_20260731_222545_results.csv")
makeStatsTable("Simulation_Results/3_Parameter_Sensitivity/Repro_Sens_Test/Reproduction_Fecundity_20260731_222545_results.csv")
makeCalculatedStats("Simulation_Results/3_Parameter_Sensitivity/Repro_Sens_Test/Reproduction_Fecundity_20260731_222545_results.csv")

#Sim 3, reproduction x disturbance regime (mixed conditions)
makeSummaryTable("Simulation_Results/3_Parameter_Sensitivity/Repro_Dist_sens_Test/Reproduction_DisturbanceRegime_20260801_123052_results.csv")
makeStatsTable("Simulation_Results/3_Parameter_Sensitivity/Repro_Dist_sens_Test/Reproduction_DisturbanceRegime_20260801_123052_results.csv")
makeCalculatedStats("Simulation_Results/3_Parameter_Sensitivity/Repro_Dist_sens_Test/Reproduction_DisturbanceRegime_20260801_123052_results.csv")
makeConditionStatsTable("Simulation_Results/3_Parameter_Sensitivity/Repro_Dist_sens_Test/Reproduction_DisturbanceRegime_20260801_123052_results.csv")

#Sim 3.1, background recruitment
makeSummaryTable("Simulation_Results/3.1_Extended_Sensitivity/BackRecruit_Sens/BackgroundRecruitment_20260802_025433_results.csv")
makeStatsTable("Simulation_Results/3.1_Extended_Sensitivity/BackRecruit_Sens/BackgroundRecruitment_20260802_025433_results.csv")
makeCalculatedStats("Simulation_Results/3.1_Extended_Sensitivity/BackRecruit_Sens/BackgroundRecruitment_20260802_025433_results.csv")

#Sim 3.1, growth rate
makeSummaryTable("Simulation_Results/3.1_Extended_Sensitivity/GrowthRate_Sens/Growth_Rate_20260801_230519_results.csv")
makeStatsTable("Simulation_Results/3.1_Extended_Sensitivity/GrowthRate_Sens/Growth_Rate_20260801_230519_results.csv")
makeCalculatedStats("Simulation_Results/3.1_Extended_Sensitivity/GrowthRate_Sens/Growth_Rate_20260801_230519_results.csv")

#Sim 3.1, maturity age
makeSummaryTable("Simulation_Results/3.1_Extended_Sensitivity/MaturityAge_Sens/MaturityAge_20260802_071441_results.csv")
makeStatsTable("Simulation_Results/3.1_Extended_Sensitivity/MaturityAge_Sens/MaturityAge_20260802_071441_results.csv")
makeCalculatedStats("Simulation_Results/3.1_Extended_Sensitivity/MaturityAge_Sens/MaturityAge_20260802_071441_results.csv")

#Sim 3.1, size competition strength
makeSummaryTable("Simulation_Results/3.1_Extended_Sensitivity/SizeComp_Sens/SizeCompetitionStrength_20260802_115434_results.csv")
makeStatsTable("Simulation_Results/3.1_Extended_Sensitivity/SizeComp_Sens/SizeCompetitionStrength_20260802_115434_results.csv")
makeCalculatedStats("Simulation_Results/3.1_Extended_Sensitivity/SizeComp_Sens/SizeCompetitionStrength_20260802_115434_results.csv")

#Sim 4, base model exploration - RPS robustness check
makeStatsTable("Simulation_Results/4_Base_model_exploration/Base_Models/Run_20260724_Initial_RPS_RawOutputs/my_run_20260724_174603_results.csv")
makeCalculatedStats("Simulation_Results/4_Base_model_exploration/RPS_Model/my_run_20260724_174603_results.csv")

#Sim 4, base model exploration - one combined design table across the 3 masters
sim4_files <- c(
  "Simulation_Results/4_Base_model_exploration/Full_Run_20260724_Classic-RPS-Linear-Neutral_n357_reef50_t1000_r30/my_run_20260724_192134_results.csv",
  "Simulation_Results/4_Base_model_exploration/Run_20260724_Initial_RPS_RawOutputs/my_run_20260724_174603_results.csv",
  "Simulation_Results/4_Base_model_exploration/Run_20260726_Remaining_Basic_Model_Exploration/AllModels_sweep_20260726_153222_results.csv"
)
makeCombinedSummary(sim4_files, out_name = "Sim4_BaseModelExploration_combined")

#Sim 4, full sweep - per-metric matrices + grouped stats + combined calculated stats
makeSweepStatsTable(
  c("Simulation_Results/4_Base_model_exploration/Base_Models/Full_Run_20260724_Classic-RPS-Linear-Neutral_n357_reef50_t1000_r30/my_run_20260724_192134_results.csv",
    "Simulation_Results/4_Base_model_exploration/Base_Models/Run_20260726_Remaining_Basic_Model_Exploration/AllModels_sweep_20260726_153222_results.csv"),
  out_name = "Sim4_FullSweep_combined")
makeSweepStatsSummary(
  c("Simulation_Results/4_Base_model_exploration/Base_Models/Full_Run_20260724_Classic-RPS-Linear-Neutral_n357_reef50_t1000_r30/my_run_20260724_192134_results.csv",
    "Simulation_Results/4_Base_model_exploration/Base_Models/Run_20260726_Remaining_Basic_Model_Exploration/AllModels_sweep_20260726_153222_results.csv"),
  out_name = "Sim4_FullSweep_combined")
makeCalculatedStats(
  c("Simulation_Results/4_Base_model_exploration/Full_Model/Run_20260726_Remaining_Basic_Model_Exploration/AllModels_sweep_20260726_153222_results.csv",
    "Simulation_Results/4_Base_model_exploration/Full_Model/Full_Run_20260724_Classic-RPS-Linear-Neutral_n357_reef50_t1000_r30/my_run_20260724_192134_results.csv"),
  out_name = "Sim4_FullModel_combined")

#Sim 5, initial placement
makeSummaryTable("Simulation_Results/5_InitialPlacement/InitialPlacement_20260802_215722_results.csv")
makeStatsTable("Simulation_Results/5_InitialPlacement/InitialPlacement_20260802_215722_results.csv")
makeCalculatedStats("Simulation_Results/5_InitialPlacement/InitialPlacement_20260802_215722_results.csv")

#Retention / anomalous-loss analyses (need run_index.rds - run Consolidate_Datasets.r first)
retention_tables <- makeRetentionTables()
loss_tables      <- makeLossTables()

}
