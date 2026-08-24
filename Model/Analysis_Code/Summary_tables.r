# ==========================================================================
#  ACTUAL SUMMARY TABLES
#  One row per scenario describing its design/settings. Reads a master results
#  CSV and writes a Word-usable table (a real, editable table - not an image):
#    - <name>_scenario_summary.csv   (open in Excel, or Insert > Table in Word)
#    - <name>_scenario_summary.doc   (HTML table; Word opens it as a real table)
#  Run makeSummaryTable(csv_path) once per master results file.
# ==========================================================================

# Shared paths + helpers (summary_out_dir, registry_path, assignScenarioIDs,
# writeWordTable(s), writeTablePair, stat_metrics/labels, decomposeCombo) all come
# from Analysis_Utils.r - defined in one place, no side effects on source.
source("Model/Analysis_Code/Analysis_Utils.r")

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

# --- Usage: one master file per call (add more calls for more files) ---------

#For sim 2 ALL DONE

#Sim 1, base RPS checks ALL DONE

#Sim 3, Parameter Sensitivity (disturbance regime) ALL DONE

#Sim 3, bias sensitivity - register + scenario table
makeSummaryTable("Simulation_Results/3_Parameter_Sensitivity/Bias_Sensitivity_Checks/Master_Sensitivity_Bias_20260731_010903_results.csv")

#Sim 3, intraspecific competition sensitivity - register + scenario table
makeSummaryTable("Simulation_Results/3_Parameter_Sensitivity/Intra_Sens_Test/Master_Sensitivity_Intraspecific_20260731_054833_results.csv")

#Sim 3 reproduction sensitivity 
makeSummaryTable("Simulation_Results/3_Parameter_Sensitivity/Repro_Sens_Test/Reproduction_Fecundity_20260731_222545_results.csv")

#Sim 3 reproduction x disturbance regime 
makeSummaryTable("Simulation_Results/3_Parameter_Sensitivity/Repro_Dist_sens_Test/Reproduction_DisturbanceRegime_20260801_123052_results.csv")

#Sim 3.1 Background recruitment
makeSummaryTable("Simulation_Results/3.1_Extended_Sensitivity/BackRecruit_Sens/BackgroundRecruitment_20260802_025433_results.csv")

#Sim 3.1 
makeSummaryTable("Simulation_Results/3.1_Extended_Sensitivity/GrowthRate_Sens/Growth_Rate_20260801_230519_results.csv")

#Sim 3.1 Maturity Age
makeSummaryTable("Simulation_Results/3.1_Extended_Sensitivity/MaturityAge_Sens/MaturityAge_20260802_071441_results.csv")

#Sim 3.1 Size competition
makeSummaryTable("Simulation_Results/3.1_Extended_Sensitivity/SizeComp_Sens/SizeCompetitionStrength_20260802_115434_results.csv")

#Sim 5 initial placement
makeSummaryTable("Simulation_Results/5_InitialPlacement/InitialPlacement_20260802_215722_results.csv")

#Sim 4, model exploration - ONE combined scenario table across the 3 master files
sim4_files <- c(
  "Simulation_Results/4_Base_model_exploration/Full_Run_20260724_Classic-RPS-Linear-Neutral_n357_reef50_t1000_r30/my_run_20260724_192134_results.csv",
  "Simulation_Results/4_Base_model_exploration/Run_20260724_Initial_RPS_RawOutputs/my_run_20260724_174603_results.csv",
  "Simulation_Results/4_Base_model_exploration/Run_20260726_Remaining_Basic_Model_Exploration/AllModels_sweep_20260726_153222_results.csv"
)
makeCombinedSummary(sim4_files, out_name = "Sim4_BaseModelExploration_combined")






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

# --- Usage: one master file per call (edit the path / n_check) ---------------


#SIM 2 Length (timesteps) test
makeStatsTable("Simulation_Results/2_Robustness_Checks/TimeConvergence(Timesteps)/Master_Robustness_TimeConvergence_20260727_110034_results.csv")

#Sim 2 reef size test
makeStatsTable("Simulation_Results/2_Robustness_Checks/ReefSize/Master_Robustness_ReefSize_20260727_132335_results.csv")

#Sim 2 number of replicates test
makeStatsTable("Simulation_Results/2_Robustness_Checks/NumReplicates/Master_Robustness_Replicates_20260727_115226_results.csv")

#Sim 3 disturbance sensitivity
makeStatsTable("Simulation_Results/3_Parameter_Sensitivity/Run_20260727_Sensitivity_DisturbanceRegime/Master_Sensitivity_DisturbanceRegime_20260727_213926_results.csv")

#Sim 3 bias
makeStatsTable("Simulation_Results/3_Parameter_Sensitivity/Bias_Sensitivity_Checks/Master_Sensitivity_Bias_20260731_010903_results.csv")

#Sim 3 intraspecific competition 
makeStatsTable("Simulation_Results/3_Parameter_Sensitivity/Intra_Sens_Test/Master_Sensitivity_Intraspecific_20260731_054833_results.csv")

#Sim 3 reproduction
makeStatsTable("Simulation_Results/3_Parameter_Sensitivity/Repro_Sens_Test/Reproduction_Fecundity_20260731_222545_results.csv")

#Sim 3 reproduction x disturbance regime (mixed conditions)
makeStatsTable("Simulation_Results/3_Parameter_Sensitivity/Repro_Dist_sens_Test/Reproduction_DisturbanceRegime_20260801_123052_results.csv")



#Sim 3.1 growth rate
makeStatsTable("Simulation_Results/3.1_Extended_Sensitivity/GrowthRate_Sens/Growth_Rate_20260801_230519_results.csv")

#Sim 3.1 Size  competition
makeStatsTable("Simulation_Results/3.1_Extended_Sensitivity/SizeComp_Sens/SizeCompetitionStrength_20260802_115434_results.csv")

#Sim 3.1 Background recruitment
makeStatsTable("Simulation_Results/3.1_Extended_Sensitivity/BackRecruit_Sens/BackgroundRecruitment_20260802_025433_results.csv")

#Sim 3.1 Maturity Age
makeStatsTable("Simulation_Results/3.1_Extended_Sensitivity/MaturityAge_Sens/MaturityAge_20260802_071441_results.csv")



#Sim 5 - initial placement
makeStatsTable("Simulation_Results/5_InitialPlacement/InitialPlacement_20260802_215722_results.csv")


#Sim 4 - RPS ROBUSTNESS check
makeStatsTable("Simulation_Results/4_Base_model_exploration/Base_Models/Run_20260724_Initial_RPS_RawOutputs/my_run_20260724_174603_results.csv")



############SIMULATION 4 ISSUES AS BIG SWEEP#########



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

# Combine the Full_Run baseline sweep + the AllModels variant sweep into one set:
makeSweepStatsTable(
  c("Simulation_Results/4_Base_model_exploration/Base_Models/Full_Run_20260724_Classic-RPS-Linear-Neutral_n357_reef50_t1000_r30/my_run_20260724_192134_results.csv",
    "Simulation_Results/4_Base_model_exploration/Base_Models/Run_20260726_Remaining_Basic_Model_Exploration/AllModels_sweep_20260726_153222_results.csv"),
  out_name = "Sim4_FullSweep_combined")


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

makeSweepStatsSummary(
  c("Simulation_Results/4_Base_model_exploration/Base_Models/Full_Run_20260724_Classic-RPS-Linear-Neutral_n357_reef50_t1000_r30/my_run_20260724_192134_results.csv",
    "Simulation_Results/4_Base_model_exploration/Base_Models/Run_20260726_Remaining_Basic_Model_Exploration/AllModels_sweep_20260726_153222_results.csv"),
  out_name = "Sim4_FullSweep_combined")


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

#Sim 3 - reproduction x disturbance regime (mixed conditions)
makeConditionStatsTable("Simulation_Results/3_Parameter_Sensitivity/Repro_Dist_sens_Test/Reproduction_DisturbanceRegime_20260801_123052_results.csv")


# ============================================================================
#  DIVERSITY RETENTION - which combinations keep diversity, replicate by replicate
# ----------------------------------------------------------------------------
#  Works from the COMBINED, de-duplicated run index (one row per replicate). Only
#  the RPS and Random networks are considered (Linear/Neutral excluded), because
#  those two collapse deterministically and would dominate any "worst" tail.
#
#  WHY REPLICATE-LEVEL, NOT THE MEAN: in these networks a single scenario's
#  replicates are often BIMODAL - some runs hold full coexistence, others collapse
#  to one species - so the scenario mean is an artefact of two modes and is
#  misleading. We therefore summarise the DISTRIBUTION across replicates:
#    * p_retain  - the PROBABILITY a replicate retains diversity (with a Wilson
#                  95% CI), the ecologically meaningful, mode-robust statistic;
#    * p_full / p_collapse - the two modes (all species kept / down to one);
#    * div_median (not just the mean) and Sarle's bimodality coefficient;
#  and we FLAG scenarios where the mean cannot be trusted:
#    * flag_bimodal        - both coexistence and collapse modes are populated;
#    * flag_mean_misleading- bimodal, or mean and median diverge markedly;
#    * flag_rescue         - >=1 replicate retains diversity while MOST fail
#                            (the rare "survivor" run a mean would erase).
#
#  Retention is measured RELATIVE to the (even) initial community using Hill
#  numbers, so 3- and 7-species runs are comparable and evenness is not ignored:
#    div_ret = exp(Shannon)_end / n_species   (proportion of initial Hill-1 kept).
#  A replicate "retains" if div_ret >= ret_thresh (default 0.5). Because retention
#  falls with richness and disturbance, n_species and disturbance are part of the
#  grouping, so combinations are compared like-for-like.
# ============================================================================
if (!requireNamespace("data.table", quietly = TRUE)) stop("please install.packages('data.table')")
suppressMessages(library(data.table))

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

retention_tables <- makeRetentionTables()


library(data.table)
s <- retention_tables$summary          # numeric summary; or: readRDS + re-run makeRetentionTables()

# 1. ROBUST maintenance — the safe region of parameter space.
#    Diversity reliably persists; the mean is trustworthy here.
reliable <- s[p_retain >= 0.8 & flag_bimodal == FALSE][order(-p_retain, -div_median)]

# 2. STOCHASTIC / alternative-stable-state maintenance — the interesting middle.
#    Same settings coexist in some runs, collapse in others (priority/founder effects).
bistable <- s[flag_bimodal == TRUE][order(-abs(p_full - p_collapse))]   # most evenly split first

# 3. RESCUE — coexistence maintained against the odds (the tail).
#    A parameter combo sitting right on the edge of a coexistence threshold.
rescue <- s[flag_rescue == TRUE][order(p_retain)]      # lowest P first = most surprising


# the single most surprising survivor (e.g. Classic_Random n7, 1-in-120)
surv   <- retention_tables$rescue[outcome == "rescue"][order(p_retain)][1]
idx    <- readRDS("…/0_Combined_Master/run_index.rds")

# the whole scenario (survivor + the collapsers), then load their trajectories
cell   <- idx[combination == surv$combination & n_species == surv$n_species &
              disturbance_on == surv$disturbance_on]      # add other matching params if needed
traj   <- loadTrajectories(cell)                          # per-timestep rows for that cell
# survivor in bold, collapsers faint — is the survivor at a stable equilibrium or just slow?


# ============================================================================
#  ANOMALOUS RICHNESS LOSS - the MIRROR of the rescue analysis
# ----------------------------------------------------------------------------
#  Same combined, de-duplicated index, but the NETWORKS INCLUDED ARE RPS, Random
#  AND Linear (only Neutral is dropped - it retains everything by design, so
#  "loss when normally retained" is meaningless for it).
#
#  Where the rescue analysis flagged rare SURVIVORS amid collapse, this flags the
#  reverse: within a scenario that NORMALLY RETAINS (p_retain >= loss_norm_thresh,
#  default 0.75), the individual replicates that nonetheless LOST diversity. These
#  are the unexpected collapses a scenario mean would smooth over - the fragile
#  runs revealing a coexistence that is reliable but not guaranteed.
#
#  metric = "richness" (default, as requested): retention is judged on the fraction
#  of SPECIES kept (rich_ret = richness / n_species), so a flagged loss is a genuine
#  drop in species number. metric = "diversity" judges on Hill-1 (exp(Shannon), so
#  a collapse of evenness also counts). Severity is reported as species dropped and
#  the fraction of the metric lost.
# ============================================================================

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

loss_tables <- makeLossTables()
