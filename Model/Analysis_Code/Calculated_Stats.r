# ==========================================================================
#  Calculated_Stats.r
#  DERIVED / post-hoc summary statistics -- measures from the "Table of
#  measurements" that are CALCULATED rather than logged directly. They are
#  computed from a master results CSV (per-checkpoint metrics + per-species
#  cover Sp1_cover..Sp7_cover + colony counts + the full time series).
#
#  Each measure is a small function of ONE replicate's time series, so adding a
#  new one is a one-liner in the `measures` list. The driver runs them per
#  replicate, then reports mean +/- SD across replicates, per scenario.
#
#  Implemented (from the metrics CSV): Simpson (D), Surviving fraction,
#    Coexistence (species-level), Stability = inverse CV (cover & richness),
#    Repeatability (= the +/-SD), dominance frequency & duration, mortality rate,
#    exclusion rate, spatial turnover, and time-to full cover / first extinction
#    / complete dominance / equilibrium / persistence, plus n small/large.
#  Added for LATER models (work now with caveats):
#    Colonisation rate - uses a logged recruits/births column if present, else a
#      proxy from gross colony gains (counts splits, not just true recruitment);
#    Recovery rate - post-disturbance, per-species cover regrowth slope; NA for
#      no-disturbance runs.
#  Still NOT possible here (need saved colony grid/IDs, not the metrics table):
#    per-colony cover, full colony size distribution, colony-level coexistence.
# ==========================================================================

# calc_out_dir (+ the shared table writers) come from Analysis_Utils.r.
source("Model/Analysis_Code/Analysis_Utils.r")
dir.create(calc_out_dir, showWarnings = FALSE, recursive = TRUE)

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

# --- Usage ------------------------------------------------------------------
# State metrics are the MEAN over the last n_report checkpoints (default 10);
# time-to metrics scan the full run. Point it at any master results CSV.

#Sim 2 - Reef size
makeCalculatedStats("Simulation_Results/2_Robustness_Checks/ReefSize/Master_Robustness_ReefSize_20260727_132335_results.csv")

#Sim 2 - Number of replicates
makeCalculatedStats("Simulation_Results/2_Robustness_Checks/NumReplicates/Master_Robustness_Replicates_20260727_115226_results.csv")


#Sim 2 - Length of sim
makeCalculatedStats("Simulation_Results/2_Robustness_Checks/TimeConvergence(Timesteps)/Master_Robustness_TimeConvergence_20260727_110034_results.csv")

#Sim 3 - Parameter sensitivity
#Sim 3 - disturbance sensitivity
makeCalculatedStats("Simulation_Results/3_Parameter_Sensitivity/Dist_Sens_Test/Master_Sensitivity_DisturbanceRegime_20260727_213926_results.csv")

#Sim 3 - bias check
makeCalculatedStats("Simulation_Results/3_Parameter_Sensitivity/Bias_Sensitivity_Checks/Master_Sensitivity_Bias_20260731_010903_results.csv")

#Sim 3 - Intraspecific competition
makeCalculatedStats("Simulation_Results/3_Parameter_Sensitivity/Intra_Sens_Test/Master_Sensitivity_Intraspecific_20260731_054833_results.csv")

#Sim 3 - reproduction
makeCalculatedStats("Simulation_Results/3_Parameter_Sensitivity/Repro_Sens_Test/Reproduction_Fecundity_20260731_222545_results.csv")

#Sim 3 - reproduction x disturbance regime (mixed conditions)
makeCalculatedStats("Simulation_Results/3_Parameter_Sensitivity/Repro_Dist_sens_Test/Reproduction_DisturbanceRegime_20260801_123052_results.csv")



#Sim 3.1 growth rate
makeCalculatedStats("Simulation_Results/3.1_Extended_Sensitivity/GrowthRate_Sens/Growth_Rate_20260801_230519_results.csv")

#Sim 3.1 Size  competition
makeCalculatedStats("Simulation_Results/3.1_Extended_Sensitivity/SizeComp_Sens/SizeCompetitionStrength_20260802_115434_results.csv")

#Sim 3.1 Background recruitment
makeCalculatedStats("Simulation_Results/3.1_Extended_Sensitivity/BackRecruit_Sens/BackgroundRecruitment_20260802_025433_results.csv")

#Sim 3.1 Maturity Age
makeCalculatedStats("Simulation_Results/3.1_Extended_Sensitivity/MaturityAge_Sens/MaturityAge_20260802_071441_results.csv")

#Sim 5 intial placement
makeCalculatedStats("Simulation_Results/5_InitialPlacement/InitialPlacement_20260802_215722_results.csv")

#Sim 4 - Full model runs
#Sim 4 only RPS CHECK
makeCalculatedStats("Simulation_Results/4_Base_model_exploration/RPS_Model/my_run_20260724_174603_results.csv")


#Sim 4 - Full model runs: both big files COMBINED into one calculated table
makeCalculatedStats(
  c("Simulation_Results/4_Base_model_exploration/Full_Model/Run_20260726_Remaining_Basic_Model_Exploration/AllModels_sweep_20260726_153222_results.csv",
    "Simulation_Results/4_Base_model_exploration/Full_Model/Full_Run_20260724_Classic-RPS-Linear-Neutral_n357_reef50_t1000_r30/my_run_20260724_192134_results.csv"),
  out_name = "Sim4_FullModel_combined")
