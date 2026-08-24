# =============================================================================
#  CONSOLIDATE DATASETS  -  one unified, de-duplication-aware view of every
#  experiment sweep, so all data can be filtered together without accidentally
#  double-counting scenarios that were re-run across separate batches.
# -----------------------------------------------------------------------------
#  This machine has limited RAM, and the pooled per-timestep data is ~2.8 GB, so
#  we DO NOT hold it all in memory. Instead we build:
#
#    (1) run_index  - ONE row per run (unique run_id): every parameter column,
#        the end-state metrics (final checkpoint), provenance (source_file), and
#        duplicate flags. Small, loads instantly, and is what you FILTER on.
#
#    (2) loadTrajectories(runs) - an on-demand loader that reads the full
#        per-timestep rows for ONLY the runs you filtered to, from their source
#        files. Never loads more than one source file at a time.
#
#  Nothing is deleted: duplicates are FLAGGED (is_duplicate / is_duplicate_loose)
#  so you decide per-analysis whether to drop them.
# =============================================================================

suppressMessages(library(data.table))

root    <- "C:/Users/dell/Documents/Research Project 2/Simulation_Results"
out_dir <- file.path(root, "0_Combined_Master")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Which files: every raw run file EXCEPT derived summaries and the
#      robustness checks (non-standard t / reef size / replicate counts) --------
all_raw <- list.files(root, pattern = "results[.]csv$", recursive = TRUE, full.names = TRUE)
all_raw <- all_raw[!grepl("(StatSum|SceSum|ExtraStatSum)", basename(all_raw))]
files   <- all_raw[!grepl("2_Robustness_Checks", all_raw)]      # experiment sweeps only
cat("Consolidating", length(files), "experiment files\n")

relpath <- function(p) sub(paste0("^", root, "/"), "", p)
# Batch date (yyyymmdd_hhmmss) parsed from the filename, so the EARLIEST run of a
# repeated simulation is treated as the canonical original when flagging copies.
fileDate <- function(p) {
  m <- regmatches(basename(p), regexpr("[0-9]{8}_[0-9]{6}", basename(p)))
  if (length(m) == 0) "00000000_000000" else m
}

# Columns that are NOT part of a run's identity (they vary within a run, or are
# provenance / outputs). Everything else is a parameter -> part of the signature.
DYNAMIC <- c("phase", "timestep",
             "species_richness", "shannon", "evenness", "total_cover", "n_colonies",
             "structural_complexity", "turnover", "deaths_since_last", "founders_alive",
             "n_small", "n_medium", "n_large", "mean_colony_size", "median_colony_size",
             "max_colony_size", paste0("Sp", 1:7, "_cover"))
PROVENANCE <- c("run_id", "experiment", "replicate", "source_file",
                "file_date", "n_checkpoints", "max_timestep")

# =============================================================================
#  BUILD THE PER-RUN INDEX
#  Read each file once (fast, ~0.5 GB peak), keep the FINAL-checkpoint row per run
#  (carries all params + end-state metrics), and tag provenance. One row per run.
# =============================================================================
buildRunIndex <- function() {
  parts <- lapply(files, function(f) {
    d <- fread(f)
    fin <- d[phase == "final"]
    # runs with no explicit "final" phase: fall back to the max-timestep row
    if (nrow(fin) == 0) fin <- d[d[, .I[timestep == max(timestep)], by = run_id]$V1]
    nchk <- d[, .(n_checkpoints = .N, max_timestep = max(timestep)), by = run_id]
    fin <- merge(fin, nchk, by = "run_id", all.x = TRUE)
    fin[, `:=`(source_file = relpath(f), file_date = fileDate(f))]
    fin[, (DYNAMIC[DYNAMIC %in% c("phase", "timestep")]) := NULL]  # drop phase/timestep (all "final")
    fin
  })
  idx <- rbindlist(parts, fill = TRUE, use.names = TRUE)

  # --- Duplicate flagging ----------------------------------------------------
  # STRICT: same simulation = identical in every PARAMETER column + seed (the
  # model is deterministic given its seed). LOOSE: same (scenario, seed) only -
  # catches copies even when files differ in which param columns they record.
  sig_cols <- setdiff(names(idx), c(DYNAMIC, PROVENANCE))
  sigOf <- function(cols) {
    m <- as.data.frame(lapply(idx[, ..cols], function(x) { x <- as.character(x); x[is.na(x)] <- "NA"; x }))
    do.call(paste, c(m, sep = ""))
  }
  setorder(idx, file_date, run_id)                 # earliest batch = canonical original
  idx[, sim_signature := sigOf(sig_cols)]
  idx[, sim_id        := .GRP,               by = sim_signature]
  idx[, dup_n         := .N,                 by = sim_signature]
  idx[, dup_rank      := seq_len(.N),        by = sim_signature]   # 1 = original, >1 = copy
  idx[, is_duplicate  := dup_rank > 1L]
  loose <- if (all(c("scenario", "seed") %in% names(idx))) sigOf(c("scenario", "seed")) else idx$sim_signature
  idx[, is_duplicate_loose := duplicated(loose)]
  idx[, sim_signature := NULL]                      # drop the long helper string

  idx[]
}

run_index <- buildRunIndex()

# --- Save the index (rds for speed, csv for portability) ---------------------
saveRDS(run_index, file.path(out_dir, "run_index.rds"))
fwrite(run_index,  file.path(out_dir, "run_index.csv"))

# --- Report ------------------------------------------------------------------
cat("\nrun_index:", nrow(run_index), "runs,", ncol(run_index), "columns\n")
cat("strict duplicates flagged (identical params + seed):", sum(run_index$is_duplicate), "\n")
cat("loose  duplicates flagged (scenario + seed):        ", sum(run_index$is_duplicate_loose), "\n")
cat("unique simulations (strict):", uniqueN(run_index$sim_id), "\n")
cat("\nruns per source file (dup = strict copies):\n")
print(run_index[, .(runs = .N, dup = sum(is_duplicate)), by = source_file][order(-runs)])
cat("\nSaved: ", file.path(out_dir, "run_index.rds"), " (+ .csv)\n", sep = "")


# =============================================================================
#  HELPERS  (source this file, then use these interactively)
# =============================================================================

# Filter the index to the runs you want. `subset` is ONE logical expression over the
# index columns (combine terms with &), e.g.
#     pickRuns(subset = combination == "Classic_RPS" & n_species == 5)
# Filter on ACTUAL PARAMETER COLUMNS, not just `scenario`: the same scenario string
# can span sweeps that differ in an un-encoded focal parameter (background_chance,
# maturity_age, ...). `drop_dups` removes strictly-identical copies (default) so a
# pooled scenario is not over-counted.
pickRuns <- function(index = run_index, subset = TRUE, drop_dups = TRUE, loose = FALSE) {
  keep <- eval(substitute(subset), envir = index, enclos = parent.frame())
  d <- index[keep]
  if (drop_dups) d <- d[if (loose) !is_duplicate_loose else !is_duplicate]
  d[]
}

# Load the FULL per-timestep trajectories for a set of runs (a subset of the
# index). Reads only each needed source file, once, keeping just those run_ids.
# `cols` optionally restricts which columns to return (run_id is always kept).
loadTrajectories <- function(runs, cols = NULL) {
  stopifnot(all(c("run_id", "source_file") %in% names(runs)))
  bysrc <- tapply(runs$run_id, runs$source_file, unique)   # run_ids grouped by source file
  out <- rbindlist(lapply(names(bysrc), function(sf) {
    d   <- fread(file.path(root, sf))
    d   <- d[run_id %in% bysrc[[sf]]]
    d[, source_file := sf]
    if (!is.null(cols)) d <- d[, intersect(unique(c("run_id", "source_file", cols)), names(d)), with = FALSE]
    d
  }), fill = TRUE, use.names = TRUE)
  out[]
}

# --- Usage examples (run after sourcing) -------------------------------------
# runs <- pickRuns(subset = combination == "Classic_RPS" & n_species == 5 & disturbance_on == FALSE)
# traj <- loadTrajectories(runs)                       # per-timestep rows for those runs
# plot(traj[run_id == runs$run_id[1], .(timestep, species_richness)], type = "l")
#
# # how many replicates does each scenario really have, after de-duplication?
# run_index[is_duplicate == FALSE, .N, by = scenario][order(-N)]
#
# # NOTE: filter on real PARAMETER columns, not `scenario` alone - the same scenario
# # string spans sweeps that differ in an un-encoded focal parameter, e.g.:
# run_index[scenario == run_index$scenario[1], .N, by = .(background_chance, maturity_age, size_beta_max, repro_base)]
