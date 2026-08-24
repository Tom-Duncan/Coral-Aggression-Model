# One de-duplication-aware index over every experiment sweep, so data can be filtered
# together without double-counting re-run scenarios. The full per-timestep data is
# large, so this builds a small per-run index (one final-checkpoint row per run_id,
# with params + end-state metrics + duplicate flags) plus loadTrajectories() to read
# full trajectories on demand. Duplicates are flagged, never deleted.

suppressMessages(library(data.table))

root    <- "Simulation_Results"
out_dir <- file.path(root, "0_Combined_Master")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# every raw run file except derived summaries and the robustness checks
all_raw <- list.files(root, pattern = "results[.]csv$", recursive = TRUE, full.names = TRUE)
all_raw <- all_raw[!grepl("(StatSum|SceSum|ExtraStatSum)", basename(all_raw))]
files   <- all_raw[!grepl("2_Robustness_Checks", all_raw)]
cat("Consolidating", length(files), "experiment files\n")

relpath <- function(p) sub(paste0("^", root, "/"), "", p)
# batch date from the filename, so the earliest run is the canonical original
fileDate <- function(p) {
  m <- regmatches(basename(p), regexpr("[0-9]{8}_[0-9]{6}", basename(p)))
  if (length(m) == 0) "00000000_000000" else m
}

# columns that are NOT part of a run's identity (dynamic / provenance / outputs);
# everything else is a parameter -> part of the signature
DYNAMIC <- c("phase", "timestep",
             "species_richness", "shannon", "evenness", "total_cover", "n_colonies",
             "structural_complexity", "turnover", "deaths_since_last", "founders_alive",
             "n_small", "n_medium", "n_large", "mean_colony_size", "median_colony_size",
             "max_colony_size", paste0("Sp", 1:7, "_cover"))
PROVENANCE <- c("run_id", "experiment", "replicate", "source_file",
                "file_date", "n_checkpoints", "max_timestep")

# Per-run index: one final-checkpoint row per run (params + end-state + provenance).
buildRunIndex <- function() {
  parts <- lapply(files, function(f) {
    d <- fread(f)
    fin <- d[phase == "final"]
    if (nrow(fin) == 0) fin <- d[d[, .I[timestep == max(timestep)], by = run_id]$V1]   # no "final": max timestep
    nchk <- d[, .(n_checkpoints = .N, max_timestep = max(timestep)), by = run_id]
    fin <- merge(fin, nchk, by = "run_id", all.x = TRUE)
    fin[, `:=`(source_file = relpath(f), file_date = fileDate(f))]
    fin[, (DYNAMIC[DYNAMIC %in% c("phase", "timestep")]) := NULL]
    fin
  })
  idx <- rbindlist(parts, fill = TRUE, use.names = TRUE)

  # duplicate flags: STRICT = identical params + seed (model is deterministic given
  # seed); LOOSE = same (scenario, seed) only.
  sig_cols <- setdiff(names(idx), c(DYNAMIC, PROVENANCE))
  sigOf <- function(cols) {
    m <- as.data.frame(lapply(idx[, ..cols], function(x) { x <- as.character(x); x[is.na(x)] <- "NA"; x }))
    do.call(paste, c(m, sep = ""))
  }
  setorder(idx, file_date, run_id)                 # earliest batch = canonical
  idx[, sim_signature := sigOf(sig_cols)]
  idx[, sim_id        := .GRP,               by = sim_signature]
  idx[, dup_n         := .N,                 by = sim_signature]
  idx[, dup_rank      := seq_len(.N),        by = sim_signature]   # 1 = original, >1 = copy
  idx[, is_duplicate  := dup_rank > 1L]
  loose <- if (all(c("scenario", "seed") %in% names(idx))) sigOf(c("scenario", "seed")) else idx$sim_signature
  idx[, is_duplicate_loose := duplicated(loose)]
  idx[, sim_signature := NULL]

  idx[]
}

run_index <- buildRunIndex()

saveRDS(run_index, file.path(out_dir, "run_index.rds"))
fwrite(run_index,  file.path(out_dir, "run_index.csv"))

cat("\nrun_index:", nrow(run_index), "runs,", ncol(run_index), "columns\n")
cat("strict duplicates flagged (identical params + seed):", sum(run_index$is_duplicate), "\n")
cat("loose  duplicates flagged (scenario + seed):        ", sum(run_index$is_duplicate_loose), "\n")
cat("unique simulations (strict):", uniqueN(run_index$sim_id), "\n")
cat("\nruns per source file (dup = strict copies):\n")
print(run_index[, .(runs = .N, dup = sum(is_duplicate)), by = source_file][order(-runs)])
cat("\nSaved: ", file.path(out_dir, "run_index.rds"), " (+ .csv)\n", sep = "")


# ---- Helpers (source this file, then use interactively) ---------------------

# Filter the index. `subset` is one logical expression over the columns; filter on
# real parameter columns (not just `scenario`). drop_dups removes strict copies.
pickRuns <- function(index = run_index, subset = TRUE, drop_dups = TRUE, loose = FALSE) {
  keep <- eval(substitute(subset), envir = index, enclos = parent.frame())
  d <- index[keep]
  if (drop_dups) d <- d[if (loose) !is_duplicate_loose else !is_duplicate]
  d[]
}

# Load full per-timestep trajectories for a set of runs (reads each source file once).
loadTrajectories <- function(runs, cols = NULL) {
  stopifnot(all(c("run_id", "source_file") %in% names(runs)))
  bysrc <- tapply(runs$run_id, runs$source_file, unique)
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
# traj <- loadTrajectories(runs)
# plot(traj[run_id == runs$run_id[1], .(timestep, species_richness)], type = "l")
# run_index[is_duplicate == FALSE, .N, by = scenario][order(-N)]
