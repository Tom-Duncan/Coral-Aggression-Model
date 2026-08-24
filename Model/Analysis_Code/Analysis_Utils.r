# Shared helpers for the results-summary scripts (definitions only, no side effects):
# the scenario-ID codebook, Word/CSV table writers, the shared metric set, and
# decomposeCombo(). Sourced by Summary_tables.r, Calculated_Stats.r and others.

# --- shared paths ------------------------------------------------------------
RESULTS_DIR     <- "Simulation_Results"
summary_out_dir <- file.path(RESULTS_DIR, "ScenarioID_Sum_Table")   # tables + ID codebook
calc_out_dir    <- summary_out_dir                                  # (Calculated_Stats.r alias)
registry_path   <- file.path(summary_out_dir, "scenario_ID_registry.csv")

# --- community-level metrics summarised across the scripts -------------------
stat_metrics <- c("species_richness", "shannon", "evenness",
                  "total_cover", "structural_complexity", "turnover")
stat_labels  <- c(species_richness = "Richness", shannon = "Shannon",
                  evenness = "Evenness", total_cover = "Total cover",
                  structural_complexity = "Structural complexity",
                  turnover = "Turnover")

# --- scenario ID codebook ----------------------------------------------------
# A permanent integer ID per run-scenario (keyed on run_id), handed out in
# first-seen order and never reused, stored in `registry`. Keeps IDs stable and
# consistent across every summary table.
assignScenarioIDs <- function(run_ids, scenarios, registry = registry_path) {
  reg <- if (file.exists(registry))
           read.csv(registry, stringsAsFactors = FALSE, colClasses = "character")
         else data.frame(ID = character(), run_id = character(),
                         Scenario = character(), stringsAsFactors = FALSE)
  new <- setdiff(run_ids, reg$run_id)                    # run-scenarios never seen before
  if (length(new) > 0) {
    start <- if (nrow(reg) > 0) max(as.integer(reg$ID)) else 0L
    reg <- rbind(reg, data.frame(
      ID       = as.character(start + seq_along(new)),
      run_id   = new,
      Scenario = scenarios[match(new, run_ids)],
      stringsAsFactors = FALSE))
    dir.create(dirname(registry), showWarnings = FALSE, recursive = TRUE)
    tryCatch(write.csv(reg, registry, row.names = FALSE),
             error = function(e)
               cat("  NOTE: could not update the ID registry (open in Excel?):\n    ", registry, "\n"))
  }
  as.integer(reg$ID[match(run_ids, reg$run_id)])
}

# --- table writers -----------------------------------------------------------
# Write ONE table as an HTML .doc that Word opens as a real, editable table.
writeWordTable <- function(tab, file, title = NULL) {
  cell <- "border:1px solid #444;padding:4px 8px;"
  th <- paste0("<th style='", cell, "background:#e8e8e8;text-align:left;'>",
               names(tab), "</th>", collapse = "")
  body <- apply(tab, 1, function(r)
    paste0("<tr>", paste0("<td style='", cell, "'>", r, "</td>", collapse = ""), "</tr>"))
  html <- c("<html><head><meta charset='utf-8'></head><body>",
            if (!is.null(title)) paste0("<h3 style='font-family:Calibri;'>", title, "</h3>"),
            "<table style='border-collapse:collapse;font-family:Calibri;font-size:11pt;'>",
            paste0("<tr>", th, "</tr>"), body, "</table></body></html>")
  writeLines(html, file)
}

# Write SEVERAL named tables into one .doc, each with its own heading.
writeWordTables <- function(named_tabs, file, title = NULL) {
  cell <- "border:1px solid #444;padding:4px 8px;"
  parts <- c("<html><head><meta charset='utf-8'></head><body>",
             if (!is.null(title)) paste0("<h2 style='font-family:Calibri;'>", title, "</h2>"))
  for (nm in names(named_tabs)) {
    tab <- named_tabs[[nm]]
    th  <- paste0("<th style='", cell, "background:#e8e8e8;text-align:left;'>",
                  names(tab), "</th>", collapse = "")
    body <- apply(tab, 1, function(r)
      paste0("<tr>", paste0("<td style='", cell, "'>", r, "</td>", collapse = ""), "</tr>"))
    parts <- c(parts, paste0("<h3 style='font-family:Calibri;'>", nm, "</h3>"),
               "<table style='border-collapse:collapse;font-family:Calibri;font-size:11pt;'>",
               paste0("<tr>", th, "</tr>"), body, "</table>")
  }
  writeLines(c(parts, "</body></html>"), file)
}

# Write a table as CSV + Word .doc; reports (not stops) if a file is locked open.
writeTablePair <- function(tab, base_name, out_dir = summary_out_dir,
                           suffix = "_scenario_summary") {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  csv_out <- file.path(out_dir, paste0(base_name, suffix, ".csv"))
  doc_out <- file.path(out_dir, paste0(base_name, suffix, ".doc"))
  safe_write <- function(path, fun)
    tryCatch({ fun(path); cat("  wrote:", path, "\n") },
             error = function(e)
               cat("  COULD NOT write:", path,
                   "\n     (is it open in Word/Excel? close it and re-run)\n"))
  safe_write(csv_out, function(p) write.csv(tab, p, row.names = FALSE))
  safe_write(doc_out, function(p) writeWordTable(tab, p, title = base_name))
}

# --- combination decomposition ----------------------------------------------
# Split a combination like "RPS_sizecompImpact" into base model + manipulation,
# plus a broad manip_type (aligned across models, for heatmaps).
decomposeCombo <- function(comb) {
  is_classic <- grepl("^Classic_", comb)
  base  <- ifelse(is_classic, sub("^Classic_", "", comb), sub("_.*$", "", comb))
  manip <- ifelse(is_classic, "none", sub("^[^_]+_", "", comb))
  type  <- ifelse(manip == "none",              "none",
           ifelse(grepl("size|Size", manip),    "size",
           ifelse(grepl("growth|Growth", manip),"growth",
           ifelse(grepl("repro|Repro", manip),  "reproduction", "other"))))
  data.frame(base_model = base, manipulation = manip, manip_type = type,
             stringsAsFactors = FALSE)
}
