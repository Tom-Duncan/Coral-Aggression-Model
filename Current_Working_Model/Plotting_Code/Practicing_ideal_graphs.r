# ==========================================================================
#  Practicing ideal graphs
#  Scratch/plotting script for the reef-model simulation results.
#  Run top-to-bottom, or source it once to load the functions in Section 3
#  and then call them interactively.
#
#  Contents
#    0. Config          - data file paths and the output folder (edit here)
#    1. Single metric   - one metric over time for one scenario
#                         (faint replicate lines + bold mean)
#    2. Exploratory     - Simon's multi-metric panels written to PNG
#    3. Convergence     - reusable time-convergence tools for robustness checks
#         3a. plotConvergence     - one scenario, mean +/-1 SD
#         3b. plotAllConvergence  - many scenarios on one plot
#         3c. plotMetricPanels    - all metrics for one model, faceted
# ==========================================================================


# ==========================================================================
# 0. Config -- all file locations in one place
# ==========================================================================
# RPS optimisation run (used by Sections 1 and 2)
rps_path <- "C:/Users/dell/Documents/Research Project 2/Simulation_Results/1_Optim_vs_Unop/Initial_Testing(RPS)/my_run_20260724_174603_results(optimised).csv"

# Timestep robustness run (used by Section 3)
robust_path <- "C:/Users/dell/Documents/Research Project 2/Simulation_Results/2_Robustness_Checks/TimeConvergence(Timesteps)/Master_Robustness_TimeConvergence_20260727_110034_results.csv"

# Where the PNGs from Section 2 are saved (created if it doesn't exist)
out_dir <- "C:/Users/dell/Documents/Research Project 2/Modelling/Current_Working_Model/Plotting_Code/Plots"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)


# ==========================================================================
# 1. Single metric over time for ONE scenario
#    Every replicate run is a faint translucent line; the mean across
#    replicates is bold on top.
# ==========================================================================
rps <- read.csv(rps_path, stringsAsFactors = FALSE)

# --- Choose what to plot ----------------------------------------------------
metric <- "evenness"          # any numeric column: evenness, species_richness, shannon, total_cover
model  <- "Classic_RPS"       # the ONE model to plot: Classic_RPS / Classic_Linear / Classic_Neutral

# --- Keep one clean scenario (its replicate runs get averaged) --------------
#  (add another condition, e.g. rps$individuals == 2, to pin it down further)
sub <- rps[rps$combination == model &
           rps$disturbance_on == FALSE &
           rps$n_species == 3, ]

# --- Draw the single plot ---------------------------------------------------
ts    <- sort(unique(sub$timestep))
faint <- adjustcolor("steelblue", alpha.f = 0.20)   # translucent = faint replicate lines

plot(NA, xlim = range(ts), ylim = range(sub[[metric]], na.rm = TRUE),
     xlab = "Timestep", ylab = metric,
     main = paste(metric, "over time -", model))
grid(col = "grey90", lty = 1)

# faint line for each replicate run (a run = one unique run_id trajectory)
for (r in unique(sub$run_id)) {
  one <- sub[sub$run_id == r, ]
  one <- one[order(one$timestep), ]
  lines(one$timestep, one[[metric]], col = faint, lwd = 1)
}

# bold mean across replicates, drawn last so it sits on top
means <- tapply(sub[[metric]], sub$timestep, mean, na.rm = TRUE)[as.character(ts)]
lines(ts, means, col = "navy", lwd = 3)

legend("topright", c("replicate runs", "mean"),
       col = c(faint, "navy"), lwd = c(1, 3), bty = "n")


# ==========================================================================
# 2. Exploratory multi-metric panels (Simon)
#    Scatter/quantile panels for every metric, saved as PNGs to out_dir.
# ==========================================================================
alldata <- read.csv(rps_path, stringsAsFactors = FALSE)

# one colour per initial number of species...
sp_col <- heat.colors(length(unique(alldata$n_species)))
names(sp_col) <- as.character(unique(alldata$n_species))

# ...and one colour per disturbance state
dcol <- c("red", "black")
names(dcol) <- c("TRUE", "FALSE")

# --- Every metric vs time, coloured by disturbance and by species count -----
png(file.path(out_dir, "allmetrics_time_full.png"), width = 1200, height = 1600, pointsize = 17)
par(mfrow = c(6, 2))
for (metrics in c("species_richness", "structural_complexity", "turnover", "shannon", "total_cover", "evenness")) {
    plot(alldata$timestep, alldata[[metrics]], col = adjustcolor(dcol[as.character(alldata$disturbance_on)], .4), pch = 20, cex = 2, main = metrics, xlab = "time", ylab = metrics)
    if (metrics == "evenness") legend("topright", legend = c("yes", "no"), col = c("red", "black"), bg = "white", pch = 20, title = "disturbance")
    plot(alldata$timestep, alldata[[metrics]], col = adjustcolor(sp_col[as.character(alldata$n_species)], .4), pch = 20, cex = 1.8, main = metrics, xlab = "time", ylab = metrics)
    # adjustcolor here creates transparency.
    # note: alldata$timestep and alldata[["timestep"]] are the same; the [[ ]] form
    # lets you pass a changing column name in a variable, which the $ form cannot.
    if (metrics == "evenness") legend("topright", legend = names(sp_col), col = sp_col, bg = "white", pch = 20, title = "initial number of species")
}
dev.off()

# --- Median + inter-quartile lines per group, one row per metric ------------
png(file.path(out_dir, "allmetrics_time_disturbance_quantiles.png"), width = 1200, height = 1600, pointsize = 17)
par(mfrow = c(6, 2))
for (metrics in c("species_richness", "structural_complexity", "turnover", "shannon", "total_cover", "evenness")) {
    print(metrics)
    # coloured by disturbance state
    plot(1, 1, type = "n", xlim = range(alldata$timestep), ylim = range(alldata[[metrics]], na.rm = TRUE))
    for (ds in unique(alldata$disturbance_on)) {
        tmp <- alldata[alldata$disturbance_on == ds, ]
        sum_met <- do.call(cbind, tapply(tmp[[metrics]], tmp$timestep, quantile, na.rm = TRUE))
        lines(sort(unique(tmp$timestep)), sum_met[3, ], lwd = 2, col = dcol[as.character(ds)])           # median
        lines(sort(unique(tmp$timestep)), sum_met[2, ], lwd = 2, lty = 3, col = dcol[as.character(ds)])  # lower quartile
        lines(sort(unique(tmp$timestep)), sum_met[4, ], lwd = 2, lty = 3, col = dcol[as.character(ds)])  # upper quartile
    }
    # coloured by initial number of species
    plot(1, 1, type = "n", xlim = range(alldata$timestep), ylim = range(alldata[[metrics]], na.rm = TRUE))
    for (ns in unique(alldata$n_species)) {
        tmp <- alldata[alldata$n_species == ns, ]
        sum_met <- do.call(cbind, tapply(tmp[[metrics]], tmp$timestep, quantile, na.rm = TRUE))
        lines(sort(unique(tmp$timestep)), sum_met[3, ], lwd = 2, col = sp_col[as.character(ns)])
        lines(sort(unique(tmp$timestep)), sum_met[2, ], lwd = 2, lty = 3, col = sp_col[as.character(ns)])
        lines(sort(unique(tmp$timestep)), sum_met[4, ], lwd = 2, lty = 3, col = sp_col[as.character(ns)])
    }
    if (metrics == "evenness") legend("topright", legend = c("yes", "no"), col = c("red", "black"), bg = "white", pch = 20, title = "disturbance")
}
dev.off()

# --- End-state interaction: metric at the final timestep, by group ----------
final <- alldata[alldata$timestep == max(alldata$timestep), ]

png(file.path(out_dir, "interaction.png"), width = 1200, height = 1200, pointsize = 17)
par(mfrow = c(2, 3))
for (metrics in c("species_richness", "structural_complexity", "turnover", "shannon", "total_cover", "evenness")) {
    boxplot(final[[metrics]] ~ final$n_species + final$disturbance_on, col = adjustcolor(sp_col, .4), pch = 20, cex = 1.8, main = metrics, xlab = "disturbance", ylab = metrics, range = 0, xaxt = "n", at = c(1, 2, 3.5, 4.5))
    if (metrics == "evenness") legend("bottomright", legend = names(sp_col), fill = adjustcolor(sp_col, .4), title = "initial number of species", bty = "n")
    axis(1, label = c("no", "yes"), at = c(1.5, 4))
}
dev.off()


# ==========================================================================
# 3. Time-convergence tools (timestep robustness checks)
# ==========================================================================
res <- read.csv(robust_path, stringsAsFactors = FALSE)

# --------------------------------------------------------------------------
# 3a. plotConvergence -- ONE scenario (model, species count, disturbance).
#     Averages the chosen metric across all replicates at each timestep and
#     draws the mean with a +/-1 SD band. Returns the mean/sd/se table.
#     show_lines = TRUE overlays each replicate's raw trajectory faintly in the
#     background (like the Sim 4 RPS plots), with the bold mean on top.
# --------------------------------------------------------------------------
plotConvergence <- function(results, model, nsp, dist_on,
                            metric = "species_richness", show_band = TRUE,
                            show_lines = FALSE) {

  # keep only this scenario's rows (all its replicates)
  keep <- results$combination == model &
          results$n_species    == nsp &
          results$disturbance_on == dist_on
  sub <- results[keep, , drop = FALSE]
  if (nrow(sub) == 0) stop("No rows for that scenario.")

  # mean and spread across replicates at each timestep
  ts     <- sort(unique(sub$timestep))
  mean_y <- tapply(sub[[metric]], sub$timestep, mean, na.rm = TRUE)[as.character(ts)]
  sd_y   <- tapply(sub[[metric]], sub$timestep, sd,   na.rm = TRUE)[as.character(ts)]
  n_rep  <- length(unique(sub$replicate))

  # y-range spans the faint raw replicate lines when shown, else the +/-1 SD band
  ylim <- if (show_lines) range(sub[[metric]], na.rm = TRUE)
          else            range(c(mean_y - sd_y, mean_y + sd_y), na.rm = TRUE)
  plot(ts, mean_y, type = "n", ylim = ylim,
       xlab = "Timestep", ylab = paste("Mean", metric),
       main = sprintf("%s  |  %d species  |  disturbance %s\n(mean of %d replicates)",
                      model, nsp, if (isTRUE(dist_on)) "ON" else "OFF", n_rep))
  # faint individual-replicate trajectories drawn FIRST (background), bold mean on top.
  # Replicates are keyed by run_id if present, else by the replicate number.
  if (show_lines) {
    rep_key <- if ("run_id" %in% names(sub)) "run_id" else "replicate"
    faint   <- adjustcolor("steelblue", alpha.f = 0.12)
    for (r in unique(sub[[rep_key]])) {
      one <- sub[sub[[rep_key]] == r, c("timestep", metric)]
      one <- one[order(one$timestep), ]
      lines(one$timestep, one[[metric]], col = faint, lwd = 1)
    }
  }
  if (show_band)                                     # +/-1 SD across replicates
    polygon(c(ts, rev(ts)), c(mean_y - sd_y, rev(mean_y + sd_y)),
            col = adjustcolor("steelblue", 0.20), border = NA)
  lines(ts, mean_y, col = "steelblue", lwd = 2.5)
  abline(v = 1000, lty = 2, col = "grey40")          # standard run length
  mtext("standard length", side = 1, at = 1000, line = 2.5, col = "grey40", cex = 0.8)

  invisible(data.frame(timestep = ts, mean = mean_y, sd = sd_y,
                       se = sd_y / sqrt(n_rep), row.names = NULL))
}

# (examples are under "Simulation 2, length of sim" below)


# --------------------------------------------------------------------------
# 3b. plotAllConvergence -- MANY groups on one plot, each a coloured line with
#     a +/-1 SD band. models/species/dist are optional filters (NULL = all).
#     `by` = the columns that define one line/colour. Default groups by
#     scenario (matrix x species x disturbance); set e.g. by = "reef_x" to
#     compare reef sizes (filter to one matrix/species/disturbance first).
# --------------------------------------------------------------------------
plotAllConvergence <- function(results, metric = "species_richness",
                               models = NULL, species = NULL, dist = NULL,
                               show_band = TRUE, show_lines = FALSE, show_legend = TRUE,
                               main = NULL,
                               by = c("combination", "n_species", "disturbance_on")) {

  if (!is.null(models))  results <- results[results$combination   %in% models,  ]
  if (!is.null(species)) results <- results[results$n_species     %in% species, ]
  if (!is.null(dist))    results <- results[results$disturbance_on %in% dist,    ]

  # one coloured line per unique combination of the `by` columns
  scen <- unique(results[, by, drop = FALSE])
  scen <- scen[do.call(order, scen[by]), , drop = FALSE]
  n    <- nrow(scen)
  cols <- hcl.colors(n, "Dark 3")            # one distinct colour per group

  ts    <- sort(unique(results$timestep))
  # row mask for each group, reused for the mean/SD and (optionally) the raw lines
  grp_keep <- lapply(seq_len(n), function(i) {
    keep <- rep(TRUE, nrow(results))
    for (b in by) keep <- keep & results[[b]] == scen[[b]][i]
    keep
  })
  stats <- lapply(seq_len(n), function(i) {
    sub <- results[grp_keep[[i]], ]
    list(mean = tapply(sub[[metric]], sub$timestep, mean, na.rm = TRUE)[as.character(ts)],
         sd   = tapply(sub[[metric]], sub$timestep, sd,   na.rm = TRUE)[as.character(ts)])
  })
  means <- lapply(stats, `[[`, "mean")
  sds   <- lapply(stats, `[[`, "sd")

  # y-range spans whatever is drawn: the faint raw replicate lines (widest) when
  # shown, else the +/-1 SD bands, else just the mean lines
  ylim <- if (show_lines)
            range(unlist(lapply(grp_keep, function(k) results[[metric]][k])), na.rm = TRUE)
          else if (show_band)
            range(unlist(Map(function(m, s) c(m - s, m + s), means, sds)), na.rm = TRUE)
          else
            range(unlist(means), na.rm = TRUE)
  # wide right margin only needed when this panel carries its own legend
  right_mar <- if (show_legend) 11 else 1.5
  op <- par(mar = c(4.5, 4.5, 3, right_mar), xpd = NA); on.exit(par(op))
  plot(NA, xlim = range(ts), ylim = ylim,
       xlab = "Timestep", ylab = paste("Mean", metric),
       main = if (is.null(main)) paste("Time convergence -", metric,
                                       "(mean of replicates)") else main)
  # faint individual-replicate trajectories per group, drawn FIRST (background),
  # matching the Sim 4 RPS style; the bold mean is drawn on top further below.
  # Replicates are keyed by run_id if present, else by the replicate number.
  if (show_lines) for (i in seq_len(n)) {
    sub     <- results[grp_keep[[i]], ]
    rep_key <- if ("run_id" %in% names(sub)) "run_id" else "replicate"
    faint   <- adjustcolor(cols[i], alpha.f = 0.10)
    for (r in unique(sub[[rep_key]])) {
      one <- sub[sub[[rep_key]] == r, c("timestep", metric)]
      one <- one[order(one$timestep), ]
      lines(one$timestep, one[[metric]], col = faint, lwd = 1)
    }
  }

  # dashed marker for standard run length; segments (not abline) so it stays
  # inside the plot region even though clipping is off (xpd = NA) for the legend
  segments(1000, ylim[1], 1000, ylim[2], lty = 2, col = "grey60")

  # +/-1 SD band per scenario, drawn first so the mean lines sit on top
  if (show_band) for (i in seq_len(n)) {
    m <- means[[i]]; s <- sds[[i]]
    ok <- !is.na(m) & !is.na(s)
    if (any(ok))
      polygon(c(ts[ok], rev(ts[ok])), c((m - s)[ok], rev((m + s)[ok])),
              col = adjustcolor(cols[i], 0.15), border = NA)
  }

  for (i in seq_len(n)) lines(ts, means[[i]], col = cols[i], lwd = 2)

  # readable label per group, built from whatever columns are in `by`
  prettify <- function(col, val) {
    if      (col == "n_species")      paste0("n", val)
    else if (col == "disturbance_on") ifelse(as.logical(val), "on", "off")
    else if (col == "reef_x")         paste0("reef ", val, "x", val)
    else                              as.character(val)
  }
  labs <- apply(scen, 1, function(row) paste(mapply(prettify, by, row), collapse = " | "))
  leg_title <- if (identical(by, c("combination", "n_species", "disturbance_on")))
                 "scenario" else paste(by, collapse = " / ")
  if (show_legend)
    legend(par("usr")[2], par("usr")[4],
           legend = c(labs, "t = 1000 (standard length)"),
           col = c(cols, "grey60"), lwd = 2, lty = c(rep(1, n), 2),
           cex = 0.75, bty = "n", xjust = 0, title = leg_title)

  # return legend ingredients so a caller can draw one shared key for many panels
  invisible(list(scen = scen, labels = labs, cols = cols))
}

# (examples are under "Simulation 2, length of sim" below)


# --------------------------------------------------------------------------
# 3c. plotMetricPanels -- ALL metrics for one model on a single figure: one
#     panel per metric, each titled with just its metric, plus a single shared
#     key and one overall title in the reserved strip at the top. Every scenario
#     (species x disturbance) of the chosen model is a coloured line.
#     Needs plotAllConvergence (3b) to be defined.
#       results : the results data frame (e.g. res)
#       models  : one model name, or a vector of names, to show in every panel
#       metrics : which metric columns to panel (defaults to the usual six)
#       title   : overall figure title above the shared key
# --------------------------------------------------------------------------
plotMetricPanels <- function(results, models,
                             metrics = c("species_richness", "structural_complexity",
                                         "turnover", "shannon", "total_cover", "evenness"),
                             show_lines = FALSE,
                             title = "Time convergence (mean of replicates)") {

  # shared-legend ingredients built straight from the data (same recipe the
  # panel function uses), so the key never depends on a function return value
  sc <- unique(results[results$combination %in% models,
                       c("combination", "n_species", "disturbance_on")])
  sc <- sc[order(sc$combination, sc$n_species, sc$disturbance_on), ]
  sc_labs <- sprintf("%s | n%d | %s", sc$combination, sc$n_species,
                     ifelse(sc$disturbance_on, "on", "off"))
  sc_cols <- hcl.colors(nrow(sc), "Dark 3")

  # panel grid + top strip for the title/key; restore the user's par on exit
  op <- par(mfrow = c(3, 2), oma = c(0, 0, 11, 0)); on.exit(par(op))
  for (m in metrics)
    plotAllConvergence(results, show_band = FALSE, show_lines = show_lines, show_legend = FALSE,
                       metric = m, models = models, main = m)

  # one overall title across the reserved top strip (outer margin, robust)
  mtext(title, outer = TRUE, line = 7.5, font = 2, cex = 1.4)

  # overlay a whole-device region JUST for the single shared key, below the title
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot(0, 0, type = "n", xlim = c(0, 1), ylim = c(0, 1), axes = FALSE, xlab = "", ylab = "")
  legend(0.5, 0.965, xjust = 0.5, yjust = 1,
         ncol = 3,               # 2 rows keeps text big without over-wide labels
         legend = c(sc_labs, "t = 1000 (standard length)"),
         col = c(sc_cols, "grey60"), lwd = 2.5,
         lty = c(rep(1, length(sc_labs)), 2),
         seg.len = 1,            # shorter line samples -> less "stretched" key
         cex = 1.05, x.intersp = 0.7, text.width = 0.22, bty = "n")
  invisible(sc)
}

# (examples are under "Simulation 2, length of sim" below)


#Simulation 2, lenght of sim (number of timesteps) -------------------------
# The timestep-length robustness run (up to t = 3000). Uses the Section 3 tools
# (plotConvergence / plotAllConvergence / plotMetricPanels); the grey dashed line
# marks t = 1000, the "standard length".
res <- read.csv("C:/Users/dell/Documents/Research Project 2/Simulation_Results/2_Robustness_Checks/TimeConvergence(Timesteps)/Master_Robustness_TimeConvergence_20260727_110034_results.csv",
                stringsAsFactors = FALSE)

# One scenario at a time (mean +/-1 SD, "standard length" marker at t = 1000):
plotConvergence(res, "Classic_Neutral", 7, FALSE)                         # richness
plotConvergence(res, "Classic_Linear",  7, TRUE, metric = "total_cover")  # other metric/scenario
plotConvergence(res, "Classic_RPS", 7, FALSE, show_lines = TRUE)          # faint replicate lines behind the mean

# Many scenarios on one plot:
plotAllConvergence(res)                                                 # every scenario
plotAllConvergence(res, models = c("Classic_RPS", "Classic_Random"))    # just two models
plotAllConvergence(res, metric = "total_cover")                         # a different metric
plotAllConvergence(res, models = "Classic_Neutral", show_band = FALSE)  # no SD band
plotAllConvergence(res, models = "Classic_RPS", show_band = FALSE,
                   show_lines = TRUE)                                   # faint replicate lines behind the means

# All metrics for one model, faceted:
plotMetricPanels(res, "Classic_Neutral")
plotMetricPanels(res, "Classic_Random")
plotMetricPanels(res, "Classic_Linear")
plotMetricPanels(res, "Classic_RPS")
plotMetricPanels(res, c("Classic_RPS", "Classic_Random"))                     # compare two models
plotMetricPanels(res, "Classic_Neutral", metrics = c("shannon", "evenness"))  # fewer metrics
plotMetricPanels(res, "Classic_RPS", show_lines = TRUE)                        # every replicate faint in the background, bold mean on top





#Sim 2, reef size -----------------------------------------------------------
reef <- read.csv("C:/Users/dell/Documents/Research Project 2/Simulation_Results/2_Robustness_Checks/ReefSize/Master_Robustness_ReefSize_20260727_132335_results.csv",
                 stringsAsFactors = FALSE)

# --------------------------------------------------------------------------
# plotReefComparison -- one model + one species + one metric on a single plot,
# MEANS ONLY (no SD band). Two visual channels:
#   colour    = reef size
#   line type = disturbance (solid = on, dashed = off)
# Call once per species for separate graphs.
# --------------------------------------------------------------------------
plotReefComparison <- function(results, model, nsp, metric = "species_richness",
                               main = NULL, show_legend = TRUE) {
  d <- results[results$combination == model & results$n_species == nsp, ]

  reefs <- sort(unique(d$reef_x))
  rcol  <- setNames(hcl.colors(length(reefs), "Dark 3"), reefs)
  ts    <- sort(unique(d$timestep))

  # mean-over-replicates series for every reef x disturbance combo
  combos <- expand.grid(reef = reefs, dist = c(TRUE, FALSE), stringsAsFactors = FALSE)
  series <- lapply(seq_len(nrow(combos)), function(i) {
    k <- d$reef_x == combos$reef[i] & d$disturbance_on == combos$dist[i]
    if (!any(k)) NULL
    else tapply(d[[metric]][k], d$timestep[k], mean, na.rm = TRUE)[as.character(ts)]
  })
  ylim <- range(unlist(series), na.rm = TRUE)

  # wide right margin only when this panel carries its own legend
  right_mar <- if (show_legend) 9 else 1.5
  op <- par(mar = c(4.5, 4.5, 3, right_mar), xpd = NA); on.exit(par(op))
  plot(NA, xlim = range(ts), ylim = ylim, xlab = "Timestep", ylab = paste("Mean", metric),
       main = if (is.null(main))
                paste0(sub("^Classic_", "", model), " | n", nsp, " - ", metric, " vs reef size")
              else main)
  segments(1000, ylim[1], 1000, ylim[2], lty = 3, col = "grey70")   # standard run length

  for (i in seq_len(nrow(combos))) {
    y <- series[[i]]; if (is.null(y)) next
    lines(ts, y, col = rcol[as.character(combos$reef[i])],
          lty = if (combos$dist[i]) 1 else 2, lwd = 2.2)
  }

  if (show_legend) {   # two small legends stacked in the right margin
    ux <- par("usr")[2]; yhi <- par("usr")[4]; yr <- yhi - par("usr")[3]
    legend(ux, yhi,             title = "reef size",   xjust = 0, bty = "n", cex = 0.8,
           legend = paste0(reefs, "x", reefs), col = rcol, lwd = 3)
    legend(ux, yhi - 0.34 * yr, title = "disturbance", xjust = 0, bty = "n", cex = 0.8,
           legend = c("on", "off"), lty = c(1, 2), lwd = 2.2)
  }
  invisible(NULL)
}

# --------------------------------------------------------------------------
# plotReefPanels -- ALL metrics for one model + species on a single figure:
# one panel per metric (each a reef-size comparison), with a single shared key
# and one overall title at the top. Call once per species -> separate figures.
# --------------------------------------------------------------------------
plotReefPanels <- function(results, model, nsp,
                           metrics = c("species_richness", "structural_complexity",
                                       "turnover", "shannon", "total_cover", "evenness"),
                           title = NULL) {
  d     <- results[results$combination == model & results$n_species == nsp, ]
  reefs <- sort(unique(d$reef_x))
  rcol  <- setNames(hcl.colors(length(reefs), "Dark 3"), reefs)

  op <- par(mfrow = c(3, 2), oma = c(0, 0, 9, 0)); on.exit(par(op))   # panels + top strip
  for (m in metrics)
    plotReefComparison(results, model, nsp, metric = m, main = m, show_legend = FALSE)

  ttl <- if (is.null(title))
           paste0(sub("^Classic_", "", model), " | n", nsp, " - metrics vs reef size")
         else title
  mtext(ttl, outer = TRUE, line = 6, font = 2, cex = 1.4)

  # one shared key across the top strip: reef colours + disturbance line type
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot(0, 0, type = "n", xlim = c(0, 1), ylim = c(0, 1), axes = FALSE, xlab = "", ylab = "")
  legend(0.5, 0.965, xjust = 0.5, yjust = 1, ncol = length(reefs) + 2, bty = "n",
         cex = 1.0, seg.len = 1.6, x.intersp = 0.7,
         legend = c(paste0("reef ", reefs, "x", reefs), "dist on", "dist off"),
         col    = c(rcol, "grey30", "grey30"),
         lwd    = 2.5, lty = c(rep(1, length(reefs)), 1, 2))
  invisible(NULL)
}

# All metrics for Classic RPS, one figure per species (colour = reef size,
# solid = disturbance on, dashed = off). Swap the model to do another.
plotReefPanels(reef, "Classic_RPS", 3)
plotReefPanels(reef, "Classic_RPS", 7)
plotReefPanels(reef, "Classic_Neutral", 3)
plotReefPanels(reef, "Classic_Neutral", 7)
plotReefPanels(reef, "Classic_Linear", 3)
plotReefPanels(reef, "Classic_Linear", 7)
plotReefPanels(reef, "Classic_Random", 3)
plotReefPanels(reef, "Classic_Random", 7)


#Sim 2, number of replicates ------------------------------------------------
# Same 16 scenarios as the timestep check, but run 100 replicates each. The
# question here is "how many replicates are enough?", so instead of a time
# series we show how each metric's estimate settles as replicates are added.
reps_data <- read.csv("C:/Users/dell/Documents/Research Project 2/Simulation_Results/2_Robustness_Checks/NumReplicates/Master_Robustness_Replicates_20260727_115226_results.csv",
                      stringsAsFactors = FALSE)

# --------------------------------------------------------------------------
# plotReplicatePanels -- for one model + species, one panel per metric showing
# the CUMULATIVE MEAN of that metric (taken at the final timestep) as replicates
# are added. Where a curve flattens, extra replicates stop changing the estimate.
#   colour = disturbance (off / on); grey dashed vertical = ref_n replicates.
# Call once per species -> separate figures.
# --------------------------------------------------------------------------
plotReplicatePanels <- function(results, model, nsp,
                                metrics = c("species_richness", "structural_complexity",
                                            "turnover", "shannon", "total_cover", "evenness"),
                                ref_n = 30, title = NULL) {
  d    <- results[results$combination == model & results$n_species == nsp, ]
  d    <- d[d$timestep == max(d$timestep), ]          # final-state value per replicate
  dcol <- c("steelblue", "firebrick")                 # off, on
  nrep <- max(d$replicate)

  op <- par(mfrow = c(3, 2), oma = c(0, 0, 9, 0), mar = c(4.5, 4.5, 3, 1.5)); on.exit(par(op))
  for (m in metrics) {
    # running mean over replicates (in replicate order) for each disturbance state,
    # skipping NA/NaN values (e.g. evenness is undefined when richness collapses to 1)
    series <- lapply(c(FALSE, TRUE), function(ds) {
      v  <- d[d$disturbance_on == ds, ]
      x  <- v[[m]][order(v$replicate)]
      ok <- !is.na(x)
      cumsum(ifelse(ok, x, 0)) / cumsum(ok)          # mean of non-NA values so far
    })
    ylim <- range(unlist(series), na.rm = TRUE)
    plot(NA, xlim = c(1, nrep), ylim = ylim,
         xlab = "Number of replicates", ylab = paste("Cumulative mean", m), main = m)
    abline(v = ref_n, lty = 3, col = "grey60")        # a typical replicate count
    for (j in seq_along(series)) lines(seq_along(series[[j]]), series[[j]], col = dcol[j], lwd = 2)
  }

  ttl <- if (is.null(title))
           paste0(sub("^Classic_", "", model), " | n", nsp, " - replicate convergence")
         else title
  mtext(ttl, outer = TRUE, line = 6, font = 2, cex = 1.4)

  # one shared key across the top strip
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot(0, 0, type = "n", xlim = c(0, 1), ylim = c(0, 1), axes = FALSE, xlab = "", ylab = "")
  legend(0.5, 0.965, xjust = 0.5, yjust = 1, ncol = 3, bty = "n", cex = 1.0, seg.len = 1.6,
         legend = c("disturbance off", "disturbance on", paste0("n = ", ref_n, " (typical)")),
         col = c("steelblue", "firebrick", "grey60"), lwd = c(2, 2, 1), lty = c(1, 1, 3))
  invisible(NULL)
}
#It is final timestep values, mean of this every time 1 replicate added
# One figure per model + species; the grey line marks a typical replicate count.
plotReplicateCI(reps_data, "Classic_RPS", 3)
plotReplicateCI(reps_data, "Classic_RPS", 7)
plotReplicateCI(reps_data, "Classic_Linear", 3)
plotReplicateCI(reps_data, "Classic_Linear", 7)
plotReplicateCI(reps_data, "Classic_Neutral", 3)
plotReplicateCI(reps_data, "Classic_Neutral", 7)
plotReplicateCI(reps_data, "Classic_Random", 7)
plotReplicateCI(reps_data, "Classic_Random", 3)



# --------------------------------------------------------------------------
# plotReplicateCI -- same running mean as plotReplicatePanels, but with a
# shaded 95% CI band (mean +/- 1.96 * SE, where SE = sd / sqrt(N)). The band
# shows the UNCERTAINTY of the estimate, and it narrows as replicates are added.
# When the band is acceptably thin, you have enough replicates - a firmer test
# than just watching the mean flatten.
#   colour = disturbance (off / on); grey dashed vertical = ref_n replicates.
# --------------------------------------------------------------------------
plotReplicateCI <- function(results, model, nsp,
                            metrics = c("species_richness", "structural_complexity",
                                        "turnover", "shannon", "total_cover", "evenness"),
                            ref_n = 30, title = NULL) {
  d    <- results[results$combination == model & results$n_species == nsp, ]
  d    <- d[d$timestep == max(d$timestep), ]          # final-state value per replicate
  dcol <- c("steelblue", "firebrick")                 # off, on
  nrep <- max(d$replicate)

  # running mean and 95% CI over the first N replicates (NA-robust)
  running <- function(x) {
    ok  <- !is.na(x)
    cnt <- cumsum(ok)                                  # valid count so far
    sx  <- cumsum(ifelse(ok, x, 0))
    sxx <- cumsum(ifelse(ok, x^2, 0))
    mean_N <- sx / cnt
    var_N  <- (sxx - sx^2 / cnt) / (cnt - 1)           # running sample variance
    var_N[cnt < 2] <- NA
    se_N   <- sqrt(var_N / cnt)
    list(mean = mean_N, lo = mean_N - 1.96 * se_N, hi = mean_N + 1.96 * se_N)
  }

  op <- par(mfrow = c(3, 2), oma = c(0, 0, 10, 0), mar = c(4.5, 4.5, 3, 1.5)); on.exit(par(op))
  for (m in metrics) {
    stats <- lapply(c(FALSE, TRUE), function(ds) {
      v <- d[d$disturbance_on == ds, ]
      running(v[[m]][order(v$replicate)])
    })
    ylim <- range(unlist(lapply(stats, function(s) c(s$lo, s$hi))), na.rm = TRUE)
    plot(NA, xlim = c(1, nrep), ylim = ylim,
         xlab = "Number of replicates", ylab = paste("Mean", m), main = m)
    abline(v = ref_n, lty = 3, col = "grey60")
    for (j in seq_along(stats)) {
      s <- stats[[j]]; N <- seq_along(s$mean); ok <- !is.na(s$lo)
      if (any(ok))
        polygon(c(N[ok], rev(N[ok])), c(s$lo[ok], rev(s$hi[ok])),
                col = adjustcolor(dcol[j], 0.18), border = NA)
      lines(N, s$mean, col = dcol[j], lwd = 2)
    }
  }

  ttl <- if (is.null(title))
           paste0(sub("^Classic_", "", model), " | n", nsp, " - replicate convergence (95% CI)")
         else title
  # title, subtitle and shared key drawn together in one whole-device overlay
  # with fixed spacing, so they cannot overlap each other on any device
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot(0, 0, type = "n", xlim = c(0, 1), ylim = c(0, 1), axes = FALSE, xlab = "", ylab = "")
  text(0.5, 0.992, ttl, font = 2, cex = 1.4)
  text(0.5, 0.956, "shaded band = 95% CI of the mean (+/-1.96 SE); it narrows as replicates are added",
       cex = 0.85, col = "grey30")
  legend(0.5, 0.928, xjust = 0.5, yjust = 1, ncol = 3, bty = "n", cex = 1.0, seg.len = 1.6,
         legend = c("disturbance off", "disturbance on", paste0("n = ", ref_n, " (typical)")),
         col = c("steelblue", "firebrick", "grey60"), lwd = c(2, 2, 1), lty = c(1, 1, 3))
  invisible(NULL)
}

# Same scenarios, now with the shrinking 95% CI band.
plotReplicateCI(reps_data, "Classic_RPS", 3)
plotReplicateCI(reps_data, "Classic_RPS", 7)

######BOOTSTRAP VERSION
# --------------------------------------------------------------------------
# plotReplicateBootstrap -- order-independent version. The running-mean plots
# walk replicates in one fixed order; here, for each candidate N we instead draw
# B random subsamples of N replicates (WITHOUT replacement) from all 100, take
# each subsample's mean, and show the median (line) and 2.5-97.5% envelope
# (band). The band = the range of means you might report if you had run only N
# replicates, free of any ordering bias. It collapses to a point at N = 100
# (subsampling all 100 = the full set).
#   colour = disturbance (off / on); grey dashed vertical = ref_n replicates.
# --------------------------------------------------------------------------
plotReplicateBootstrap <- function(results, model, nsp,
                                   metrics = c("species_richness", "structural_complexity",
                                               "turnover", "shannon", "total_cover", "evenness"),
                                   B = 200, ref_n = 30, title = NULL, seed = 1) {
  d    <- results[results$combination == model & results$n_species == nsp, ]
  d    <- d[d$timestep == max(d$timestep), ]          # final-state value per replicate
  dcol <- c("steelblue", "firebrick")                 # off, on
  nrep <- max(d$replicate)
  Ns   <- seq_len(nrep)
  set.seed(seed)                                       # reproducible resampling

  # for a vector of final-state values, return an nrep x 3 matrix of
  # (2.5%, 50%, 97.5%) quantiles of the subsample mean at each N
  envelope <- function(x) {
    x <- x[!is.na(x)]; n <- length(x)
    t(vapply(Ns, function(N) {
      if (N > n) return(c(NA, NA, NA))
      means <- replicate(B, mean(x[sample.int(n, N)]))   # subsample without replacement
      quantile(means, c(0.025, 0.5, 0.975), names = FALSE)
    }, numeric(3)))
  }

  op <- par(mfrow = c(3, 2), oma = c(0, 0, 10, 0), mar = c(4.5, 4.5, 3, 1.5)); on.exit(par(op))
  for (m in metrics) {
    envs <- lapply(c(FALSE, TRUE), function(ds) envelope(d[d$disturbance_on == ds, m]))
    ylim <- range(unlist(envs), na.rm = TRUE)
    plot(NA, xlim = c(1, nrep), ylim = ylim,
         xlab = "Number of replicates", ylab = paste("Mean", m), main = m)
    abline(v = ref_n, lty = 3, col = "grey60")
    for (j in seq_along(envs)) {
      e <- envs[[j]]; ok <- !is.na(e[, 1])
      if (any(ok))
        polygon(c(Ns[ok], rev(Ns[ok])), c(e[ok, 1], rev(e[ok, 3])),
                col = adjustcolor(dcol[j], 0.18), border = NA)
      lines(Ns, e[, 2], col = dcol[j], lwd = 2)
    }
  }

  ttl <- if (is.null(title))
           paste0(sub("^Classic_", "", model), " | n", nsp, " - replicate subsampling")
         else title
  # title, subtitle and shared key drawn together in one whole-device overlay
  # with fixed spacing, so they cannot overlap each other on any device
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot(0, 0, type = "n", xlim = c(0, 1), ylim = c(0, 1), axes = FALSE, xlab = "", ylab = "")
  text(0.5, 0.992, ttl, font = 2, cex = 1.4)
  text(0.5, 0.956, paste0("band = 2.5-97.5% of the mean over ", B,
                          " random subsamples of N replicates (order-independent)"),
       cex = 0.85, col = "grey30")
  legend(0.5, 0.928, xjust = 0.5, yjust = 1, ncol = 3, bty = "n", cex = 1.0, seg.len = 1.6,
         legend = c("disturbance off", "disturbance on", paste0("n = ", ref_n, " (typical)")),
         col = c("steelblue", "firebrick", "grey60"), lwd = c(2, 2, 1), lty = c(1, 1, 3))
  invisible(NULL)
}

# Order-independent subsampling version.
plotReplicateBootstrap(reps_data, "Classic_RPS", 3)
plotReplicateBootstrap(reps_data, "Classic_RPS", 7)




#SIM 3
#Parameter sensitivity
#Disturbance regime differences
#Can do ANOVA to test end state importance of disturbance regime on final state metrics, but first need to read in the data
parares <- read.csv("C:/Users/dell/Documents/Research Project 2/Simulation_Results/3_Parameter_Sensitivity/Run_20260727_Sensitivity_DisturbanceRegime/Master_Sensitivity_DisturbanceRegime_20260727_213926_results.csv",
                        stringsAsFactors = FALSE)

# --------------------------------------------------------------------------
# plotRegimePanels -- for one matrix + one species, all disturbance regimes on
# each metric panel (mean over replicates, means only). A regime = frequency x
# size; the no-disturbance run is a black reference.
#   colour    = frequency (rarely / often / very_often)
#   line type = size (small dotted / medium dashed / large solid)
#   black     = no disturbance (baseline)
# --------------------------------------------------------------------------
#   smooth    = rolling-mean window in timesteps to de-noise the lines (0/1 = off)
plotRegimePanels <- function(results, model, nsp,
                             metrics = c("species_richness", "structural_complexity",
                                         "turnover", "shannon", "total_cover", "evenness"),
                             title = NULL, smooth = 25) {
  d  <- results[results$combination == model & results$n_species == nsp, ]
  ts <- sort(unique(d$timestep))
  # colours/line types from the regime levels actually present (ordered), so the
  # key never shows levels this file doesn't have
  freq_ord <- c("rarely", "often", "very_often"); size_ord <- c("small", "medium", "large")
  freqs <- freq_ord[freq_ord %in% d$dist_freq];      sizes <- size_ord[size_ord %in% d$dist_size]
  fcol <- setNames(c("forestgreen", "darkorange", "firebrick")[match(freqs, freq_ord)], freqs)
  slty <- setNames(c(3, 2, 1)[match(sizes, size_ord)], sizes)   # small dotted .. large solid

  # style each scenario (regime) from its freq/size; no-disturbance = black
  style_of <- function(on, fr, sz)
    if (!isTRUE(on)) list(col = "black", lty = 1, lwd = 2.6)
    else             list(col = fcol[[fr]], lty = slty[[sz]], lwd = 1.8)

  # centred rolling mean (window shrinks at the ends, NA-robust); off if w <= 1
  roll <- function(y, w) {
    if (w <= 1) return(y)
    h <- floor(w / 2)
    vapply(seq_along(y), function(i)
      mean(y[max(1, i - h):min(length(y), i + h)], na.rm = TRUE), numeric(1))
  }

  scen <- unique(d$scenario)
  op <- par(mfrow = c(3, 2), oma = c(0, 0, 11, 0), mar = c(4.5, 4.5, 3, 1.5)); on.exit(par(op))
  for (m in metrics) {
    series <- lapply(scen, function(s) {
      v <- d[d$scenario == s, ]
      list(y  = roll(tapply(v[[m]], v$timestep, mean, na.rm = TRUE)[as.character(ts)], smooth),
           st = style_of(v$disturbance_on[1], v$dist_freq[1], v$dist_size[1]))
    })
    series <- series[order(vapply(series, function(s) s$st$col == "black", logical(1)))]  # baseline last
    ylim <- range(unlist(lapply(series, `[[`, "y")), na.rm = TRUE)
    plot(NA, xlim = range(ts), ylim = ylim, xlab = "Timestep", ylab = paste("Mean", m), main = m)
    for (s in series) lines(ts, s$y, col = s$st$col, lty = s$st$lty, lwd = s$st$lwd)
  }

  ttl <- if (is.null(title))
           paste0(sub("^Classic_", "", model), " | n", nsp, " - disturbance regimes")
         else title

  # title + two keys (frequency colours, size line types) in one overlay
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot(0, 0, type = "n", xlim = c(0, 1), ylim = c(0, 1), axes = FALSE, xlab = "", ylab = "")
  text(0.5, 0.992, ttl, font = 2, cex = 1.4)
  legend(0.34, 0.958, title = "frequency (colour)", horiz = TRUE, xjust = 0.5, yjust = 1,
         bty = "n", cex = 0.85, seg.len = 1.6,
         legend = c("none", gsub("_", " ", freqs)),
         col = c("black", fcol), lwd = c(2.6, rep(1.8, length(freqs))), lty = 1)
  legend(0.75, 0.958, title = "size (line type)", horiz = TRUE, xjust = 0.5, yjust = 1,
         bty = "n", cex = 0.85, seg.len = 1.6,
         legend = sizes, col = "grey30", lwd = 1.8, lty = slty)
  invisible(NULL)
}

# One figure per matrix + species: all disturbance regimes overlaid.
plotRegimePanels(parares, "Classic_RPS", 3)
plotRegimePanels(parares, "Classic_RPS", 7)
plotRegimePanels(parares, "Classic_Linear", 3)
plotRegimePanels(parares, "Classic_Linear", 7)
plotRegimePanels(parares, "Classic_Neutral", 3)
plotRegimePanels(parares, "Classic_Neutral", 7)
plotRegimePanels(parares, "Classic_Random", 3)
plotRegimePanels(parares, "Classic_Random", 7)

#Sim 4, test of if RPS dyanmics work and are robust
rpsres <- read.csv("C:/Users/dell/Documents/Research Project 2/Simulation_Results/4_Base_model_exploration/Base_Models/Run_20260724_Initial_RPS_RawOutputs/my_run_20260724_174603_results.csv",
                   stringsAsFactors = FALSE)

# small centred rolling mean, shared by the RPS plots below (off if w <= 1)
rpsRoll <- function(y, w) {
  if (w <= 1) return(y)
  h <- floor(w / 2)
  vapply(seq_along(y), function(i)
    mean(y[max(1, i - h):min(length(y), i + h)], na.rm = TRUE), numeric(1))
}

# --------------------------------------------------------------------------
# plotSpeciesCoexistence -- the headline "RPS is sustained" plot: mean cover of
# each species through time (mean over replicates) for ONE scenario. If the
# species coexist, every line stays well above zero rather than one taking over.
#   nsp     : number of species (only Sp1..Sp<nsp> are drawn)
#   indiv   : founders per species; dist_on : disturbance TRUE/FALSE
#   smooth  : rolling-mean window in timesteps (0/1 = off)
# --------------------------------------------------------------------------
plotSpeciesCoexistence <- function(results, nsp, indiv, dist_on, smooth = 10,
                                   show_reps = TRUE, main = NULL, ylim = NULL) {
  d  <- results[results$n_species == nsp & results$individuals == indiv &
                results$disturbance_on == dist_on, ]
  ts <- sort(unique(d$timestep))
  spcols <- paste0("Sp", seq_len(nsp), "_cover")
  cols   <- hcl.colors(nsp, "Dark 3")

  means <- lapply(spcols, function(cc)
    rpsRoll(tapply(d[[cc]], d$timestep, mean, na.rm = TRUE)[as.character(ts)], smooth))

  # y-range spans the faint replicate lines too, not just the means (unless the
  # caller passes a fixed ylim, e.g. to share one scale across stacked panels)
  if (is.null(ylim))
    ylim <- if (show_reps) range(0, unlist(d[spcols]), na.rm = TRUE)
            else           range(0, unlist(means), na.rm = TRUE)

  op <- par(mar = c(4.5, 4.5, 3, 1.5)); on.exit(par(op))
  plot(NA, xlim = range(ts), ylim = ylim, xlab = "Timestep", ylab = "Cover (% of reef)",
       main = if (is.null(main))
                sprintf("RPS coexistence | n%d | %d founders/sp | disturbance %s",
                        nsp, indiv, if (isTRUE(dist_on)) "on" else "off")
              else main)
  grid(col = "grey92")

  # faint individual-replicate trajectories per species, drawn first (background)
  if (show_reps) for (i in seq_len(nsp)) {
    faint <- adjustcolor(cols[i], alpha.f = 0.08)
    for (r in unique(d$run_id)) {
      one <- d[d$run_id == r, c("timestep", spcols[i])]
      one <- one[order(one$timestep), ]
      lines(one$timestep, one[[spcols[i]]], col = faint, lwd = 1)
    }
  }

  # bold mean cover per species on top
  for (i in seq_len(nsp)) lines(ts, means[[i]], col = cols[i], lwd = 2.6)
  legend("topright", legend = paste0("Sp", seq_len(nsp)), col = cols, lwd = 2.6,
         bty = "n", title = "species", cex = 0.9)
  invisible(NULL)
}

# --------------------------------------------------------------------------
# plotCoexistencePair -- disturbance OFF (top) and ON (bottom) stacked in one
# figure, sharing a common y-axis so the two are directly comparable.
# --------------------------------------------------------------------------
plotCoexistencePair <- function(results, nsp, indiv, smooth = 10, show_reps = TRUE) {
  d      <- results[results$n_species == nsp & results$individuals == indiv, ]
  spcols <- paste0("Sp", seq_len(nsp), "_cover")

  # shared y-range across BOTH disturbance states
  ylim <- if (show_reps)
            range(0, unlist(d[spcols]), na.rm = TRUE)
          else
            range(0, unlist(lapply(c(FALSE, TRUE), function(ds) {
              v <- d[d$disturbance_on == ds, ]
              lapply(spcols, function(cc) tapply(v[[cc]], v$timestep, mean, na.rm = TRUE))
            })), na.rm = TRUE)

  op <- par(mfrow = c(2, 1), oma = c(0, 0, 3, 0)); on.exit(par(op))
  plotSpeciesCoexistence(results, nsp, indiv, dist_on = FALSE, smooth = smooth,
                         show_reps = show_reps, ylim = ylim, main = "Disturbance OFF")
  plotSpeciesCoexistence(results, nsp, indiv, dist_on = TRUE,  smooth = smooth,
                         show_reps = show_reps, ylim = ylim, main = "Disturbance ON")
  mtext(sprintf("RPS coexistence | n%d | %d founders/sp", nsp, indiv),
        outer = TRUE, line = 0.7, font = 2, cex = 1.3)
  invisible(NULL)
}

# --------------------------------------------------------------------------
# plotRPSmetrics -- the "other metrics" view: 6 metric panels over time for one
# (n_species, founders) config, disturbance off vs on as two lines. Shows the
# summary metrics hold up (richness near N, evenness/shannon high, cover stable).
# --------------------------------------------------------------------------
plotRPSmetrics <- function(results, nsp, indiv, smooth = 10, show_band = TRUE,
                           metrics = c("species_richness", "structural_complexity",
                                       "turnover", "shannon", "total_cover", "evenness"),
                           title = NULL) {
  d    <- results[results$n_species == nsp & results$individuals == indiv, ]
  ts   <- sort(unique(d$timestep))
  dcol <- c("steelblue", "firebrick")                 # off, on

  op <- par(mfrow = c(3, 2), oma = c(0, 0, 9, 0), mar = c(4.5, 4.5, 3, 1.5)); on.exit(par(op))
  for (m in metrics) {
    stats <- lapply(c(FALSE, TRUE), function(ds) {
      v  <- d[d$disturbance_on == ds, ]
      list(mu = rpsRoll(tapply(v[[m]], v$timestep, mean, na.rm = TRUE)[as.character(ts)], smooth),
           sd = rpsRoll(tapply(v[[m]], v$timestep, sd,   na.rm = TRUE)[as.character(ts)], smooth))
    })
    # y-range spans the +/-1 SD bands when shown, else just the means
    ylim <- if (show_band)
              range(unlist(lapply(stats, function(s) c(s$mu - s$sd, s$mu + s$sd))), na.rm = TRUE)
            else range(unlist(lapply(stats, `[[`, "mu")), na.rm = TRUE)
    plot(NA, xlim = range(ts), ylim = ylim, xlab = "Timestep", ylab = paste("Mean", m), main = m)
    # +/-1 SD band per disturbance state, drawn first so the means sit on top
    if (show_band) for (j in seq_along(stats)) {
      s <- stats[[j]]; ok <- !is.na(s$mu) & !is.na(s$sd)
      if (any(ok))
        polygon(c(ts[ok], rev(ts[ok])), c((s$mu - s$sd)[ok], rev((s$mu + s$sd)[ok])),
                col = adjustcolor(dcol[j], 0.15), border = NA)
    }
    for (j in seq_along(stats)) lines(ts, stats[[j]]$mu, col = dcol[j], lwd = 2)
  }
  ttl <- if (is.null(title))
           sprintf("RPS | n%d | %d founders/sp - metrics over time", nsp, indiv) else title

  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot(0, 0, type = "n", xlim = c(0, 1), ylim = c(0, 1), axes = FALSE, xlab = "", ylab = "")
  text(0.5, 0.99, ttl, font = 2, cex = 1.4)
  legend(0.5, 0.955, xjust = 0.5, yjust = 1, horiz = TRUE, bty = "n", cex = 1.0, seg.len = 1.8,
         legend = c("disturbance off", "disturbance on"), col = dcol, lwd = 2)
  invisible(NULL)
}

# Headline coexistence plot (one scenario), or both disturbance states stacked.
#For one
plotSpeciesCoexistence(rpsres, nsp = 3, indiv = 2, dist_on = FALSE)

#For two with disturbance on/off -- one coexistence pair per (species, founders) combo
plotCoexistencePair(rpsres, nsp = 3, indiv = 2)
plotCoexistencePair(rpsres, nsp = 3, indiv = 4)
plotCoexistencePair(rpsres, nsp = 5, indiv = 2)
plotCoexistencePair(rpsres, nsp = 5, indiv = 4)


#Metrics -- one RPS metric panel per (species, founders) combo
plotRPSmetrics(rpsres, nsp = 3, indiv = 2)
plotRPSmetrics(rpsres, nsp = 3, indiv = 4)
plotRPSmetrics(rpsres, nsp = 5, indiv = 2)
plotRPSmetrics(rpsres, nsp = 5, indiv = 4)



#SIM 3 - Parameter sensitivity: metric panels, one line per parameter level ----
# Generalised from plotReefPanels: colour = a gradient parameter (any numeric
# column, e.g. "bias" or "intraspecific"), line type = disturbance (solid on /
# dashed off), one panel per metric. One figure per model + species.
plotParamComparison <- function(results, model, nsp, param, metric = "species_richness",
                                main = NULL, show_legend = TRUE) {
  d    <- results[results$combination == model & results$n_species == nsp, ]
  levs <- sort(unique(d[[param]]))
  pcol <- setNames(hcl.colors(length(levs), "viridis"), levs)   # sequential = ordered gradient
  ts   <- sort(unique(d$timestep))

  combos <- expand.grid(lev = levs, dist = c(TRUE, FALSE), stringsAsFactors = FALSE)
  series <- lapply(seq_len(nrow(combos)), function(i) {
    k <- d[[param]] == combos$lev[i] & d$disturbance_on == combos$dist[i]
    if (!any(k)) NULL else tapply(d[[metric]][k], d$timestep[k], mean, na.rm = TRUE)[as.character(ts)]
  })
  ylim <- range(unlist(series), na.rm = TRUE)

  right_mar <- if (show_legend) 9 else 1.5
  op <- par(mar = c(4.5, 4.5, 3, right_mar), xpd = NA); on.exit(par(op))
  plot(NA, xlim = range(ts), ylim = ylim, xlab = "Timestep", ylab = paste("Mean", metric),
       main = if (is.null(main)) paste0(sub("^Classic_", "", model), " | n", nsp, " - ", metric, " vs ", param) else main)
  for (i in seq_len(nrow(combos))) {
    y <- series[[i]]; if (is.null(y)) next
    lines(ts, y, col = pcol[as.character(combos$lev[i])], lty = if (combos$dist[i]) 1 else 2, lwd = 2.2)
  }
  if (show_legend) {
    ux <- par("usr")[2]; yhi <- par("usr")[4]; yr <- yhi - par("usr")[3]
    legend(ux, yhi,             title = param,        xjust = 0, bty = "n", cex = 0.8, legend = levs, col = pcol, lwd = 3)
    legend(ux, yhi - 0.34 * yr, title = "disturbance", xjust = 0, bty = "n", cex = 0.8, legend = c("on", "off"), lty = c(1, 2), lwd = 2.2)
  }
  invisible(NULL)
}

plotParamPanels <- function(results, model, nsp, param,
                            metrics = c("species_richness", "structural_complexity",
                                        "turnover", "shannon", "total_cover", "evenness"),
                            title = NULL) {
  d    <- results[results$combination == model & results$n_species == nsp, ]
  levs <- sort(unique(d[[param]]))
  pcol <- setNames(hcl.colors(length(levs), "viridis"), levs)

  op <- par(mfrow = c(3, 2), oma = c(0, 0, 9, 0)); on.exit(par(op))
  for (m in metrics)
    plotParamComparison(results, model, nsp, param, metric = m, main = m, show_legend = FALSE)

  ttl <- if (is.null(title)) paste0(sub("^Classic_", "", model), " | n", nsp, " - metrics vs ", param) else title
  mtext(ttl, outer = TRUE, line = 6, font = 2, cex = 1.4)

  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot(0, 0, type = "n", xlim = c(0, 1), ylim = c(0, 1), axes = FALSE, xlab = "", ylab = "")
  legend(0.5, 0.965, xjust = 0.5, yjust = 1, ncol = length(levs) + 2, bty = "n",
         cex = 1.0, seg.len = 1.6, x.intersp = 0.7,
         legend = c(paste0(param, " ", levs), "dist on", "dist off"),
         col    = c(pcol, "grey30", "grey30"),
         lwd    = 2.5, lty = c(rep(1, length(levs)), 1, 2))
  invisible(NULL)
}

#SIM 3 - BIAS IMPACT (colour = bias level, solid/dashed = disturbance)
bias <- read.csv("C:/Users/dell/Documents/Research Project 2/Simulation_Results/3_Parameter_Sensitivity/Bias_Sensitivity_Checks/Master_Sensitivity_Bias_20260731_010903_results.csv",
                 stringsAsFactors = FALSE)
plotParamPanels(bias, "Classic_RPS",     3, "bias")
plotParamPanels(bias, "Classic_RPS",     7, "bias")
plotParamPanels(bias, "Classic_Linear",  3, "bias")
plotParamPanels(bias, "Classic_Linear",  7, "bias")
plotParamPanels(bias, "Classic_Neutral", 3, "bias")
plotParamPanels(bias, "Classic_Neutral", 7, "bias")
plotParamPanels(bias, "Classic_Random",  3, "bias")
plotParamPanels(bias, "Classic_Random",  7, "bias")



#SIM 3 - INTRASPECIFIC COMP IMPACT (colour = intraspecific level, 0..1)
# NOTE: intraspecific is a 5-level gradient (0, 0.25, 0.5, 0.75, 1), not on/off.
intra <- read.csv("C:/Users/dell/Documents/Research Project 2/Simulation_Results/3_Parameter_Sensitivity/Intra_Sens_Test/Master_Sensitivity_Intraspecific_20260731_054833_results.csv",
                  stringsAsFactors = FALSE)
plotParamPanels(intra, "Classic_RPS",     3, "intraspecific")
plotParamPanels(intra, "Classic_RPS",     7, "intraspecific")
plotParamPanels(intra, "Classic_Linear",  3, "intraspecific")
plotParamPanels(intra, "Classic_Linear",  7, "intraspecific")
plotParamPanels(intra, "Classic_Neutral", 3, "intraspecific")
plotParamPanels(intra, "Classic_Neutral", 7, "intraspecific")
plotParamPanels(intra, "Classic_Random",  3, "intraspecific")
plotParamPanels(intra, "Classic_Random",  7, "intraspecific")



#Sim 3 - Reproduction fecundity (colour = fecundity level = repro_base: 0.02..0.8)
# 9 matrix-variants (base model x reproduction mode); fecundity is the `repro_base`
# column. Runs go to t = 2500. Same key: colour = fecundity, solid = disturbance
# on, dashed = off. One figure per matrix-variant + species (change 7 -> 3 for n3).
repro <- read.csv("C:/Users/dell/Documents/Research Project 2/Simulation_Results/3_Parameter_Sensitivity/Repro_Sens_Test/Reproduction_Fecundity_20260731_222545_results.csv",
                  stringsAsFactors = FALSE)

plotParamPanels(repro, "RPS_reproduction",         3, "repro_base")
plotParamPanels(repro, "RPS_reproductionEven",     3, "repro_base")
plotParamPanels(repro, "RPS_reproduction",         7, "repro_base")
plotParamPanels(repro, "RPS_reproductionEven",     7, "repro_base")

plotParamPanels(repro, "Linear_reproductionEven",     3, "repro_base")
plotParamPanels(repro, "Linear_reproductionNormal",   3, "repro_base")
plotParamPanels(repro, "Linear_reproductionOpposite", 3, "repro_base")
plotParamPanels(repro, "Linear_reproductionEven",     7, "repro_base")
plotParamPanels(repro, "Linear_reproductionNormal",   7, "repro_base")
plotParamPanels(repro, "Linear_reproductionOpposite", 7, "repro_base")

plotParamPanels(repro, "Random_reproduction",      3, "repro_base")
plotParamPanels(repro, "Random_reproductionEven",  3, "repro_base")
plotParamPanels(repro, "Random_reproduction",      7, "repro_base")
plotParamPanels(repro, "Random_reproductionEven",  7, "repro_base")

plotParamPanels(repro, "Neutral_reproduction",     7, "repro_base")
plotParamPanels(repro, "Neutral_reproductionEven", 7, "repro_base")
plotParamPanels(repro, "Neutral_reproduction",     3, "repro_base")
plotParamPanels(repro, "Neutral_reproductionEven", 3, "repro_base")




#Sim 3 - reproduction x disturbance regime combo ----------------------------
# Same figures as the disturbance-regime sim (colour = frequency, line type =
# size, black = no disturbance). Reproduction is encoded in the combination:
# Classic_* = repro OFF, *_reproductionEven = repro ON -- so each matrix gives an
# OFF figure and an ON figure to compare. (2 freq x 2 size here, not 3 x 3.)
reprodist <- read.csv("C:/Users/dell/Documents/Research Project 2/Simulation_Results/3_Parameter_Sensitivity/Repro_Dist_sens_Test/Reproduction_DisturbanceRegime_20260801_123052_results.csv",
                      stringsAsFactors = FALSE)



plotRegimePanels(reprodist, "RPS_reproductionEven",    3) 
plotRegimePanels(reprodist, "RPS_reproductionEven",    7)  

 # RPS   | reproduction on
plotRegimePanels(reprodist, "Classic_Linear",          3)   # Linear| off
plotRegimePanels(reprodist, "Linear_reproductionEven", 3)
plotRegimePanels(reprodist, "Classic_Linear",          7)   # Linear| off
plotRegimePanels(reprodist, "Linear_reproductionEven", 7)  
 # Linear| on
plotRegimePanels(reprodist, "Classic_Neutral",         3)
plotRegimePanels(reprodist, "Neutral_reproductionEven",3)
plotRegimePanels(reprodist, "Classic_Neutral",         7)
plotRegimePanels(reprodist, "Neutral_reproductionEven",7)

plotRegimePanels(reprodist, "Classic_Random",          3)
plotRegimePanels(reprodist, "Random_reproductionEven", 3)
plotRegimePanels(reprodist, "Classic_Random",          7)
plotRegimePanels(reprodist, "Random_reproductionEven", 7)



#Sim 3 - GROWTH RATE IMPACT (colour = growth rate, solid/dashed = disturbance) ---
# growth is a 3-level gradient (1 = slow, 3 = medium, 5 = fast, uniform across all
# species). Metric panels per matrix x species count, one coloured line per growth
# level; same style as the bias / intraspecific / fecundity sweeps above.
growth <- read.csv("C:/Users/dell/Documents/Research Project 2/Simulation_Results/3.1_Extended_Sensitivity/GrowthRate_Sens/Growth_Rate_20260801_230519_results.csv",
                   stringsAsFactors = FALSE)

plotParamPanels(growth, "Classic_RPS",     3, "growth")
plotParamPanels(growth, "Classic_RPS",     7, "growth")
plotParamPanels(growth, "Classic_Linear",  3, "growth")
plotParamPanels(growth, "Classic_Linear",  7, "growth")
plotParamPanels(growth, "Classic_Neutral", 3, "growth")
plotParamPanels(growth, "Classic_Neutral", 7, "growth")
plotParamPanels(growth, "Classic_Random",  3, "growth")
plotParamPanels(growth, "Classic_Random",  7, "growth")



#Sim 3 - BACKGROUND RECRUITMENT IMPACT (colour = background_chance, solid/dashed =
# disturbance) --- external larval supply per species per pulse: 0 (closed), 0.5,
# 1. Metric panels per matrix x species count; same style as the sweeps above.
back <- read.csv("C:/Users/dell/Documents/Research Project 2/Simulation_Results/3.1_Extended_Sensitivity/BackRecruit_Sens/BackgroundRecruitment_20260802_025433_results.csv",
                 stringsAsFactors = FALSE)

plotParamPanels(back, "Classic_RPS",     3, "background_chance")
plotParamPanels(back, "Classic_RPS",     7, "background_chance")
plotParamPanels(back, "Classic_Linear",  3, "background_chance")
plotParamPanels(back, "Classic_Linear",  7, "background_chance")
plotParamPanels(back, "Classic_Neutral", 3, "background_chance")
plotParamPanels(back, "Classic_Neutral", 7, "background_chance")
plotParamPanels(back, "Classic_Random",  3, "background_chance")
plotParamPanels(back, "Classic_Random",  7, "background_chance")



#Sim 3 - MATURITY AGE IMPACT (colour = maturity_age, solid/dashed = disturbance) --
# timesteps a colony must be alive before it can reproduce: 5, 30, 150. Only
# meaningful when reproduction is ON, so the sweep uses the *_reproductionEven
# variants (one per matrix). Metric panels per matrix x species count.
mat <- read.csv("C:/Users/dell/Documents/Research Project 2/Simulation_Results/3.1_Extended_Sensitivity/MaturityAge_Sens/MaturityAge_20260802_071441_results.csv",
                stringsAsFactors = FALSE)

plotParamPanels(mat, "RPS_reproductionEven",     3, "maturity_age")
plotParamPanels(mat, "RPS_reproductionEven",     7, "maturity_age")
plotParamPanels(mat, "Linear_reproductionEven",  3, "maturity_age")
plotParamPanels(mat, "Linear_reproductionEven",  7, "maturity_age")
plotParamPanels(mat, "Neutral_reproductionEven", 3, "maturity_age")
plotParamPanels(mat, "Neutral_reproductionEven", 7, "maturity_age")
plotParamPanels(mat, "Random_reproductionEven",  3, "maturity_age")
plotParamPanels(mat, "Random_reproductionEven",  7, "maturity_age")



#Sim 3 - SIZE-COMPETITION STRENGTH IMPACT (colour = size_beta_max, solid/dashed =
# disturbance) --- strength of the colony-size advantage in contests (log-ratio):
# 0 (size-independent), 0.3, 0.6, 0.9. Uses the size-variant combinations (RPS /
# Random sizecompImpact; Linear SizeImpact Normal & Reverse). Panels per matrix x
# species count.
sizecomp <- read.csv("C:/Users/dell/Documents/Research Project 2/Simulation_Results/3.1_Extended_Sensitivity/SizeComp_Sens/SizeCompetitionStrength_20260802_115434_results.csv",
                     stringsAsFactors = FALSE)

plotParamPanels(sizecomp, "RPS_sizecompImpact",       3, "size_beta_max")
plotParamPanels(sizecomp, "RPS_sizecompImpact",       7, "size_beta_max")
plotParamPanels(sizecomp, "Linear_SizeImpactNormal",  3, "size_beta_max")
plotParamPanels(sizecomp, "Linear_SizeImpactNormal",  7, "size_beta_max")
plotParamPanels(sizecomp, "Linear_SizeImpactReverse", 3, "size_beta_max")
plotParamPanels(sizecomp, "Linear_SizeImpactReverse", 7, "size_beta_max")
plotParamPanels(sizecomp, "Random_sizecompImpact",    3, "size_beta_max")
plotParamPanels(sizecomp, "Random_sizecompImpact",    7, "size_beta_max")












###########SIM 4 ISSUES AS BIG SWEEP#######




# ==========================================================================
# Simulation 4, remaining model robustness -- BIG factorial sweeps
#   The AllModels sweep is ~1.2M rows / 445 MB. Golden rule: read the giant file
#   ONCE, boil it down to small per-scenario summaries, save those as .rds, and
#   do ALL plotting from the summaries (never re-parse the CSV).
# ==========================================================================
sweep_paths <- c(
  full  = "C:/Users/dell/Documents/Research Project 2/Simulation_Results/4_Base_model_exploration/Base_Models/Full_Run_20260724_Classic-RPS-Linear-Neutral_n357_reef50_t1000_r30/my_run_20260724_192134_results.csv",
  sweep = "C:/Users/dell/Documents/Research Project 2/Simulation_Results/4_Base_model_exploration/Base_Models/Run_20260726_Remaining_Basic_Model_Exploration/AllModels_sweep_20260726_153222_results.csv"
)
sweep_rds_dir <- "C:/Users/dell/Documents/Research Project 2/Modelling/Current_Working_Model/Plotting_Code/Sweep_Summaries"
dir.create(sweep_rds_dir, showWarnings = FALSE, recursive = TRUE)

sweep_metrics <- c("species_richness", "structural_complexity", "turnover",
                   "shannon", "total_cover", "evenness")

# decomposeCombo() (base model + manipulation type) comes from Analysis_Utils.r,
# shared with the summary/stats scripts.
source("C:/Users/dell/Documents/Research Project 2/Modelling/Current_Working_Model/Plotting_Code/Analysis_Utils.r")

# Read a sweep ONCE (fread if data.table is installed, else read.csv), decompose
# the model names, and summarise to:
#   traj     = mean over replicates at each (scenario, timestep)   [396 x 1000]
#   endstate = mean over the last `tail_frac` of the run, per scenario  [396]
# Both are saved to out_dir as <name>_traj.rds / <name>_endstate.rds.
loadSweep <- function(path, name, out_dir = sweep_rds_dir,
                      metrics = sweep_metrics, tail_frac = 0.1) {
  reader <- if (requireNamespace("data.table", quietly = TRUE))
              function(p) as.data.frame(data.table::fread(p))
            else function(p) read.csv(p, stringsAsFactors = FALSE)
  cat("reading", name, "(this is the slow bit for the big file)...\n")
  d <- reader(path)
  cat("  ", nrow(d), "rows,", length(unique(d$scenario)), "scenarios\n")

  d <- cbind(d, decomposeCombo(d$combination))
  facs <- c("scenario", "base_model", "manipulation", "manip_type", "combination",
            "n_species", "individuals", "disturbance_on", "reproduction_on")
  info <- unique(d[, facs])

  traj <- aggregate(d[metrics], by = list(scenario = d$scenario, timestep = d$timestep),
                    FUN = mean, na.rm = TRUE)
  traj <- merge(traj, info, by = "scenario")

  tmax <- max(d$timestep)
  win  <- d[d$timestep >= tmax * (1 - tail_frac), ]           # late window
  endstate <- aggregate(win[metrics], by = list(scenario = win$scenario),
                        FUN = mean, na.rm = TRUE)

  # initial (first timestep) mean per scenario, so plots can show CHANGE (end - init)
  d0   <- d[d$timestep == min(d$timestep), ]
  init <- aggregate(d0[metrics], by = list(scenario = d0$scenario), FUN = mean, na.rm = TRUE)
  names(init)[-1] <- paste0(names(init)[-1], "_init")

  endstate <- merge(merge(endstate, init, by = "scenario"), info, by = "scenario")

  # per-REPLICATE late-window mean (keeps the run-to-run spread for boxplots),
  # carrying the initial value so relative/change modes work per replicate
  repvals <- aggregate(win[metrics], by = list(scenario = win$scenario, replicate = win$replicate),
                       FUN = mean, na.rm = TRUE)
  repvals <- merge(merge(repvals, init, by = "scenario"), info, by = "scenario")

  saveRDS(traj,     file.path(out_dir, paste0(name, "_traj.rds")))
  saveRDS(endstate, file.path(out_dir, paste0(name, "_endstate.rds")))
  saveRDS(repvals,  file.path(out_dir, paste0(name, "_repvals.rds")))
  cat("  saved", name, "_traj/_endstate/_repvals.rds\n")
  invisible(list(traj = traj, endstate = endstate, repvals = repvals))
}

# --- Readable plots from the small end-state summary ------------------------
# `mode`: "relative" colours by end / initial -- the PROPORTION of the starting
#          value retained, so it is comparable across species counts (3->3 and
#          7->7 both = 1 = yellow; 3->1 and 7->2.3 both ~ 0.33 = dark).
#         "change" colours by (end - initial); "absolute" by the raw end value.
# One disturbance state per call (all manipulations/models shown for that slice).

# value each cell/point is coloured by, and its label, per mode
sweepValue <- function(d, metric, mode) {
  init <- d[[paste0(metric, "_init")]]
  if (mode != "absolute" && is.null(init))
    stop("no '", metric, "_init' column - this summary is stale; re-run loadSweep(...) ",
         "to regenerate it (or use mode = 'absolute').")
  switch(mode,
    relative = ifelse(abs(init) < 1e-9, NA, d[[metric]] / init),   # proportion retained
    change   = d[[metric]] - init,
    absolute = d[[metric]],
    stop("mode must be 'relative', 'change' or 'absolute'"))
}
sweepLab <- function(metric, mode) switch(mode,
  relative = paste("proportion of initial", metric, "retained (end / initial)"),
  change   = paste("change in", metric, "(end - initial)"),
  absolute = paste("end-state", metric))

# Boxplots of the per-replicate spread: one box per scenario (manipulation),
# dodged by species count, within each base model. Grey mean +/-1 SD bars over
# each box. Takes the `repvals` summary (per-replicate late-window means).
plotSweepEndstate <- function(repvals, metric = "species_richness",
                              mode = "relative", disturbance = FALSE, sd_bars = TRUE) {
  library(ggplot2)
  d <- repvals[repvals$disturbance_on == disturbance, ]
  d$.val <- sweepValue(d, metric, mode)
  ref <- if (mode == "relative") geom_hline(yintercept = 1, colour = "grey80")
         else if (mode == "change") geom_hline(yintercept = 0, colour = "grey80") else NULL
  ylab  <- sweepLab(metric, mode)
  dodge <- position_dodge(width = 0.8)
  p <- ggplot(d, aes(x = manipulation, y = .val, fill = factor(n_species))) +
    ref +
    geom_boxplot(position = dodge, outlier.size = 0.4, linewidth = 0.3) +
    facet_grid(individuals ~ base_model, scales = "free_x", labeller = label_both) +
    labs(title = sprintf("%s  (disturbance %s)", ylab, if (disturbance) "on" else "off"),
         x = "scenario (manipulation)", y = ylab, fill = "n species",
         caption = paste("Boxes collapse to a median line where replicate variance is zero",
                         "(e.g. Neutral retains richness exactly, Linear pins to 1);",
                         "the coloured point marks the per-species mean.")) +
    theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1),
                       plot.caption = element_text(hjust = 0, colour = "grey40"))
  if (sd_bars)
    p <- p + stat_summary(aes(group = factor(n_species)),
             fun.data = function(x) data.frame(y = mean(x), ymin = mean(x) - sd(x), ymax = mean(x) + sd(x)),
             geom = "errorbar", position = dodge, width = 0.6, colour = "grey25", linewidth = 0.35)
  # A species-coloured mean point, so EVERY species stays visible even when its box
  # collapses to a flat line (zero replicate variance -> geom_boxplot draws no filled
  # rectangle, so its fill colour vanishes). This is why the 5-species (green) box
  # looks "missing" for Neutral (richness retained exactly) and Linear (pinned to 1):
  # the data is there, the box just has no height. The point restores it. The default
  # colour scale matches the fill hue, so the marker colour tracks the box colour.
  p <- p + stat_summary(aes(group = factor(n_species), colour = factor(n_species)),
           fun = mean, geom = "point", position = dodge, size = 1.6,
           stroke = 0, show.legend = FALSE)
  print(p); invisible(p)
}

# Heatmap overview: base model x manipulation TYPE, fill = proportion retained
# (or change/absolute), faceted by species x founders. Cells average sub-variants.
plotSweepHeatmap <- function(endstate, metric = "species_richness",
                             mode = "relative", disturbance = FALSE) {
  library(ggplot2)
  d <- endstate[endstate$disturbance_on == disturbance, ]
  d$.val <- sweepValue(d, metric, mode)
  agg <- aggregate(.val ~ base_model + manip_type + n_species + individuals,
                   data = d, FUN = mean, na.rm = TRUE)
  agg$manip_type <- factor(agg$manip_type, c("none", "size", "growth", "reproduction", "other"))
  # for "relative", pin the scale to 0..1 so yellow always = fully retained
  fillscale <- if (mode == "relative")
                 scale_fill_viridis_c(limits = c(0, 1), oob = scales::squish)
               else scale_fill_viridis_c()
  p <- ggplot(agg, aes(x = manip_type, y = base_model, fill = .val)) +
    geom_tile(colour = "white") +
    facet_grid(individuals ~ n_species, labeller = label_both) +
    fillscale +
    labs(title = sprintf("%s  (disturbance %s)",
                         tools::toTitleCase(sweepLab(metric, mode)),
                         if (disturbance) "on" else "off"),
         x = "manipulation type", y = "base model",
         fill = if (mode == "relative") "retained" else if (mode == "change") "change" else metric) +
    theme_minimal() + theme(axis.text.x = element_text(angle = 30, hjust = 1))
  print(p); invisible(p)
}

# --- Usage ------------------------------------------------------------------
# Read + summarise each file ONCE (the sweep read takes a few minutes; writes .rds):
full_sum  <- loadSweep(sweep_paths["full"],  "full")
sweep_sum <- loadSweep(sweep_paths["sweep"], "sweep")

# Afterwards, reload the tiny summaries instantly instead of re-reading the CSV:
# endstate <- readRDS(file.path(sweep_rds_dir, "sweep_endstate.rds"))

# Readable overviews from the end-state summary. Default mode = "relative":
# proportion of the STARTING value retained, so a model that keeps 3->3 and one
# that keeps 7->7 are both yellow; the darker a cell, the bigger the relative
# collapse. One call per disturbance state.
plotSweepHeatmap(sweep_sum$endstate,  "species_richness")                     # relative, disturbance off
plotSweepHeatmap(sweep_sum$endstate,  "species_richness", disturbance = TRUE) # relative, disturbance on
plotSweepHeatmap(sweep_sum$endstate,  "shannon", mode = "absolute")          # raw end value instead

# Boxplots of the per-replicate spread (from repvals): one box per scenario,
# dodged by species, faceted by founders x base model, with mean +/-1 SD bars.
# richness: relative (proportion of species retained) is meaningful
plotSweepEndstate(sweep_sum$repvals, "species_richness")
plotSweepEndstate(sweep_sum$repvals, "species_richness", disturbance = TRUE)

# all other metrics: absolute end value (relative is undefined/misleading for these)
plotSweepEndstate(sweep_sum$repvals, "evenness",              mode = "absolute")
plotSweepEndstate(sweep_sum$repvals, "evenness",              mode = "absolute", disturbance = TRUE)

plotSweepEndstate(sweep_sum$repvals, "shannon",              mode = "absolute")
plotSweepEndstate(sweep_sum$repvals, "shannon",              mode = "absolute", disturbance = TRUE)

plotSweepEndstate(sweep_sum$repvals, "structural_complexity", mode = "absolute")
plotSweepEndstate(sweep_sum$repvals, "structural_complexity", mode = "absolute", disturbance = TRUE)

plotSweepEndstate(sweep_sum$repvals, "total_cover",          mode = "absolute")
plotSweepEndstate(sweep_sum$repvals, "total_cover",          mode = "absolute", disturbance = TRUE)

plotSweepEndstate(sweep_sum$repvals, "turnover",             mode = "absolute")
plotSweepEndstate(sweep_sum$repvals, "turnover",             mode = "absolute", disturbance = TRUE)




#Heatmap plots -- every metric (absolute) x both disturbance states
plotSweepHeatmap(sweep_sum$endstate, "species_richness")
plotSweepHeatmap(sweep_sum$endstate, "species_richness", disturbance = TRUE)

plotSweepHeatmap(sweep_sum$endstate, "evenness",              mode = "absolute")
plotSweepHeatmap(sweep_sum$endstate, "evenness",              mode = "absolute", disturbance = TRUE)

plotSweepHeatmap(sweep_sum$endstate, "shannon",              mode = "absolute")
plotSweepHeatmap(sweep_sum$endstate, "shannon",              mode = "absolute", disturbance = TRUE)

plotSweepHeatmap(sweep_sum$endstate, "structural_complexity", mode = "absolute")
plotSweepHeatmap(sweep_sum$endstate, "structural_complexity", mode = "absolute", disturbance = TRUE)

plotSweepHeatmap(sweep_sum$endstate, "total_cover",          mode = "absolute")
plotSweepHeatmap(sweep_sum$endstate, "total_cover",          mode = "absolute", disturbance = TRUE)

plotSweepHeatmap(sweep_sum$endstate, "turnover",             mode = "absolute")
plotSweepHeatmap(sweep_sum$endstate, "turnover",             mode = "absolute", disturbance = TRUE)
