# ============================================================================
#  Fecundity dose-response figures  (base R, no packages)
#    #2  Linear competition-colonization trade-off (richness vs fecundity)
#    #4  Colony-size structure vs fecundity (the recruitment mechanism)
# ============================================================================

# --- point this at your fecundity results file --------------------------------
RESULTS <- "Reproduction_Fecundity_20260731_222545_results.rds"
res <- readRDS(RESULTS)

se <- function(x) { x <- x[is.finite(x)]; if (length(x) < 2) return(NA); sd(x)/sqrt(length(x)) }
FEC <- sort(unique(res$repro_base))

# ============================================================================
#  #2  Linear trade-off: Normal vs Even vs Opposite, richness vs fecundity
#      faceted (2x2) by species x disturbance
# ============================================================================
tr <- subset(res, phase == "final" &
             combination %in% c("Linear_reproductionNormal",
                                "Linear_reproductionEven",
                                "Linear_reproductionOpposite"))
lev  <- c("Linear_reproductionNormal", "Linear_reproductionEven", "Linear_reproductionOpposite")
labs <- c("Normal (dominant breeds fastest)", "Even", "Opposite (subordinate breeds fastest)")
cols <- c("#d73027", "#7f7f7f", "#1a9850"); names(cols) <- lev

png("fig2_linear_tradeoff.png", width = 1100, height = 850)
op <- par(mfrow = c(2, 2), mar = c(4.2, 4.2, 2.6, 1), oma = c(0, 0, 2.4, 0))
ylim <- c(0.8, max(tr$n_species))
for (nsp in c(3, 7)) for (dst in c(FALSE, TRUE)) {
  sub <- subset(tr, n_species == nsp & disturbance_on == dst)
  plot(NA, xlim = range(FEC), ylim = ylim, log = "x",
       xlab = "Base fecundity  (repro_base)", ylab = "Final species richness",
       main = sprintf("%d species  |  disturbance %s", nsp, ifelse(dst, "ON", "OFF")))
  for (mo in lev) {
    s <- subset(sub, combination == mo)
    m  <- tapply(s$species_richness, s$repro_base, mean)
    er <- tapply(s$species_richness, s$repro_base, se)
    x  <- as.numeric(names(m))
    arrows(x, m - er, x, m + er, angle = 90, code = 3, length = 0.03, col = cols[mo])
    lines(x, m, col = cols[mo], lwd = 2.5); points(x, m, col = cols[mo], pch = 19)
  }
  if (nsp == 3 && !dst)
    legend("topright", legend = labs, col = cols, lwd = 2.5, pch = 19, cex = 0.75, bty = "n")
}
mtext("Linear hierarchy: does the competition-colonization trade-off rescue coexistence?",
      outer = TRUE, font = 2, cex = 1.05)
par(op); dev.off(); cat("saved fig2_linear_tradeoff.png\n")

# ============================================================================
#  #4a  Mean colony size vs fecundity (all 9 models), faceted by species x disturbance
# ============================================================================
fin  <- subset(res, phase == "final")
mods <- sort(unique(fin$combination))
mcol <- hcl.colors(length(mods), "Dark 3"); names(mcol) <- mods

png("fig4a_mean_colony_size.png", width = 1100, height = 850)
op <- par(mfrow = c(2, 2), mar = c(4.2, 4.2, 2.6, 1), oma = c(0, 0, 2.4, 0))
ylim <- range(tapply(fin$mean_colony_size, list(fin$combination, fin$repro_base,
                     fin$n_species, fin$disturbance_on), mean, na.rm = TRUE), na.rm = TRUE)
for (nsp in c(3, 7)) for (dst in c(FALSE, TRUE)) {
  sub <- subset(fin, n_species == nsp & disturbance_on == dst)
  plot(NA, xlim = range(FEC), ylim = ylim, log = "x",
       xlab = "Base fecundity", ylab = "Mean colony size (cells)",
       main = sprintf("%d species  |  disturbance %s", nsp, ifelse(dst, "ON", "OFF")))
  for (mo in mods) {
    s <- subset(sub, combination == mo)
    m <- tapply(s$mean_colony_size, s$repro_base, mean, na.rm = TRUE)
    x <- as.numeric(names(m)); lines(x, m, col = mcol[mo], lwd = 2); points(x, m, col = mcol[mo], pch = 19)
  }
  if (nsp == 3 && !dst)
    legend("topright", legend = mods, col = mcol, lwd = 2, cex = 0.55, bty = "n")
}
mtext("Colony size shrinks as fecundity rises (more recruits, smaller colonies)",
      outer = TRUE, font = 2, cex = 1.05)
par(op); dev.off(); cat("saved fig4a_mean_colony_size.png\n")

# ============================================================================
#  #4b  Colony size-class distribution vs fecundity (stacked proportions)
#       averaged over the 9 models; one panel per species x disturbance
# ============================================================================
png("fig4b_size_distribution.png", width = 1100, height = 850)
op <- par(mfrow = c(2, 2), mar = c(4.2, 4.2, 2.6, 1), oma = c(0, 0, 2.4, 0))
sccol <- c("#fee08b", "#fc8d59", "#b30000")
for (nsp in c(3, 7)) for (dst in c(FALSE, TRUE)) {
  sub <- subset(fin, n_species == nsp & disturbance_on == dst)
  a  <- aggregate(cbind(n_small, n_medium, n_large) ~ repro_base, sub, mean, na.rm = TRUE)
  pr <- as.matrix(a[, c("n_small", "n_medium", "n_large")]); pr <- pr / rowSums(pr)
  barplot(t(pr), names.arg = a$repro_base, col = sccol, ylim = c(0, 1),
          xlab = "Base fecundity", ylab = "Proportion of colonies",
          main = sprintf("%d species  |  disturbance %s", nsp, ifelse(dst, "ON", "OFF")))
  if (nsp == 3 && !dst)
    legend("topright", legend = c("small", "medium", "large"), fill = sccol, cex = 0.8, bty = "n")
}
mtext("Community shifts toward small colonies as fecundity rises",
      outer = TRUE, font = 2, cex = 1.05)
par(op); dev.off(); cat("saved fig4b_size_distribution.png\n")
