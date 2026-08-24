# Colony $size (a cell count) -> % of the reef. Used by reproduction scaling and
# cover metrics. (Competition/growth size effects use cell counts directly.)
coralCoverPercent <- function(coral, reef) {
  100 * coral$size / (nrow(reef) * ncol(reef))
}
