# Size helper for coral species
#
# Colony $size is a CELL COUNT (the number of cells occupied). This is the single
# place cells are turned into a % of the reef, used by reproduction scaling and the
# reported cover metrics. (Competition and growth size effects work on cell counts
# directly via a log ratio, so they no longer need % cover or size categories.)
coralCoverPercent <- function(coral, reef) {
  100 * coral$size / (nrow(reef) * ncol(reef))
}
