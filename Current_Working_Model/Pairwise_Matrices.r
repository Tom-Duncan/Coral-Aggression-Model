# =============================================================================
#  Pairwise_Matrices.r  -  curated LIBRARY of named pairwise overgrowth matrices.
# -----------------------------------------------------------------------------
#  Every simulation uses ONE pairwise overgrowth matrix M, where
#     M[i, j] = P(species i overgrows species j)   (rows = attacker, cols = occupant)
#  and the DIAGONAL holds the intraspecific (same-species) overgrowth probability.
#
#  This file assigns each matrix a STABLE id (PairMat1, PairMat2, ...) - the
#  `pairmat_id` recorded in the results table. The full matrix + its metadata are
#  kept here so a run's network is always recoverable from its id, never lost to a
#  side-file. Add new matrices by appending to PAIRWISE_MATRICES; ids never change.
#
#  The actual matrices are BUILT with the existing, tested generators in
#  Trait_Combinations.r (rpsMatrix / hierarchyMatrix / neutralMatrix via
#  makeTraitMatrix) so the RPS/linear/neutral logic lives in exactly one place.
# =============================================================================

# Needs the matrix builders (makeTraitMatrix, rpsMatrix, ...). Source them if this
# file is used on its own; harmless if they are already loaded.
if (!exists("makeTraitMatrix")) {
  source(file.path("Current_Working_Model", "Trait_Combinations.r"))
}


# --- Shared parameters for the matrices defined below -----------------------
#  bias  = the stronger species' win probability in a decided pair (0.5-1).
#  intra = the same-species (diagonal) overgrowth probability.
#  These are FIXED into each matrix here, so a PairMat id denotes one concrete
#  network. Change them per matrix in makePairMat() if a matrix needs its own.
PAIRMAT_BIAS  <- 0.9
PAIRMAT_INTRA <- 0.5


# --- Build one library entry ------------------------------------------------
#  Returns a list: the id, the dynamics label, species count, the bias/intra used,
#  free-text notes, and the built matrix. `dynamics` is one of the makeTraitMatrix
#  types: "rps" (intransitive cycle, odd n), "hierarchy"/"linear" (transitive rank),
#  or "neutral" (no overgrowth).
makePairMat <- function(id, dynamics, n, bias = PAIRMAT_BIAS,
                        intra = PAIRMAT_INTRA, notes = "") {
  M <- makeTraitMatrix(tolower(dynamics), n, bias, intra)
  list(id = id, dynamics = dynamics, n_species = n,
       bias = bias, intra = intra, notes = notes, matrix = M)
}


# =============================================================================
#  THE LIBRARY  -  append new matrices here; existing ids must never change.
# =============================================================================
PAIRWISE_MATRICES <- list(

  # ---- RPS dynamics (intransitive rock-paper-scissors cycle) ----------------
  PairMat1 = makePairMat("PairMat1", "RPS", 3,
                         notes = "RPS cycle, 3 species (classic rock-paper-scissors)"),
  PairMat2 = makePairMat("PairMat2", "RPS", 5,
                         notes = "RPS cycle, 5 species (each beats the next 2)"),
  PairMat3 = makePairMat("PairMat3", "RPS", 7,
                         notes = "RPS cycle, 7 species (each beats the next 3)")
)


# =============================================================================
#  ACCESSORS
# =============================================================================

# Look up one matrix entry by id (stops with a clear error on an unknown id).
getPairMatrix <- function(id) {
  if (!id %in% names(PAIRWISE_MATRICES)) {
    stop("Unknown pairwise matrix id: ", id,
         ". Known ids: ", paste(names(PAIRWISE_MATRICES), collapse = ", "))
  }
  PAIRWISE_MATRICES[[id]]
}

# A one-row-per-matrix summary of the whole library (no matrices, just metadata).
listPairMatrices <- function() {
  do.call(rbind, lapply(PAIRWISE_MATRICES, function(p) {
    data.frame(id = p$id, dynamics = p$dynamics, n_species = p$n_species,
               bias = p$bias, intra = p$intra, notes = p$notes,
               stringsAsFactors = FALSE)
  }))
}

# Pretty-print one matrix (its metadata then the M grid) for eyeballing.
printPairMatrix <- function(id) {
  p <- getPairMatrix(id)
  cat(sprintf("%s | %s | %d species | bias %.2f | intra %.2f\n",
              p$id, p$dynamics, p$n_species, p$bias, p$intra))
  if (nzchar(p$notes)) cat(" ", p$notes, "\n", sep = "")
  print(round(p$matrix, 3))
  invisible(p)
}
