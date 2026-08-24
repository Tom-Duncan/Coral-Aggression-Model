# Library of named pairwise overgrowth matrices M, where
#   M[i, j] = P(species i overgrows species j); diagonal = intraspecific prob.
# Each matrix has a stable id (PairMat1, ...) recorded in results, so a run's
# network is recoverable from its id. Matrices are built by the generators in
# Trait_Combinations.r. Append new matrices; ids never change.

if (!exists("makeTraitMatrix")) {
  source(file.path("Model", "Trait_Combinations.r"))
}


# bias = stronger species' win probability (0.5-1); intra = diagonal probability.
PAIRMAT_BIAS  <- 0.9
PAIRMAT_INTRA <- 0.5


# One library entry. dynamics: "rps" / "hierarchy"(="linear") / "neutral".
makePairMat <- function(id, dynamics, n, bias = PAIRMAT_BIAS,
                        intra = PAIRMAT_INTRA, notes = "") {
  M <- makeTraitMatrix(tolower(dynamics), n, bias, intra)
  list(id = id, dynamics = dynamics, n_species = n,
       bias = bias, intra = intra, notes = notes, matrix = M)
}


# The library - append here; existing ids must never change.
PAIRWISE_MATRICES <- list(

  PairMat1 = makePairMat("PairMat1", "RPS", 3,
                         notes = "RPS cycle, 3 species (classic rock-paper-scissors)"),
  PairMat2 = makePairMat("PairMat2", "RPS", 5,
                         notes = "RPS cycle, 5 species (each beats the next 2)"),
  PairMat3 = makePairMat("PairMat3", "RPS", 7,
                         notes = "RPS cycle, 7 species (each beats the next 3)")
)


# ---- Accessors --------------------------------------------------------------

# Look up one entry by id.
getPairMatrix <- function(id) {
  if (!id %in% names(PAIRWISE_MATRICES)) {
    stop("Unknown pairwise matrix id: ", id,
         ". Known ids: ", paste(names(PAIRWISE_MATRICES), collapse = ", "))
  }
  PAIRWISE_MATRICES[[id]]
}

# One-row-per-matrix metadata summary.
listPairMatrices <- function() {
  do.call(rbind, lapply(PAIRWISE_MATRICES, function(p) {
    data.frame(id = p$id, dynamics = p$dynamics, n_species = p$n_species,
               bias = p$bias, intra = p$intra, notes = p$notes,
               stringsAsFactors = FALSE)
  }))
}

# Print one matrix's metadata + grid.
printPairMatrix <- function(id) {
  p <- getPairMatrix(id)
  cat(sprintf("%s | %s | %d species | bias %.2f | intra %.2f\n",
              p$id, p$dynamics, p$n_species, p$bias, p$intra))
  if (nzchar(p$notes)) cat(" ", p$notes, "\n", sep = "")
  print(round(p$matrix, 3))
  invisible(p)
}
