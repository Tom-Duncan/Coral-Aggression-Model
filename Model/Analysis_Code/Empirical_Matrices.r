# =============================================================================
#  EMPIRICAL INTERACTION MATRICES - overgrowth-probability matrices built from
#  published/observed coral-aggression data, for use as model communities.
#  (The Connell network in the paper uses the Heron Island data referred to here.)
#
#  1. LOGAN 1984 - Bermuda (Atlantic), 17 species, aquarium + field aggression
#     trials (Logan A. (1984) Interspecific aggression in hermatypic corals from
#     Bermuda. Coral Reefs 3:131-138, Table 1 p.133).
#  2. DAI 1990   - Nanwan Bay, Taiwan (Indo-Pacific), 19 scleractinians, field
#     observations (Dai C.F. (1990) Interspecific competition in Taiwanese
#     corals with special reference to interactions between alcyonaceans and
#     scleractinians. Mar Ecol Prog Ser 60:291-297, Table 1). The 9 alcyonacean
#     (soft-coral) taxa in Dai's table are EXCLUDED (different functional guild).
#
#  Convention (matches the engine): M[i, j] = P(row species i overgrows column
#  species j GIVEN a decisive contest); complementary, M[i,j] + M[j,i] = 1;
#  diagonal = intraspecific (0.5 default).
#
#  Excel exports (kept as SEPARATE workbooks, one per study):
#    writeLoganExcel()  ->  Interaction_Matrices/Logan_1984_interaction_matrix.xlsx
#    writeDaiExcel()    ->  Interaction_Matrices/Dai_1990_interaction_matrix.xlsx
# =============================================================================

interaction_xlsx_dir <- file.path("Model", "Interaction_Matrices")   # relative to the repo root


# =============================================================================
#  1a. LOGAN 1984 - GRADED matrix (primary)
#  Cell values are Laplace-smoothed win fractions among decisive contests,
#    p = (B + 1) / (B + E + 2)
#  where B = conclusive results in the consensus direction and E = results
#  inconsistent with the consensus (Table 1's corner counts). Per pair the
#  better-replicated of the aquarium / field cell was used; pairs without
#  consensus (dots / '?' in Logan) are 0.50. Direction-only cells whose counts
#  were not legible carry the modal-cell default 0.75 (see transcription notes).
# =============================================================================

logan_species <- c("Ds","Dl","Ma","Mc","Pa","Pp","Ff","Md","Mm","Oc",
                   "Is","Me","Sm","Dk","Sy","Si","Ag")

#Full binomials, same order (Bermudian hermatypic corals, Logan's Table 1)
logan_names <- c(Ds = "Diploria strigosa",        Dl = "Diploria labyrinthiformis",
                 Ma = "Montastrea annularis",      Mc = "Montastrea cavernosa",
                 Pa = "Porites astreoides",        Pp = "Porites porites",
                 Ff = "Favia fragum",              Md = "Madracis decactis",
                 Mm = "Madracis mirabilis",        Oc = "Oculina spp.",
                 Is = "Isophyllia sinuosa",        Me = "Meandrina meandrites",
                 Sm = "Stephanocoenia michelinii", Dk = "Dichocoenia stokesi",
                 Sy = "Scolymia cubensis",         Si = "Siderastrea spp.",
                 Ag = "Agaricia fragilis")

M_logan_graded <- rbind(
  Ds = c(.50,.14,.06,.08,.98,.80,.85,.75,.80,.75,.25,.25,.71,.25,.50,.80,.75),
  Dl = c(.86,.50,.08,.25,.98,.80,.75,.80,.88,.75,.25,.33,.75,.33,.25,.75,.80),
  Ma = c(.94,.92,.50,.39,.98,.75,.75,.80,.75,.75,.25,.33,.71,.25,.25,.67,.75),
  Mc = c(.92,.75,.61,.50,.94,.75,.67,.67,.67,.67,.33,.33,.83,.50,.50,.67,.67),
  Pa = c(.02,.02,.02,.06,.50,.86,.75,.75,.75,.50,.13,.25,.78,.25,.25,.75,.75),
  Pp = c(.20,.20,.25,.25,.14,.50,.33,.67,.33,.17,.17,.17,.67,.25,.50,.75,.75),
  Ff = c(.15,.25,.25,.33,.25,.67,.50,.80,.50,.75,.12,.50,.25,.25,.17,.75,.67),
  Md = c(.25,.20,.20,.33,.25,.33,.20,.50,.25,.67,.17,.33,.33,.25,.20,.75,.86),
  Mm = c(.20,.12,.25,.33,.25,.67,.50,.75,.50,.50,.09,.33,.50,.25,.20,.75,.75),
  Oc = c(.25,.25,.25,.33,.50,.83,.25,.33,.50,.50,.17,.25,.50,.50,.33,.50,.75),
  Is = c(.75,.75,.75,.67,.87,.83,.88,.83,.91,.83,.50,.50,.80,.50,.75,.67,.83),
  Me = c(.75,.67,.67,.67,.75,.83,.50,.67,.67,.75,.50,.50,.75,.75,.50,.67,.50),
  Sm = c(.29,.25,.29,.17,.22,.33,.75,.67,.50,.50,.20,.25,.50,.25,.25,.67,.50),
  Dk = c(.75,.67,.75,.50,.75,.75,.75,.75,.75,.50,.50,.25,.75,.50,.50,.80,.80),
  Sy = c(.50,.75,.75,.50,.75,.50,.83,.80,.80,.67,.25,.50,.75,.50,.50,.80,.80),
  Si = c(.20,.25,.33,.33,.25,.25,.25,.25,.25,.50,.33,.33,.33,.20,.20,.50,.25),
  Ag = c(.25,.20,.25,.33,.25,.25,.33,.14,.25,.25,.17,.50,.50,.20,.20,.75,.50))
colnames(M_logan_graded) <- logan_species


# =============================================================================
#  1b. LOGAN 1984 - CATEGORICAL matrix (alternative)
#  Same consensus directions, but fixed bias: winner 0.9 / loser 0.1, no
#  consensus 0.5 (the model's standard bias convention). Use when you want the
#  hierarchy STRUCTURE without Logan's win-rate magnitudes.
# =============================================================================
M_logan_categorical <- rbind(
  Ds = c(.5,.1,.1,.1,.9,.9,.9,.9,.9,.9,.1,.1,.9,.1,.5,.9,.9),
  Dl = c(.9,.5,.1,.1,.9,.9,.9,.9,.9,.9,.1,.1,.9,.1,.1,.9,.9),
  Ma = c(.9,.9,.5,.1,.9,.9,.9,.9,.9,.9,.1,.1,.9,.1,.1,.9,.9),
  Mc = c(.9,.9,.9,.5,.9,.9,.9,.9,.9,.9,.1,.1,.9,.5,.5,.9,.9),
  Pa = c(.1,.1,.1,.1,.5,.9,.9,.9,.9,.5,.1,.1,.9,.1,.1,.9,.9),
  Pp = c(.1,.1,.1,.1,.1,.5,.1,.9,.1,.1,.1,.1,.9,.1,.5,.9,.9),
  Ff = c(.1,.1,.1,.1,.1,.9,.5,.9,.5,.9,.1,.5,.1,.1,.1,.9,.9),
  Md = c(.1,.1,.1,.1,.1,.1,.1,.5,.1,.9,.1,.1,.1,.1,.1,.9,.9),
  Mm = c(.1,.1,.1,.1,.1,.9,.5,.9,.5,.5,.1,.1,.5,.1,.1,.9,.9),
  Oc = c(.1,.1,.1,.1,.5,.9,.1,.1,.5,.5,.1,.1,.5,.5,.1,.5,.9),
  Is = c(.9,.9,.9,.9,.9,.9,.9,.9,.9,.9,.5,.5,.9,.5,.9,.9,.9),
  Me = c(.9,.9,.9,.9,.9,.9,.5,.9,.9,.9,.5,.5,.9,.9,.5,.9,.5),
  Sm = c(.1,.1,.1,.1,.1,.1,.9,.9,.5,.5,.1,.1,.5,.1,.1,.9,.5),
  Dk = c(.9,.9,.9,.5,.9,.9,.9,.9,.9,.5,.5,.1,.9,.5,.5,.9,.9),
  Sy = c(.5,.9,.9,.5,.9,.5,.9,.9,.9,.9,.1,.5,.9,.5,.5,.9,.9),
  Si = c(.1,.1,.1,.1,.1,.1,.1,.1,.1,.5,.1,.1,.1,.1,.1,.5,.1),
  Ag = c(.1,.1,.1,.1,.1,.1,.1,.1,.1,.1,.1,.5,.5,.1,.1,.9,.5))
colnames(M_logan_categorical) <- logan_species

#Transcription caveats (both Logan matrices):
# - Aquarium half of Table 1 primary; field used where better replicated.
# - Field-inversion pairs (aquarium direction kept; Logan: ~30% of pairs invert
#   partially/completely in field): Ma-Me, Ma-Dk, Mc-Me, Pa-Ff, Pa-Mm, Sm-Sy.
# - Cells transcribed from a page scan; low-data cells (rare species block)
#   carry defaults - verify against a zoomed original before publication use.


# =============================================================================
#  2. DAI 1990 - CATEGORICAL matrix (Indo-Pacific; 19 scleractinians)
# -----------------------------------------------------------------------------
#  Dai's Table 1 records each pair's outcome as SYMBOLS only (arrow toward the
#  winner; solid = direct aggression, dotted = overgrowth, wavy = allelopathy;
#  double line = stand-off; blank = never observed), each based on >= 2 repeated
#  observations - NO counts. Probabilities are therefore ASSIGNED categorically:
#     consistent winner                     -> 0.9 / 0.1
#     winner with stand-offs also observed  -> 0.7 / 0.3
#     stand-off only                        -> 0.5   (observed tie)
#     blank                                 -> 0.5   (ND: no data - flagged)
#  The flags matrix (dai_flags) distinguishes observed ties ("standoff"), no
#  data ("ND"), and uncertain transcription reads ("check").
# =============================================================================

dai_species <- c("ACHY","MOVR","MOUN","PCVR","SEHY","STPI","PRAU","PRLI","PRLU",
                 "PSSP","FASP","FTAB","PTLA","MEAM","ECAS","MYEL","GAFA","ACEC","SYRA")

dai_names <- c(ACHY = "Acropora hyacinthus",    MOVR = "Montipora verrucosa",
               MOUN = "Montipora undata",       PCVR = "Pocillopora verrucosa",
               SEHY = "Seriatopora hystrix",    STPI = "Stylophora pistillata",
               PRAU = "Porites australiensis",  PRLI = "Porites lichen",
               PRLU = "Porites lutea",          PSSP = "Pachyseris speciosa",
               FASP = "Favia speciosa",         FTAB = "Favites abdita",
               PTLA = "Platygyra lamellina",    MEAM = "Merulina ampliata",
               ECAS = "Echinophyllia aspera",   MYEL = "Mycedium elephantotus",
               GAFA = "Galaxea fascicularis",   ACEC = "Acanthastrea echinata",
               SYRA = "Symphyllia radians")

M_dai_categorical <- rbind(
  ACHY = c(.5,.1,.9,.1,.1,.1,.9,.9,.9,.5,.1,.1,.3,.1,.1,.1,.1,.1,.1),
  MOVR = c(.9,.5,.9,.5,.9,.9,.5,.5,.1,.5,.1,.1,.1,.1,.1,.1,.1,.1,.3),
  MOUN = c(.1,.1,.5,.5,.5,.5,.9,.9,.9,.5,.5,.5,.1,.1,.5,.1,.5,.5,.1),
  PCVR = c(.9,.5,.5,.5,.5,.5,.1,.5,.5,.5,.5,.5,.1,.1,.5,.1,.1,.1,.1),
  SEHY = c(.9,.1,.5,.5,.5,.5,.5,.5,.5,.5,.5,.5,.1,.1,.5,.1,.1,.1,.1),
  STPI = c(.9,.1,.5,.5,.5,.5,.5,.5,.5,.5,.5,.1,.1,.1,.5,.1,.1,.1,.3),
  PRAU = c(.1,.5,.1,.9,.5,.5,.5,.5,.5,.5,.5,.5,.5,.1,.1,.1,.1,.1,.3),
  PRLI = c(.1,.5,.1,.5,.5,.5,.5,.5,.5,.5,.5,.5,.5,.1,.1,.1,.1,.1,.1),
  PRLU = c(.1,.9,.1,.5,.5,.5,.5,.5,.5,.5,.5,.5,.5,.1,.5,.1,.1,.1,.1),
  PSSP = c(.5,.5,.5,.5,.5,.5,.5,.5,.5,.5,.5,.5,.5,.5,.5,.1,.5,.5,.5),
  FASP = c(.9,.9,.5,.5,.5,.5,.5,.5,.5,.5,.5,.5,.5,.5,.5,.1,.1,.5,.3),
  FTAB = c(.9,.9,.5,.5,.5,.9,.5,.5,.5,.5,.5,.5,.5,.5,.5,.1,.1,.5,.1),
  PTLA = c(.7,.9,.9,.9,.9,.9,.5,.5,.5,.5,.5,.5,.5,.5,.5,.1,.1,.5,.3),
  MEAM = c(.9,.9,.9,.9,.9,.9,.9,.9,.9,.5,.5,.5,.5,.5,.9,.5,.5,.5,.5),
  ECAS = c(.9,.9,.5,.5,.5,.5,.9,.9,.5,.5,.5,.5,.5,.1,.5,.5,.5,.5,.5),
  MYEL = c(.9,.9,.9,.9,.9,.9,.9,.9,.9,.9,.9,.9,.9,.5,.5,.5,.5,.5,.5),
  GAFA = c(.9,.9,.5,.9,.9,.9,.9,.9,.9,.5,.9,.9,.9,.5,.5,.5,.5,.5,.5),
  ACEC = c(.9,.9,.5,.9,.9,.9,.9,.9,.9,.5,.5,.5,.5,.5,.5,.5,.5,.5,.5),
  SYRA = c(.9,.7,.9,.9,.9,.7,.7,.9,.9,.5,.7,.9,.7,.5,.5,.5,.5,.5,.5))
colnames(M_dai_categorical) <- dai_species

#Flags per pair (symmetric): "" observed winner; "standoff" observed tie;
#"mixed" winner + stand-offs (0.7 cells); "ND" never observed; "check" my
#transcription uncertain - verify against a zoomed original.
dai_flags <- matrix("", 19, 19, dimnames = list(dai_species, dai_species))
.df <- function(i, j, v) { dai_flags[i, j] <<- v; dai_flags[j, i] <<- v }
for (p in list(c(1,10), c(2,7), c(2,8), c(4,5), c(4,6), c(4,12), c(5,8), c(5,9),
               c(5,12), c(6,11), c(7,13), c(15,19), c(16,19), c(17,18), c(18,19)))
  .df(p[1], p[2], "standoff")
for (p in list(c(13,1), c(19,2), c(19,6), c(19,7), c(19,11), c(19,13)))
  .df(p[1], p[2], "mixed")
for (p in list(c(1,5), c(1,6), c(4,8), c(4,11), c(11,14), c(14,16), c(14,17),
               c(14,19), c(15,16), c(15,17), c(15,18), c(16,17), c(17,19)))
  .df(p[1], p[2], "check")
for (i in 1:18) for (j in (i + 1):19)
  if (M_dai_categorical[i, j] == .5 && dai_flags[i, j] == "") .df(i, j, "ND")
rm(.df)

# --- The unique largest COMPLETE submatrix (no ND, no uncertain cells) --------
# Found by exact max-clique search on the data-presence graph: 7 species, every
# one of the 21 pairs observed with a confident read. Exactly matches the
# model's Sp1..Sp7 output limit.
dai7_species <- c("ACHY", "MOVR", "MOUN", "PRAU", "PTLA", "MYEL", "SYRA")
M_dai7 <- M_dai_categorical[dai7_species, dai7_species]


# =============================================================================
#  EXCEL EXPORTS - one self-documenting workbook per study (NOT combined)
# =============================================================================

#matrix -> data.frame with the species codes as the first column
.matSheet <- function(M) data.frame(Species = rownames(M), as.data.frame(M),
                                    check.names = FALSE)
.keySheet <- function(nm) data.frame(Code = names(nm), Species = unname(nm))

writeLoganExcel <- function(out_dir = interaction_xlsx_dir) {
  if (!requireNamespace("writexl", quietly = TRUE)) stop("install.packages('writexl')")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  notes <- data.frame(Note = c(
    "Source: Logan A. (1984) Interspecific aggression in hermatypic corals from Bermuda. Coral Reefs 3:131-138, Table 1 (p.133). Region: Bermuda (Atlantic).",
    "Convention: cell [row, column] = probability the ROW species overgrows the COLUMN species, given a decisive contest. M[i,j] + M[j,i] = 1; diagonal 0.5.",
    "Sheet 'Graded': Laplace-smoothed win fractions among decisive contests, p = (B+1)/(B+E+2); B = conclusive results in the consensus direction, E = inconsistent results. Better-replicated of aquarium/field cell used per pair. No-consensus pairs = 0.50; count-illegible direction-only cells = 0.75 default.",
    "Sheet 'Categorical': same directions with fixed bias (winner 0.9 / loser 0.1 / no consensus 0.5).",
    "Field-inversion pairs (aquarium direction kept): Ma-Me, Ma-Dk, Mc-Me, Pa-Ff, Pa-Mm, Sm-Sy.",
    "Caveat: transcribed from a page scan; verify low-data cells against a zoomed original before publication use."))
  path <- file.path(out_dir, "Logan_1984_interaction_matrix.xlsx")
  writexl::write_xlsx(list(Graded      = .matSheet(M_logan_graded),
                           Categorical = .matSheet(M_logan_categorical),
                           Species_key = .keySheet(logan_names),
                           Notes       = notes), path = path)
  cat("Wrote:", path, "\n")
  invisible(path)
}

writeDaiExcel <- function(out_dir = interaction_xlsx_dir) {
  if (!requireNamespace("writexl", quietly = TRUE)) stop("install.packages('writexl')")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  notes <- data.frame(Note = c(
    "Source: Dai C.F. (1990) Interspecific competition in Taiwanese corals with special reference to interactions between alcyonaceans and scleractinians. Mar Ecol Prog Ser 60:291-297, Table 1. Region: Nanwan Bay, Taiwan (Indo-Pacific).",
    "Scleractinians only (19 taxa); Dai's 9 alcyonacean (soft coral) taxa excluded as a different functional guild.",
    "Convention: cell [row, column] = probability the ROW species overgrows the COLUMN species, given a decisive contest. M[i,j] + M[j,i] = 1; diagonal 0.5.",
    "Dai's table records outcomes as symbols only (arrow toward the winner; solid = direct aggression, dotted = overgrowth, wavy = allelopathy; double line = stand-off; blank = never observed), each based on >=2 repeated observations. No counts exist, so probabilities are ASSIGNED categorically, not calculated:",
    "  consistent winner = 0.9/0.1; winner with stand-offs also observed = 0.7/0.3; stand-off = 0.5; blank = 0.5 flagged ND.",
    "Sheet 'Flags': standoff = observed tie; mixed = winner + stand-offs (0.7 cells); ND = pair never observed (no data - NOT an observed tie); check = transcription uncertain, verify against a zoomed original.",
    "Sheet 'Complete_7sp': the unique largest complete submatrix (max-clique search) - 7 species, all 21 pairs observed with confident reads: ACHY MOVR MOUN PRAU PTLA MYEL SYRA.",
    "Assigned magnitudes (0.9/0.7) are modelling conventions - state in methods and sensitivity-test."))
  path <- file.path(out_dir, "Dai_1990_interaction_matrix.xlsx")
  writexl::write_xlsx(list(Categorical  = .matSheet(M_dai_categorical),
                           Complete_7sp = .matSheet(M_dai7),
                           Flags        = .matSheet(dai_flags),
                           Species_key  = .keySheet(dai_names),
                           Notes        = notes), path = path)
  cat("Wrote:", path, "\n")
  invisible(path)
}

# --- Usage -------------------------------------------------------------------
# writeLoganExcel()    # Bermuda workbook  (Graded, Categorical, Species_key, Notes)
# writeDaiExcel()      # Taiwan workbook   (Categorical, Complete_7sp, Flags, Species_key, Notes)
# M_dai7               # engine-ready 7x7 Indo-Pacific matrix for buildManualTraits()


# =============================================================================
#  3. OBSERVED CONTEST COMPILATION (Precoda et al. 2017 source data)
# -----------------------------------------------------------------------------
#  coralInteractions_combined.csv = the digitised observed-contest records the
#  Precoda model was trained/tested on: one row per observed interaction,
#  outcome coded from species1's perspective (1 win / -1 loss / 0 draw), with
#  the source reference number and location. 117 species, 16 published sources
#  (Logan 1984 = ref 6, Dai 1990 = ref 7, Connell's Heron Island = ref 1, ...).
#
#  buildObservedPairs() computes per-pair WIN PROBABILITIES from these counts:
#    * reference 17 (Ferriz-Dominguez & Horta-Puga, Veracruz) EXCLUDED
#    * only pairs with MORE THAN `min_obs - 1` recorded interactions kept
#      (default min_obs = 3, i.e. > 2 instances, per analysis decision)
#    * p(A beats B) = (wins_A + 1) / (wins_A + wins_B + 2)  - Laplace-smoothed,
#      decisive contests only (draws reported but excluded from the
#      denominator), the same convention as the Logan graded matrix.
#  observedMatrix() then assembles any species subset into an engine-ready
#  matrix from those filtered pairs (unfiltered/unobserved pairs -> NA).
# =============================================================================

#The observed records live inside the original download zip (training + test CSVs);
#they are combined on read, so no intermediate file needs to be kept on disk.
#Point this at your local copy of the Precoda et al. (2017) contest-data download.
observed_interactions_zip <- "coralcompetition-inputOutputCode.zip"

loadObservedInteractions <- function(zip = observed_interactions_zip) {
  tr <- read.csv(unz(zip, "coralInteractions-Training.csv"), stringsAsFactors = FALSE)
  te <- read.csv(unz(zip, "coralInteractions-Test.csv"),     stringsAsFactors = FALSE)
  tr$set <- "training"; te$set <- "test"
  rbind(tr, te)
}

#Per-pair win-probability table from the observed records.
#  min_obs      : minimum recorded interactions per pair (3 = "more than 2")
#  exclude_refs : source reference numbers to drop (17 = Veracruz study)
buildObservedPairs <- function(d = loadObservedInteractions(),
                               min_obs = 3, exclude_refs = 17) {
  d <- d[!(d$reference %in% exclude_refs), ]

  #unordered pair key (alphabetical), with outcome re-signed to that order
  a <- pmin(d$species1, d$species2); b <- pmax(d$species1, d$species2)
  out <- ifelse(d$species1 == a, d$outcome, -d$outcome)

  key <- paste(a, b, sep = " | ")
  winsA <- tapply(out ==  1, key, sum)
  winsB <- tapply(out == -1, key, sum)
  draws <- tapply(out ==  0, key, sum)
  n     <- tapply(out, key, length)

  pairs <- data.frame(
    speciesA = sub(" \\|.*$", "", names(n)),
    speciesB = sub("^.*\\| ", "", names(n)),
    n = as.integer(n), winsA = as.integer(winsA),
    winsB = as.integer(winsB), draws = as.integer(draws), row.names = NULL)
  pairs <- pairs[pairs$n >= min_obs, ]

  #Laplace-smoothed win probability among decisive contests (A's perspective)
  pairs$pA <- (pairs$winsA + 1) / (pairs$winsA + pairs$winsB + 2)
  pairs[order(-pairs$n), ]
}

#Assemble an engine-ready matrix for a species subset from the filtered pairs.
#Pairs absent from the filtered table -> NA (set na_fill = 0.5 to run anyway).
observedMatrix <- function(species, pairs = buildObservedPairs(), na_fill = NA_real_) {
  n <- length(species)
  M <- matrix(na_fill, n, n, dimnames = list(species, species))
  sub <- pairs[pairs$speciesA %in% species & pairs$speciesB %in% species, ]
  for (r in seq_len(nrow(sub))) {
    M[sub$speciesA[r], sub$speciesB[r]] <- sub$pA[r]
    M[sub$speciesB[r], sub$speciesA[r]] <- 1 - sub$pA[r]
  }
  diag(M) <- 0.5
  M
}

# --- Usage -------------------------------------------------------------------
# obs <- buildObservedPairs()                     # filtered pair table (n > 2, no ref 17)
# head(obs, 20)                                   # best-replicated pairs first
# M <- observedMatrix(c("Galaxea fascicularis", "Porites lutea", ...), obs)


# =============================================================================
#  4. HERON ISLAND (GBR, Indo-Pacific) - three site tables combined, 10 species
# -----------------------------------------------------------------------------
#  Source: three site-level tables of long-term Heron Island interaction records
#  (Connell-programme reef-flat data; supplied as page images), cells =
#  Wins/Losses/Ties from the row species' perspective, combined across the three
#  sites. The 10 species below are the LARGEST subset with every pair observed
#  (exact max-clique; the alternative variant swaps A. aspera for A. formosa,
#  but A. formosa's totals failed the printed-checksum validation).
#  M = Laplace-smoothed win fractions among decisive contests (ties excluded):
#  p = (wins_A + 1) / (wins_A + wins_B + 2). heron_n = interactions per pair
#  (wins + losses + ties); heron_ties = ties per pair (large here: reef-flat
#  stand-offs are common).
#  CAUTION: these records are the source behind the observed compilation's
#  Heron Island block (ref 1) - never pool the two (double-counting).
#  Transcription validated against the tables' printed per-species totals
#  (20 of 52 exact; residual count uncertainty mainly in table-3 cells).
# =============================================================================

heron_species <- c("Aasp","Acun","Adig","Ahum","Ahya","Amil","Anas","Arob","Pdam","Por")
heron_names <- c(Aasp = "Acropora aspera",    Acun = "Acropora cuneata",
                 Adig = "Acropora digitifera", Ahum = "Acropora humilis",
                 Ahya = "Acropora hyacinthus", Amil = "Acropora millepora",
                 Anas = "Acropora nasuta",     Arob = "Acropora robusta",
                 Pdam = "Pocillopora damicornis", Por = "Porites spp.")

M_heron10 <- rbind(
  Aasp = c(.50,.60,.50,.83,.50,.50,.75,.60,.69,.88),
  Acun = c(.40,.50,.33,.14,.14,.15,.20,.20,.67,.88),
  Adig = c(.50,.67,.50,.60,.36,.67,.33,.40,.90,.96),
  Ahum = c(.17,.86,.40,.50,.09,.12,.67,.20,.80,.86),
  Ahya = c(.50,.86,.64,.91,.50,.70,.50,.67,.92,.93),
  Amil = c(.50,.85,.33,.88,.30,.50,.50,.25,.71,.92),
  Anas = c(.25,.80,.67,.33,.50,.50,.50,.40,.89,.89),
  Arob = c(.40,.80,.60,.80,.33,.75,.60,.50,.82,.91),
  Pdam = c(.31,.33,.10,.20,.08,.29,.11,.18,.50,.90),
  Por  = c(.12,.12,.04,.14,.07,.08,.11,.09,.10,.50))
colnames(M_heron10) <- heron_species

heron_n <- rbind(
  Aasp = c(0,4,1,8,8,10,6,6,13,11),
  Acun = c(4,0,15,7,13,19,16,9,6,14),
  Adig = c(1,15,0,4,16,10,8,10,23,48),
  Ahum = c(8,7,4,0,9,10,2,3,5,10),
  Ahya = c(8,13,16,9,0,11,9,1,26,12),
  Amil = c(10,19,10,10,11,0,1,7,10,16),
  Anas = c(6,16,8,2,9,1,0,6,11,21),
  Arob = c(6,9,10,3,1,7,6,0,15,14),
  Pdam = c(13,6,23,5,26,10,11,15,0,14),
  Por  = c(11,14,48,10,12,16,21,14,14,0))
colnames(heron_n) <- heron_species

heron_ties <- rbind(
  Aasp = c(0,1,1,4,6,4,4,3,2,5),
  Acun = c(1,0,11,2,8,8,8,6,2,8),
  Adig = c(1,11,0,1,7,6,4,7,15,26),
  Ahum = c(4,2,1,0,0,4,1,0,2,5),
  Ahya = c(6,8,7,0,0,3,3,0,4,0),
  Amil = c(4,8,6,4,3,0,1,5,5,6),
  Anas = c(4,8,4,1,3,1,0,3,4,14),
  Arob = c(3,6,7,0,0,5,3,0,6,5),
  Pdam = c(2,2,15,2,4,5,4,6,0,6),
  Por  = c(5,8,26,5,0,6,14,5,6,0))
colnames(heron_ties) <- heron_species

writeHeronExcel <- function(out_dir = interaction_xlsx_dir) {
  if (!requireNamespace("writexl", quietly = TRUE)) stop("install.packages('writexl')")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  notes <- data.frame(Note = c(
    "Source: three site-level Heron Island (GBR, Indo-Pacific) interaction tables (Connell-programme long-term reef-flat records), combined across sites.",
    "Convention: cell [row, column] = probability the ROW species overgrows/wins against the COLUMN species, given a decisive contest. M[i,j] + M[j,i] = 1; diagonal 0.5.",
    "Probabilities are MEASURED: Laplace-smoothed win fractions among decisive contests, p = (wins_A + 1)/(wins_A + wins_B + 2); ties excluded from the denominator (see Ties sheet - stand-offs are frequent on the reef flat).",
    "These 10 species are the unique largest subset with every pair observed (exact max-clique; alternative variant swaps A. aspera for A. formosa, rejected because A. formosa failed checksum validation).",
    "Sample sizes per pair in the Sample_sizes sheet (wins + losses + ties; range 1-48).",
    "CAUTION: same source as the observed compilation's Heron Island block (ref 1) - never pool the two datasets (double-counting).",
    "Transcription validated against printed per-species totals: 20 of 52 species-totals exact; residual count uncertainty concentrated in table-3-derived cells."))
  path <- file.path(out_dir, "Heron_interaction_matrix.xlsx")
  writexl::write_xlsx(list(Matrix_10sp  = .matSheet(M_heron10),
                           Sample_sizes = .matSheet(heron_n),
                           Ties         = .matSheet(heron_ties),
                           Species_key  = .keySheet(heron_names),
                           Notes        = notes), path = path)
  cat("Wrote:", path, "\n")
  invisible(path)
}


# =============================================================================
#  5. "PRECODA" OBSERVED COMPILATION, DE-DUPLICATED - the observed-contest
#  records EXCLUDING the studies already covered by their own matrices here:
#  ref 1 (Heron Island -> section 4), ref 6 (Logan -> section 1), ref 7 (Dai ->
#  section 2), ref 17 (Veracruz, set aside). With the standard n > 2 filter the
#  largest complete matrix in what remains is the OKINAWA PORITES set
#  (Rinkevich & Sakai 2001): 5 congeneric species, every pair measured with
#  n = 7-21 contests - a strong, fully counted linear hierarchy.
# =============================================================================

precoda_obs_excluded_refs <- c(1, 6, 7, 17)
porites5 <- c("Porites australiensis", "Porites cylindrica", "Porites lobata",
              "Porites lutea", "Porites rus")

writePrecodaObservedExcel <- function(out_dir = interaction_xlsx_dir, min_obs = 3) {
  if (!requireNamespace("writexl", quietly = TRUE)) stop("install.packages('writexl')")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  obs <- buildObservedPairs(min_obs = min_obs, exclude_refs = precoda_obs_excluded_refs)

  #the complete 5-species Porites block (verified largest complete subset)
  sub <- obs[obs$speciesA %in% porites5 & obs$speciesB %in% porites5, ]
  if (nrow(sub) < choose(length(porites5), 2))
    stop("Porites-5 block incomplete under these settings - re-check filters.")
  M <- observedMatrix(porites5, obs)
  Nn <- matrix(NA_integer_, 5, 5, dimnames = list(porites5, porites5))
  for (r in seq_len(nrow(sub))) {
    Nn[sub$speciesA[r], sub$speciesB[r]] <- sub$n[r]
    Nn[sub$speciesB[r], sub$speciesA[r]] <- sub$n[r]
  }
  diag(Nn) <- 0L

  notes <- data.frame(Note = c(
    "Source: observed coral-contest compilation (source data of Precoda et al. 2017, Am Nat 190:420-429): individual win/lose/draw records from published studies.",
    paste("EXCLUDED references:", paste(precoda_obs_excluded_refs, collapse = ", "),
          "- Heron Island (ref 1), Logan 1984 (ref 6) and Dai 1990 (ref 7) have their own dedicated matrices in this collection; ref 17 (Veracruz) set aside."),
    paste("Filter: only species pairs with more than", min_obs - 1, "recorded interactions."),
    "Convention: cell [row, column] = probability the ROW species beats the COLUMN species, given a decisive contest; Laplace-smoothed win fractions, p = (wins_A + 1)/(wins_A + wins_B + 2); draws excluded from the denominator. M[i,j] + M[j,i] = 1; diagonal 0.5.",
    "Sheet 'Complete_5sp': the largest species subset with EVERY pair observed after filtering - the Okinawa Porites set (Rinkevich & Sakai 2001, Zoology 104:91-97; Indo-Pacific): a strongly linear congeneric hierarchy, n = 7-21 contests per pair.",
    "Sheet 'All_filtered_pairs': every pair passing the filter (with counts), for building other subsets via observedMatrix()."))
  path <- file.path(out_dir, "PrecodaObserved_interaction_matrix.xlsx")
  writexl::write_xlsx(list(Complete_5sp       = .matSheet(round(M, 3)),
                           Sample_sizes       = .matSheet(Nn),
                           All_filtered_pairs = obs,
                           Notes              = notes), path = path)
  cat("Wrote:", path, "\n")
  invisible(path)
}


# =============================================================================
#  WRITE ALL FOUR WORKBOOKS (one per study; never combined)
# =============================================================================
writeAllInteractionExcels <- function() {
  writeLoganExcel()            # 1. Logan 1984      - Bermuda, 17 sp, measured (graded)
  writeDaiExcel()              # 2. Dai 1990        - Taiwan, 19 sp categorical + complete 7 sp
  writePrecodaObservedExcel()  # 3. Observed compilation (de-duplicated) - Okinawa Porites 5 sp, measured
  writeHeronExcel()            # 4. Heron Island    - GBR, 10 sp, measured
}


# =============================================================================
#  GROWTH FORM -> GROWTH-RATE TRAIT  (a togglable per-species trait)
# -----------------------------------------------------------------------------
#  Per-species growth traits derived from documented colony morphology (growth
#  forms + depths from coraltraits.org / Corals of the World / WoRMS / IUCN).
#  Three tiers on the model's 1-5 scale, compressed to 2/3/4 (a 2:1 fast:slow
#  ratio). A pilot showed the maximal 1/5 spread SWAMPS the interaction matrix
#  (growth alone decided outcomes); 2/3/4 keeps the matrix the primary driver
#  while growth form still matters as a secondary axis.
#      Fast (4)  = branching / arborescent / digitate / corymbose / tabular
#      Intermed (3) = submassive / foliose / laminar / columnar
#      Slow (2)  = massive / encrusting / solitary (free-living)
#  These are OPTIONAL: the final runs toggle growth variation on/off (off =
#  uniform growth 3, so the interaction matrix is the sole species asymmetry).
# =============================================================================

GROWTH_TIER <- c(fast = 4L, intermediate = 3L, slow = 2L)

#Logan (Bermuda) - morphology + growth trait, keyed by species code
logan_morph <- c(
  Ds = "massive",  Dl = "massive",  Ma = "massive",  Mc = "massive",
  Pa = "massive-encrusting", Pp = "branching", Ff = "massive", Md = "encrusting",
  Mm = "branching", Oc = "branching", Is = "massive", Me = "submassive",
  Sm = "massive-encrusting", Dk = "massive", Sy = "solitary/free-living",
  Si = "massive", Ag = "laminar/foliose")
logan_growth <- c(
  Ds = 2L, Dl = 2L, Ma = 2L, Mc = 2L, Pa = 2L, Pp = 4L, Ff = 2L, Md = 2L,
  Mm = 4L, Oc = 4L, Is = 2L, Me = 3L, Sm = 2L, Dk = 2L, Sy = 2L, Si = 2L, Ag = 3L)

#Logan reduced to 10 species (to match Heron's richness), removing redundant
#near-duplicates - one of each pair sharing genus/morphology/growth: dropped
#Ds (~Dl brain), Ma (~Mc boulder), Dk (~Me meandrinid), Sm (~Pa massive-encr.),
#Mm (branching; Madracis kept via encrusting Md), Ff & Si (redundant weak
#massives). Keeps 8 morphologies, the full competitive gradient (0.31-0.74),
#and every genus unique except Porites (massive Pa vs branching Pp - distinct
#forms). Logan's full matrix is complete, so all 45 pairs are observed.
logan10_species <- c("Is","Me","Mc","Sy","Dl","Pa","Oc","Pp","Md","Ag")
M_logan10 <- M_logan_graded[logan10_species, logan10_species]

#Logan 7- and 3-species subsets (same reduce-redundancy method applied further).
#Logan-7 drops Mc (redundant strong massive ~ Is/Dl), Pp (branching ~ Oc), Md
#(encrusting ~ Pa): keeps 6 morphologies, competition 0.28-0.74.
#Logan-3 = the 3 most distinct: massive/submassive/branching, growth 2/3/4,
#genera Isophyllia/Meandrina/Oculina - two strong slow competitors vs one weak
#fast colonizer (a competition-colonization contrast).
logan7_species <- c("Is","Me","Sy","Dl","Oc","Pa","Ag")
logan3_species <- c("Is","Me","Oc")
M_logan7 <- M_logan_graded[logan7_species, logan7_species]
M_logan3 <- M_logan_graded[logan3_species, logan3_species]

#Heron 7- and 3-species subsets.
#Heron-7 drops Ahum (digitate ~ Adig), Anas (corymbose ~ Amil), Aasp (branching
#~ Arob): keeps 7 distinct morphologies across 4 genera, competition 0.14-0.71.
#Heron-3 = tabular/submassive-encrusting/massive, growth 4/3/2, genera Acropora/
#Isopora/Porites - a clean strong->mid->weak gradient (0.71/0.36/0.14).
heron7_species <- c("Ahya","Arob","Adig","Amil","Acun","Pdam","Por")
heron3_species <- c("Ahya","Acun","Por")
M_heron7 <- M_heron10[heron7_species, heron7_species]
M_heron3 <- M_heron10[heron3_species, heron3_species]

#Heron Island (GBR) - morphology + growth trait, keyed by species code
heron_morph <- c(
  Aasp = "branching", Acun = "submassive-encrusting", Adig = "digitate",
  Ahum = "digitate", Ahya = "tabular", Amil = "corymbose", Anas = "corymbose",
  Arob = "robust branching", Pdam = "branching", Por = "massive")
heron_growth <- c(
  Aasp = 4L, Acun = 3L, Adig = 4L, Ahum = 4L, Ahya = 4L, Amil = 4L,
  Anas = 4L, Arob = 4L, Pdam = 4L, Por = 2L)

#Look up the per-species growth vector for a named community (else uniform 3).
#Vectors are keyed by species code, so a subset (e.g. logan10) aligns by name.
communityGrowth <- function(community) {
  switch(as.character(community),
         logan = , logan10 = , logan7 = , logan3 = logan_growth,
         heron = , heron10 = , heron7 = , heron3 = heron_growth,
         mia = , mia10 = , mia7 = , mia3 = if (exists("mia_growth")) mia_growth else 3L,
         3L)
}

#Readable morphology/growth table for a community, for checking.
growthFormTable <- function(community = c("logan", "logan10", "logan7", "logan3",
                                          "heron", "heron7", "heron3")) {
  community <- match.arg(community)
  codes <- switch(community, logan = logan_species, logan10 = logan10_species,
                  logan7 = logan7_species, logan3 = logan3_species,
                  heron7 = heron7_species, heron3 = heron3_species, heron_species)
  logan <- grepl("^logan", community)
  nm    <- if (logan) logan_names  else heron_names
  mo    <- if (logan) logan_morph  else heron_morph
  gr    <- if (logan) logan_growth else heron_growth
  tier  <- c(`2` = "slow", `3` = "intermediate", `4` = "fast")
  data.frame(code = codes, species = unname(nm[codes]),
             morphology = unname(mo[codes]),
             growth_rate = unname(gr[codes]),
             tier = unname(tier[as.character(gr[codes])]),
             row.names = NULL, stringsAsFactors = FALSE)
}
