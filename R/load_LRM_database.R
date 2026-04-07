# ==============================================================================
# FILE:    load_LRM_database.R
# PURPOSE: Load a ligand-receptor pair database, either from a built-in
#          shorthand ("fantom5", "omnipath") or from a user-supplied
#          data.frame. Always returns a standardised data.frame with at least
#          two columns: "ligand" and "receptor".
#
# FUNCTIONS:
#   load_LRM_database()  — main function
#
# DEPS: data.table
# ==============================================================================


# ==============================================================================
#' Load a Ligand-Receptor Pair Database
#'
#' @description
#' Returns a standardised ligand-receptor pair table for use as the
#' \code{LRM.db} argument of \code{create_NICHESObject()} and
#' \code{compute_CellToCell()}. Two curated databases can be loaded by name;
#' alternatively, any user-supplied \code{data.frame} is accepted and column
#' names are remapped to the canonical \code{"ligand"} / \code{"receptor"}
#' names.
#'
#' @section Built-in databases:
#' \describe{
#'   \item{\code{"fantom5"}}{
#'     Ramilowski et al. (2015) FANTOM5 ligand-receptor pairs, bundled within
#'     the package for human, mouse, rat, and pig. No internet access or
#'     additional packages required. Species-specific datasets were derived
#'     using the methodology in Raredon et al. 2019
#'     (DOI: 10.1126/sciadv.aaw3851).
#'   }
#'   \item{\code{"omnipath"}}{
#'     OmniPath ligand-receptor extra interactions (Türei et al. 2021).
#'     Requires internet access and the \pkg{OmnipathR} package.
#'     Install with: \code{BiocManager::install("OmnipathR")}.
#'     Supported species: \code{"human"}, \code{"mouse"}, \code{"rat"}.
#'   }
#' }
#'
#' @section Complex ligands/receptors:
#' Both built-in databases and user-supplied tables may contain gene complexes.
#' Complexes are represented as underscore-joined gene symbols, e.g.
#' \code{"COL1A1_COL1A2"} for the collagen I heterodimer. The scoring
#' functions in NICHESv2 handle complexes automatically.
#'
#' @param db            Character string (\code{"fantom5"} or \code{"omnipath"})
#'   OR a \code{data.frame} / \code{data.table} with at least two columns
#'   containing ligand and receptor gene symbols.
#' @param species       Species for built-in databases. One of \code{"human"}
#'   (default), \code{"mouse"}, \code{"rat"}, or \code{"pig"}
#'   (FANTOM5 only — OmniPath does not support pig). Ignored silently when
#'   \code{db} is a \code{data.frame}.
#' @param ligand.col    Column name for ligand symbols when \code{db} is a
#'   \code{data.frame}. Default: \code{"ligand"}.
#' @param receptor.col  Column name for receptor symbols when \code{db} is a
#'   \code{data.frame}. Default: \code{"receptor"}.
#' @param verbose       Print progress messages. Default \code{TRUE}.
#'
#' @return A \code{data.frame} with columns \code{"ligand"}, \code{"receptor"},
#'   and any additional annotation columns present in the source. Rows with
#'   \code{NA} in ligand or receptor are silently removed.
#'
#' @examples
#' \dontrun{
#' # Built-in databases — species must match your count matrix gene symbols
#' lrm <- load_LRM_database("fantom5", species = "human")
#' lrm <- load_LRM_database("fantom5", species = "rat")
#' lrm <- load_LRM_database("omnipath", species = "mouse")
#'
#' # User-supplied table (species argument ignored)
#' my.db <- data.frame(
#'   Ligand   = c("MIF", "VEGFA", "Cxcl12_Pf4"),
#'   Receptor = c("CD74", "KDR",  "Ackr3")
#' )
#' lrm <- load_LRM_database(my.db, ligand.col = "Ligand", receptor.col = "Receptor")
#' }
#'
#' @export
load_LRM_database <- function(db,
                               species      = "human",
                               ligand.col   = "ligand",
                               receptor.col = "receptor",
                               verbose      = TRUE) {

  # ── Dispatch on type ────────────────────────────────────────────────────────
  if (is.character(db) && length(db) == 1L) {
    db.name <- tolower(db)

    if (db.name == "omnipath") {
      result <- .load_omnipath(species, verbose)
    } else if (db.name == "fantom5") {
      result <- .load_fantom5(species, verbose)
    } else {
      stop(sprintf(
        "[load_LRM_database] Unknown built-in database '%s'.\n  Supported shorthands: 'omnipath', 'fantom5'.\n  To use a custom table, supply a data.frame directly.",
        db
      ))
    }

  } else if (is.data.frame(db) || is.data.table(db)) {
    result <- .load_custom(db, ligand.col, receptor.col, verbose)

  } else {
    stop(
      "[load_LRM_database] 'db' must be either a character string ",
      "('omnipath' or 'fantom5') or a data.frame."
    )
  }

  # ── Final cleaning: remove NA rows, reset rownames ─────────────────────────
  result <- result[!is.na(result$ligand) & !is.na(result$receptor), , drop = FALSE]
  rownames(result) <- NULL

  result
}


# ==============================================================================
# INTERNAL LOADERS
# ==============================================================================

#' @noRd
.load_fantom5 <- function(species, verbose) {

  valid.species <- c("human", "mouse", "rat", "pig")
  if (!species %in% valid.species)
    stop(sprintf(
      "[load_LRM_database] Invalid species '%s' for FANTOM5.\n  Supported species: %s.",
      species, paste(valid.species, collapse = ", ")
    ))

  data.name <- switch(species,
    human = "ncomms8866_human",
    mouse = "ncomms8866_mouse",
    rat   = "ncomms8866_rat",
    pig   = "ncomms8866_pig"
  )

  if (verbose)
    message(sprintf(
      "[load_LRM_database] Loading FANTOM5 (%s) from bundled package data ...",
      species
    ))

  env <- new.env(parent = emptyenv())
  utils::data(list = data.name, package = "NICHESv2", envir = env)
  fantom <- env[[data.name]]

  result <- data.frame(
    ligand   = fantom$Ligand.ApprovedSymbol,
    receptor = fantom$Receptor.ApprovedSymbol,
    stringsAsFactors = FALSE
  )
  result <- unique(result)

  if (verbose)
    message(sprintf(
      "[load_LRM_database] FANTOM5 (%s): %d unique LR pairs loaded.",
      species, nrow(result)
    ))

  result
}


#' @noRd
.load_omnipath <- function(species, verbose) {

  valid.species  <- c("human", "mouse", "rat")
  organism.codes <- c(human = 9606L, mouse = 10090L, rat = 10116L)

  if (!species %in% valid.species)
    stop(sprintf(
      "[load_LRM_database] Invalid species '%s' for OmniPath.\n  Supported species: %s.",
      species, paste(valid.species, collapse = ", ")
    ))

  if (!requireNamespace("OmnipathR", quietly = TRUE))
    stop(
      "[load_LRM_database] The 'OmnipathR' package is required to load the ",
      "OmniPath database.\n  Install it with: BiocManager::install('OmnipathR')"
    )

  organism <- organism.codes[[species]]

  if (verbose)
    message(sprintf(
      "[load_LRM_database] Fetching OmniPath ligand-receptor interactions (species = '%s', organism = %d) ...",
      species, organism
    ))

  raw <- OmnipathR::import_ligrecextra_interactions(organism = organism)

  result <- data.frame(
    ligand   = raw$source_genesymbol,
    receptor = raw$target_genesymbol,
    stringsAsFactors = FALSE
  )

  extra.cols <- intersect(
    c("category_intercell_source", "category_intercell_target",
      "database", "curation_effort"),
    colnames(raw)
  )
  if (length(extra.cols) > 0L) result <- cbind(result, raw[, extra.cols, drop = FALSE])

  result <- unique(result)

  if (verbose)
    message(sprintf(
      "[load_LRM_database] OmniPath (%s): %d unique LR pairs loaded.",
      species, nrow(result)
    ))

  result
}


#' @noRd
.load_custom <- function(db, ligand.col, receptor.col, verbose) {

  db <- as.data.frame(db)

  missing.cols <- setdiff(c(ligand.col, receptor.col), colnames(db))
  if (length(missing.cols) > 0L)
    stop(sprintf(
      "[load_LRM_database] Column(s) not found in supplied data.frame: %s\n  Use 'ligand.col' and 'receptor.col' to specify the correct column names.",
      paste(missing.cols, collapse = ", ")
    ))

  # Rename to canonical names; carry any extra annotation columns through
  other.cols <- setdiff(colnames(db), c(ligand.col, receptor.col))
  result <- data.frame(
    ligand   = db[[ligand.col]],
    receptor = db[[receptor.col]],
    db[, other.cols, drop = FALSE],
    stringsAsFactors = FALSE
  )

  if (verbose)
    message(sprintf("[load_LRM_database] Custom table: %d LR pairs loaded.", nrow(result)))

  result
}
