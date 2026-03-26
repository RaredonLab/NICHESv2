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
#' \code{compute_CellToCell_Spatial()}. Two curated databases can be loaded
#' by name; alternatively, any user-supplied \code{data.frame} is accepted and
#' column names are remapped to the canonical \code{"ligand"} /
#' \code{"receptor"} names.
#'
#' @section Built-in databases:
#' \describe{
#'   \item{\code{"fantom5"}}{
#'     Ramilowski et al. (2015) FANTOM5 human ligand-receptor pairs.
#'     Requires internet access on first use; downloaded from the OmniPath
#'     REST API (\url{https://omnipathdb.org}) via the \pkg{OmnipathR} package.
#'     Install with: \code{BiocManager::install("OmnipathR")}.
#'   }
#'   \item{\code{"omnipath"}}{
#'     OmniPath intercellular communication network (Türei et al. 2021).
#'     Also requires \pkg{OmnipathR}.
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
#' @param ligand.col    Column name for ligand symbols when \code{db} is a
#'   \code{data.frame}. Default: \code{"ligand"}.
#' @param receptor.col  Column name for receptor symbols when \code{db} is a
#'   \code{data.frame}. Default: \code{"receptor"}.
#' @param verbose       Print a summary message after loading. Default \code{TRUE}.
#'
#' @return A \code{data.frame} with columns \code{"ligand"}, \code{"receptor"},
#'   and any additional annotation columns present in the source. Rows with
#'   \code{NA} in ligand or receptor are silently removed.
#'
#' @examples
#' \dontrun{
#' # Built-in databases (require OmnipathR)
#' lrm <- load_LRM_database("omnipath")
#' lrm <- load_LRM_database("fantom5")
#'
#' # User-supplied table
#' my.db <- data.frame(
#'   Ligand   = c("MIF", "VEGFA", "CXCL12_PF4"),
#'   Receptor = c("CD74", "KDR",  "ACKR3")
#' )
#' lrm <- load_LRM_database(my.db, ligand.col = "Ligand", receptor.col = "Receptor")
#' }
#'
#' @export
load_LRM_database <- function(db,
                               ligand.col   = "ligand",
                               receptor.col = "receptor",
                               verbose      = TRUE) {

  # ── Dispatch on type ────────────────────────────────────────────────────────
  if (is.character(db) && length(db) == 1L) {
    db.name <- tolower(db)

    if (db.name == "omnipath") {
      result <- .load_omnipath(verbose)
    } else if (db.name == "fantom5") {
      result <- .load_fantom5(verbose)
    } else {
      stop(sprintf(
        "[load_LRM_database] Unknown built-in database '%s'.
  Supported shorthands: 'omnipath', 'fantom5'.
  To use a custom table, supply a data.frame directly.",
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

  # ── Final cleaning: remove NA rows, reset rownames ────────────────────────
  result <- result[!is.na(result$ligand) & !is.na(result$receptor), , drop = FALSE]
  rownames(result) <- NULL

  result
}


# ==============================================================================
# INTERNAL LOADERS
# ==============================================================================

#' @noRd
.load_omnipath <- function(verbose) {

  if (!requireNamespace("OmnipathR", quietly = TRUE))
    stop(
      "[load_LRM_database] The 'OmnipathR' package is required to load the ",
      "OmniPath database.\n  Install it with: BiocManager::install('OmnipathR')"
    )

  if (verbose) message("[load_LRM_database] Fetching OmniPath intercell network ...")

  raw <- OmnipathR::import_intercell_network(
    transmitter_categories = c("ligand"),
    receiver_categories    = c("receptor")
  )

  result <- data.frame(
    ligand   = raw$source_genesymbol,
    receptor = raw$target_genesymbol,
    stringsAsFactors = FALSE
  )

  # Carry useful annotation columns if present
  extra.cols <- intersect(
    c("category_intercell_source", "category_intercell_target",
      "database", "curation_effort"),
    colnames(raw)
  )
  if (length(extra.cols) > 0L) result <- cbind(result, raw[, extra.cols, drop = FALSE])

  result <- unique(result)

  if (verbose)
    message(sprintf("[load_LRM_database] OmniPath: %d unique LR pairs loaded.", nrow(result)))

  result
}


#' @noRd
.load_fantom5 <- function(verbose) {

  if (!requireNamespace("OmnipathR", quietly = TRUE))
    stop(
      "[load_LRM_database] The 'OmnipathR' package is required to load the ",
      "FANTOM5 database via OmniPath.\n  Install it with: BiocManager::install('OmnipathR')"
    )

  if (verbose) message("[load_LRM_database] Fetching FANTOM5 LR pairs via OmniPath ...")

  raw <- OmnipathR::import_intercell_network(
    transmitter_categories = c("ligand"),
    receiver_categories    = c("receptor"),
    datasets               = c("fantom5")
  )

  result <- data.frame(
    ligand   = raw$source_genesymbol,
    receptor = raw$target_genesymbol,
    stringsAsFactors = FALSE
  )

  extra.cols <- intersect(
    c("category_intercell_source", "category_intercell_target", "curation_effort"),
    colnames(raw)
  )
  if (length(extra.cols) > 0L) result <- cbind(result, raw[, extra.cols, drop = FALSE])

  result <- unique(result)

  if (verbose)
    message(sprintf("[load_LRM_database] FANTOM5: %d unique LR pairs loaded.", nrow(result)))

  result
}


#' @noRd
.load_custom <- function(db, ligand.col, receptor.col, verbose) {

  db <- as.data.frame(db)

  missing.cols <- setdiff(c(ligand.col, receptor.col), colnames(db))
  if (length(missing.cols) > 0L)
    stop(sprintf(
      "[load_LRM_database] Column(s) not found in supplied data.frame: %s
  Use 'ligand.col' and 'receptor.col' to specify the correct column names.",
      paste(missing.cols, collapse = ", ")
    ))

  # Rename to canonical names
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
