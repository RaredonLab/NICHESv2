# load_LRM_database.R -----------------------------------------------------
# Load a ligand-receptor pair database, either from a built-in
# shorthand ("connectomedb2025", "fantom5", "omnipath") or from a user-supplied
# data.frame. Always returns a standardised data.frame with at least
# two columns: "ligand" and "receptor".
#
# FUNCTIONS:
#   load_LRM_database()  -- main function
#
# DEPS: data.table


# ==============================================================================
#' Load a Ligand-Receptor Pair Database
#'
#' @description
#' Returns a standardized ligand-receptor pair table for use as the
#' \code{LRM.db} argument of \code{create_NICHESObject()} and
#' \code{compute_CellToCell()}. Two curated databases can be loaded by name;
#' alternatively, any user-supplied \code{data.frame} is accepted and column
#' names are remapped to the canonical \code{"ligand"} / \code{"receptor"}
#' names.
#'
#' @section Built-in databases:
#' \describe{
#'   \item{\code{"connectomedb2025"}}{
#'     ConnectomeDB2025 ligand-receptor pairs (default), bundled within the
#'     package for human, mouse, rat, and pig. No internet access or additional
#'     packages required. Includes experimentally verified and orthology-inferred
#'     interactions with evidence labels and subcellular location annotations.
#'   }
#'   \item{\code{"fantom5"}}{
#'     Ramilowski et al. (2015) FANTOM5 ligand-receptor pairs, bundled within
#'     the package for human, mouse, rat, and pig. No internet access or
#'     additional packages required. Species-specific datasets were derived
#'     using the methodology in Raredon et al. 2019
#'     (DOI: 10.1126/sciadv.aaw3851).
#'   }
#'   \item{\code{"omnipath"}}{
#'     OmniPath ligand-receptor extra interactions (Turei et al. 2021).
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
#' @param db            Character string (\code{"connectomedb2025"} (default),
#'   \code{"fantom5"}, or \code{"omnipath"}) OR a \code{data.frame} /
#'   \code{data.table} with at least two columns containing ligand and receptor
#'   gene symbols.
#' @param species       Species for built-in databases. One of \code{"human"}
#'   (default), \code{"mouse"}, \code{"rat"}, or \code{"pig"}
#'   (ConnectomeDB2025 and FANTOM5 only; OmniPath does not support pig).
#'   Ignored silently when \code{db} is a \code{data.frame}.
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
#' # ConnectomeDB2025 is the default -- no db argument needed.
#' lrm.human <- load_LRM_database(species = "human")
#' head(lrm.human)
#'
#' # FANTOM5 -- bundled, always available, no internet.
#' lrm.f5  <- load_LRM_database("fantom5", species = "human")
#' lrm.rat <- load_LRM_database("fantom5", species = "rat")
#'
#' # User-supplied table -- species argument is ignored. Column names are
#' # remapped to canonical "ligand" / "receptor".
#' my.db <- data.frame(
#'   Ligand   = c("MIF", "VEGFA", "Cxcl12_Pf4"),
#'   Receptor = c("CD74", "KDR",  "Ackr3")
#' )
#' lrm.user <- load_LRM_database(my.db,
#'                               ligand.col   = "Ligand",
#'                               receptor.col = "Receptor")
#'
#' # OmniPath needs internet access, so this runs in interactive sessions only.
#' if (interactive()) {
#'   lrm.op <- load_LRM_database("omnipath", species = "mouse")
#' }
#'
#' @export
load_LRM_database <- function(db           = "connectomedb2025",
                              species      = "human",
                              ligand.col   = "ligand",
                              receptor.col = "receptor",
                              verbose      = TRUE) {

  # -- Dispatch on type --------------------------------------------------------
  if (is.character(db) && length(db) == 1L) {
    db.name <- tolower(db)

    if (db.name == "connectomedb2025") {
      result <- .load_connectomedb2025(species, verbose)
    } else if (db.name == "omnipath") {
      result <- .load_omnipath(species, verbose)
    } else if (db.name == "fantom5") {
      result <- .load_fantom5(species, verbose)
    } else {
      stop(sprintf(
        "[load_LRM_database] Unknown built-in database '%s'.\n  Supported shorthands: 'connectomedb2025', 'omnipath', 'fantom5'.\n  To use a custom table, supply a data.frame directly.",
        db
      ))
    }
    
  } else if (is.data.frame(db)) {
    result <- .load_custom(db, ligand.col, receptor.col, verbose)
    
  } else {
    stop(
      "[load_LRM_database] 'db' must be either a character string ",
      "('connectomedb2025', 'omnipath', or 'fantom5') or a data.frame."
    )
  }

  # -- Final cleaning: remove NA rows, reset rownames -------------------------
  result <- result[!is.na(result$ligand) & !is.na(result$receptor), , drop = FALSE]
  rownames(result) <- NULL

  result
}


## INTERNAL LOADERS ----

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
.load_connectomedb2025 <- function(species, verbose) {
  
  valid.species <- c("human", "mouse", "rat", "pig")
  if (!species %in% valid.species)
    stop(sprintf(
      "[load_LRM_database] Invalid species '%s' for ConnectomeDB2025.\n  Supported species: %s.",
      species, paste(valid.species, collapse = ", ")
    ))
  
  data.name <- paste0("connectomedb2025_", species)
  
  if (verbose)
    message(sprintf(
      "[load_LRM_database] Loading ConnectomeDB2025 (%s) from bundled package data ...",
      species
    ))
  
  env <- new.env(parent = emptyenv())
  utils::data(list = data.name, package = "NICHESv2", envir = env)
  cdb <- env[[data.name]]
  cdb <- unique(cdb)
  
  if (!all(c("ligand", "receptor") %in% colnames(cdb)))
    stop(sprintf(
      "[load_LRM_database] Bundled 'connectomedb2025_%s' is missing 'ligand' or 'receptor' columns.\n  Re-run data-raw/build_connectomedb2025.R to rebuild the bundled data.",
      species
    ))
  
  if (verbose)
    message(sprintf(
      "[load_LRM_database] ConnectomeDB2025 (%s): %d unique LR pairs loaded.",
      species, nrow(cdb)
    ))
  
  cdb
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

  # NOTE (ported from NICHESv1 PR #72 / issue #70): newer OmnipathR (>= ~3.17,
  # incl. 4.1.0) resolves the `organism` argument through OmnipathR::ncbi_taxid(),
  # which builds an organism-name table by downloading and coalescing species
  # lists from Ensembl, OMA and UniProt. When one of those ancillary services is
  # unreachable (e.g. omabrowser.org returning HTTP 502), that step aborts with
  # "Can't combine `..1` and `..3`" - even though the OmniPath interaction server
  # itself is fine. Since we have already mapped the species to its NCBI taxon id
  # above, we do not need that translation. The helper below first tries the
  # normal OmnipathR path (preserving its caching when the services are healthy)
  # and, if that fails, falls back to querying the OmniPath REST API directly with
  # the known taxon id, which bypasses the organism-name lookup.
  raw <- .niches_fetch_omnipath_ligrec(organism, verbose)

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


## OMNIPATH RESILIENT FETCH (ported from NICHESv1 PR #72 / issue #70) ----

#' Resiliently fetch OmniPath ligand-receptor interactions
#'
#' Internal helper for \code{.load_omnipath()} that works around a known upstream
#' OmnipathR failure (NICHESv1 issue #70). Newer OmnipathR versions resolve the
#' \code{organism} argument via an organism-name table built from Ensembl, OMA and
#' UniProt downloads; when one of those services is down the call aborts before any
#' interactions are retrieved. This helper first attempts the standard OmnipathR
#' path, and on failure falls back to a direct OmniPath REST query using the
#' already-known NCBI taxon id, which does not require the organism-name lookup.
#'
#' @param organism Integer NCBI taxonomy id (9606 human, 10090 mouse, 10116 rat).
#' @param verbose  Print progress/provenance messages.
#'
#' @return A data.frame of ligand-receptor interactions containing at least the
#'   \code{source_genesymbol} and \code{target_genesymbol} columns.
#'
#' @noRd
.niches_fetch_omnipath_ligrec <- function(organism, verbose = TRUE) {

  # 1) Preferred path: standard OmnipathR call (uses its cache when services are
  #    up). The organism-name download noise is muted; if it fails we recover.
  result <- tryCatch(
    .niches_quiet_omnipath(OmnipathR::import_ligrecextra_interactions(organism = organism)),
    error = function(e) e
  )
  if (!inherits(result, "error")) {
    if (verbose) {
      ver <- tryCatch(as.character(utils::packageVersion("OmnipathR")),
                      error = function(e) NA_character_)
      message(sprintf(
        "[load_LRM_database] Loaded OmniPath via OmnipathR %s (retrieved %s): %d interactions.",
        ver, Sys.Date(), nrow(result)))
    }
    return(result)
  }

  # 2) Fallback: query the OmniPath REST API directly with the known taxon id.
  #    This avoids OmnipathR's organism-name translation (the step that fails).
  fallback <- tryCatch(.niches_omnipath_rest_ligrec(organism), error = function(e) e)
  if (!inherits(fallback, "error") && nrow(fallback) > 0L) {
    if (verbose) {
      ver <- .niches_omnipath_server_version()
      message(sprintf(
        "[load_LRM_database] Loaded OmniPath via REST fallback (OmniPath server %s, retrieved %s): %d interactions.",
        if (is.na(ver)) "REST API" else ver, Sys.Date(), nrow(fallback)))
    }
    return(fallback)
  }

  # 3) Both paths failed: raise an informative error (this stays loud on purpose).
  stop(
    "[load_LRM_database] Unable to load the OmniPath ligand-receptor database.\n",
    "  OmnipathR error : ", conditionMessage(result), "\n",
    "  REST fallback   : ",
    if (inherits(fallback, "error")) conditionMessage(fallback) else "returned no interactions", "\n",
    "  The OmniPath interaction server (omnipathdb.org) may be unreachable from\n",
    "  this machine (check network/proxy access to https://omnipathdb.org).\n",
    "  See NICHESv1 issue #70 for background.",
    call. = FALSE
  )
}


#' Run an OmnipathR expression with its console logging muted
#'
#' OmnipathR reports failed ancillary downloads (e.g. the OMA/Ensembl species
#' lists) through its own \code{logger}-based console appender, which is not
#' captured by \code{suppressWarnings}/\code{suppressMessages}. Since NICHES
#' handles those failures itself (issue #70), this helper temporarily raises the
#' OmnipathR console log threshold so the expected noise is hidden, restoring the
#' previous level afterwards. Genuine errors still propagate as R conditions.
#'
#' @param expr An expression calling OmnipathR (evaluated lazily inside here).
#'
#' @noRd
.niches_quiet_omnipath <- function(expr) {
  old <- getOption("omnipathr.console_loglevel", default = "success")
  try(OmnipathR::omnipath_set_console_loglevel("fatal"), silent = TRUE)
  on.exit(try(OmnipathR::omnipath_set_console_loglevel(old), silent = TRUE), add = TRUE)
  suppressWarnings(suppressMessages(expr))
}


#' Best-effort OmniPath web server version
#'
#' Reads the OmniPath service banner (\code{https://omnipathdb.org/about}) and
#' extracts the server version for provenance reporting. Returns \code{NA} if the
#' banner cannot be read or parsed - version reporting must never block a load.
#'
#' @return Character server version (e.g. "0.1.5"), or \code{NA_character_}.
#'
#' @noRd
.niches_omnipath_server_version <- function() {
  tryCatch({
    txt <- paste(readLines("https://omnipathdb.org/about", warn = FALSE, n = 3), collapse = " ")
    m <- regmatches(txt, regexec("omnipath-server[[:space:]]+([0-9]+(\\.[0-9]+)*)", txt))[[1]]
    if (length(m) >= 2 && nzchar(m[2])) m[2] else NA_character_
  }, error = function(e) NA_character_)
}


#' Direct OmniPath REST query for ligrecextra interactions
#'
#' Downloads the \code{ligrecextra} ligand-receptor dataset straight from the
#' OmniPath web service for a given NCBI taxon id, bypassing OmnipathR's
#' organism-name translation. The web service performs ortholog translation
#' server-side, so mouse (10090) and rat (10116) return organism-appropriate
#' gene symbols, matching what OmnipathR would otherwise produce.
#'
#' @param organism Integer NCBI taxonomy id.
#'
#' @return A data.frame with (at least) \code{source_genesymbol} and
#'   \code{target_genesymbol} columns.
#'
#' @noRd
.niches_omnipath_rest_ligrec <- function(organism) {

  url <- sprintf(
    "https://omnipathdb.org/interactions?datasets=ligrecextra&organisms=%s&genesymbols=yes",
    as.integer(organism)
  )

  tmp <- tempfile(fileext = ".tsv")
  on.exit(unlink(tmp), add = TRUE)
  utils::download.file(url, destfile = tmp, quiet = TRUE, mode = "wb")

  df <- utils::read.delim(
    tmp,
    sep = "\t",
    header = TRUE,
    quote = "",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  if (!all(c("source_genesymbol", "target_genesymbol") %in% colnames(df)))
    stop("OmniPath REST response did not contain the expected gene symbol columns.")

  # Drop rows without gene symbols on either side (mirrors usable interactions).
  df <- df[nzchar(df$source_genesymbol) & nzchar(df$target_genesymbol) &
             !is.na(df$source_genesymbol) & !is.na(df$target_genesymbol), , drop = FALSE]

  df
}
