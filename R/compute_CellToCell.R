# ==============================================================================
# FILE:    compute_CellToCell.R
# PURPOSE: Universal cell-to-cell ligand-receptor signal scoring engine.
#          Accepts any edge list produced by create_Edgelist_Spatial(),
#          create_Edgelist_PseudoBulk(), or create_Edgelist_SamplingByType().
#          Supports single-gene and multi-gene complex ligands / receptors
#          (underscore-delimited, e.g. "COL1A1_COL1A2"). Returns an
#          EdgeLRResult S3 object suitable for assembly into a NICHESObject
#          via create_NICHESObject().
#
# FUNCTIONS:
#   .sparse_emin()        (internal) element-wise min of sparse matrices
#   compute_CellToCell()  main function
#   print.EdgeLRResult    S3 print method
#   as_matrix()           S3 generic + method to reconstruct dgCMatrix
#
# OUTPUT STRUCTURE (class "EdgeLRResult"):
#   $edges     — data.table, one row per non-zero (edge × LRM)
#   $cell.meta — data.table, one row per cell/type (with x, y if available)
#   $lr.meta   — data.table, one row per valid LRM pair
#   $params    — list, run provenance
#
# DEPS: Matrix (>= 1.5), data.table (>= 1.14), parallel
# ==============================================================================


# ==============================================================================
# INTERNAL: Element-wise minimum of two non-negative sparse matrices
# ==============================================================================

#' @noRd
.sparse_emin <- function(A, B) {

  stopifnot(identical(dim(A), dim(B)))

  At <- as(A, "TsparseMatrix")
  Bt <- as(B, "TsparseMatrix")

  dtA <- data.table(i = At@i, j = At@j, vA = At@x, key = c("i", "j"))
  dtB <- data.table(i = Bt@i, j = Bt@j, vB = Bt@x, key = c("i", "j"))

  common <- dtA[dtB, nomatch = 0L]

  if (nrow(common) == 0L)
    return(Matrix::sparseMatrix(
      i = integer(0), j = integer(0), x = numeric(0),
      dims = dim(A), dimnames = dimnames(A)
    ))

  Matrix::sparseMatrix(
    i        = common$i + 1L,
    j        = common$j + 1L,
    x        = pmin(common$vA, common$vB),
    dims     = dim(A),
    dimnames = dimnames(A)
  )
}


# ==============================================================================
#' Compute Cell-to-Cell Ligand-Receptor Signals
#'
#' @description
#' Universal scoring engine for all NICHESv2 edge modes. For each edge and each
#' valid LR pair, computes a signal score. The edge list may come from any of
#' the three edgelist constructors — spatial, pseudobulk, or type-sampled — and
#' the count matrix is expected to match accordingly (cell barcodes as column
#' names for spatial and sampled modes; cell type names for pseudobulk mode).
#'
#' Supports \strong{multi-gene complexes}: a ligand or receptor represented as
#' \code{"gene1_gene2"} is split on \code{"_"} and all component expressions
#' are combined via the chosen method (product or minimum) \emph{across all
#' ligand and receptor components together}.
#'
#' Built-in LRM databases can be loaded by name:
#' \itemize{
#'   \item \code{"fantom5"}  — Ramilowski et al. 2015 (requires \pkg{OmnipathR})
#'   \item \code{"omnipath"} — OmniPath intercell network (requires \pkg{OmnipathR})
#' }
#'
#' @section Complex scoring example:
#' LRM \code{"COL1A1_COL1A2 → ITGB1"} with \code{method = "product"}:
#' \deqn{score = \text{expr}_{send}(COL1A1) \times \text{expr}_{send}(COL1A2) \times \text{expr}_{recv}(ITGB1)}
#' With \code{method = "minimum"}:
#' \deqn{score = \min(\text{expr}_{send}(COL1A1),\ \text{expr}_{send}(COL1A2),\ \text{expr}_{recv}(ITGB1))}
#'
#' @param count.mtx    Gene × cell (or gene × type for pseudobulk) expression
#'   matrix (sparse or dense). Row names = gene symbols; column names = cell
#'   barcodes or cell type labels matching \code{edge.list$from} / \code{$to}.
#' @param edge.list    \code{data.table} from \code{create_Edgelist_Spatial()},
#'   \code{create_Edgelist_PseudoBulk()}, or
#'   \code{create_Edgelist_SamplingByType()}; or any data.frame with \code{from}
#'   and \code{to} columns whose values match \code{colnames(count.mtx)}.
#' @param LRM.db       LR pair table (data.frame) \strong{or} a character string:
#'   \code{"fantom5"} or \code{"omnipath"} (loaded via \code{load_LRM_database()}).
#' @param ligand.col   Column name for ligand symbols in \code{LRM.db}.
#'   Default: \code{"ligand"}.
#' @param receptor.col Column name for receptor symbols in \code{LRM.db}.
#'   Default: \code{"receptor"}.
#' @param method       Scoring method. \code{"product"} (default) or
#'   \code{"minimum"}. See Complex scoring section.
#' @param n.cores      CPU cores for parallel LRM-chunk processing. Default 1.
#'   Uses socket clusters (\code{parallel::makeCluster}) — works on Windows,
#'   macOS, and Linux without OS restrictions.
#'
#' @return An \code{EdgeLRResult} object with slots:
#' \describe{
#'   \item{\code{$edges}}{data.table: \code{from}, \code{to}, \code{ligand},
#'     \code{receptor}, \code{LRM} (\code{"ligand|receptor"}), \code{x.from},
#'     \code{y.from}, \code{x.to}, \code{y.to} (spatial only), \code{score}.
#'     One row per non-zero (edge × LRM). Keyed on \code{(from, to, LRM)}.}
#'   \item{\code{$cell.meta}}{data.table: \code{cell}, \code{x}, \code{y}
#'     (when coordinates are embedded in \code{edge.list}). Keyed on
#'     \code{cell}.}
#'   \item{\code{$lr.meta}}{data.table: \code{ligand}, \code{receptor},
#'     \code{LRM}, plus any extra columns from \code{LRM.db}. Keyed on
#'     \code{(ligand, receptor)}.}
#'   \item{\code{$params}}{Named list: \code{method}, \code{n.cells},
#'     \code{n.edges}, \code{n.lr}, \code{n.nonzero}, \code{created.at}.}
#' }
#'
#' @seealso \code{\link{create_Edgelist_Spatial}},
#'   \code{\link{create_Edgelist_PseudoBulk}},
#'   \code{\link{create_Edgelist_SamplingByType}},
#'   \code{\link{create_NICHESObject}},
#'   \code{\link{load_LRM_database}}
#'
#' @export
compute_CellToCell <- function(count.mtx,
                                edge.list,
                                LRM.db,
                                ligand.col   = "ligand",
                                receptor.col = "receptor",
                                method       = c("product", "minimum"),
                                n.cores      = 1L) {

  method  <- match.arg(method)
  n.cores <- as.integer(n.cores)

  # ── Resolve LRM.db (accept name string or data.frame) ─────────────────────
  if (is.character(LRM.db) && length(LRM.db) == 1L) {
    LRM.db <- load_LRM_database(LRM.db, verbose = TRUE)
    ligand.col   <- "ligand"
    receptor.col <- "receptor"
  }

  if (!is.data.frame(LRM.db) && !is.data.table(LRM.db))
    stop("[compute_CellToCell] 'LRM.db' must be a data.frame or a ",
         "character string ('fantom5' or 'omnipath').")

  # ════════════════════════════════════════════════════════════════════════════
  # SECTION 1 ── Validate & filter LRM pairs
  #
  # A pair is valid only when ALL gene components of BOTH ligand and receptor
  # are present as row names in count.mtx.
  # ════════════════════════════════════════════════════════════════════════════

  ligands   <- as.character(LRM.db[[ligand.col]])
  receptors <- as.character(LRM.db[[receptor.col]])

  if (is.null(ligands) || is.null(receptors))
    stop("[compute_CellToCell] 'ligand.col' or 'receptor.col' not found in LRM.db.")

  genes.in.mtx <- rownames(count.mtx)

  valid.lr <- mapply(
    function(lig, rec) .genes_available(lig, genes.in.mtx) && .genes_available(rec, genes.in.mtx),
    ligands, receptors
  )

  n.dropped <- sum(!valid.lr)
  if (n.dropped > 0L)
    message(sprintf(
      "[compute_CellToCell] Dropping %d LR pair(s): gene(s) absent from count matrix.",
      n.dropped
    ))

  LRM.valid  <- LRM.db[valid.lr, , drop = FALSE]
  ligands    <- ligands[valid.lr]
  receptors  <- receptors[valid.lr]
  n.lr       <- length(ligands)

  if (n.lr == 0L)
    stop("[compute_CellToCell] No valid LR pairs remain after gene filtering.")

  # ════════════════════════════════════════════════════════════════════════════
  # SECTION 2 ── Validate & filter edges
  # ════════════════════════════════════════════════════════════════════════════

  cells.in.mtx   <- colnames(count.mtx)
  valid.edges    <- edge.list$from %in% cells.in.mtx & edge.list$to %in% cells.in.mtx
  n.edge.dropped <- sum(!valid.edges)

  if (n.edge.dropped > 0L)
    message(sprintf(
      "[compute_CellToCell] Dropping %d edge(s): identifiers absent from count matrix.",
      n.edge.dropped
    ))

  edge.list     <- edge.list[valid.edges, ]
  sending.cells <- edge.list$from
  receiv.cells  <- edge.list$to
  n.edges       <- length(sending.cells)

  if (n.edges == 0L)
    stop("[compute_CellToCell] No valid edges remain after identifier filtering.")

  message(sprintf(
    "[compute_CellToCell] %d edges × %d LR pairs | method = '%s' | cores = %d",
    n.edges, n.lr, method, n.cores
  ))

  # ════════════════════════════════════════════════════════════════════════════
  # SECTION 3 ── Parallel computation over LRM chunks (socket cluster)
  #
  # LRM pairs are split into n.cores chunks; each worker processes one chunk,
  # calling .compute_lrm_score() per pair and returning non-zero (edge, LRM,
  # score) rows for its chunk.
  #
  # Internal scoring helpers are exported from the NICHESv2 namespace via
  # asNamespace(), which works both with devtools::load_all() and installed
  # builds. on.exit() guarantees stopCluster() even if an error occurs mid-run.
  # data.table threads are throttled to 1 per worker to prevent OpenMP
  # contention across socket-cluster workers.
  # ════════════════════════════════════════════════════════════════════════════

  lr.idx.chunks <- if (n.cores == 1L) {
    list(seq_len(n.lr))
  } else {
    split(seq_len(n.lr),
          cut(seq_len(n.lr), breaks = n.cores, labels = FALSE))
  }

  .process_chunk <- function(chunk.idx) {

    chunk.results <- vector("list", length(chunk.idx))

    for (k in seq_along(chunk.idx)) {
      i <- chunk.idx[k]

      score.vec <- .compute_lrm_score(
        lig.str       = ligands[i],
        rec.str       = receptors[i],
        count.mtx     = count.mtx,
        sending.cells = sending.cells,
        receiv.cells  = receiv.cells,
        method        = method
      )

      nonzero <- which(score.vec != 0)
      if (length(nonzero) == 0L) next

      chunk.results[[k]] <- data.table(
        from     = sending.cells[nonzero],
        to       = receiv.cells[nonzero],
        ligand   = ligands[i],
        receptor = receptors[i],
        score    = score.vec[nonzero]
      )
    }

    filtered <- Filter(Negate(is.null), chunk.results)
    if (length(filtered) == 0L) NULL else rbindlist(filtered)
  }

  if (n.cores > 1L) {
    cl <- parallel::makeCluster(n.cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterEvalQ(cl, {
      library(data.table)
      library(Matrix)
      data.table::setDTthreads(1)
    })
    parallel::clusterExport(cl,
      varlist = c("ligands", "receptors", "count.mtx",
                  "sending.cells", "receiv.cells", "method"),
      envir = environment())
    parallel::clusterExport(cl,
      varlist = c(".compute_lrm_score", ".parse_complex", ".genes_available"),
      envir = asNamespace("NICHESv2"))
    chunk.results <- parallel::parLapply(cl, lr.idx.chunks, .process_chunk)
  } else {
    chunk.results <- lapply(lr.idx.chunks, .process_chunk)
  }

  # ════════════════════════════════════════════════════════════════════════════
  # SECTION 4 ── Assemble edge data.table
  # ════════════════════════════════════════════════════════════════════════════

  raw.chunks <- Filter(Negate(is.null), chunk.results)
  rm(chunk.results)

  has.coords <- all(c("x.from", "y.from", "x.to", "y.to") %in% colnames(edge.list))

  # When every score is zero, every chunk worker returns NULL and raw.chunks
  # is empty. rbindlist on an empty list produces a zero-column data.table,
  # so column operations like LRM := paste(ligand, ...) would crash because
  # the columns do not exist. We build the correct zero-row schema explicitly.
  if (length(raw.chunks) == 0L) {
    coord.cols <- if (has.coords)
      list(x.from = numeric(0), y.from = numeric(0),
           x.to   = numeric(0), y.to   = numeric(0))
    else list()
    edges <- do.call(data.table, c(
      list(from = character(0), to = character(0),
           ligand = character(0), receptor = character(0),
           LRM = character(0), score = numeric(0)),
      coord.cols
    ))
  } else {
    edges <- rbindlist(raw.chunks)

    # LRM column: "ligand|receptor"
    edges[, LRM := paste(ligand, receptor, sep = "|")]

    # Attach sending/receiving coordinates if present in edge.list (spatial only)
    if (has.coords) {
      coord.lut <- unique(as.data.table(edge.list)[, .(from, to, x.from, y.from, x.to, y.to)])
      setkeyv(coord.lut, c("from", "to"))
      edges <- coord.lut[edges, on = c("from", "to")]
    }
  }

  setkeyv(edges, c("from", "to", "LRM"))
  n.nonzero <- nrow(edges)

  message(sprintf("[compute_CellToCell] Done. %d non-zero edge signals.", n.nonzero))

  # ════════════════════════════════════════════════════════════════════════════
  # SECTION 5 ── Build cell.meta and lr.meta slots
  # ════════════════════════════════════════════════════════════════════════════

  all.cells <- unique(c(sending.cells, receiv.cells))
  cell.meta <- data.table(cell = all.cells)

  if (has.coords) {
    coords.from  <- unique(as.data.table(edge.list)[, .(cell = from, x = x.from, y = y.from)])
    coords.to    <- unique(as.data.table(edge.list)[, .(cell = to,   x = x.to,   y = y.to)])
    coord.lookup <- unique(rbindlist(list(coords.from, coords.to)), by = "cell")
    cell.meta    <- coord.lookup[cell.meta, on = "cell"]
  }
  setkeyv(cell.meta, "cell")

  # lr.meta: one row per valid LRM pair, standardise column names
  lr.meta <- as.data.table(LRM.valid)
  setnames(lr.meta,
           old = c(ligand.col, receptor.col),
           new = c("ligand", "receptor"),
           skip_absent = FALSE)
  lr.meta <- unique(lr.meta, by = c("ligand", "receptor"))
  lr.meta[, LRM := paste(ligand, receptor, sep = "|")]
  setkeyv(lr.meta, c("ligand", "receptor"))

  # ════════════════════════════════════════════════════════════════════════════
  # SECTION 6 ── Assemble and return EdgeLRResult
  # ════════════════════════════════════════════════════════════════════════════

  structure(
    list(
      edges     = edges,
      cell.meta = cell.meta,
      lr.meta   = lr.meta,
      params    = list(
        method     = method,
        n.cells    = length(all.cells),
        n.edges    = n.edges,
        n.lr       = n.lr,
        n.nonzero  = n.nonzero,
        created.at = Sys.time()
      )
    ),
    class = "EdgeLRResult"
  )
}


# ==============================================================================
# S3 METHODS FOR EdgeLRResult
# ==============================================================================

#' @export
print.EdgeLRResult <- function(x, ...) {
  p <- x$params
  cat("EdgeLRResult\n")
  cat(sprintf("  Cells        : %d\n", p$n.cells))
  cat(sprintf("  Edges        : %d\n", p$n.edges))
  cat(sprintf("  LR pairs     : %d\n", p$n.lr))
  cat(sprintf("  Non-zero     : %d\n", p$n.nonzero))
  cat(sprintf("  Method       : %s\n", p$method))
  cat(sprintf("  Created      : %s\n", format(p$created.at, "%Y-%m-%d %H:%M:%S")))
  cat("\n  Slots: $edges, $cell.meta, $lr.meta, $params\n")
  invisible(x)
}


#' Reconstruct a dgCMatrix (edges × LR pairs) from an EdgeLRResult on demand.
#'
#' The matrix is NOT stored in the object — it is recomputed each call.
#' Use only when a downstream function requires a matrix; for most operations
#' work directly with \code{$edges}.
#'
#' @param x   An \code{EdgeLRResult} object.
#' @param ... Unused.
#' @return A \code{dgCMatrix} with rows = unique edges, columns = LRM pairs.
#' @export
as_matrix <- function(x, ...) UseMethod("as_matrix")

#' @export
as_matrix.EdgeLRResult <- function(x, ...) {
  sig        <- x$edges
  edge.ids   <- paste(sig$from, sig$to, sep = "|")
  edge.lvls  <- unique(edge.ids)
  lr.lvls    <- unique(sig$LRM)

  mat <- Matrix::sparseMatrix(
    i        = match(edge.ids, edge.lvls),
    j        = match(sig$LRM,  lr.lvls),
    x        = sig$score,
    dims     = c(length(edge.lvls), length(lr.lvls)),
    dimnames = list(edge.lvls, lr.lvls)
  )
  as(mat, "CsparseMatrix")
}
