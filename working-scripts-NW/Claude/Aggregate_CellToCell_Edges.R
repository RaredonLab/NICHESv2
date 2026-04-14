# ==============================================================================
# FILE:    Aggregate_CellToCell_Edges.R
# PURPOSE: Aggregate cell-to-cell edge signals per cell from an EdgeLRResult
#          object. Computes sending and receiving statistics (sum, mean, max)
#          across all neighbor edges for every cell × LR pair combination.
#
# FUNCTIONS:
#   Aggregate_CellToCell_Edges()  — main function
#
# DEPS: data.table (>= 1.14)
# ==============================================================================

suppressPackageStartupMessages(library(data.table))


# ==============================================================================
#' Aggregate_CellToCell_Edges
#'
#' @title Aggregate Cell-to-Cell Edge Signals Per Cell
#'
#' @description
#' Aggregates the \code{$edges} slot of an \code{EdgeLRResult} object to
#' produce per-cell sending and receiving statistics. All three summary
#' statistics (sum, mean, max) are always computed and stored together, so
#' you can retrieve whichever you need without re-running the function.
#'
#' @section Output structure:
#' A named list with two data.tables:
#' \describe{
#'   \item{\code{$sending}}{
#'     One row per \strong{cell × LRM combination} — every cell in
#'     \code{ctc$cell.meta} appears for every LR pair in \code{ctc$lr.meta}.
#'     For cell A with neighbors B, C, D, the value for LRM "VEGFA.KDR" is
#'     the aggregation of edge(A→B) + edge(A→C) + edge(A→D) (when
#'     \code{stat = "sum"}). Cells with no signal for a given LRM get 0.
#'     Total rows = n.cells × n.lr.
#'     Columns: \code{cell}, \code{LR}, \code{ligand}, \code{receptor},
#'     \code{n.edges}, \code{edge.sum}, \code{edge.mean}, \code{edge.max},
#'     \code{edge.stat}, plus all columns from \code{ctc$cell.meta}.
#'   }
#'   \item{\code{$receiving}}{
#'     Same structure — for cell A receiving from neighbors B, C, D, the
#'     value for each LRM is the aggregation of edge(B→A) + edge(C→A)
#'     + edge(D→A). Total rows = n.cells × n.lr.
#'   }
#' }
#'
#' @param ctc     An \code{EdgeLRResult} object from
#'   \code{compute_CellToCell_Spatial()}.
#' @param stat    Character. Determines the \code{edge.stat} alias column.
#'   One of \code{"sum"} (default), \code{"mean"}, or \code{"max"}.
#'   All three statistics are always stored regardless.
#'
#' @return A named list with \code{$sending} and \code{$receiving} data.tables,
#'   both keyed on \code{(cell, LR)}, each with \code{n.cells × n.lr} rows.
# ==============================================================================

Aggregate_CellToCell_Edges <- function(ctc,
                                       stat = c("sum", "mean", "max")) {

  if (!inherits(ctc, "EdgeLRResult"))
    stop("[Aggregate_CellToCell_Edges] 'ctc' must be an EdgeLRResult object.")

  # Guard: $edges slot must exist. Old objects built with computeCellToCellSpatial()
  # stored this as $signals — re-run compute_CellToCell_Spatial() to update.
  if (is.null(ctc$edges))
    stop("[Aggregate_CellToCell_Edges] 'ctc$edges' is NULL.\n",
         "  Your ctc object was likely built with an older version of the function\n",
         "  which stored signals in '$signals' instead of '$edges'.\n",
         "  Please re-run compute_CellToCell_Spatial() to generate an updated object.")

  # Guard: $cell.meta and $lr.meta must exist (new dot-style slot names)
  if (is.null(ctc$cell.meta))
    stop("[Aggregate_CellToCell_Edges] 'ctc$cell.meta' is NULL.\n",
         "  Please re-run compute_CellToCell_Spatial() to generate an updated object.")

  if (is.null(ctc$lr.meta))
    stop("[Aggregate_CellToCell_Edges] 'ctc$lr.meta' is NULL.\n",
         "  Please re-run compute_CellToCell_Spatial() to generate an updated object.")

  stat <- match.arg(stat)
  sig  <- ctc$edges

  # ════════════════════════════════════════════════════════════════════════════
  # SECTION 1 ── Sparse aggregation by (cell × LRM)
  #
  # For each sending cell and each LR pair, sum/mean/max the signal across
  # all neighbor edges. Result is sparse — only (cell, LRM) combinations
  # with at least one non-zero signal edge are present.
  # ════════════════════════════════════════════════════════════════════════════

  sending.sparse <- sig[,
    .(n.edges   = .N,
      edge.sum  = sum(edge),
      edge.mean = mean(edge),
      edge.max  = max(edge)),
    by = .(cell = from, LR)
  ]

  receiving.sparse <- sig[,
    .(n.edges   = .N,
      edge.sum  = sum(edge),
      edge.mean = mean(edge),
      edge.max  = max(edge)),
    by = .(cell = to, LR)
  ]

  # ════════════════════════════════════════════════════════════════════════════
  # SECTION 2 ── Build full cell × LRM grid
  #
  # Cross-join all cells (from cell_meta) with all LR pairs (from lr_meta)
  # to get the complete grid. Then left-join the sparse aggregation onto it
  # so every (cell, LRM) row exists, with zeros for silent combinations.
  # ════════════════════════════════════════════════════════════════════════════

  # All cells (authoritative list from cell.meta)
  all.cells <- ctc$cell.meta[, .(cell)]

  # All LR pairs — use LR column from lr.meta, keep ligand/receptor for output
  all.lr <- ctc$lr.meta[, .(LR, ligand, receptor)]

  # Cross join: every cell paired with every LR pair.
  # merge(..., by = NULL) is the clearest way to express a Cartesian product.
  full.grid <- setDT(merge(
    as.data.frame(all.cells),
    as.data.frame(all.lr),
    by = NULL
  ))

  # ── Left-join sparse aggregation onto full grid ───────────────────────────
  sending   <- sending.sparse[  full.grid, on = c("cell", "LR")]
  receiving <- receiving.sparse[full.grid, on = c("cell", "LR")]

  # Fill NAs (silent cell × LRM combinations) with zeros
  fill.cols <- c("n.edges", "edge.sum", "edge.mean", "edge.max")
  for (col in fill.cols) {
    set(sending,   which(is.na(sending[[col]])),   col, 0)
    set(receiving, which(is.na(receiving[[col]])), col, 0)
  }

  # ════════════════════════════════════════════════════════════════════════════
  # SECTION 3 ── edge.stat convenience alias
  # ════════════════════════════════════════════════════════════════════════════

  stat.col <- paste0("edge.", stat)
  sending[,   edge.stat := .SD[[stat.col]], .SDcols = stat.col]
  receiving[, edge.stat := .SD[[stat.col]], .SDcols = stat.col]

  # ════════════════════════════════════════════════════════════════════════════
  # SECTION 4 ── Join cell metadata
  #
  # Left-join cell.meta (x, y, cell_type, etc.) onto both tables.
  # The aggregated table is the left side so no rows are dropped.
  # ════════════════════════════════════════════════════════════════════════════

  meta.cols <- setdiff(colnames(ctc$cell.meta), "cell")

  if (length(meta.cols) > 0L) {
    sending   <- sending[  ctc$cell.meta, on = "cell", nomatch = NA]
    receiving <- receiving[ctc$cell.meta, on = "cell", nomatch = NA]
  }

  # ════════════════════════════════════════════════════════════════════════════
  # SECTION 5 ── Sort and set keys
  # ════════════════════════════════════════════════════════════════════════════

  setorder(sending,   cell, LR)
  setorder(receiving, cell, LR)

  setkeyv(sending,   c("cell", "LR"))
  setkeyv(receiving, c("cell", "LR"))

  # ════════════════════════════════════════════════════════════════════════════
  # SECTION 6 ── Report
  # ════════════════════════════════════════════════════════════════════════════

  n.cells <- nrow(all.cells)
  n.lr    <- nrow(all.lr)

  message(sprintf(
    "[Aggregate_CellToCell_Edges] %d cells x %d LR pairs = %d rows each | stat = '%s'",
    n.cells, n.lr, n.cells * n.lr, stat
  ))

  return(list(
    sending   = sending,
    receiving = receiving
  ))
}
