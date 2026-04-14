# ==============================================================================
# FILE:    Aggregate_CellToCell_Edges.R
# PURPOSE: Aggregate per-cell sending and receiving LRM signals from a
#          NICHESObject. Produces full cell × LRM grids (silent combinations
#          zero-filled) for both the sending and receiving perspectives.
#
# FUNCTIONS:
#   Aggregate_CellToCell_Edges()  — main function
#
# DEPS: data.table (>= 1.14)
# ==============================================================================


# ==============================================================================
#' Aggregate Cell-to-Cell Edge Signals Per Cell
#'
#' @description
#' Collapses the \code{$edge.data} slot of a \code{NICHESObject} into
#' per-cell sending and receiving statistics. All three summary statistics
#' (sum, mean, max) are computed and stored in a single pass; the
#' \code{stat} argument determines the \code{score.stat} alias column.
#'
#' The output is a full cell × LRM grid — every cell in \code{$cell.list}
#' paired with every LRM pair in \code{$LRM.list} — with zero-fill for
#' silent combinations.
#'
#' @section Sending vs receiving:
#' \describe{
#'   \item{\strong{Sending}}{For cell A with neighbours B, C, D:
#'     the value for LRM "VEGFA|KDR" is the aggregation of
#'     score(A→B) + score(A→C) + score(A→D).}
#'   \item{\strong{Receiving}}{For cell A receiving from neighbours B, C, D:
#'     the value is the aggregation of score(B→A) + score(C→A) + score(D→A).}
#' }
#'
#' @section Output columns:
#' Both \code{$sending} and \code{$receiving} data.tables share the same
#' column schema:
#' \code{cell}, \code{LRM}, \code{ligand}, \code{receptor},
#' \code{n.edges}, \code{score.sum}, \code{score.mean}, \code{score.max},
#' \code{score.stat} (alias for the chosen \code{stat}),
#' plus all columns from \code{$cell.meta}.
#' Keyed on \code{(cell, LRM)}.
#'
#' @param niches.obj A \code{NICHESObject} from \code{create_NICHESObject()}.
#' @param stat       Character. Alias for \code{score.stat}. One of
#'   \code{"sum"} (default), \code{"mean"}, or \code{"max"}. All three are
#'   always computed and stored regardless of this choice.
#'
#' @return A named list with two elements:
#' \describe{
#'   \item{\code{$sending}}{data.table, n.cells × n.LRM rows, keyed on
#'     \code{(cell, LRM)}.}
#'   \item{\code{$receiving}}{Same structure for receiving perspective.}
#' }
#'
#' @examples
#' \dontrun{
#' agg <- Aggregate_CellToCell_Edges(niches.obj, stat = "mean")
#' head(agg$sending)
#' head(agg$receiving)
#' }
#'
#' @export
Aggregate_CellToCell_Edges <- function(niches.obj,
                                        stat = c("sum", "mean", "max")) {

  if (!is_NICHESObject(niches.obj))
    stop("[Aggregate_CellToCell_Edges] 'niches.obj' must be a NICHESObject.")

  stat <- match.arg(stat)

  ed        <- niches.obj$edge.data       # edge, LRM, score, score.norm
  edge.list <- niches.obj$edge.list       # edge, SendingCell, ReceivingCell
  cell.meta <- niches.obj$cell.meta       # cell, x, y, ...
  LRM.meta  <- niches.obj$LRM.meta        # LRM, ligand, receptor, ...

  # ════════════════════════════════════════════════════════════════════════════
  # SECTION 1 ── Attach SendingCell / ReceivingCell to edge.data
  #
  # We use match() rather than a data.table join because ed is keyed on the
  # composite key (edge, LRM). Passing a composite-keyed table as the i
  # argument of another data.table's [ operator confuses data.table into
  # falling back to [.data.frame, which cannot handle a data.table subscript.
  # match() on the edge character column is unambiguous and equally fast here.
  # ════════════════════════════════════════════════════════════════════════════

  ed.with.cells <- copy(ed)
  idx <- match(ed.with.cells$edge, edge.list$edge)
  ed.with.cells[, SendingCell   := edge.list$SendingCell[idx]]
  ed.with.cells[, ReceivingCell := edge.list$ReceivingCell[idx]]

  # ════════════════════════════════════════════════════════════════════════════
  # SECTION 2 ── Sparse aggregation by (cell × LRM) for send and receive
  # ════════════════════════════════════════════════════════════════════════════

  sending.sparse <- ed.with.cells[,
    .(n.edges    = .N,
      score.sum  = sum(score),
      score.mean = mean(score),
      score.max  = max(score)),
    by = .(cell = SendingCell, LRM)
  ]

  receiving.sparse <- ed.with.cells[,
    .(n.edges    = .N,
      score.sum  = sum(score),
      score.mean = mean(score),
      score.max  = max(score)),
    by = .(cell = ReceivingCell, LRM)
  ]

  # ════════════════════════════════════════════════════════════════════════════
  # SECTION 3 ── Build full cell × LRM grid via cross join
  # ════════════════════════════════════════════════════════════════════════════

  all.cells <- data.table(cell = niches.obj$cell.list)
  all.lrm   <- LRM.meta[, .(LRM, ligand, receptor)]

  # Cartesian product (cross join)
  full.grid <- setDT(merge(
    as.data.frame(all.cells),
    as.data.frame(all.lrm),
    by = NULL
  ))

  # ── Left-join sparse aggregations onto full grid ─────────────────────────
  sending   <- sending.sparse[  full.grid, on = c("cell", "LRM")]
  receiving <- receiving.sparse[full.grid, on = c("cell", "LRM")]

  # ════════════════════════════════════════════════════════════════════════════
  # SECTION 4 ── Zero-fill silent (cell × LRM) combinations
  # ════════════════════════════════════════════════════════════════════════════

  fill.cols <- c("n.edges", "score.sum", "score.mean", "score.max")
  for (col in fill.cols) {
    set(sending,   which(is.na(sending[[col]])),   col, 0)
    set(receiving, which(is.na(receiving[[col]])), col, 0)
  }

  # ════════════════════════════════════════════════════════════════════════════
  # SECTION 5 ── score.stat convenience alias
  # ════════════════════════════════════════════════════════════════════════════

  stat.col <- paste0("score.", stat)
  sending[,   score.stat := .SD[[stat.col]], .SDcols = stat.col]
  receiving[, score.stat := .SD[[stat.col]], .SDcols = stat.col]

  # ════════════════════════════════════════════════════════════════════════════
  # SECTION 6 ── Join cell metadata onto both tables
  # ════════════════════════════════════════════════════════════════════════════

  meta.cols <- setdiff(colnames(cell.meta), "cell")

  if (length(meta.cols) > 0L) {
    sending   <- sending[  cell.meta, on = "cell", nomatch = NA]
    receiving <- receiving[cell.meta, on = "cell", nomatch = NA]
  }

  # ════════════════════════════════════════════════════════════════════════════
  # SECTION 7 ── Sort, key, and report
  # ════════════════════════════════════════════════════════════════════════════

  setorder(sending,   cell, LRM)
  setorder(receiving, cell, LRM)
  setkeyv(sending,   c("cell", "LRM"))
  setkeyv(receiving, c("cell", "LRM"))

  n.cells <- nrow(all.cells)
  n.lrm   <- nrow(all.lrm)

  message(sprintf(
    "[Aggregate_CellToCell_Edges] %d cells × %d LRM pairs = %d rows each | stat = '%s'",
    n.cells, n.lrm, n.cells * n.lrm, stat
  ))

  list(sending = sending, receiving = receiving)
}
