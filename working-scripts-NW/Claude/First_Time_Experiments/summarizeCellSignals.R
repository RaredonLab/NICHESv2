# ==============================================================================
# FILE:    summarizeCellSignals.R
# PURPOSE: Summarize ligand-receptor signals per cell from an EdgeLRResult
#          object. Computes sending and receiving statistics (sum, mean, max)
#          for every cell, and optionally breaks them down by LR pair.
#
# FUNCTIONS:
#   summarizeCellSignals()  — main function
#
# DEPS: data.table (>= 1.14)
# ==============================================================================

suppressPackageStartupMessages(library(data.table))


# ==============================================================================
#' summarizeCellSignals
#'
#' @title Summarize Ligand-Receptor Signals Per Cell
#'
#' @description
#' Aggregates the \code{$signals} slot of an \code{EdgeLRResult} object to
#' produce per-cell sending and receiving statistics. All three summary
#' statistics (sum, mean, max) are always computed and stored together, so
#' you can retrieve whichever you need without re-running the function.
#'
#' @section Output structure:
#' A named list with two data.tables:
#' \describe{
#'   \item{\code{$sending}}{
#'     One row per \strong{cell × LRM combination} — every cell in
#'     \code{ctc$cell_meta} appears for every LR pair in \code{ctc$lr_meta}.
#'     For cell A with neighbors B, C, D, the value for LRM "VEGFA→KDR" is
#'     the aggregation of signal(A→B) + signal(A→C) + signal(A→D) (when
#'     \code{stat = "sum"}). Cells with no signal for a given LRM get 0.
#'     Total rows = n_cells × n_lr_pairs.
#'     Columns: \code{cell}, \code{ligand}, \code{receptor}, \code{n_edges},
#'     \code{signal_sum}, \code{signal_mean}, \code{signal_max},
#'     \code{signal_stat}, plus all columns from \code{ctc$cell_meta}.
#'   }
#'   \item{\code{$receiving}}{
#'     Same structure — for cell A receiving from neighbors B, C, D, the
#'     value for each LRM is the aggregation of signal(B→A) + signal(C→A)
#'     + signal(D→A). Total rows = n_cells × n_lr_pairs.
#'   }
#' }
#'
#' @param ctc     An \code{EdgeLRResult} object from
#'   \code{computeCellToCellSpatial()}.
#' @param stat    Character. Determines the \code{signal_stat} alias column.
#'   One of \code{"sum"} (default), \code{"mean"}, or \code{"max"}.
#'   All three statistics are always stored regardless.
#'
#' @return A named list with \code{$sending} and \code{$receiving} data.tables,
#'   both keyed on \code{(cell, ligand, receptor)}, each with
#'   \code{n_cells × n_lr_pairs} rows.
# ==============================================================================

summarizeCellSignals <- function(ctc,
                                 stat  = c("sum", "mean", "max")) {

  if (!inherits(ctc, "EdgeLRResult"))
    stop("[summarizeCellSignals] 'ctc' must be an EdgeLRResult object.")

  stat <- match.arg(stat)
  sig  <- ctc$signals

  # ════════════════════════════════════════════════════════════════════════════
  # SECTION 1 ── Sparse aggregation by (cell × LRM)
  #
  # For each sending cell and each LR pair, sum/mean/max the signal across
  # all neighbor edges. Result is sparse — only (cell, LRM) combinations
  # with at least one non-zero signal edge are present.
  # ════════════════════════════════════════════════════════════════════════════

  sending_sparse <- sig[,
    .(n_edges     = .N,
      signal_sum  = sum(signal),
      signal_mean = mean(signal),
      signal_max  = max(signal)),
    by = .(cell = from, ligand, receptor)
  ]

  receiving_sparse <- sig[,
    .(n_edges     = .N,
      signal_sum  = sum(signal),
      signal_mean = mean(signal),
      signal_max  = max(signal)),
    by = .(cell = to, ligand, receptor)
  ]

  # ════════════════════════════════════════════════════════════════════════════
  # SECTION 2 ── Build full cell × LRM grid
  #
  # Cross-join all cells (from cell_meta) with all LR pairs (from lr_meta)
  # to get the complete grid. Then left-join the sparse aggregation onto it
  # so every (cell, LRM) row exists, with zeros for silent combinations.
  # ════════════════════════════════════════════════════════════════════════════

  # All cells (authoritative list)
  all_cells <- ctc$cell_meta[, .(cell)]

  # All LR pairs (canonical ligand/receptor columns only)
  all_lr <- ctc$lr_meta[, .(ligand, receptor)]

  # Cross join: every cell paired with every LR pair.
  # merge(..., by = NULL) is the clearest way to express a Cartesian product.
  full_grid <- setDT(merge(
    as.data.frame(all_cells),
    as.data.frame(all_lr),
    by = NULL
  ))

  # ── Left-join sparse aggregation onto full grid ───────────────────────────
  sending   <- sending_sparse[  full_grid, on = c("cell", "ligand", "receptor")]
  receiving <- receiving_sparse[full_grid, on = c("cell", "ligand", "receptor")]

  # Fill NAs (silent cell × LRM combinations) with zeros
  fill_cols <- c("n_edges", "signal_sum", "signal_mean", "signal_max")
  for (col in fill_cols) {
    set(sending,   which(is.na(sending[[col]])),   col, 0)
    set(receiving, which(is.na(receiving[[col]])), col, 0)
  }

  # ════════════════════════════════════════════════════════════════════════════
  # SECTION 3 ── signal_stat convenience alias
  # ════════════════════════════════════════════════════════════════════════════

  stat_col <- paste0("signal_", stat)
  sending[,   signal_stat := .SD[[stat_col]], .SDcols = stat_col]
  receiving[, signal_stat := .SD[[stat_col]], .SDcols = stat_col]

  # ════════════════════════════════════════════════════════════════════════════
  # SECTION 4 ── Join cell metadata
  #
  # Left-join cell_meta (x, y, cell_type, etc.) onto both tables.
  # The signal table is the left side so no rows are dropped.
  # ════════════════════════════════════════════════════════════════════════════

  meta_cols <- setdiff(colnames(ctc$cell_meta), "cell")

  if (length(meta_cols) > 0L) {
    sending   <- sending[  ctc$cell_meta, on = "cell", nomatch = NA]
    receiving <- receiving[ctc$cell_meta, on = "cell", nomatch = NA]
  }

  # ════════════════════════════════════════════════════════════════════════════
  # SECTION 5 ── Sort and set keys
  # ════════════════════════════════════════════════════════════════════════════

  setorder(sending,   cell, ligand, receptor)
  setorder(receiving, cell, ligand, receptor)

  setkeyv(sending,   c("cell", "ligand", "receptor"))
  setkeyv(receiving, c("cell", "ligand", "receptor"))

  # ════════════════════════════════════════════════════════════════════════════
  # SECTION 6 ── Report
  # ════════════════════════════════════════════════════════════════════════════

  n_cells <- nrow(all_cells)
  n_lr    <- nrow(all_lr)

  message(sprintf(
    "[summarizeCellSignals] %d cells x %d LR pairs = %d rows each | stat = '%s'",
    n_cells, n_lr, n_cells * n_lr, stat
  ))

  return(list(
    sending   = sending,
    receiving = receiving
  ))
}
