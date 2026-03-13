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
#'     One row per (cell x LR pair) combination where the cell is a sender.
#'     Columns: \code{cell}, \code{ligand}, \code{receptor}, \code{n_edges},
#'     \code{signal_sum}, \code{signal_mean}, \code{signal_max},
#'     \code{signal_stat} (alias for the chosen stat), plus all columns
#'     from \code{ctc$cell_meta} (e.g. cell type, sample, cluster).
#'   }
#'   \item{\code{$receiving}}{
#'     Same structure for receiving cells.
#'   }
#' }
#' Cell metadata columns are joined automatically from \code{ctc$cell_meta}.
#' If \code{cell_meta} only contains the bare \code{cell} column (not yet
#' enriched), the output is still correct — no extra columns are added.
#'
#' @section Retrieving a specific statistic:
#' \preformatted{
#' smry <- summarizeCellSignals(ctc)
#'
#' # Get sum of signals sent by each cell (default-like)
#' smry$sending[, .(cell, signal_sum)]
#'
#' # Get mean signal received, sorted
#' smry$receiving[order(-signal_mean), .(cell, signal_mean)]
#'
#' # cell_meta is already joined — use directly:
#' smry$sending[, .(mean_sum = mean(signal_sum)), by = cell_type]
#' }
#'
#' @param ctc     An \code{EdgeLRResult} object from
#'   \code{computeCellToCellSpatial()}.
#' @param stat    Character. Determines the \code{signal_stat} alias column
#'   for convenient sorting. One of \code{"sum"} (default), \code{"mean"},
#'   or \code{"max"}. All three statistics are always stored regardless.
#'
#' @return A named list with \code{$sending} and \code{$receiving} data.tables,
#'   both keyed on \code{(cell, ligand, receptor)}.
# ==============================================================================

summarizeCellSignals <- function(ctc,
                                 stat  = c("sum", "mean", "max")) {

  if (!inherits(ctc, "EdgeLRResult"))
    stop("[summarizeCellSignals] 'ctc' must be an EdgeLRResult object.")

  stat <- match.arg(stat)
  sig  <- ctc$signals

  # ════════════════════════════════════════════════════════════════════════════
  # SECTION 1 ── Compute per (cell x LR pair) summaries
  # ════════════════════════════════════════════════════════════════════════════

  sending <- sig[,
    .(n_edges     = .N,
      signal_sum  = sum(signal),
      signal_mean = mean(signal),
      signal_max  = max(signal)),
    by = .(cell = from, ligand, receptor)
  ]

  receiving <- sig[,
    .(n_edges     = .N,
      signal_sum  = sum(signal),
      signal_mean = mean(signal),
      signal_max  = max(signal)),
    by = .(cell = to, ligand, receptor)
  ]

  # ════════════════════════════════════════════════════════════════════════════
  # SECTION 2 ── Add a 'signal_stat' convenience column for the chosen stat
  #
  # This is just an alias column — it does not replace sum/mean/max.
  # It lets downstream code use a consistent column name regardless of
  # which statistic was chosen, e.g. smry$sending[order(-signal_stat)].
  # ════════════════════════════════════════════════════════════════════════════

  stat_col <- paste0("signal_", stat)
  sending[,   signal_stat := .SD[[stat_col]], .SDcols = stat_col]
  receiving[, signal_stat := .SD[[stat_col]], .SDcols = stat_col]

  # ════════════════════════════════════════════════════════════════════════════
  # SECTION 3 ── Sort and set keys
  # ════════════════════════════════════════════════════════════════════════════

  setorder(sending,   -signal_stat)
  setorder(receiving, -signal_stat)

  setkeyv(sending,   c("cell", "ligand", "receptor"))
  setkeyv(receiving, c("cell", "ligand", "receptor"))

  # ════════════════════════════════════════════════════════════════════════════
  # SECTION 4 ── Join cell metadata
  #
  # Left-join ctc$cell_meta onto both tables on 'cell'.
  # Rows with no metadata match are kept (all.x = TRUE) so cells that
  # were not enriched yet still appear in the output.
  # If cell_meta has only the bare 'cell' column, this is a no-op.
  # ════════════════════════════════════════════════════════════════════════════

  meta_cols <- setdiff(colnames(ctc$cell_meta), "cell")

  if (length(meta_cols) > 0L) {
    sending   <- ctc$cell_meta[sending,   on = "cell"]
    receiving <- ctc$cell_meta[receiving, on = "cell"]
  }

  # Restore key after join (join drops it)
  setkeyv(sending,   c("cell", "ligand", "receptor"))
  setkeyv(receiving, c("cell", "ligand", "receptor"))

  # ════════════════════════════════════════════════════════════════════════════
  # SECTION 5 ── Report
  # ════════════════════════════════════════════════════════════════════════════

  message(sprintf(
    "[summarizeCellSignals] Sending: %d (cell x LR) rows | Receiving: %d (cell x LR) rows | stat = '%s'",
    nrow(sending), nrow(receiving), stat
  ))

  return(list(
    sending   = sending,
    receiving = receiving
  ))
}
