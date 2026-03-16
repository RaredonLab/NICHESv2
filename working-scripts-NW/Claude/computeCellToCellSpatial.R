# ==============================================================================
# FILE:    computeCellToCellSpatial.R
# PURPOSE: Compute cell-to-cell ligand-receptor signals for spatial
#          transcriptomics data. Returns an EdgeLRResult object — a named list
#          with four slots optimized for metadata enrichment, flexible
#          aggregation, and memory efficiency.
#
# FUNCTIONS:
#   .sparse_emin()              — internal: element-wise min of sparse matrices
#   computeCellToCellSpatial()  — main function
#   as_matrix.EdgeLRResult()    — reconstruct dgCMatrix from signals on demand
#
# OUTPUT STRUCTURE (class "EdgeLRResult"):
#   $signals   — data.table: one row per non-zero (edge × LR) signal
#   $cell_meta — data.table: one row per unique cell (ready for metadata join)
#   $lr_meta   — data.table: one row per valid LR pair
#   $params    — list: provenance and run parameters
#
# DEPS: Matrix (>= 1.5), data.table (>= 1.14), parallel
# ==============================================================================

suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
  library(parallel)
})


# ==============================================================================
# INTERNAL HELPER — .sparse_emin
#
# Element-wise minimum of two non-negative sparse matrices A and B.
# Non-zero pattern of min(A,B) equals intersection of non-zero patterns of A
# and B (since min(a,b)=0 iff a=0 or b=0 for non-negatives).
# Computed via inner join on triplet coordinates — no densification.
#
# @param A, B  Non-negative sparse matrices of identical dimensions.
# @return      Sparse matrix; [i,j] = min(A[i,j], B[i,j]).
# ==============================================================================

.sparse_emin <- function(A, B) {

  stopifnot(identical(dim(A), dim(B)))

  At <- as(A, "TsparseMatrix")
  Bt <- as(B, "TsparseMatrix")

  dtA <- data.table(i = At@i, j = At@j, vA = At@x, key = c("i", "j"))
  dtB <- data.table(i = Bt@i, j = Bt@j, vB = Bt@x, key = c("i", "j"))

  common <- dtA[dtB, nomatch = 0L]

  if (nrow(common) == 0L) {
    return(Matrix::sparseMatrix(
      i = integer(0), j = integer(0), x = numeric(0),
      dims = dim(A), dimnames = dimnames(A)
    ))
  }

  Matrix::sparseMatrix(
    i        = common$i + 1L,
    j        = common$j + 1L,
    x        = pmin(common$vA, common$vB),
    dims     = dim(A),
    dimnames = dimnames(A)
  )
}


# ==============================================================================
#' computeCellToCellSpatial
#'
#' @title Compute Cell-to-Cell Ligand-Receptor Signals (Spatial)
#'
#' @description
#' For each spatial edge (sending cell -> receiving cell) and each valid
#' ligand-receptor (LR) pair, computes a signal score. Returns an
#' \code{EdgeLRResult} object — a named list designed for downstream metadata
#' enrichment and flexible aggregation.
#'
#' @section Output structure (\code{EdgeLRResult}):
#' \describe{
#'   \item{\code{$signals}}{
#'     \code{data.table} with columns \code{from}, \code{to}, \code{ligand},
#'     \code{receptor}, \code{signal}. Contains only non-zero entries
#'     (implicit zeros omitted, equivalent to sparse storage).
#'     Keyed on \code{(from, to, ligand, receptor)} for fast subsetting.
#'   }
#'   \item{\code{$cell_meta}}{
#'     \code{data.table} with one row per unique cell appearing in any edge.
#'     Contains \code{cell} plus \code{x} and \code{y} if the edge list was
#'     produced by \code{createEdgelistSpatial()} (which embeds coordinates
#'     automatically). Add further metadata columns via a keyed join:
#'     \code{obj$cell_meta[your_meta, on = "cell"]}.
#'   }
#'   \item{\code{$lr_meta}}{
#'     \code{data.table} with one row per valid LR pair. Carries \code{ligand}
#'     and \code{receptor} columns, plus any additional columns present in
#'     \code{LRM.list} (pathway class, database score, etc.).
#'   }
#'   \item{\code{$params}}{
#'     Named list with run provenance: \code{method}, \code{n_cells},
#'     \code{n_edges}, \code{n_lr}, \code{n_nonzero}, \code{sparsity},
#'     \code{created_at}.
#'   }
#' }
#'
#' @section Aggregation examples:
#' \preformatted{
#' # Total signal sent by each cell
#' obj$signals[, .(total_sent = sum(signal)), by = from]
#'
#' # After adding cell_type to cell_meta:
#' obj$signals[
#'   obj$cell_meta, on = .(from = cell)
#' ][, .(mean_signal = mean(signal)), by = sending_cell_type]
#'
#' # Sum signal per LR pair
#' obj$signals[, .(total = sum(signal)), by = .(ligand, receptor)]
#' }
#'
#' @param count.mtx    Gene x cell expression matrix (sparse or dense).
#'   Row names = gene symbols; column names = cell barcodes.
#' @param edge.list    Data frame or data.table with \code{from} and \code{to}
#'   columns (cell barcodes). Self-loops allowed.
#' @param LRM.list     Ligand-receptor pair table (e.g., connectomeDB2025).
#'   Extra columns beyond ligand/receptor are carried into \code{$lr_meta}.
#' @param ligand.col   Column name for ligand symbols in \code{LRM.list}.
#'   Default: \code{"Ligand.ApprovedSymbol"}.
#' @param receptor.col Column name for receptor symbols in \code{LRM.list}.
#'   Default: \code{"Receptor.ApprovedSymbol"}.
#' @param method       \code{"product"} (default, matches NICHES) or
#'   \code{"minimum"}. See \code{compute_edge_LR_matrix} for details.
#' @param n_cores      CPU cores for parallel LR-chunk computation.
#'   Default = 1. Fork-based on Unix/macOS; serial fallback on Windows.
#'
#' @return An \code{EdgeLRResult} object (named list). See Output structure.
# ==============================================================================

computeCellToCellSpatial <- function(count.mtx,
                                     edge.list,
                                     LRM.list,
                                     ligand.col   = "Ligand.ApprovedSymbol",
                                     receptor.col = "Receptor.ApprovedSymbol",
                                     method       = c("product", "minimum"),
                                     n_cores      = 1L) {

  method  <- match.arg(method)
  n_cores <- as.integer(n_cores)

  # ════════════════════════════════════════════════════════════════════════════
  # SECTION 1 ── Validate & filter LR pairs against available genes
  # ════════════════════════════════════════════════════════════════════════════

  ligands   <- LRM.list[[ligand.col]]
  receptors <- LRM.list[[receptor.col]]

  if (is.null(ligands) || is.null(receptors))
    stop("[computeCellToCellSpatial] ligand.col or receptor.col not found in LRM.list.")

  genes_in_mtx <- rownames(count.mtx)
  valid_lr     <- ligands %in% genes_in_mtx & receptors %in% genes_in_mtx

  n_dropped <- sum(!valid_lr)
  if (n_dropped > 0)
    message(sprintf(
      "[computeCellToCellSpatial] Dropping %d LR pairs: genes absent from count matrix.",
      n_dropped
    ))

  LRM.valid  <- LRM.list[valid_lr, , drop = FALSE]   # keep full row for lr_meta
  ligands    <- ligands[valid_lr]
  receptors  <- receptors[valid_lr]
  n_lr       <- length(ligands)

  if (n_lr == 0L)
    stop("[computeCellToCellSpatial] No valid LR pairs remain after gene filtering.")

  # ════════════════════════════════════════════════════════════════════════════
  # SECTION 2 ── Validate & filter edges against available cell barcodes
  # ════════════════════════════════════════════════════════════════════════════

  cells_in_mtx   <- colnames(count.mtx)
  valid_edges    <- edge.list$from %in% cells_in_mtx & edge.list$to %in% cells_in_mtx
  n_edge_dropped <- sum(!valid_edges)

  if (n_edge_dropped > 0)
    message(sprintf(
      "[computeCellToCellSpatial] Dropping %d edges: barcodes absent from count matrix.",
      n_edge_dropped
    ))

  edge.list     <- edge.list[valid_edges, ]
  sending_cells <- edge.list$from
  receiv_cells  <- edge.list$to
  n_edges       <- length(sending_cells)

  if (n_edges == 0L)
    stop("[computeCellToCellSpatial] No valid edges remain after barcode filtering.")

  message(sprintf(
    "[computeCellToCellSpatial] Computing %d edges x %d LR pairs | method = '%s' | cores = %d",
    n_edges, n_lr, method, n_cores
  ))

  # ════════════════════════════════════════════════════════════════════════════
  # SECTION 3 ── Parallel computation over LR pair chunks
  #
  # Same strategy as compute_edge_LR_matrix: split LR pairs into n_cores
  # chunks, each worker computes its slice of the signal matrix.
  # Workers return a sparse chunk_LR x n_edges matrix, transposed to
  # n_edges x chunk_LR before returning.
  # ════════════════════════════════════════════════════════════════════════════

  lr_idx_chunks <- if (n_cores == 1L) {
    list(seq_len(n_lr))
  } else {
    split(
      seq_len(n_lr),
      cut(seq_len(n_lr), breaks = n_cores, labels = FALSE)
    )
  }

  .process_chunk <- function(chunk_idx) {
    lig_sub <- count.mtx[ligands[chunk_idx],   sending_cells, drop = FALSE]
    rec_sub <- count.mtx[receptors[chunk_idx], receiv_cells,  drop = FALSE]

    chunk_signal <- if (method == "product") {
      lig_sub * rec_sub
    } else {
      .sparse_emin(lig_sub, rec_sub)
    }

    t(chunk_signal)   # n_edges x chunk_LR
  }

  if (n_cores > 1L) {
    if (.Platform$OS.type == "unix") {
      chunk_results <- parallel::mclapply(lr_idx_chunks, .process_chunk,
                                          mc.cores = n_cores)
    } else {
      message("[computeCellToCellSpatial] Fork parallelism unavailable on Windows. Running serially.")
      chunk_results <- lapply(lr_idx_chunks, .process_chunk)
    }
  } else {
    chunk_results <- lapply(lr_idx_chunks, .process_chunk)
  }

  # ════════════════════════════════════════════════════════════════════════════
  # SECTION 4 ── Convert sparse matrix chunks to long-format data.table
  #
  # Instead of cbind-ing chunks into a full n_edges x n_lr matrix and then
  # converting, we convert each chunk to triplets immediately and rbind them.
  # This avoids ever materializing the full matrix in memory.
  #
  # Each chunk produces a TsparseMatrix; we extract (i, j, x) triplets and
  # map them back to (from, to, ligand, receptor, signal) using the chunk's
  # LR index offset.
  # ════════════════════════════════════════════════════════════════════════════

  signals_list <- vector("list", length(chunk_results))

  for (k in seq_along(chunk_results)) {
    chunk_idx <- lr_idx_chunks[[k]]          # which LR pairs this chunk covers
    cmat      <- as(chunk_results[[k]], "TsparseMatrix")   # n_edges x chunk_LR

    if (length(cmat@x) == 0L) {
      signals_list[[k]] <- NULL
      next
    }

    # i = edge index (0-based), j = within-chunk LR index (0-based)
    lr_global_idx <- chunk_idx[cmat@j + 1L]   # map back to global LR indices

    signals_list[[k]] <- data.table(
      from     = sending_cells[cmat@i + 1L],
      to       = receiv_cells[cmat@i + 1L],
      ligand   = ligands[lr_global_idx],
      receptor = receptors[lr_global_idx],
      signal   = cmat@x
    )
  }

  # Free chunk matrices before rbind to reduce peak memory
  rm(chunk_results)

  signals <- rbindlist(signals_list, use.names = TRUE)
  rm(signals_list)

  # Key for fast subsetting and joining
  setkeyv(signals, c("from", "to", "ligand", "receptor"))

  n_nonzero <- nrow(signals)
  sparsity  <- 1 - n_nonzero / (n_edges * n_lr)

  message(sprintf(
    "[computeCellToCellSpatial] Done. %d non-zero signals | Sparsity: %.1f%%",
    n_nonzero, sparsity * 100
  ))

  # ════════════════════════════════════════════════════════════════════════════
  # SECTION 5 ── Build cell_meta and lr_meta slots
  # ════════════════════════════════════════════════════════════════════════════

  # cell_meta: one row per unique cell appearing in any edge.
  # If edge.list carries x_from/y_from/x_to/y_to columns (as produced by
  # createEdgelistSpatial), extract x/y here so plotting functions need
  # nothing beyond the ctc object itself.
  all_cells <- unique(c(sending_cells, receiv_cells))
  cell_meta <- data.table(cell = all_cells)

  has_coords <- all(c("x_from", "y_from", "x_to", "y_to") %in% colnames(edge.list))
  if (has_coords) {
    # Build a cell → (x, y) lookup from both the from-side and to-side columns,
    # then deduplicate (both sides carry the same coordinate for a given cell).
    coords_from <- unique(edge.list[, .(cell = from, x = x_from, y = y_from)])
    coords_to   <- unique(edge.list[, .(cell = to,   x = x_to,   y = y_to)])
    coord_lookup <- unique(rbindlist(list(coords_from, coords_to)), by = "cell")
    cell_meta <- coord_lookup[cell_meta, on = "cell"]   # left join: all cells kept
  }

  setkeyv(cell_meta, "cell")

  # lr_meta: one row per valid LR pair, carrying all columns from LRM.list
  # Rename user-specified columns to canonical names for consistency
  lr_meta <- as.data.table(LRM.valid)
  setnames(lr_meta, old = c(ligand.col, receptor.col),
                    new = c("ligand",    "receptor"),
                    skip_absent = FALSE)
  # Deduplicate (LRM.list may have duplicate pairs)
  lr_meta <- unique(lr_meta, by = c("ligand", "receptor"))
  setkeyv(lr_meta, c("ligand", "receptor"))

  # ════════════════════════════════════════════════════════════════════════════
  # SECTION 6 ── Assemble and return EdgeLRResult object
  # ════════════════════════════════════════════════════════════════════════════

  result <- structure(
    list(
      signals   = signals,
      cell_meta = cell_meta,
      lr_meta   = lr_meta,
      params    = list(
        method     = method,
        n_cells    = length(all_cells),
        n_edges    = n_edges,
        n_lr       = n_lr,
        n_nonzero  = n_nonzero,
        sparsity   = sparsity,
        created_at = Sys.time()
      )
    ),
    class = "EdgeLRResult"
  )

  return(result)
}


# ==============================================================================
# S3 METHODS FOR EdgeLRResult
# ==============================================================================

#' Print method: concise summary
print.EdgeLRResult <- function(x, ...) {
  p <- x$params
  cat("EdgeLRResult\n")
  cat(sprintf("  Cells    : %d\n",   p$n_cells))
  cat(sprintf("  Edges    : %d\n",   p$n_edges))
  cat(sprintf("  LR pairs : %d\n",   p$n_lr))
  cat(sprintf("  Signals  : %d non-zero (sparsity %.1f%%)\n",
              p$n_nonzero, p$sparsity * 100))
  cat(sprintf("  Method   : %s\n",   p$method))
  cat(sprintf("  Created  : %s\n",   format(p$created_at, "%Y-%m-%d %H:%M:%S")))
  cat("\n  Slots: $signals, $cell_meta, $lr_meta, $params\n")
  cat("  Cell metadata columns:", paste(colnames(x$cell_meta), collapse = ", "), "\n")
  invisible(x)
}

#' Reconstruct a dgCMatrix (edges x LR) from signals on demand.
#'
#' Use this when a downstream function specifically requires a matrix.
#' The matrix is NOT stored in the object — it is recomputed each call.
#'
#' @param x     An EdgeLRResult object.
#' @param ...   Unused.
#' @return      A dgCMatrix of dimensions n_edges x n_lr_pairs.
as_matrix <- function(x, ...) UseMethod("as_matrix")

as_matrix.EdgeLRResult <- function(x, ...) {

  sig <- x$signals

  # Build factor levels to define row and column indices
  edge_ids <- paste(sig$from, sig$to, sep = "->")
  lr_ids   <- paste(sig$ligand, sig$receptor, sep = "_")

  edge_levels <- unique(edge_ids)
  lr_levels   <- unique(lr_ids)

  i <- match(edge_ids, edge_levels)
  j <- match(lr_ids,   lr_levels)

  mat <- Matrix::sparseMatrix(
    i        = i,
    j        = j,
    x        = sig$signal,
    dims     = c(length(edge_levels), length(lr_levels)),
    dimnames = list(edge_levels, lr_levels)
  )

  as(mat, "CsparseMatrix")
}
