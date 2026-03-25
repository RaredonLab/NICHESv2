# ==============================================================================
# FILE:    compute_edge_LR_matrix.R
# PURPOSE: Compute a cell-edge × ligand-receptor signal matrix from spatial
#          single-cell data, with optional parallelization and two signal
#          aggregation methods (product or minimum).
#
# FUNCTIONS:
#   .sparse_emin()            — internal helper: element-wise min of two
#                               non-negative sparse matrices via triplet join
#   compute_edge_LR_matrix()  — main function
#
# DEPS: Matrix (>= 1.5), data.table (>= 1.14), parallel
# AUTHOR:  [Your Name]
# DATE:    2026-03-12
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
#
# Key property exploited: for non-negative values,
#   min(a, b) == 0   iff   a == 0  OR  b == 0
# Therefore the non-zero PATTERN of min(A, B) equals the INTERSECTION of the
# non-zero patterns of A and B — identical to A * B.  We compute this via an
# inner join on the (row, col) triplet coordinates, which avoids ever
# constructing a dense intermediate matrix.
#
# @param A, B  Non-negative sparse matrices of identical dimensions.
# @return      Sparse matrix of same dimensions; [i,j] = min(A[i,j], B[i,j]).
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
#' compute_edge_LR_matrix
#'
#' @title Compute an Edge × Ligand-Receptor Signal Matrix
#'
#' @description
#' For each spatial edge (sending cell → receiving cell) and each valid
#' ligand-receptor (LR) pair, computes a signal score and returns the full
#' edge × LR sparse matrix.
#'
#' Two scoring methods are supported:
#' \itemize{
#'   \item \code{"product"} (default) — \eqn{score = ligand\_expr \times receptor\_expr}.
#'         Matches NICHES \code{CellToCellSpatial} scoring. Score is zero if
#'         either gene is unexpressed; shrinks sub-unit values further.
#'   \item \code{"minimum"} — \eqn{score = \min(ligand\_expr, receptor\_expr)}.
#'         Also zero when either gene is unexpressed, but does not further
#'         penalize sub-unit expression values. Preferable for normalized data.
#' }
#'
#' @param count.mtx    Gene × cell expression matrix (sparse or dense).
#'   Row names must be gene symbols; column names must be cell barcodes.
#' @param edge.list    Data frame or data.table with columns \code{from} and
#'   \code{to} (cell barcodes). Autocrine self-loops (from == to) are allowed.
#' @param LRM.list     Ligand-receptor pair table (e.g., connectomeDB2025).
#' @param ligand.col   Column name for ligand gene symbols in \code{LRM.list}.
#'   Default: \code{"Ligand.ApprovedSymbol"}.
#' @param receptor.col Column name for receptor gene symbols in \code{LRM.list}.
#'   Default: \code{"Receptor.ApprovedSymbol"}.
#' @param method       Scoring method: \code{"product"} (default) or
#'   \code{"minimum"}. See Description.
#' @param n_cores      Number of CPU cores for parallel computation of LR
#'   chunks. Default = 1 (serial). Fork-based on Unix/macOS; falls back to
#'   serial on Windows with a message.
#'
#' @return A \code{dgCMatrix} of dimensions \code{n_edges × n_LR_pairs}.
#' \itemize{
#'   \item Row names: \code{"sendingCell->receivingCell"}.
#'   \item Column names: \code{"LigandGene_ReceptorGene"}.
#' }
#' Carries attributes \code{method}, \code{n_edges}, \code{n_lr}.
#'
#' @details
#' \strong{Parallelization:}
#' LR pairs are split into \code{n_cores} chunks. Each worker independently
#' extracts its ligand and receptor sub-matrices from \code{count.mtx} and
#' computes its portion of the result. Results are \code{cbind}'d at the end.
#' On Unix/macOS, \code{mclapply} uses fork-based shared memory so
#' \code{count.mtx} is not physically copied per worker. Windows falls back
#' to serial execution.
#'
#' \strong{Sparsity note:}
#' Both methods produce the same non-zero pattern: a result entry is non-zero
#' only when both the ligand (in the sending cell) and receptor (in the
#' receiving cell) are expressed. The \code{"minimum"} method uses an inner
#' join on triplet coordinates to avoid densification.
# ==============================================================================

compute_edge_LR_matrix <- function(count.mtx,
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
    stop("[compute_edge_LR_matrix] ligand.col or receptor.col not found in LRM.list.")

  genes_in_mtx <- rownames(count.mtx)
  valid_lr     <- ligands %in% genes_in_mtx & receptors %in% genes_in_mtx

  n_dropped <- sum(!valid_lr)
  if (n_dropped > 0)
    message(sprintf(
      "[compute_edge_LR_matrix] Dropping %d LR pairs: genes absent from count matrix.",
      n_dropped
    ))

  ligands   <- ligands[valid_lr]
  receptors <- receptors[valid_lr]
  lr_names  <- paste(ligands, receptors, sep = "_")
  n_lr      <- length(ligands)

  if (n_lr == 0L)
    stop("[compute_edge_LR_matrix] No valid LR pairs remain after gene filtering.")

  # ════════════════════════════════════════════════════════════════════════════
  # SECTION 2 ── Validate & filter edges against available cell barcodes
  # ════════════════════════════════════════════════════════════════════════════

  cells_in_mtx   <- colnames(count.mtx)
  valid_edges    <- edge.list$from %in% cells_in_mtx & edge.list$to %in% cells_in_mtx
  n_edge_dropped <- sum(!valid_edges)

  if (n_edge_dropped > 0)
    message(sprintf(
      "[compute_edge_LR_matrix] Dropping %d edges: barcodes absent from count matrix.",
      n_edge_dropped
    ))

  edge.list     <- edge.list[valid_edges, ]
  sending_cells <- edge.list$from
  receiv_cells  <- edge.list$to
  n_edges       <- length(sending_cells)

  if (n_edges == 0L)
    stop("[compute_edge_LR_matrix] No valid edges remain after barcode filtering.")

  message(sprintf(
    "[compute_edge_LR_matrix] Computing %d edges × %d LR pairs | method = '%s' | cores = %d",
    n_edges, n_lr, method, n_cores
  ))

  # ════════════════════════════════════════════════════════════════════════════
  # SECTION 3 ── Parallel computation over LR pair chunks
  #
  # LR pairs are split into n_cores chunks. Each worker:
  #   a) Extracts its ligand rows from count.mtx for sending cells
  #   b) Extracts its receptor rows from count.mtx for receiving cells
  #   c) Computes score (product or min) for its chunk → chunk_LR × n_edges
  #   d) Transposes to n_edges × chunk_LR
  # Results are cbind'd into the final n_edges × n_lr matrix.
  #
  # On Unix/macOS, mclapply fork-shares count.mtx (copy-on-write),
  # so no physical per-worker duplication occurs.
  # Falls back to serial on Windows.
  # ════════════════════════════════════════════════════════════════════════════

  # cut() requires breaks >= 2 to define intervals, so when n_cores == 1
  # we skip it and treat all LR pairs as a single chunk.
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

    t(chunk_signal)   # n_edges × chunk_LR
  }

  if (n_cores > 1L) {
    if (.Platform$OS.type == "unix") {
      chunk_results <- parallel::mclapply(lr_idx_chunks, .process_chunk,
                                          mc.cores = n_cores)
    } else {
      message("[compute_edge_LR_matrix] Fork parallelism unavailable on Windows. Running serially.")
      chunk_results <- lapply(lr_idx_chunks, .process_chunk)
    }
  } else {
    chunk_results <- lapply(lr_idx_chunks, .process_chunk)
  }

  # ════════════════════════════════════════════════════════════════════════════
  # SECTION 4 ── Assemble final matrix and attach metadata
  # ════════════════════════════════════════════════════════════════════════════

  edge_LR_matrix <- do.call(cbind, chunk_results)
  edge_LR_matrix <- as(edge_LR_matrix, "CsparseMatrix")   # canonical dgCMatrix

  rownames(edge_LR_matrix) <- paste(sending_cells, receiv_cells, sep = "->")
  colnames(edge_LR_matrix) <- lr_names

  sparsity <- 1 - nnzero(edge_LR_matrix) / prod(dim(edge_LR_matrix))
  message(sprintf(
    "[compute_edge_LR_matrix] Done. Output: %d × %d | Sparsity: %.1f%%",
    nrow(edge_LR_matrix), ncol(edge_LR_matrix), sparsity * 100
  ))

  attr(edge_LR_matrix, "method")  <- method
  attr(edge_LR_matrix, "n_edges") <- n_edges
  attr(edge_LR_matrix, "n_lr")    <- n_lr

  return(edge_LR_matrix)
}
