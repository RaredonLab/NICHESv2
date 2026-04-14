# ==============================================================================
# FILE:    compute_CellToCell_Spatial.R
# PURPOSE: Compute cell-to-cell ligand-receptor signals for spatial
#          transcriptomics data. Returns an EdgeLRResult object — a named list
#          with four slots optimized for metadata enrichment, flexible
#          aggregation, and memory efficiency.
#
# FUNCTIONS:
#   .sparse_emin()              — internal: element-wise min of sparse matrices
#   compute_CellToCell_Spatial()  — main function
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
#' compute_CellToCell_Spatial
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
#'   \item{\code{$edges}}{
#'     \code{data.table} with columns \code{from}, \code{to}, \code{ligand},
#'     \code{receptor}, \code{LR} (\code{ligand.receptor} concatenated with \code{"."}),
#'     \code{x.from}, \code{y.from}, \code{x.to}, \code{y.to} (sending/receiving
#'     coordinates, present when edge list was built by \code{create_Edgelist_Spatial()}),
#'     and \code{edge} (signal value). Contains only non-zero entries
#'     (implicit zeros omitted, equivalent to sparse storage).
#'     Keyed on \code{(from, to, LR)} for fast subsetting.
#'   }
#'   \item{\code{$cell_meta}}{
#'     \code{data.table} with one row per unique cell appearing in any edge.
#'     Contains \code{cell} plus \code{x} and \code{y} if the edge list was
#'     produced by \code{createEdgelistSpatial()} (which embeds coordinates
#'     automatically). Add further metadata columns via a keyed join:
#'     \code{obj$cell_meta[your_meta, on = "cell"]}.
#'   }
#'   \item{\code{$lr_meta}}{
#'     \code{data.table} with one row per valid LR pair. Carries \code{ligand},
#'     \code{receptor}, and \code{LR} (\code{"ligand.receptor"}) columns,
#'     plus any additional columns present in
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
#' obj$edges[, .(total_sent = sum(edge)), by = from]
#'
#' # After adding cell_type to cell_meta:
#' obj$edges[
#'   obj$cell.meta, on = .(from = cell)
#' ][, .(mean_signal = mean(edge)), by = sending_cell_type]
#'
#' # Sum signal per LR pair (using LR column)
#' obj$edges[, .(total = sum(edge)), by = LR]
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
#' @param n.cores      CPU cores for parallel LR-chunk computation.
#'   Default = 1. Fork-based on Unix/macOS; serial fallback on Windows.
#'
#' @return An \code{EdgeLRResult} object (named list). See Output structure.
# ==============================================================================

compute_CellToCell_Spatial <- function(count.mtx,
                                       edge.list,
                                       LRM.list,
                                       ligand.col   = "Ligand.ApprovedSymbol",
                                       receptor.col = "Receptor.ApprovedSymbol",
                                       method       = c("product", "minimum"),
                                       n.cores      = 1L) {

  method  <- match.arg(method)
  n.cores <- as.integer(n.cores)

  # ════════════════════════════════════════════════════════════════════════════
  # SECTION 1 ── Validate & filter LR pairs against available genes
  # ════════════════════════════════════════════════════════════════════════════

  ligands   <- LRM.list[[ligand.col]]
  receptors <- LRM.list[[receptor.col]]

  if (is.null(ligands) || is.null(receptors))
    stop("[compute_CellToCell_Spatial] ligand.col or receptor.col not found in LRM.list.")

  genes.in.mtx <- rownames(count.mtx)
  valid.lr     <- ligands %in% genes.in.mtx & receptors %in% genes.in.mtx

  n.dropped <- sum(!valid.lr)
  if (n.dropped > 0)
    message(sprintf(
      "[compute_CellToCell_Spatial] Dropping %d LR pairs: genes absent from count matrix.",
      n.dropped
    ))

  LRM.valid  <- LRM.list[valid.lr, , drop = FALSE]   # keep full row for lr.meta
  ligands    <- ligands[valid.lr]
  receptors  <- receptors[valid.lr]
  n.lr       <- length(ligands)

  if (n.lr == 0L)
    stop("[compute_CellToCell_Spatial] No valid LR pairs remain after gene filtering.")

  # ════════════════════════════════════════════════════════════════════════════
  # SECTION 2 ── Validate & filter edges against available cell barcodes
  # ════════════════════════════════════════════════════════════════════════════

  cells.in.mtx   <- colnames(count.mtx)
  valid.edges    <- edge.list$from %in% cells.in.mtx & edge.list$to %in% cells.in.mtx
  n.edge.dropped <- sum(!valid.edges)

  if (n.edge.dropped > 0)
    message(sprintf(
      "[compute_CellToCell_Spatial] Dropping %d edges: barcodes absent from count matrix.",
      n.edge.dropped
    ))

  edge.list     <- edge.list[valid.edges, ]
  sending.cells <- edge.list$from
  receiv.cells  <- edge.list$to
  n.edges       <- length(sending.cells)

  if (n.edges == 0L)
    stop("[compute_CellToCell_Spatial] No valid edges remain after barcode filtering.")

  message(sprintf(
    "[compute_CellToCell_Spatial] Computing %d edges x %d LR pairs | method = '%s' | cores = %d",
    n.edges, n.lr, method, n.cores
  ))

  # ════════════════════════════════════════════════════════════════════════════
  # SECTION 3 ── Parallel computation over LR pair chunks
  #
  # Same strategy as compute_edge_LR_matrix: split LR pairs into n.cores
  # chunks, each worker computes its slice of the edge signal matrix.
  # Workers return a sparse chunk.LR x n.edges matrix, transposed to
  # n.edges x chunk.LR before returning.
  # ════════════════════════════════════════════════════════════════════════════

  lr.idx.chunks <- if (n.cores == 1L) {
    list(seq_len(n.lr))
  } else {
    split(
      seq_len(n.lr),
      cut(seq_len(n.lr), breaks = n.cores, labels = FALSE)
    )
  }

  .process_chunk <- function(chunk_idx) {
    lig_sub <- count.mtx[ligands[chunk_idx],   sending.cells, drop = FALSE]
    rec_sub <- count.mtx[receptors[chunk_idx], receiv.cells,  drop = FALSE]

    chunk_signal <- if (method == "product") {
      lig_sub * rec_sub
    } else {
      .sparse_emin(lig_sub, rec_sub)
    }

    t(chunk_signal)   # n.edges x chunk.LR
  }

  if (n.cores > 1L) {
    if (.Platform$OS.type == "unix") {
      chunk.results <- parallel::mclapply(lr.idx.chunks, .process_chunk,
                                          mc.cores = n.cores)
    } else {
      message("[compute_CellToCell_Spatial] Fork parallelism unavailable on Windows. Running serially.")
      chunk.results <- lapply(lr.idx.chunks, .process_chunk)
    }
  } else {
    chunk.results <- lapply(lr.idx.chunks, .process_chunk)
  }

  # ════════════════════════════════════════════════════════════════════════════
  # SECTION 4 ── Convert sparse matrix chunks to long-format data.table
  #
  # Instead of cbind-ing chunks into a full n.edges x n.lr matrix and then
  # converting, we convert each chunk to triplets immediately and rbind them.
  # This avoids ever materializing the full matrix in memory.
  #
  # Each chunk produces a TsparseMatrix; we extract (i, j, x) triplets and
  # map them back to (from, to, ligand, receptor, edge) using the chunk's
  # LR index offset.
  # ════════════════════════════════════════════════════════════════════════════

  edge.list.tmp <- vector("list", length(chunk.results))

  for (k in seq_along(chunk.results)) {
    chunk.idx <- lr.idx.chunks[[k]]          # which LR pairs this chunk covers
    cmat      <- as(chunk.results[[k]], "TsparseMatrix")   # n.edges x chunk.LR

    if (length(cmat@x) == 0L) {
      edge.list.tmp[[k]] <- NULL
      next
    }

    # i = edge index (0-based), j = within-chunk LR index (0-based)
    lr.global.idx <- chunk.idx[cmat@j + 1L]   # map back to global LR indices

    edge.list.tmp[[k]] <- data.table(
      from     = sending.cells[cmat@i + 1L],
      to       = receiv.cells[cmat@i + 1L],
      ligand   = ligands[lr.global.idx],
      receptor = receptors[lr.global.idx],
      edge     = cmat@x
    )
  }

  # Free chunk matrices before rbind to reduce peak memory
  rm(chunk.results)

  edge <- rbindlist(edge.list.tmp, use.names = TRUE)
  rm(edge.list.tmp)

  # Add LR convenience column: "ligand.receptor"
  edge[, LR := paste(ligand, receptor, sep = ".")]

  # ── Attach sending/receiving coordinates to edge table ───────────────────────
  # edge.list carries x.from, y.from, x.to, y.to when produced by
  # create_Edgelist_Spatial(). Join them onto edge so every edge row
  # carries the full spatial context of its edge.
  has.coords <- all(c("x.from", "y.from", "x.to", "y.to") %in% colnames(edge.list))
  if (has.coords) {
    edge.coords <- unique(as.data.table(edge.list)[, .(from, to, x.from, y.from, x.to, y.to)])
    setkeyv(edge.coords, c("from", "to"))
    edge <- edge.coords[edge, on = c("from", "to")]
  }

  # Key for fast subsetting and joining
  setkeyv(edge, c("from", "to", "LR"))

  n.nonzero <- nrow(edge)

  message(sprintf(
    "[compute_CellToCell_Spatial] Done. %d non-zero edge signals.",
    n.nonzero
  ))

  # ════════════════════════════════════════════════════════════════════════════
  # SECTION 5 ── Build cell.meta and lr.meta slots
  # ════════════════════════════════════════════════════════════════════════════

  # cell.meta: one row per unique cell appearing in any edge.
  # Coordinates are extracted automatically from edge.list columns x.from, y.from,
  # x.to, y.to — these are embedded by create_Edgelist_Spatial().
  all.cells <- unique(c(sending.cells, receiv.cells))
  cell.meta <- data.table(cell = all.cells)

  if (all(c("x.from", "y.from", "x.to", "y.to") %in% colnames(edge.list))) {
    coords.from  <- unique(as.data.table(edge.list)[, .(cell = from, x = x.from, y = y.from)])
    coords.to    <- unique(as.data.table(edge.list)[, .(cell = to,   x = x.to,   y = y.to)])
    coord.lookup <- unique(rbindlist(list(coords.from, coords.to)), by = "cell")
    cell.meta <- coord.lookup[cell.meta, on = "cell"]
  } else {
    message("[compute_CellToCell_Spatial] Coordinates not found in edge.list. ",
            "Use create_Edgelist_Spatial() to embed them automatically.")
  }

  setkeyv(cell.meta, "cell")

  # lr.meta: one row per valid LR pair, carrying all columns from LRM.list
  # Rename user-specified columns to canonical names for consistency
  lr.meta <- as.data.table(LRM.valid)
  setnames(lr.meta, old = c(ligand.col, receptor.col),
                    new = c("ligand",    "receptor"),
                    skip_absent = FALSE)
  # Deduplicate (LRM.list may have duplicate pairs)
  lr.meta <- unique(lr.meta, by = c("ligand", "receptor"))
  # Add LR convenience column: "ligand.receptor"
  lr.meta[, LR := paste(ligand, receptor, sep = ".")]
  setkeyv(lr.meta, c("ligand", "receptor"))

  # ════════════════════════════════════════════════════════════════════════════
  # SECTION 6 ── Assemble and return EdgeLRResult object
  # ════════════════════════════════════════════════════════════════════════════

  result <- structure(
    list(
      edges     = edge,
      cell.meta = cell.meta,
      lr.meta   = lr.meta,
      params    = list(
        method     = method,
        n.cells    = length(all.cells),
        n.edges    = n.edges,
        n.lr       = n.lr,
        n.nonzero  = n.nonzero,
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
  cat(sprintf("  Cells    : %d\n",   p$n.cells))
  cat(sprintf("  Edges    : %d\n",   p$n.edges))
  cat(sprintf("  LR pairs : %d\n",   p$n.lr))
  cat(sprintf("  Edge signals: %d non-zero\n", p$n.nonzero))
  cat(sprintf("  Method   : %s\n",   p$method))
  cat(sprintf("  Created  : %s\n",   format(p$created_at, "%Y-%m-%d %H:%M:%S")))
  cat("\n  Slots: $edges, $cell.meta, $lr.meta, $params\n")
  cat("  Cell metadata columns:", paste(colnames(x$cell.meta), collapse = ", "), "\n")
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

  sig <- x$edgess

  # Build factor levels to define row and column indices
  edge_ids <- paste(sig$from, sig$to, sep = "->")
  lr_ids   <- sig$LR

  edge_levels <- unique(edge_ids)
  lr_levels   <- unique(lr_ids)

  i <- match(edge_ids, edge_levels)
  j <- match(lr_ids,   lr_levels)

  mat <- Matrix::sparseMatrix(
    i        = i,
    j        = j,
    x        = sig$edge,
    dims     = c(length(edge_levels), length(lr_levels)),
    dimnames = list(edge_levels, lr_levels)
  )

  as(mat, "CsparseMatrix")
}
