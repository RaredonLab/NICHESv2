# ==============================================================================
# FILE:    smooth_expression.R
# PURPOSE: PCA + kNN smoothing of a count matrix for sampled-mode workflows.
# ==============================================================================

#' Smooth a count matrix by averaging k nearest neighbors in PCA space
#'
#' Replaces the columns of \code{count.mtx} corresponding to \code{barcodes}
#' with the average of their k nearest neighbors in PCA space (computed
#' across ALL cells in \code{count.mtx}). Non-target columns are returned
#' unchanged.
#'
#' Workflow:
#' \enumerate{
#'   \item Library-size normalize + log1p (temporary, used only for PCA).
#'   \item Truncated SVD via \code{irlba::irlba} on all cells.
#'   \item For each barcode in \code{barcodes}, find k nearest neighbors in
#'         PC space via \code{RANN::nn2}.
#'   \item Average the RAW counts (from the original \code{count.mtx}) across
#'         the k + 1 cells (the barcode itself plus its k neighbors).
#'   \item Replace the corresponding columns in \code{count.mtx} with the
#'         smoothed values.
#' }
#'
#' The PCA embedding is computed internally and immediately discarded.
#'
#' Note on \code{irlba} call: \code{fastpath = FALSE} and no \code{center}
#' argument. Omitting centering avoids ALL irlba/Matrix version compatibility
#' issues with \code{as_cholmod_sparse} (removed in Matrix >= 1.6-2). For kNN
#' neighbor finding, Euclidean distances are identical with or without
#' centering, so kNN results are unaffected.
#'
#' @param count.mtx Gene x cell sparse or dense matrix.
#' @param barcodes  Character vector. Subset of \code{colnames(count.mtx)}
#'   whose columns will be smoothed.
#' @param k         Integer. Number of nearest neighbors (excluding the cell
#'   itself). Smoothed value averages \code{k + 1} cells in total.
#' @param n.pcs     Integer. Number of principal components for the neighbor
#'   search. Capped internally at \code{min(n.pcs, ncol - 1, nrow - 1)}.
#' @return A matrix of the same dimensions and type as \code{count.mtx} with
#'   the targeted columns replaced by smoothed raw counts.
#' @export
smooth_expression <- function(count.mtx, barcodes, k = 9L, n.pcs = 30L) {
  
  n.cells      <- ncol(count.mtx)
  n.genes      <- nrow(count.mtx)
  n.pcs.actual <- min(as.integer(n.pcs),
                      n.cells - 1L,
                      n.genes - 1L)
  
  if (n.pcs.actual < 1L)
    stop("[smooth_expression] Not enough cells or genes to compute PCA.")
  
  lib.sizes             <- Matrix::colSums(count.mtx)
  lib.sizes[lib.sizes == 0L] <- 1L
  med.lib               <- stats::median(lib.sizes)
  norm.mtx              <- Matrix::t(Matrix::t(count.mtx) * (med.lib / lib.sizes))
  norm.mtx              <- log1p(norm.mtx)
  
  svd.res   <- irlba::irlba(A  = Matrix::t(norm.mtx),
                            nv = n.pcs.actual, fastpath = FALSE)
  embedding <- svd.res$u %*% diag(svd.res$d, n.pcs.actual, n.pcs.actual)
  rownames(embedding) <- colnames(count.mtx)
  rm(norm.mtx, svd.res)
  
  barcode.idx <- match(barcodes, colnames(count.mtx))
  query.emb   <- embedding[barcode.idx, , drop = FALSE]
  nn.idx      <- RANN::nn2(data  = embedding,
                           query = query.emb,
                           k     = k + 1L)$nn.idx
  rm(embedding)
  
  smoothed.mat <- vapply(
    seq_along(barcodes),
    function(b) {
      Matrix::rowMeans(count.mtx[, nn.idx[b, ], drop = FALSE])
    },
    numeric(n.genes)
  )
  
  for (b in seq_along(barcodes)) {
    count.mtx[, barcode.idx[b]] <- smoothed.mat[, b]
  }
  
  count.mtx
}
