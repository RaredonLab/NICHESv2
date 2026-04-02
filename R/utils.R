# ==============================================================================
# FILE:    utils.R
# PURPOSE: Internal helper functions shared across NICHESv2.
#          Not exported. Prefixed with "." to signal private scope.
# ==============================================================================


# ------------------------------------------------------------------------------
# Complex gene-set parsing
# ------------------------------------------------------------------------------

#' Split a potentially complex gene string (e.g. "gene1_gene2") into parts.
#' Returns a character vector of length >= 1.
#' @noRd
.parse_complex <- function(gene.str) {
  strsplit(gene.str, "_", fixed = TRUE)[[1L]]
}

#' Check that every component of a (possibly complex) gene string is present
#' in a reference gene set. Returns TRUE/FALSE.
#' @noRd
.genes_available <- function(gene.str, gene.set) {
  all(.parse_complex(gene.str) %in% gene.set)
}

#' Compute the LRM score vector (length = n.edges) for one ligand-receptor pair.
#'
#' Supports complexes: "gene1_gene2" → each component expression extracted from
#' the appropriate cell set (sending for ligand, receiving for receptor), then
#' combined via product or element-wise minimum across ALL components.
#'
#' @param lig.str      Character. Ligand gene or complex, e.g. "VEGFA" or "CXCL12_PF4".
#' @param rec.str      Character. Receptor gene or complex, e.g. "KDR" or "ACKR3_CXCR4".
#' @param count.mtx    Gene x cell matrix (sparse or dense).
#' @param sending.cells Character vector, length = n.edges — sending cell barcodes.
#' @param receiv.cells  Character vector, length = n.edges — receiving cell barcodes.
#' @param method        "product" or "minimum".
#' @return Numeric vector of length n.edges (may contain zeros).
#' @noRd
.compute_lrm_score <- function(lig.str, rec.str,
                                count.mtx,
                                sending.cells, receiv.cells,
                                method) {

  lig.parts <- .parse_complex(lig.str)
  rec.parts <- .parse_complex(rec.str)

  # Extract expression vectors: each has length = n.edges
  lig.vecs <- lapply(lig.parts, function(g)
    as.numeric(count.mtx[g, sending.cells, drop = TRUE]))

  rec.vecs <- lapply(rec.parts, function(g)
    as.numeric(count.mtx[g, receiv.cells, drop = TRUE]))

  all.vecs <- c(lig.vecs, rec.vecs)

  if (method == "product") {
    Reduce("*", all.vecs)
  } else {
    Reduce(pmin, all.vecs)
  }
}


# ------------------------------------------------------------------------------
# Object identity checks
# ------------------------------------------------------------------------------

#' @noRd
.is_NICHESObject <- function(x) inherits(x, "NICHESObject")


# ------------------------------------------------------------------------------
# Categorical column detection (used for Joint.* columns in edge.meta)
# ------------------------------------------------------------------------------

#' Return the names of columns in a data.table that are character or factor.
#' Excludes the key columns supplied in \code{exclude}.
#' @noRd
.categorical_cols <- function(dt, exclude = character(0)) {
  candidates <- setdiff(colnames(dt), exclude)
  candidates[vapply(candidates, function(col) {
    v <- dt[[col]]
    is.character(v) || is.factor(v)
  }, logical(1L))]
}


# ------------------------------------------------------------------------------
# Normalisation helpers
# ------------------------------------------------------------------------------

#' Proportional normalisation within group.
#' Divides each value by the group total. Returns 0 when the group total is 0.
#' @noRd
.prop_norm <- function(x) {
  s <- sum(x, na.rm = TRUE)
  if (s == 0) rep(0, length(x)) else x / s
}


# ------------------------------------------------------------------------------
# Neighbourhood smoothing for sampled mode
# ------------------------------------------------------------------------------

#' Replace sampled cell expression values with the average of their
#' k-nearest-neighbour hood in PCA space.
#'
#' Workflow:
#'   1. Library-size normalise + log1p (temporary — used only to compute PCA).
#'   2. Truncated PCA via irlba on all cells in count.mtx.
#'   3. For each barcode in \code{barcodes}, find k nearest neighbours in PC
#'      space using RANN::nn2.
#'   4. Average the RAW counts (from the original count.mtx) across the
#'      k + 1 cells (the barcode itself plus its k neighbours).
#'   5. Replace the corresponding columns in count.mtx with the smoothed
#'      values.  All non-sampled columns are returned unchanged.
#'
#' The PCA embedding is computed internally and immediately discarded — it is
#' never stored or returned.
#'
#' @param count.mtx  Gene × cell sparse or dense matrix.
#'   Row names = gene symbols; column names = cell barcodes.
#' @param barcodes   Character vector. Unique barcodes whose columns will be
#'   smoothed. Must be a subset of \code{colnames(count.mtx)}.
#' @param k          Integer. Number of nearest neighbours (not counting the
#'   cell itself). The smoothed value averages k + 1 cells in total.
#' @param n.pcs      Integer. Number of principal components to compute for the
#'   neighbour search. Capped internally at
#'   \code{min(n.pcs, ncol(count.mtx) - 1, nrow(count.mtx) - 1)}.
#' @return A matrix of the same dimensions and type as \code{count.mtx} with
#'   the sampled barcode columns replaced by their neighbourhood-averaged
#'   raw counts.
#' @noRd
.smooth_expression <- function(count.mtx, barcodes, k, n.pcs) {

  n.cells      <- ncol(count.mtx)
  n.genes      <- nrow(count.mtx)
  n.pcs.actual <- min(as.integer(n.pcs),
                      n.cells - 1L,
                      n.genes - 1L)

  if (n.pcs.actual < 1L)
    stop("[.smooth_expression] Not enough cells or genes to compute PCA.")

  # ── Step 1: library-size normalise + log1p (for PCA only) ─────────────────
  lib.sizes             <- Matrix::colSums(count.mtx)
  lib.sizes[lib.sizes == 0L] <- 1L          # guard against empty cells
  med.lib               <- stats::median(lib.sizes)
  norm.mtx              <- Matrix::t(Matrix::t(count.mtx) * (med.lib / lib.sizes))
  norm.mtx              <- log1p(norm.mtx)  # genes x cells, sparse

  # ── Step 2: truncated SVD via irlba — no explicit centering ─────────────
  # We deliberately omit center= here. For kNN neighbour finding (the only
  # purpose of this embedding), Euclidean distances between points are
  # IDENTICAL with or without centering: subtracting the same mean vector
  # from every row shifts all points by the same offset, leaving pairwise
  # distances unchanged. kNN results are therefore equivalent.
  #
  # Omitting centering avoids ALL irlba/Matrix version compatibility issues:
  # as_cholmod_sparse is only called when irlba handles centering internally,
  # via both its fast C path and Matrix's own sparse %*% operator on some
  # Matrix versions. With no center argument, irlba performs a plain sparse
  # SVD with no CHOLMOD dependency, working on any irlba >= 2.3.0 and any
  # Matrix version.
  svd.res   <- irlba::irlba(A  = Matrix::t(norm.mtx),
                              nv = n.pcs.actual, fastpath = FALSE)
  # u (cells x PCs) weighted by singular values = cell embeddings in PC space
  embedding <- svd.res$u %*% diag(svd.res$d, n.pcs.actual, n.pcs.actual)
  rownames(embedding) <- colnames(count.mtx)
  rm(norm.mtx, svd.res)                      # discard immediately

  # ── Step 3: kNN in PC space ───────────────────────────────────────────────
  # Search within the full embedding; query only the sampled barcodes.
  # nn2 returns k + 1 columns (column 1 is the query cell itself).
  barcode.idx <- match(barcodes, colnames(count.mtx))
  query.emb   <- embedding[barcode.idx, , drop = FALSE]
  nn.idx      <- RANN::nn2(data  = embedding,
                            query = query.emb,
                            k     = k + 1L)$nn.idx   # length(barcodes) x (k+1)
  rm(embedding)

  # ── Step 4 & 5: average raw counts over neighbourhood; replace columns ────
  # Build a genes x length(barcodes) dense matrix of smoothed values.
  # vapply over barcodes keeps peak memory low relative to densifying count.mtx.
  smoothed.mat <- vapply(
    seq_along(barcodes),
    function(b) {
      Matrix::rowMeans(count.mtx[, nn.idx[b, ], drop = FALSE])
    },
    numeric(n.genes)
  )
  # smoothed.mat is genes x length(barcodes), dense

  # Replace sampled columns in-place.  Column assignment on a sparseMatrix
  # coerces the dense vector to a sparse column automatically.
  for (b in seq_along(barcodes)) {
    count.mtx[, barcode.idx[b]] <- smoothed.mat[, b]
  }

  count.mtx
}
