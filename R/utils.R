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


# ------------------------------------------------------------------------------
# Null-coalescing operator
# ------------------------------------------------------------------------------

#' @noRd
`%||%` <- function(x, y) if (is.null(x)) y else x


# ==============================================================================
# .merge_NICHESObject_list()
# Unexported internal worker. Merges a list of NICHESObjects into one.
# Called by both merge_NICHESObjects() (public API) and
# .create_NICHESObject_by_sample() (sample.col split/merge path).
# ==============================================================================

#' @noRd
.merge_NICHESObject_list <- function(obj.list, verbose = TRUE) {

  # ----------------------------------------------------------
  # 0. Length guard
  # ----------------------------------------------------------
  if (length(obj.list) == 0L)
    stop("obj.list is empty — nothing to merge.", call. = FALSE)

  if (length(obj.list) == 1L) {
    if (verbose) message("Only one object supplied — returning it unchanged.")
    return(obj.list[[1L]])
  }

  # ----------------------------------------------------------
  # 1. All inputs must be NICHESObjects
  # ----------------------------------------------------------
  not.niches <- !vapply(obj.list, is_NICHESObject, logical(1L))
  if (any(not.niches)) {
    bad.idx <- which(not.niches)
    stop(sprintf(
      "Element(s) %s of obj.list are not NICHESObjects.",
      paste(bad.idx, collapse = ", ")
    ), call. = FALSE)
  }

  # ----------------------------------------------------------
  # 2. Mode compatibility — hard error on mismatch
  # ----------------------------------------------------------
  modes <- vapply(obj.list,
                  function(x) attr(x, "mode") %||% NA_character_,
                  character(1L))
  if (length(unique(modes)) > 1L)
    stop(sprintf(
      paste0("Cannot merge NICHESObjects with different modes: %s.\n",
             "All objects must share the same mode (\"spatial\" or \"sampled\")."),
      paste(unique(modes), collapse = ", ")
    ), call. = FALSE)

  # unname() prevents a named vector propagating into attr(obj, "mode")
  merged.mode <- unname(modes[1L])

  # ----------------------------------------------------------
  # 3. Barcode collision check — hard error
  # ----------------------------------------------------------
  all.barcodes <- unlist(lapply(obj.list, function(x) x$cell.list),
                         use.names = FALSE)
  dup.mask <- duplicated(all.barcodes)
  if (any(dup.mask)) {
    dup.barcodes <- unique(all.barcodes[dup.mask])
    n.shown      <- min(length(dup.barcodes), 10L)
    stop(sprintf(
      paste0(
        "Duplicate cell barcodes detected across objects (%d total).\n",
        "First %d: %s\n",
        "Tip: prepend a sample-specific prefix to barcodes before calling\n",
        "create_NICHESObject() (e.g. paste0(sample.id, \"_\", barcode))."
      ),
      length(dup.barcodes), n.shown,
      paste(head(dup.barcodes, n.shown), collapse = ", ")
    ), call. = FALSE)
  }

  # ----------------------------------------------------------
  # 4. $cell.list — concatenate (already validated unique)
  # ----------------------------------------------------------
  merged.cell.list <- all.barcodes

  # ----------------------------------------------------------
  # 5. $cell.data — rbindlist, re-key
  # ----------------------------------------------------------
  merged.cell.data <- rbindlist(lapply(obj.list, function(x) x$cell.data),
                                use.names = TRUE, fill = FALSE)
  if (nrow(merged.cell.data) > 0L)
    setkeyv(merged.cell.data, c("cell", "gene"))

  # ----------------------------------------------------------
  # 6. $cell.meta — rbindlist with fill, re-key
  # ----------------------------------------------------------
  merged.cell.meta <- rbindlist(lapply(obj.list, function(x) x$cell.meta),
                                use.names = TRUE, fill = TRUE)
  if (nrow(merged.cell.meta) > 0L)
    setkeyv(merged.cell.meta, "cell")

  # ----------------------------------------------------------
  # 7. $edge.list — rbindlist, re-key
  # ----------------------------------------------------------
  merged.edge.list <- rbindlist(lapply(obj.list, function(x) x$edge.list),
                                use.names = TRUE, fill = FALSE)
  if (nrow(merged.edge.list) > 0L)
    setkeyv(merged.edge.list, "edge")

  # ----------------------------------------------------------
  # 8. $edge.data — rbindlist, re-key
  # ----------------------------------------------------------
  merged.edge.data <- rbindlist(lapply(obj.list, function(x) x$edge.data),
                                use.names = TRUE, fill = FALSE)
  if (nrow(merged.edge.data) > 0L)
    setkeyv(merged.edge.data, c("edge", "LRM"))

  # ----------------------------------------------------------
  # 9. $edge.meta — rbindlist with fill, re-key
  # ----------------------------------------------------------
  merged.edge.meta <- rbindlist(lapply(obj.list, function(x) x$edge.meta),
                                use.names = TRUE, fill = TRUE)
  if (nrow(merged.edge.meta) > 0L)
    setkeyv(merged.edge.meta, "edge")

  # ----------------------------------------------------------
  # 10. $gene.list — union; warn if lists differ
  # ----------------------------------------------------------
  gene.lists       <- lapply(obj.list, function(x) x$gene.list)
  merged.gene.list <- Reduce(union, gene.lists)
  if (verbose) {
    all.same <- all(vapply(gene.lists[-1L],
                           function(g) identical(sort(g), sort(gene.lists[[1L]])),
                           logical(1L)))
    if (!all.same)
      message(sprintf(
        paste0("Gene lists differ across objects. Merged $gene.list contains ",
               "%d genes (union). Missing genes are implicitly zero in ",
               "$cell.data (sparse)."),
        length(merged.gene.list)
      ))
  }

  # ----------------------------------------------------------
  # 11. $gene.meta — rbindlist + deduplicate by gene (keep first)
  # ----------------------------------------------------------
  merged.gene.meta <- unique(
    rbindlist(lapply(obj.list, function(x) x$gene.meta),
              use.names = TRUE, fill = TRUE),
    by = "gene"
  )
  if (nrow(merged.gene.meta) > 0L) setkeyv(merged.gene.meta, "gene")

  # ----------------------------------------------------------
  # 12. $LRM.list — union; warn if lists differ
  # ----------------------------------------------------------
  lrm.lists       <- lapply(obj.list, function(x) x$LRM.list)
  merged.LRM.list <- Reduce(union, lrm.lists)
  if (verbose) {
    all.same <- all(vapply(lrm.lists[-1L],
                           function(l) identical(sort(l), sort(lrm.lists[[1L]])),
                           logical(1L)))
    if (!all.same)
      message(sprintf(
        "LRM lists differ across objects. Merged $LRM.list contains %d LRM pairs (union).",
        length(merged.LRM.list)
      ))
  }

  # ----------------------------------------------------------
  # 13. $LRM.meta — rbindlist + deduplicate by LRM (keep first)
  # ----------------------------------------------------------
  merged.LRM.meta <- unique(
    rbindlist(lapply(obj.list, function(x) x$LRM.meta),
              use.names = TRUE, fill = TRUE),
    by = "LRM"
  )
  if (nrow(merged.LRM.meta) > 0L) setkeyv(merged.LRM.meta, "LRM")

  # ----------------------------------------------------------
  # 14. $segmentation — sf::rbind if ALL non-NULL; warn if partial
  # ----------------------------------------------------------
  seg.list   <- lapply(obj.list, function(x) x$segmentation)
  has.seg    <- !vapply(seg.list, is.null, logical(1L))
  merged.seg <- NULL

  if (all(has.seg)) {
    if (!requireNamespace("sf", quietly = TRUE)) {
      warning(
        "All objects have $segmentation but package 'sf' is not installed. ",
        "$segmentation set to NULL in merged object.",
        call. = FALSE
      )
    } else {
      merged.seg <- do.call(rbind, seg.list)
    }
  } else if (any(has.seg)) {
    warning(
      sprintf(
        "Only %d of %d objects have $segmentation populated. $segmentation set to NULL in merged object.",
        sum(has.seg), length(obj.list)
      ),
      call. = FALSE
    )
  }

  # ----------------------------------------------------------
  # 15. sample.col attribute — propagate if consistent
  # ----------------------------------------------------------
  sample.cols <- vapply(obj.list,
                        function(x) attr(x, "sample.col") %||% NA_character_,
                        character(1L))
  merged.sample.col <- NULL
  if (!all(is.na(sample.cols))) {
    non.na <- unique(sample.cols[!is.na(sample.cols)])
    if (length(non.na) == 1L) {
      merged.sample.col <- non.na
    } else if (verbose) {
      warning(
        "Objects have inconsistent 'sample.col' attributes — dropping from merged object.",
        call. = FALSE
      )
    }
  }

  # ----------------------------------------------------------
  # 16. Assemble via .new_NICHESObject()
  # ----------------------------------------------------------
  merged.obj <- .new_NICHESObject(
    cell.list              = merged.cell.list,
    cell.data              = merged.cell.data,
    cell.meta              = merged.cell.meta,
    edge.list              = merged.edge.list,
    edge.data              = merged.edge.data,
    edge.meta              = merged.edge.meta,
    neighborhood.cell.list = NULL,
    neighborhood.edge.list = NULL,
    neighborhood.meta      = NULL,
    gene.list              = merged.gene.list,
    gene.meta              = merged.gene.meta,
    LRM.list               = merged.LRM.list,
    LRM.meta               = merged.LRM.meta,
    segmentation           = merged.seg
  )

  attr(merged.obj, "mode") <- merged.mode
  if (!is.null(merged.sample.col))
    attr(merged.obj, "sample.col") <- merged.sample.col

  if (verbose)
    message(sprintf(
      "Merged %d NICHESObjects: %d cells, %d edges, %d LRM pairs.",
      length(obj.list),
      length(merged.cell.list),
      nrow(merged.edge.list),
      length(merged.LRM.list)
    ))

  merged.obj
}
