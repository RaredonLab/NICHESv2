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
