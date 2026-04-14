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
# Null-coalescing operator
# ------------------------------------------------------------------------------

#' @noRd
`%||%` <- function(x, y) if (is.null(x)) y else x


