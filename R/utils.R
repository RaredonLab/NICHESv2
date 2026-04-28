# utils.R -----------------------------------------------------------------
# Internal helper functions shared across NICHESv2.
# Not exported. Prefixed with "." to signal private scope.


## Complex gene-set parsing ----

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


## Categorical column detection (used for Joint.* columns in edge.meta) ----

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


## Normalisation helpers ----

#' Proportional normalisation within group.
#' Divides each value by the group total. Returns 0 when the group total is 0.
#' @noRd
.prop_norm <- function(x) {
  s <- sum(x, na.rm = TRUE)
  if (s == 0) rep(0, length(x)) else x / s
}



## Null-coalescing operator ----

#' @noRd
`%||%` <- function(x, y) if (is.null(x)) y else x


## Suppress R CMD check NOTEs for NSE column names ----
# These names are used as bare symbols inside data.table j-expressions and
# dplyr-style verbs throughout the package. They are not missing globals;
# declaring them here is the CRAN-recommended pattern.
utils::globalVariables(c(
  # data.table j-expression alias (.(col1, col2) == list(col1, col2))
  ".",
  
  # edge-table columns
  "edge", "from", "to",
  "SendingCell", "ReceivingCell",
  "SendingType", "ReceivingType",
  
  # LRM / score columns
  "LRM", "ligand", "receptor",
  "score", "score.norm", "score.stat",
  
  # cell / spatial columns
  "cell", "x", "y",
  "x.from", "y.from", "x.to", "y.to",
  
  # neighborhood columns
  "neighbor", "radius",
  "mean.neighbors", "sd.neighbors",
  
  # summary stat columns (used as column *names*, not functions)
  "median", "sd",
  
  # plot internals
  ".plot.val",
  
  # add_cell_centric_edge_agg -- cell-centric aggregation columns
  "count.out", "count.in",
  "score.out.sum",   "score.in.sum",
  "score.out.mean",  "score.in.mean",
  "score.out.cmean", "score.in.cmean",
  
  # add_neighborhood_centric_edge_agg -- neighborhood aggregation columns
  "n.edges",
  "count.total",
  "direction",
  "score.total.sum",   "score.total.mean",   "score.total.cmean",
  "count.focal.in",    "count.focal.out",    "count.cross",
  "score.focal.in.sum",    "score.focal.in.mean",    "score.focal.in.cmean",
  "score.focal.out.sum",   "score.focal.out.mean",   "score.focal.out.cmean",
  "score.cross.sum",       "score.cross.mean",        "score.cross.cmean",
  
  # add_neighborhood_composition -- composition columns
  "cell.type", "cell.type.polished",
  
  # add_cell_centric_summaries -- derived summary columns
  # (score.total.sum/mean/cmean and count.total already declared above)
  "score.sum.inout.ratio",      "score.mean.inout.ratio",      "score.cmean.inout.ratio",
  "score.sum.inout.log2ratio",  "score.mean.inout.log2ratio",  "score.cmean.inout.log2ratio",
  "asymmetry",
  
  # add_neighborhood_centric_summaries -- derived summary columns
  "score.focal.sum.inout.ratio",      "score.focal.mean.inout.ratio",      "score.focal.cmean.inout.ratio",
  "score.focal.sum.inout.log2ratio",  "score.focal.mean.inout.log2ratio",  "score.focal.cmean.inout.log2ratio",
  "count.focal.asymmetry", "count.cross.fraction", "count.focal.fraction"
))