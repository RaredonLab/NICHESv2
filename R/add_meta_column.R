# add_meta_column.R -----------------------------------------------------------
# Append or overwrite a single metadata column in $cell.meta on a
# NICHESObject, with optional $edge.meta rebuild.
#
# FUNCTIONS:
#   add_meta_column()  -- exported
#
# DEPS: data.table (>= 1.14)


# ==============================================================================
#' Append or Overwrite a Single Column in \code{$cell.meta}
#'
#' @description
#' Adds one new column to \code{niches.obj$cell.meta} from a named vector,
#' or overwrites an existing column of the same name. All other columns in
#' \code{$cell.meta} are preserved unchanged.
#'
#' Use \code{\link{add_cell_meta}} instead when you want to replace the entire
#' metadata table.
#'
#' \code{col.values} must be a named vector whose names are cell barcodes. It
#' must supply a value for every cell in \code{niches.obj$cell.list}. Extra
#' names (cells not in \code{$cell.list}) are silently ignored.
#'
#' If \code{col.name} already exists in \code{$cell.meta}, a message is
#' emitted and the column is overwritten.
#'
#' @param niches.obj        A \code{NICHESObject}.
#' @param col.name          Character scalar. Name of the column to add or
#'   overwrite.
#' @param col.values        Named vector. Names must be cell barcodes covering
#'   every cell in \code{niches.obj$cell.list}. The vector type determines the
#'   column type (character, numeric, integer, logical, factor).
#' @param rebuild.edge.meta Logical. If \code{TRUE} (default), regenerates
#'   \code{niches.obj$edge.meta} via \code{.build_edge_meta()} so the new
#'   column appears as \code{Sending.<col.name>} / \code{Receiving.<col.name>}
#'   (and \code{Joint.<col.name>} for categorical types) on every edge.
#' @param verbose           Logical. If \code{TRUE} (default), emits a message
#'   when overwriting an existing column and when aggregation slots may be
#'   stale.
#'
#' @return Updated \code{NICHESObject}.
#'
#' @seealso \code{\link{add_cell_meta}} to replace the full metadata table;
#'   \code{\link{add_neighborhood_composition}} and
#'   \code{\link{aggregate_NICHESObject}} to recompute stale aggregations.
#'
#' @examples
#' \donttest{
#' data(NICHESv2_example)
#'
#' # Add a new numeric column
#' scores <- setNames(
#'   runif(length(NICHESv2_example$cell.list)),
#'   NICHESv2_example$cell.list
#' )
#' obj.updated <- add_meta_column(NICHESv2_example,
#'                                col.name   = "quality.score",
#'                                col.values = scores)
#' head(obj.updated$cell.meta[, c("cell", "quality.score")])
#'
#' # Overwrite an existing categorical column
#' new.types <- setNames(
#'   sample(c("TypeA", "TypeB"), length(NICHESv2_example$cell.list),
#'          replace = TRUE),
#'   NICHESv2_example$cell.list
#' )
#' obj.updated <- add_meta_column(NICHESv2_example,
#'                                col.name   = "celltype",
#'                                col.values = new.types)
#' }
#'
#' @export
#' @import data.table
add_meta_column <- function(niches.obj,
                             col.name,
                             col.values,
                             rebuild.edge.meta = TRUE,
                             verbose           = TRUE) {

  if (!inherits(niches.obj, "NICHESObject"))
    stop("[add_meta_column] 'niches.obj' must be a NICHESObject.")

  if (!is.character(col.name) || length(col.name) != 1L || is.na(col.name))
    stop("[add_meta_column] 'col.name' must be a single non-NA character string.")

  if (!is.atomic(col.values))
    stop(
      "[add_meta_column] 'col.values' must be a named atomic vector with cell ",
      "barcodes as names."
    )
  if (is.null(names(col.values)))
    stop(
      "[add_meta_column] 'col.values' must be a named atomic vector with cell ",
      "barcodes as names."
    )

  ## Section 1: cell coverage validation -------------------------------------

  expected      <- niches.obj$cell.list
  provided      <- names(col.values)
  missing.cells <- setdiff(expected, provided)

  if (length(missing.cells) > 0L)
    stop(sprintf(
      "[add_meta_column] 'col.values' is missing %d cell(s) present in $cell.list.\n  First missing (up to 3): %s",
      length(missing.cells),
      paste(head(missing.cells, 3L), collapse = ", ")
    ))

  ## Section 2: overwrite notification ----------------------------------------

  if (col.name %in% colnames(niches.obj$cell.meta) && verbose)
    message(sprintf(
      "[add_meta_column] Column \"%s\" already exists in $cell.meta and will be overwritten.",
      col.name
    ))

  ## Section 3: assign column -------------------------------------------------
  #
  # Index col.values by the cell column of $cell.meta (which is already aligned
  # to $cell.list order and keyed). Extra names in col.values are ignored by
  # named-vector indexing.

  niches.obj$cell.meta <- copy(niches.obj$cell.meta)
  niches.obj$cell.meta[, (col.name) := col.values[niches.obj$cell.meta$cell]]

  ## Section 4: staleness notification ----------------------------------------

  if (!is.null(niches.obj$aggregations$neighborhood.composition) && verbose)
    message(sprintf(
      "[add_meta_column] $aggregations$neighborhood.composition may be stale if \"%s\" was the 'cell.type.col' used to compute it. Call add_neighborhood_composition() or aggregate_NICHESObject() to update it.",
      col.name
    ))

  ## Section 5: rebuild edge.meta ---------------------------------------------

  if (rebuild.edge.meta)
    niches.obj$edge.meta <- .build_edge_meta(niches.obj$edge.list,
                                             niches.obj$cell.meta)
  niches.obj
}
