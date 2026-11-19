# add_cell_meta.R -------------------------------------------------------------
# Replace $cell.meta on a NICHESObject with a new metadata table and
# (by default) rebuild $edge.meta so that new columns propagate to edges.
#
# FUNCTIONS:
#   add_cell_meta()  -- exported
#
# DEPS: data.table (>= 1.14)


# ==============================================================================
#' Replace \code{$cell.meta} on a NICHESObject
#'
#' @description
#' Fully replaces \code{niches.obj$cell.meta} with \code{new.meta} and
#' (by default) rebuilds \code{$edge.meta} so the new columns propagate as
#' \code{Sending.*} / \code{Receiving.*} / \code{Joint.*} columns on every
#' edge.
#'
#' The replacement is complete: any column present in the existing
#' \code{$cell.meta} but absent from \code{new.meta} is dropped. The user
#' is responsible for supplying all columns they wish to retain, including
#' spatial coordinates (\code{x}, \code{y}).
#'
#' \code{new.meta} must cover every cell in \code{niches.obj$cell.list}.
#' Extra rows (cells not in \code{$cell.list}) are silently dropped with a
#' message when \code{verbose = TRUE}.
#'
#' @param niches.obj       A \code{NICHESObject}.
#' @param new.meta         A \code{data.frame} or \code{data.table} with either
#'   a \code{"cell"} column or rownames matching barcodes in
#'   \code{niches.obj$cell.list}.
#' @param rebuild.edge.meta Logical. If \code{TRUE} (default), regenerates
#'   \code{niches.obj$edge.meta} via \code{.build_edge_meta()} so new metadata
#'   columns appear on each edge. Set to \code{FALSE} when calling repeatedly
#'   in batch -- call once manually at the end.
#' @param verbose          Logical. If \code{TRUE} (default), emits messages
#'   for dropped extra cells and stale aggregation slots.
#'
#' @return Updated \code{NICHESObject}.
#'
#' @seealso \code{\link{add_meta_column}} to append or overwrite a single
#'   column without replacing the full table;
#'   \code{\link{add_neighborhood_composition}} and
#'   \code{\link{aggregate_NICHESObject}} to recompute stale aggregations.
#'
#' @examples
#' data(NICHESv2_example)
#'
#' # Build a complete replacement metadata table
#' new.meta <- data.frame(
#'   cell     = NICHESv2_example$cell.list,
#'   x        = NICHESv2_example$cell.meta$x,
#'   y        = NICHESv2_example$cell.meta$y,
#'   celltype = NICHESv2_example$cell.meta$celltype,
#'   region   = sample(c("core", "edge"),
#'                     length(NICHESv2_example$cell.list),
#'                     replace = TRUE),
#'   stringsAsFactors = FALSE
#' )
#' obj.updated <- add_cell_meta(NICHESv2_example, new.meta)
#'
#' # New column propagates to edges
#' head(obj.updated$edge.meta[, c("SendingCell", "ReceivingCell",
#'                                "Joint.region")])
#'
#' @export
#' @import data.table
add_cell_meta <- function(niches.obj,
                           new.meta,
                           rebuild.edge.meta = TRUE,
                           verbose           = TRUE) {

  if (!inherits(niches.obj, "NICHESObject"))
    stop("[add_cell_meta] 'niches.obj' must be a NICHESObject.")

  if (!is.data.frame(new.meta))
    stop("[add_cell_meta] 'new.meta' must be a data.frame or data.table.")

  ## Section 1: resolve cell column ------------------------------------------
  #
  # Check for rowname-style input BEFORE as.data.table() -- coercion drops
  # rownames. Default data.frame rownames are "1", "2", ... (character); real
  # barcodes will not match that sequence.

  if (!"cell" %in% colnames(new.meta)) {
    rn <- rownames(new.meta)
    has.real.rownames <- !is.null(rn) &&
      !identical(rn, as.character(seq_len(nrow(new.meta))))
    if (has.real.rownames) {
      new.meta        <- as.data.table(new.meta)
      new.meta[, cell := rn]
    } else {
      stop(
        "[add_cell_meta] 'new.meta' must have a 'cell' column or ",
        "rownames matching cell barcodes in $cell.list."
      )
    }
  } else {
    new.meta <- as.data.table(new.meta)
  }

  ## Section 2: cell coverage validation -------------------------------------

  expected     <- niches.obj$cell.list
  provided     <- new.meta$cell
  missing.cells <- setdiff(expected, provided)

  if (length(missing.cells) > 0L)
    stop(sprintf(
      "[add_cell_meta] 'new.meta' is missing %d cell(s) present in $cell.list.\n  First missing (up to 3): %s",
      length(missing.cells),
      paste(head(missing.cells, 3L), collapse = ", ")
    ))

  extra.cells <- setdiff(provided, expected)
  if (length(extra.cells) > 0L && verbose)
    message(sprintf(
      "[add_cell_meta] 'new.meta' contains %d cell(s) not in $cell.list -- these rows will be dropped.",
      length(extra.cells)
    ))

  ## Section 3: full replacement ----------------------------------------------
  #
  # Subset to cells in cell.list (drops extras), align to cell.list order,
  # re-key on "cell".

  idx              <- match(expected, new.meta$cell)
  new.meta         <- new.meta[idx]
  setkeyv(new.meta, "cell")
  niches.obj$cell.meta <- new.meta

  ## Section 4: staleness notification ----------------------------------------

  if (!is.null(niches.obj$aggregations$neighborhood.composition) && verbose)
    message(
      "[add_cell_meta] $aggregations$neighborhood.composition is now stale ",
      "because $cell.meta was replaced. Call add_neighborhood_composition() ",
      "or aggregate_NICHESObject() to update it."
    )

  ## Section 5: rebuild edge.meta ---------------------------------------------

  if (rebuild.edge.meta)
    niches.obj$edge.meta <- .build_edge_meta(niches.obj$edge.list,
                                             niches.obj$cell.meta)
  niches.obj
}
