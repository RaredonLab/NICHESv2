# aggregate_NICHESObject.R ---------------------------------------------------
# Convenience wrapper that runs the full aggregation pipeline on a
# NICHESObject in a single call. Calls the following functions in order:
#
#   1. add_neighborhood_edge_agg()     -- four-category edge aggregation
#   2. add_neighborhood_composition()  -- cell-type composition
#      (skipped if cell.type.col = NULL or mode = "nonspatial")
#
# Both steps always fully recompute and replace their slots. There is no
# incremental update mode -- call the function when you want a fresh result.
#
# DEPS: (none beyond sub-functions)


## AGGREGATE NICHESOBJECT ----

#' Run the Full Aggregation Pipeline on a NICHESObject
#'
#' A convenience wrapper that calls all aggregation functions in the correct
#' dependency order, populating the \code{$aggregations} slot in a single
#' call. Both steps always fully recompute and replace their slots -- there
#' is no incremental update mode. Individual functions can still be called
#' directly for finer control.
#'
#' @section Pipeline steps:
#' \enumerate{
#'   \item \code{\link{add_neighborhood_edge_agg}}: for each focal cell and
#'     each LRM pair, aggregates scores across all neighborhood edges using
#'     a four-way classification (in, out, cross, autocrine). All three
#'     score methods (sum, mean, cmean) and the in/out ratio are always
#'     computed.
#'   \item \code{\link{add_neighborhood_composition}}: computes cell-type
#'     composition of each cell's local neighborhood from \code{$edge.list}
#'     and \code{$cell.meta}. \emph{Skipped if \code{cell.type.col = NULL}
#'     or if the object mode is \code{"nonspatial"}.}
#' }
#'
#' @section Neighborhood slots requirement:
#' Step 1 requires the three neighborhood slots to be populated.
#' \code{\link{create_NICHESObject}} populates them automatically using
#' \code{edge.filter.mode = "both_in_neighborhood"} and
#' \code{remove.autocrine = FALSE}, so this requirement is met in the
#' standard workflow.
#'
#' @param niches.obj    A \code{NICHESObject} with neighborhood slots
#'   populated.
#' @param cell.type.col Character or \code{NULL}. Name of the column in
#'   \code{$cell.meta} holding cell type annotations, passed to
#'   \code{\link{add_neighborhood_composition}}. When \code{NULL} (default),
#'   step 2 is skipped.
#' @param n.cores       Positive integer. Socket-cluster workers passed to
#'   \code{\link{add_neighborhood_edge_agg}}. Default \code{1L}.
#' @param verbose       Logical. Emit progress messages. Default \code{TRUE}.
#'
#' @return The input \code{niches.obj} with \code{$aggregations} populated:
#' \describe{
#'   \item{\code{$neighborhood.edge.agg}}{Always populated (four-category
#'     edge aggregation keyed on \code{(cell, LRM)}, with all three score
#'     methods and in/out ratio).}
#'   \item{\code{$neighborhood.composition}}{Populated when
#'     \code{cell.type.col} is provided and \code{mode = "spatial"}.}
#' }
#'
#' @seealso
#'   \code{\link{add_neighborhood_edge_agg}},
#'   \code{\link{add_neighborhood_composition}}
#'
#' @examples
#' \donttest{
#' data(NICHESv2_example)
#'
#' # Edge aggregation only
#' obj <- aggregate_NICHESObject(NICHESv2_example, n.cores = 1L)
#' names(obj$aggregations)
#'
#' # With cell-type composition
#' obj <- aggregate_NICHESObject(
#'   NICHESv2_example,
#'   cell.type.col = "celltype",
#'   n.cores       = 1L
#' )
#' }
#'
#' @export
aggregate_NICHESObject <- function(
    niches.obj,
    cell.type.col = NULL,
    n.cores       = 1L,
    verbose       = TRUE
) {

  ## Input validation ----------------------------------------------------------
  if (!inherits(niches.obj, "NICHESObject"))
    stop("[aggregate_NICHESObject] 'niches.obj' must be a NICHESObject.",
         call. = FALSE)

  n.cores <- as.integer(n.cores)
  if (is.na(n.cores) || n.cores < 1L)
    stop("[aggregate_NICHESObject] 'n.cores' must be a positive integer.",
         call. = FALSE)

  if (verbose)
    message("[aggregate_NICHESObject] Starting aggregation pipeline ...")

  ## Step 1: Neighborhood edge aggregation ------------------------------------
  if (verbose)
    message("[aggregate_NICHESObject] Step 1/2: neighborhood edge aggregation.")

  niches.obj <- add_neighborhood_edge_agg(
    niches.obj,
    n.cores = n.cores,
    verbose = verbose
  )

  ## Step 2: Neighborhood composition -----------------------------------------
  if (is.null(cell.type.col)) {
    if (verbose)
      message(
        "[aggregate_NICHESObject] Step 2/2: SKIPPED -- 'cell.type.col' not ",
        "provided. Supply cell.type.col to compute neighborhood composition."
      )
  } else if (identical(attr(niches.obj, "mode"), "nonspatial")) {
    if (verbose)
      message(
        "[aggregate_NICHESObject] Step 2/2: SKIPPED -- neighborhood ",
        "composition is not supported for mode = 'nonspatial'."
      )
  } else {
    if (verbose)
      message("[aggregate_NICHESObject] Step 2/2: neighborhood composition.")
    niches.obj <- add_neighborhood_composition(
      niches.obj,
      cell.type.col = cell.type.col,
      include.self  = TRUE,
      normalize     = TRUE,
      verbose       = verbose
    )
  }

  if (verbose)
    message(sprintf(
      "[aggregate_NICHESObject] Pipeline complete. Slots populated: %s.",
      paste(names(niches.obj$aggregations), collapse = ", ")
    ))

  niches.obj
}
