# aggregate_NICHESObject.R ---------------------------------------------------
# Convenience wrapper that runs the full aggregation and summary pipeline on
# a NICHESObject in a single call. Calls the following functions in order:
#
#   1. add_cell_centric_edge_agg()         -- all three methods
#   2. add_cell_centric_summaries()        -- all summaries
#   3. add_neighborhood_centric_edge_agg() -- all methods; optional by.direction
#      (skipped with a warning if neighborhood slots are NULL)
#   4. add_neighborhood_centric_summaries()
#      (skipped if step 3 skipped or by.direction = FALSE)
#   5. add_neighborhood_composition()
#      (skipped if cell.type.col = NULL or mode = "sampled")
#
# Each sub-function is called with overwrite = overwrite and n.cores = n.cores.
# Verbose output from sub-functions provides per-step progress.
#
# DEPS: (none beyond sub-functions)


## AGGREGATE NICHESOBJECT ----

#' Run the Full Aggregation and Summary Pipeline on a NICHESObject
#'
#' A convenience wrapper that calls all aggregation and summary functions in
#' the correct dependency order, populating the entire
#' \code{$aggregations} slot in a single call. Individual functions can still
#' be called directly for finer control (e.g.\ to add a single summary type
#' or to operate on a specific neighborhood key).
#'
#' @section Pipeline steps:
#' \enumerate{
#'   \item \code{\link{add_cell_centric_edge_agg}}: aggregates scores across
#'     all edges directly linking each cell, using all three methods
#'     (\code{"sum"}, \code{"mean"}, \code{"cmean"}).
#'   \item \code{\link{add_cell_centric_summaries}}: appends all eleven derived
#'     summary columns to the cell-centric aggregation table.
#'   \item \code{\link{add_neighborhood_centric_edge_agg}}: aggregates scores
#'     across all neighborhood edges for each cell, using all three methods.
#'     When \code{by.direction = TRUE} (default), also computes the focal-in,
#'     focal-out, and cross directional decomposition. \emph{Skipped with a
#'     warning if the neighborhood slots are \code{NULL}.}
#'   \item \code{\link{add_neighborhood_centric_summaries}}: appends all nine
#'     derived summary columns to the neighborhood-centric aggregation entry.
#'     \emph{Skipped if step 3 was skipped or if \code{by.direction = FALSE}.}
#'   \item \code{\link{add_neighborhood_composition}}: computes cell-type
#'     composition of each cell's local neighborhood. \emph{Skipped if
#'     \code{cell.type.col = NULL} or if the object mode is
#'     \code{"sampled"}.}
#' }
#'
#' @section Neighborhood slots requirement:
#' Steps 3 and 4 require the three neighborhood slots to be populated.
#' \code{\link{create_NICHESObject}} populates them automatically (default
#' \code{edge.filter.mode = "both_in_neighborhood"}), so this requirement is
#' met in the standard workflow. If neighborhood slots are \code{NULL} (e.g.\
#' after subsetting the object), steps 3 and 4 are skipped with an
#' informational message rather than a hard error.
#'
#' @section Incremental use:
#' When \code{overwrite = FALSE} (default), each sub-function skips columns
#' that are already present and only appends genuinely new results. This means
#' \code{aggregate_NICHESObject} can be re-called safely after manually running
#' individual aggregation steps -- it will fill in any missing pieces without
#' recomputing existing ones.
#'
#' @param niches.obj    A \code{NICHESObject}.
#' @param cell.type.col Character or \code{NULL}. Name of the column in
#'   \code{$cell.meta} holding cell type annotations, passed to
#'   \code{\link{add_neighborhood_composition}}. When \code{NULL} (default),
#'   step 5 is skipped.
#' @param by.direction  Logical. When \code{TRUE} (default), computes
#'   directional decomposition (focal-in, focal-out, cross) in
#'   \code{\link{add_neighborhood_centric_edge_agg}} and runs
#'   \code{\link{add_neighborhood_centric_summaries}}. When \code{FALSE},
#'   only total-level neighborhood aggregation is computed and summaries are
#'   skipped.
#' @param include.self  Logical. Passed to
#'   \code{\link{add_neighborhood_composition}}. Whether to count the focal
#'   cell itself as one of its own neighbors. Default \code{TRUE}.
#' @param normalize     Logical. Passed to
#'   \code{\link{add_neighborhood_composition}}. Whether to also compute
#'   \code{.neighbor.prop} proportion columns. Default \code{TRUE}.
#' @param overwrite     Logical. When \code{FALSE} (default), each sub-function
#'   skips columns already present. When \code{TRUE}, all steps are recomputed
#'   from scratch.
#' @param n.cores       Positive integer. Number of parallel worker processes
#'   passed to the aggregation functions. Default \code{1L}.
#' @param verbose       Logical. Emit progress messages from each sub-function.
#'   Default \code{TRUE}.
#'
#' @return The input \code{niches.obj} with \code{$aggregations} populated:
#' \describe{
#'   \item{\code{$cell.centric.edge.agg}}{Always populated.}
#'   \item{\code{$neighborhood.centric.edge.agg}}{Populated when neighborhood
#'     slots are non-\code{NULL}.}
#'   \item{\code{$neighborhood.composition}}{Populated when
#'     \code{cell.type.col} is provided and \code{mode = "spatial"}.}
#' }
#'
#' @seealso
#'   \code{\link{add_cell_centric_edge_agg}},
#'   \code{\link{add_cell_centric_summaries}},
#'   \code{\link{add_neighborhood_centric_edge_agg}},
#'   \code{\link{add_neighborhood_centric_summaries}},
#'   \code{\link{add_neighborhood_composition}}
#'
#' @examples
#' \donttest{
#' data(NICHESv2_example)
#'
#' # Full pipeline with cell type composition
#' obj <- aggregate_NICHESObject(
#'   NICHESv2_example,
#'   cell.type.col = "celltype",
#'   n.cores       = 1L
#' )
#' names(obj$aggregations)
#'
#' # Without directional decomposition (faster; skips summaries step)
#' obj2 <- aggregate_NICHESObject(
#'   NICHESv2_example,
#'   by.direction = FALSE,
#'   n.cores      = 1L
#' )
#'
#' # Incremental: add composition to an already-aggregated object
#' obj3 <- aggregate_NICHESObject(
#'   obj,
#'   cell.type.col = "celltype",
#'   overwrite     = FALSE
#' )
#' }
#'
#' @export
aggregate_NICHESObject <- function(
    niches.obj,
    cell.type.col = NULL,
    by.direction  = TRUE,
    include.self  = TRUE,
    normalize     = TRUE,
    overwrite     = FALSE,
    n.cores       = 1L,
    verbose       = TRUE
) {

  ## Input validation ----------------------------------------------------------
  if (!inherits(niches.obj, "NICHESObject"))
    stop("[aggregate_NICHESObject] 'niches.obj' must be a NICHESObject.",
         call. = FALSE)

  if (!is.logical(by.direction) || length(by.direction) != 1L)
    stop("[aggregate_NICHESObject] 'by.direction' must be a single logical.",
         call. = FALSE)

  if (!is.logical(overwrite) || length(overwrite) != 1L)
    stop("[aggregate_NICHESObject] 'overwrite' must be a single logical.",
         call. = FALSE)

  n.cores <- as.integer(n.cores)
  if (is.na(n.cores) || n.cores < 1L)
    stop("[aggregate_NICHESObject] 'n.cores' must be a positive integer.",
         call. = FALSE)

  if (verbose)
    message("[aggregate_NICHESObject] Starting aggregation pipeline ...")

  ## Step 1: Cell-centric edge aggregation ------------------------------------
  if (verbose) message("[aggregate_NICHESObject] Step 1/5: cell-centric edge aggregation.")
  niches.obj <- add_cell_centric_edge_agg(
    niches.obj,
    methods   = c("sum", "mean", "cmean"),
    overwrite = overwrite,
    n.cores   = n.cores,
    verbose   = verbose
  )

  ## Step 2: Cell-centric summaries -------------------------------------------
  if (verbose) message("[aggregate_NICHESObject] Step 2/5: cell-centric summaries.")
  niches.obj <- add_cell_centric_summaries(
    niches.obj,
    overwrite = overwrite,
    verbose   = verbose
  )

  ## Steps 3-4: Neighborhood-centric aggregation and summaries ----------------
  has.neighborhoods <- !is.null(niches.obj$neighborhood.cell.list)

  if (!has.neighborhoods) {
    if (verbose)
      message(
        "[aggregate_NICHESObject] Steps 3-4: SKIPPED -- neighborhood slots are ",
        "NULL. Run add_neighborhoods() first if neighborhood-centric ",
        "aggregation is needed."
      )
  } else {

    ## Step 3: Neighborhood-centric edge aggregation --------------------------
    if (verbose)
      message(sprintf(
        "[aggregate_NICHESObject] Step 3/5: neighborhood-centric edge aggregation ",
        "(by.direction = %s).", by.direction
      ))
    niches.obj <- add_neighborhood_centric_edge_agg(
      niches.obj,
      methods      = c("sum", "mean", "cmean"),
      by.direction = by.direction,
      overwrite    = overwrite,
      n.cores      = n.cores,
      verbose      = verbose
    )

    ## Step 4: Neighborhood-centric summaries ---------------------------------
    if (!by.direction) {
      if (verbose)
        message(
          "[aggregate_NICHESObject] Step 4/5: SKIPPED -- neighborhood-centric ",
          "summaries require by.direction = TRUE."
        )
    } else {
      if (verbose)
        message("[aggregate_NICHESObject] Step 4/5: neighborhood-centric summaries.")
      niches.obj <- add_neighborhood_centric_summaries(
        niches.obj,
        overwrite = overwrite,
        verbose   = verbose
      )
    }
  }

  ## Step 5: Neighborhood composition -----------------------------------------
  if (is.null(cell.type.col)) {
    if (verbose)
      message(
        "[aggregate_NICHESObject] Step 5/5: SKIPPED -- 'cell.type.col' not ",
        "provided. Supply cell.type.col to compute neighborhood composition."
      )
  } else if (identical(attr(niches.obj, "mode"), "sampled")) {
    if (verbose)
      message(
        "[aggregate_NICHESObject] Step 5/5: SKIPPED -- neighborhood composition ",
        "is not supported for mode = 'sampled'."
      )
  } else {
    if (verbose)
      message("[aggregate_NICHESObject] Step 5/5: neighborhood composition.")
    niches.obj <- add_neighborhood_composition(
      niches.obj,
      cell.type.col = cell.type.col,
      include.self  = include.self,
      normalize     = normalize,
      verbose       = verbose
    )
  }

  if (verbose)
    message(sprintf(
      "[aggregate_NICHESObject] Pipeline complete. Aggregation slots populated: %s.",
      paste(names(niches.obj$aggregations), collapse = ", ")
    ))

  niches.obj
}
