#' Add Neighborhood Slots to a NICHESObject
#'
#' Extracts neighborhoods from the edge list of a \code{NICHESObject} and
#' populates the three neighborhood slots: \code{$neighborhood.cell.list},
#' \code{$neighborhood.edge.list}, and \code{$neighborhood.meta}.
#'
#' This function calls \code{\link{extract_neighborhoods}} internally.
#' If the neighborhood slots are already populated, they are overwritten
#' with a message.
#'
#' @param niches.obj A \code{NICHESObject} with a non-empty \code{$edge.list}.
#' @param edge.filter.mode Character. Passed to
#'   \code{\link{extract_neighborhoods}}. Default
#'   \code{"both_in_neighborhood"}.
#' @param remove.autocrine Logical. Passed to
#'   \code{\link{extract_neighborhoods}}. Default \code{FALSE}.
#' @param n.cores Integer. Number of cores for parallel processing. Default 1.
#' @param verbose Logical. Print progress messages. Default \code{TRUE}.
#'
#' @return The input \code{NICHESObject} with the three neighborhood slots
#'   populated. Attribute \code{"edge.filter.mode"} is set on the object.
#'
#' @examples
#' \donttest{
#' data(NICHESv2_example)
#'
#' # Default neighborhood: edges where BOTH endpoints lie in the same
#' # neighborhood as the central cell
#' obj.nb <- add_neighborhoods(NICHESv2_example)
#'
#' length(obj.nb$neighborhood.cell.list)   # one entry per cell
#' head(obj.nb$neighborhood.meta)          # n.neighbors and n.edges per cell
#'
#' # Inspect one cell's neighborhood
#' first.cell <- names(obj.nb$neighborhood.cell.list)[1]
#' obj.nb$neighborhood.cell.list[[first.cell]]
#' obj.nb$neighborhood.edge.list[[first.cell]]
#' }
#'
#' @export
add_neighborhoods <- function(niches.obj,
                              edge.filter.mode = "both_in_neighborhood",
                              remove.autocrine = FALSE,
                              n.cores = 1,
                              verbose = TRUE) {

  if (!inherits(niches.obj, "NICHESObject")) {
    stop("niches.obj must be a NICHESObject")
  }
  if (is.null(niches.obj$edge.list) || nrow(niches.obj$edge.list) == 0) {
    stop("niches.obj$edge.list is empty; cannot extract neighborhoods")
  }

  # Warn if overwriting
  if (!is.null(niches.obj$neighborhood.cell.list)) {
    if (verbose) message("Overwriting existing neighborhood slots")
  }

  # Extract neighborhoods
  result <- extract_neighborhoods(
    edge.list        = niches.obj$edge.list,
    edge.filter.mode = edge.filter.mode,
    remove.autocrine = remove.autocrine,
    n.cores          = n.cores,
    verbose          = verbose
  )

  # Write to slots
  niches.obj$neighborhood.cell.list <- result$neighborhood.cell.list
  niches.obj$neighborhood.edge.list <- result$neighborhood.edge.list
  niches.obj$neighborhood.meta      <- result$neighborhood.meta

  # Store parameter on the object
  attr(niches.obj, "edge.filter.mode") <- attr(result, "edge.filter.mode")

  return(niches.obj)
}
