# add_neighborhood_composition.R ---------------------------------------------
# For each cell, counts how many neighbors of each cell type exist within
# its local edge-defined neighborhood.
#
# "Neighbor" is defined purely from $edge.list: cell B is a neighbor of cell A
# if any directed edge connects them (in either direction). The count for each
# cell type is zero-filled against the full dataset type universe so that the
# resulting table always has the same columns regardless of a cell's actual
# neighborhood.
#
# Spatial mode only: in "nonspatial" mode edges are statistically paired rather
# than spatially proximate, so composition would reflect sampling design rather
# than local tissue structure.
#
# Sources: $edge.list, $cell.meta  -- no neighborhood slot dependency.
#
# DEPS: data.table


## NEIGHBORHOOD COMPOSITION ----

#' Add Neighborhood Cell-Type Composition to a NICHESObject
#'
#' For each cell, counts how many neighbors of each cell type fall within its
#' local neighborhood. Neighborhood membership is derived directly from
#' \code{$edge.list}: cell B is considered a neighbor of cell A if any directed
#' edge connects them (regardless of direction). Counts are zero-filled against
#' the full dataset's cell type universe so that every cell has a value for
#' every type.
#'
#' The result is stored as a single data.table at
#' \code{$aggregations$neighborhood.composition}, replacing any previously
#' computed version unconditionally.
#'
#' @section Column naming:
#' Cell type labels are sanitized for use as column names via
#' \code{gsub("[^A-Za-z0-9.]", ".", type)}, which replaces spaces, slashes,
#' parentheses, and other special characters with dots. The resulting columns
#' follow the pattern \code{<PolishedType>.neighbor.count} (and optionally
#' \code{<PolishedType>.neighbor.prop}). Columns are ordered alphabetically by
#' polished type name. For example, cell type \code{"CD8+ T cell"} becomes
#' \code{"CD8..T.cell.neighbor.count"}.
#'
#' @section Spatial mode only:
#' This function raises a hard error when
#' \code{attr(niches.obj, "mode") == "nonspatial"}. In nonspatial mode, edges are
#' statistically paired rather than spatially proximate, so neighborhood
#' composition would reflect sampling design rather than local tissue
#' structure. This is biologically meaningless and is therefore not supported.
#'
#' @section Double-counting prevention:
#' Because spatial graphs are typically bidirectional (edge A→B and B→A both
#' exist for proximate cells), the neighbor set for each focal cell is formed
#' as the \emph{union} of all cells appearing on the other end of any edge,
#' deduplicated before counting. Each physical neighbor is counted exactly once
#' regardless of how many directed edges connect them.
#'
#' @param niches.obj    A \code{NICHESObject} with \code{mode = "spatial"}.
#' @param cell.type.col Character. Name of the column in \code{$cell.meta}
#'   holding cell type annotations. Hard error if \code{NULL} or not found.
#' @param include.self  Logical. Whether to count the focal cell itself as one
#'   of its own neighbors. Default \code{TRUE}. When \code{TRUE}, cell A
#'   contributes 1 to its own type's count regardless of whether an autocrine
#'   self-loop edge exists.
#' @param normalize     Logical. When \code{TRUE} (default), also produces
#'   \code{<PolishedType>.neighbor.prop} columns: each count divided by the
#'   total number of neighbors (including self when \code{include.self = TRUE}).
#'   Proportions sum to 1 per cell; cells with zero neighbors receive 0.
#' @param verbose       Logical. Emit progress messages. Default \code{TRUE}.
#'
#' @return The input \code{niches.obj} with
#'   \code{$aggregations$neighborhood.composition} populated or replaced.
#'
#' @seealso \code{\link{add_neighborhood_edge_agg}},
#'   \code{\link{add_neighborhoods}},
#'   \code{\link{aggregate_NICHESObject}}
#'
#' @examples
#' \donttest{
#' data(NICHESv2_example)
#'
#' obj <- add_neighborhood_composition(
#'   NICHESv2_example,
#'   cell.type.col = "celltype"
#' )
#' head(obj$aggregations$neighborhood.composition)
#'
#' # Without self-counting and without proportions
#' obj2 <- add_neighborhood_composition(
#'   NICHESv2_example,
#'   cell.type.col = "celltype",
#'   include.self  = FALSE,
#'   normalize     = FALSE
#' )
#' }
#'
#' @export
add_neighborhood_composition <- function(
    niches.obj,
    cell.type.col = NULL,
    include.self  = TRUE,
    normalize     = TRUE,
    verbose       = TRUE
) {

  ## Input validation ----------------------------------------------------------
  if (!inherits(niches.obj, "NICHESObject"))
    stop("[add_neighborhood_composition] 'niches.obj' must be a NICHESObject.",
         call. = FALSE)

  ## Spatial-mode guard --------------------------------------------------------
  if (identical(attr(niches.obj, "mode"), "nonspatial"))
    stop(
      "[add_neighborhood_composition] This function requires mode = 'spatial'. ",
      "In nonspatial mode, edges are statistically paired rather than spatially ",
      "proximate, so neighborhood composition would reflect sampling design ",
      "rather than local tissue structure and is not supported.",
      call. = FALSE
    )

  ## cell.type.col guard -------------------------------------------------------
  if (is.null(cell.type.col))
    stop(
      "[add_neighborhood_composition] 'cell.type.col' must be provided. ",
      "Specify the name of the cell type annotation column in $cell.meta.",
      call. = FALSE
    )
  if (!cell.type.col %in% colnames(niches.obj$cell.meta))
    stop(sprintf(
      "[add_neighborhood_composition] 'cell.type.col' \"%s\" not found in $cell.meta. Available columns: %s",
      cell.type.col,
      paste(setdiff(colnames(niches.obj$cell.meta), "cell"), collapse = ", ")
    ), call. = FALSE)

  if (!is.logical(include.self) || length(include.self) != 1L)
    stop("[add_neighborhood_composition] 'include.self' must be a single logical.",
         call. = FALSE)

  if (!is.logical(normalize) || length(normalize) != 1L)
    stop("[add_neighborhood_composition] 'normalize' must be a single logical.",
         call. = FALSE)

  ## Build bidirectional (cell, neighbor) pairs --------------------------------
  # For each edge A→B, cell B is a neighbor of A and A is a neighbor of B.
  # Both directions are included, then deduplicated so each physical neighbor
  # is counted exactly once per focal cell.
  el    <- niches.obj$edge.list
  pairs <- unique(rbindlist(list(
    data.table(cell = el$SendingCell,   neighbor = el$ReceivingCell),
    data.table(cell = el$ReceivingCell, neighbor = el$SendingCell)
  )))

  ## Handle self-inclusion -----------------------------------------------------
  # When include.self = TRUE, add every cell as its own neighbor.
  # unique() ensures cells with existing autocrine self-loop edges are not
  # double-counted.
  if (include.self) {
    self.pairs <- data.table(
      cell     = niches.obj$cell.list,
      neighbor = niches.obj$cell.list
    )
    pairs <- unique(rbindlist(list(pairs, self.pairs)))
  } else {
    # Remove self-loops that may have arisen from autocrine edges
    pairs <- pairs[cell != neighbor]
  }

  ## Join cell type annotation onto each neighbor -----------------------------
  # match-based lookup avoids the composite-key join gotcha.
  cell.meta <- niches.obj$cell.meta
  nbr.idx   <- match(pairs$neighbor, cell.meta$cell)
  pairs[, cell.type := cell.meta[[cell.type.col]][nbr.idx]]

  ## Build polished type name mapping -----------------------------------------
  # gsub() replaces any character that is not alphanumeric or "." with ".".
  # sort() ensures stable alphabetical column order across calls.
  all.types      <- sort(unique(cell.meta[[cell.type.col]]))
  polished.types <- gsub("[^A-Za-z0-9.]", ".", all.types)
  type.map       <- setNames(polished.types, all.types)
  pairs[, cell.type.polished := type.map[cell.type]]

  if (verbose)
    message(sprintf(
      "[add_neighborhood_composition] Computing composition across %d cells | %d cell types ...",
      length(niches.obj$cell.list), length(all.types)
    ))

  ## Count per (cell, polished type) ------------------------------------------
  count.dt <- pairs[, .N, by = .(cell, cell.type.polished)]

  ## Pivot to wide format ------------------------------------------------------
  comp.dt <- data.table::dcast(
    count.dt,
    cell ~ cell.type.polished,
    value.var = "N",
    fill      = 0L
  )

  ## Rename columns to .neighbor.count convention -----------------------------
  type.cols.present <- setdiff(colnames(comp.dt), "cell")
  setnames(comp.dt, type.cols.present,
           paste0(type.cols.present, ".neighbor.count"))

  ## Zero-fill types absent from dcast output ---------------------------------
  # dcast omits levels that never appear across any cell; ensure every type
  # has a column present.
  expected.count.cols <- paste0(polished.types, ".neighbor.count")
  for (col in setdiff(expected.count.cols, colnames(comp.dt)))
    comp.dt[, (col) := 0L]

  ## Proportion columns (optional) --------------------------------------------
  if (normalize) {
    count.cols <- paste0(polished.types, ".neighbor.count")
    prop.cols  <- paste0(polished.types, ".neighbor.prop")
    row.totals <- rowSums(comp.dt[, count.cols, with = FALSE])
    for (i in seq_along(polished.types)) {
      prop.val               <- comp.dt[[count.cols[i]]] / row.totals
      prop.val[is.nan(prop.val)] <- 0  # guard: cells with zero neighbors
      comp.dt[, (prop.cols[i]) := prop.val]
    }
  }

  ## Ensure all cells are represented -----------------------------------------
  # Cells with no edges at all will be absent from comp.dt after dcast.
  # Left-join against the full cell list and zero-fill NAs.
  all.cells.dt <- data.table(cell = niches.obj$cell.list)
  comp.dt      <- merge(all.cells.dt, comp.dt, by = "cell", all.x = TRUE)

  count.na.cols <- grep("\\.neighbor\\.count$", colnames(comp.dt), value = TRUE)
  prop.na.cols  <- grep("\\.neighbor\\.prop$",  colnames(comp.dt), value = TRUE)
  for (col in count.na.cols) {
    na.rows <- is.na(comp.dt[[col]])
    if (any(na.rows)) comp.dt[na.rows, (col) := 0L]
  }
  for (col in prop.na.cols) {
    na.rows <- is.na(comp.dt[[col]])
    if (any(na.rows)) comp.dt[na.rows, (col) := 0]
  }

  ## Enforce column order: cell | counts (alpha) | props (alpha) -------------
  count.col.order <- paste0(polished.types, ".neighbor.count")
  prop.col.order  <- paste0(polished.types, ".neighbor.prop")
  col.order <- c(
    "cell",
    intersect(count.col.order, colnames(comp.dt)),
    intersect(prop.col.order,  colnames(comp.dt))
  )
  setcolorder(comp.dt, col.order)
  setkey(comp.dt, cell)

  ## Store result (always overwrites) -----------------------------------------
  niches.obj$aggregations$neighborhood.composition <- comp.dt

  if (verbose)
    message(sprintf(
      "[add_neighborhood_composition] Done. %d cells | %d cell types | columns: %s.",
      nrow(comp.dt),
      length(all.types),
      if (normalize) "count + prop" else "count"
    ))

  niches.obj
}
