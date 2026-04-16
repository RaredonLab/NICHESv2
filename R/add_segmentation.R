# add_segmentation.R ------------------------------------------------------
# Attach a pre-built sf segmentation object to a NICHESObject.
#
# FUNCTIONS:
#   add_segmentation()  -- main user-facing function
#
# DEPS: sf (Suggests -- checked at runtime)


# ==============================================================================
#' Attach Cell Segmentation Polygons to a NICHESObject
#'
#' @description
#' Stores a pre-built \code{sf} object of cell segmentation polygons in the
#' \code{$segmentation} slot of a \code{NICHESObject}.
#'
#' NICHESv2 does not extract segmentation data from Seurat objects -- that step
#' is platform-specific and belongs in the user's own preprocessing script.
#' Once you have built an \code{sf} object (using \pkg{sf}), call this function
#' to attach it so that downstream plotting utilities can overlay cell
#' boundaries.
#'
#' The \code{sf} object must have a \code{cell} column whose values are cell
#' barcodes, and every cell in \code{niches.obj$cell.list} must have a
#' corresponding polygon -- since both the count matrix and the segmentation
#' are extracted from the same source object, missing polygons indicate a
#' barcode mismatch that must be resolved upstream. An error is raised if any
#' cell is absent. Extra polygons for barcodes not in the object are permitted
#' with a warning.
#'
#' @param niches.obj A \code{NICHESObject}.
#' @param seg.sf     An \code{sf} object with at least a \code{cell} column
#'   containing cell barcodes and an \code{sf} geometry column of cell
#'   polygons. Barcodes must match those in \code{niches.obj$cell.list}.
#'
#' @return The input \code{NICHESObject} with \code{$segmentation} populated.
#'
#' @seealso \code{\link{NICHESv2_example}} for an example object whose
#'   \code{$segmentation} slot is \code{NULL} and ready to receive an
#'   attached \code{sf} object; \code{\link{merge_NICHESObjects}} which
#'   uses \code{sf::rbind()} to combine \code{$segmentation} slots across
#'   merged objects.
#'
#' @examples
#' \dontrun{
#' # NICHESv2 expects 'seg.sf' to be an sf data.frame with:
#' #   * a 'cell' column matching niches.obj$cell.list
#' #   * an sf geometry column of polygons (one per cell)
#' # The exact construction is platform-specific (Visium HD, Xenium, MERFISH,
#' # CosMx, etc.). The skeleton below shows the required shape:
#'
#' library(sf)
#' data(NICHESv2_example)
#'
#' # Replace 'polygon.list' with polygons from your segmentation pipeline:
#' polygon.list <- lapply(NICHESv2_example$cell.list, function(cell) {
#'   # ... your platform-specific polygon extraction ...
#'   st_polygon(list(matrix(c(0, 0, 1, 0, 1, 1, 0, 1, 0, 0),
#'                          ncol = 2, byrow = TRUE)))
#' })
#'
#' seg.sf <- st_sf(
#'   cell     = NICHESv2_example$cell.list,
#'   geometry = st_sfc(polygon.list)
#' )
#'
#' niches.obj <- add_segmentation(NICHESv2_example, seg.sf)
#' }
#'
#' @export
add_segmentation <- function(niches.obj, seg.sf) {

  ## SECTION 1: Input validation ----

  if (!is_NICHESObject(niches.obj))
    stop("[add_segmentation] 'niches.obj' must be a NICHESObject.")

  if (!requireNamespace("sf", quietly = TRUE))
    stop(
      "[add_segmentation] Package 'sf' is required but is not installed.\n",
      "  Install it with: install.packages('sf')"
    )

  if (!inherits(seg.sf, "sf"))
    stop("[add_segmentation] 'seg.sf' must be an sf object.")

  if (!"cell" %in% colnames(seg.sf))
    stop("[add_segmentation] 'seg.sf' must have a 'cell' column of cell barcodes.")

  ## SECTION 2: Coverage check ----
  #
  # Every cell in the object must have a polygon -- segmentation and count
  # matrix always come from the same source, so missing polygons indicate a
  # barcode mismatch that must be fixed before attaching.
  # Extra barcodes in seg.sf (not in the object) are less dangerous -- they are
  # unused -- so these produce a warning rather than an error.

  obj.cells <- niches.obj$cell.list
  seg.cells <- seg.sf$cell

  n.missing <- length(setdiff(obj.cells, seg.cells))
  if (n.missing > 0L)
    stop(sprintf(
      "[add_segmentation] %d cell(s) in the NICHESObject have no polygon in 'seg.sf'.\n  Check that cell barcodes match between your segmentation and count matrix.",
      n.missing
    ))

  n.extra <- length(setdiff(seg.cells, obj.cells))
  if (n.extra > 0L)
    warning(sprintf(
      "[add_segmentation] %d cell(s) in 'seg.sf' are not in the NICHESObject and will be kept but unused.",
      n.extra
    ))

  ## SECTION 3: Attach and return ----

  niches.obj$segmentation <- seg.sf

  message(sprintf(
    "[add_segmentation] Segmentation attached: %d polygon(s) for %d cell(s).",
    nrow(seg.sf), length(obj.cells)
  ))

  niches.obj
}
