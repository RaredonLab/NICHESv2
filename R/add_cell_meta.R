#' Add or Update Columns in \code{$cell.meta}
#'
#' Joins new metadata onto an existing NICHESObject's \code{$cell.meta} table
#' and (by default) rebuilds \code{$edge.meta} so the new columns are
#' propagated as \code{Sending.*} / \code{Receiving.*} / \code{Joint.*}
#' columns on every edge.
#'
#' @param niches.obj A \code{NICHESObject}.
#' @param new.meta   A \code{data.frame} or \code{data.table} with a
#'   \code{"cell"} column matching barcodes in \code{niches.obj$cell.list}.
#' @param rebuild.edge.meta Logical. If \code{TRUE} (default), regenerates
#'   \code{niches.obj$edge.meta} via the internal \code{.build_edge_meta()}
#'   helper so the new metadata columns appear on each edge. Set to
#'   \code{FALSE} to skip the rebuild (useful when adding many columns in
#'   batch -- call once manually at the end).
#' @return Updated \code{NICHESObject}.
#'
#' @examples
#' \donttest{
#' data(NICHESv2_example)
#'
#' # Tag each cell with a synthetic region label and propagate to edges
#' new.meta <- data.frame(
#'   cell   = NICHESv2_example$cell.list,
#'   region = sample(c("core", "edge"),
#'                   length(NICHESv2_example$cell.list),
#'                   replace = TRUE)
#' )
#' obj.tagged <- add_cell_meta(NICHESv2_example, new.meta)
#'
#' # The categorical column appears as Sending.*/Receiving.*/Joint.* on edges
#' head(obj.tagged$edge.meta[, c("SendingCell", "ReceivingCell",
#'                               "Sending.region", "Receiving.region",
#'                               "Joint.region")])
#' }
#'
#' @export
#' @import data.table
add_cell_meta <- function(niches.obj, new.meta, rebuild.edge.meta = TRUE) {
  if (!is_NICHESObject(niches.obj))
    stop("[add_cell_meta] 'niches.obj' must be a NICHESObject.")

  new.meta <- as.data.table(new.meta)
  if (!"cell" %in% colnames(new.meta))
    stop("[add_cell_meta] 'new.meta' must have a 'cell' column.")

  niches.obj$cell.meta <- new.meta[niches.obj$cell.meta, on = "cell"]

  if (rebuild.edge.meta)
    niches.obj$edge.meta <- .build_edge_meta(niches.obj$edge.list,
                                             niches.obj$cell.meta)
  niches.obj
}
