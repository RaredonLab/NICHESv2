#' Merge multiple NICHESObjects into one
#'
#' Combines a list of NICHESObjects — or any number passed as individual
#' arguments — into a single NICHESObject by concatenating all per-cell and
#' per-edge data. This is the complement to the \code{sample.col} splitting
#' that \code{\link{create_NICHESObject}} performs: it can be used to merge
#' objects that were produced independently and then need to be analysed
#' together.
#'
#' @param ... Either (a) multiple NICHESObjects passed directly as positional
#'   arguments, or (b) a single named or unnamed list of NICHESObjects, or
#'   (c) an explicit \code{list(a, b, c)} call.
#' @param verbose Logical. Whether to emit progress messages and warnings.
#'   Default \code{TRUE}.
#'
#' @details
#' \strong{Mode compatibility.} All objects must share the same
#' \code{attr(., "mode")} (\code{"spatial"} or \code{"sampled"}). A mismatch
#' is a hard error because the slot schemas differ (spatial objects carry
#' \code{x}/\code{y} in \code{$cell.meta}; edge semantics differ).
#'
#' \strong{Barcode uniqueness.} Cell barcodes must be globally unique across
#' all objects. If the same barcode string appears in two or more objects the
#' function stops with a clear error. Prepend a sample-specific prefix before
#' constructing the objects (e.g. \code{paste0(sample.id, "_", barcode)}) to
#' avoid collisions — Seurat does this automatically via \code{RenameCells()}.
#'
#' \strong{Gene / LRM union.} When gene or LRM lists differ across objects,
#' their union is taken. Genes absent in one object are implicitly zero in
#' \code{$cell.data} (which is sparse — zero rows are never stored).
#'
#' \strong{Slot-by-slot strategy:}
#' \itemize{
#'   \item \code{$cell.list}, \code{$cell.data}, \code{$edge.list},
#'         \code{$edge.data} — simple concatenation.
#'   \item \code{$cell.meta}, \code{$edge.meta} — concatenated with
#'         \code{fill = TRUE}; extra columns in some objects become \code{NA}
#'         in others.
#'   \item \code{$gene.list}, \code{$LRM.list} — union.
#'   \item \code{$gene.meta}, \code{$LRM.meta} — concatenated, deduplicated
#'         by key (first occurrence kept).
#'   \item \code{$segmentation} — merged via \code{sf::rbind()} when all
#'         objects have it; \code{NULL} with a warning when only some do.
#'   \item Neighborhood slots — always \code{NULL} (reserved).
#' }
#'
#' @return A single \code{NICHESObject}.
#'
#' @seealso \code{\link{create_NICHESObject}} for the \code{sample.col}
#'   argument that automates splitting and merging within one call.
#'
#' @examples
#' \dontrun{
#' # Three individual objects
#' merged <- merge_NICHESObjects(obj.s1, obj.s2, obj.s3)
#'
#' # A pre-built list
#' obj.list <- list(obj.s1, obj.s2, obj.s3)
#' merged <- merge_NICHESObjects(obj.list)
#'
#' # Explicit list() call
#' merged <- merge_NICHESObjects(list(obj.s1, obj.s2, obj.s3))
#' }
#'
#' @export
merge_NICHESObjects <- function(..., verbose = TRUE) {

  args <- list(...)

  # Disambiguate: single list of NICHESObjects vs. multiple direct arguments
  if (length(args) == 1L &&
      is.list(args[[1L]]) &&
      !is_NICHESObject(args[[1L]])) {
    obj.list <- args[[1L]]
  } else {
    obj.list <- args
  }

  if (!is.list(obj.list)) {
    stop("Could not interpret arguments as a list of NICHESObjects.",
         call. = FALSE)
  }

  .merge_NICHESObject_list(obj.list, verbose = verbose)
}
