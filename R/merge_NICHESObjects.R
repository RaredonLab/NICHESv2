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
#' @import data.table
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

  # ----------------------------------------------------------
  # 0. Length guard
  # ----------------------------------------------------------
  if (length(obj.list) == 0L)
    stop("obj.list is empty — nothing to merge.", call. = FALSE)

  if (length(obj.list) == 1L) {
    if (verbose) message("Only one object supplied — returning it unchanged.")
    return(obj.list[[1L]])
  }

  # ----------------------------------------------------------
  # 1. All inputs must be NICHESObjects
  # ----------------------------------------------------------
  not.niches <- !vapply(obj.list, is_NICHESObject, logical(1L))
  if (any(not.niches)) {
    bad.idx <- which(not.niches)
    stop(sprintf(
      "Element(s) %s of obj.list are not NICHESObjects.",
      paste(bad.idx, collapse = ", ")
    ), call. = FALSE)
  }

  # ----------------------------------------------------------
  # 2. Mode compatibility — hard error on mismatch
  # ----------------------------------------------------------
  modes <- vapply(obj.list,
                  function(x) attr(x, "mode") %||% NA_character_,
                  character(1L))
  if (length(unique(modes)) > 1L)
    stop(sprintf(
      paste0("Cannot merge NICHESObjects with different modes: %s.\n",
             "All objects must share the same mode (\"spatial\" or \"sampled\")."),
      paste(unique(modes), collapse = ", ")
    ), call. = FALSE)

  # unname() prevents a named vector propagating into attr(obj, "mode")
  merged.mode <- unname(modes[1L])

  # ----------------------------------------------------------
  # 3. Barcode collision check — hard error
  # ----------------------------------------------------------
  all.barcodes <- unlist(lapply(obj.list, function(x) x$cell.list),
                         use.names = FALSE)
  dup.mask <- duplicated(all.barcodes)
  if (any(dup.mask)) {
    dup.barcodes <- unique(all.barcodes[dup.mask])
    n.shown      <- min(length(dup.barcodes), 10L)
    stop(sprintf(
      paste0(
        "Duplicate cell barcodes detected across objects (%d total).\n",
        "First %d: %s\n",
        "Tip: prepend a sample-specific prefix to barcodes before calling\n",
        "create_NICHESObject() (e.g. paste0(sample.id, \"_\", barcode))."
      ),
      length(dup.barcodes), n.shown,
      paste(head(dup.barcodes, n.shown), collapse = ", ")
    ), call. = FALSE)
  }

  # ----------------------------------------------------------
  # 4. $cell.list — concatenate (already validated unique)
  # ----------------------------------------------------------
  merged.cell.list <- all.barcodes

  # ----------------------------------------------------------
  # 5. $cell.data — rbindlist, re-key
  # ----------------------------------------------------------
  merged.cell.data <- rbindlist(lapply(obj.list, function(x) x$cell.data),
                                use.names = TRUE, fill = FALSE)
  if (nrow(merged.cell.data) > 0L)
    setkeyv(merged.cell.data, c("cell", "gene"))

  # ----------------------------------------------------------
  # 6. $cell.meta — rbindlist with fill, re-key
  # ----------------------------------------------------------
  merged.cell.meta <- rbindlist(lapply(obj.list, function(x) x$cell.meta),
                                use.names = TRUE, fill = TRUE)
  if (nrow(merged.cell.meta) > 0L)
    setkeyv(merged.cell.meta, "cell")

  # ----------------------------------------------------------
  # 7. $edge.list — rbindlist, re-key
  # ----------------------------------------------------------
  merged.edge.list <- rbindlist(lapply(obj.list, function(x) x$edge.list),
                                use.names = TRUE, fill = FALSE)
  if (nrow(merged.edge.list) > 0L)
    setkeyv(merged.edge.list, "edge")

  # ----------------------------------------------------------
  # 8. $edge.data — rbindlist, re-key
  # ----------------------------------------------------------
  merged.edge.data <- rbindlist(lapply(obj.list, function(x) x$edge.data),
                                use.names = TRUE, fill = FALSE)
  if (nrow(merged.edge.data) > 0L)
    setkeyv(merged.edge.data, c("edge", "LRM"))

  # ----------------------------------------------------------
  # 9. $edge.meta — rbindlist with fill, re-key
  # ----------------------------------------------------------
  merged.edge.meta <- rbindlist(lapply(obj.list, function(x) x$edge.meta),
                                use.names = TRUE, fill = TRUE)
  if (nrow(merged.edge.meta) > 0L)
    setkeyv(merged.edge.meta, "edge")

  # ----------------------------------------------------------
  # 10. $gene.list — union; warn if lists differ
  # ----------------------------------------------------------
  gene.lists       <- lapply(obj.list, function(x) x$gene.list)
  merged.gene.list <- Reduce(union, gene.lists)
  if (verbose) {
    all.same <- all(vapply(gene.lists[-1L],
                           function(g) identical(sort(g), sort(gene.lists[[1L]])),
                           logical(1L)))
    if (!all.same)
      message(sprintf(
        paste0("Gene lists differ across objects. Merged $gene.list contains ",
               "%d genes (union). Missing genes are implicitly zero in ",
               "$cell.data (sparse)."),
        length(merged.gene.list)
      ))
  }

  # ----------------------------------------------------------
  # 11. $gene.meta — rbindlist + deduplicate by gene (keep first)
  # ----------------------------------------------------------
  merged.gene.meta <- unique(
    rbindlist(lapply(obj.list, function(x) x$gene.meta),
              use.names = TRUE, fill = TRUE),
    by = "gene"
  )
  if (nrow(merged.gene.meta) > 0L) setkeyv(merged.gene.meta, "gene")

  # ----------------------------------------------------------
  # 12. $LRM.list — union; warn if lists differ
  # ----------------------------------------------------------
  lrm.lists       <- lapply(obj.list, function(x) x$LRM.list)
  merged.LRM.list <- Reduce(union, lrm.lists)
  if (verbose) {
    all.same <- all(vapply(lrm.lists[-1L],
                           function(l) identical(sort(l), sort(lrm.lists[[1L]])),
                           logical(1L)))
    if (!all.same)
      message(sprintf(
        "LRM lists differ across objects. Merged $LRM.list contains %d LRM pairs (union).",
        length(merged.LRM.list)
      ))
  }

  # ----------------------------------------------------------
  # 13. $LRM.meta — rbindlist + deduplicate by LRM (keep first)
  # ----------------------------------------------------------
  merged.LRM.meta <- unique(
    rbindlist(lapply(obj.list, function(x) x$LRM.meta),
              use.names = TRUE, fill = TRUE),
    by = "LRM"
  )
  if (nrow(merged.LRM.meta) > 0L) setkeyv(merged.LRM.meta, "LRM")

  # ----------------------------------------------------------
  # 14. $segmentation — sf::rbind if ALL non-NULL; warn if partial
  # ----------------------------------------------------------
  seg.list   <- lapply(obj.list, function(x) x$segmentation)
  has.seg    <- !vapply(seg.list, is.null, logical(1L))
  merged.seg <- NULL

  if (all(has.seg)) {
    if (!requireNamespace("sf", quietly = TRUE)) {
      warning(
        "All objects have $segmentation but package 'sf' is not installed. ",
        "$segmentation set to NULL in merged object.",
        call. = FALSE
      )
    } else {
      merged.seg <- do.call(rbind, seg.list)
    }
  } else if (any(has.seg)) {
    warning(
      sprintf(
        "Only %d of %d objects have $segmentation populated. $segmentation set to NULL in merged object.",
        sum(has.seg), length(obj.list)
      ),
      call. = FALSE
    )
  }

  # ----------------------------------------------------------
  # 15. sample.col attribute — propagate if consistent
  # ----------------------------------------------------------
  sample.cols <- vapply(obj.list,
                        function(x) attr(x, "sample.col") %||% NA_character_,
                        character(1L))
  merged.sample.col <- NULL
  if (!all(is.na(sample.cols))) {
    non.na <- unique(sample.cols[!is.na(sample.cols)])
    if (length(non.na) == 1L) {
      merged.sample.col <- non.na
    } else if (verbose) {
      warning(
        "Objects have inconsistent 'sample.col' attributes — dropping from merged object.",
        call. = FALSE
      )
    }
  }

  # ----------------------------------------------------------
  # 16. Assemble via .new_NICHESObject()
  # ----------------------------------------------------------
  merged.obj <- .new_NICHESObject(
    cell.list              = merged.cell.list,
    cell.data              = merged.cell.data,
    cell.meta              = merged.cell.meta,
    edge.list              = merged.edge.list,
    edge.data              = merged.edge.data,
    edge.meta              = merged.edge.meta,
    neighborhood.cell.list = NULL,
    neighborhood.edge.list = NULL,
    neighborhood.meta      = NULL,
    gene.list              = merged.gene.list,
    gene.meta              = merged.gene.meta,
    LRM.list               = merged.LRM.list,
    LRM.meta               = merged.LRM.meta,
    segmentation           = merged.seg
  )

  attr(merged.obj, "mode") <- merged.mode
  if (!is.null(merged.sample.col))
    attr(merged.obj, "sample.col") <- merged.sample.col

  if (verbose)
    message(sprintf(
      "Merged %d NICHESObjects: %d cells, %d edges, %d LRM pairs.",
      length(obj.list),
      length(merged.cell.list),
      nrow(merged.edge.list),
      length(merged.LRM.list)
    ))

  merged.obj
}
