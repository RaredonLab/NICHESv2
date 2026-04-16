# NICHESObject_class.R ----------------------------------------------------
# S3 class definition for NICHESObject -- the central data structure
# of the NICHESv2 package. Defines the 14-slot specification and all
# S3 methods: print, summary, dim, length, [, and is_NICHESObject.
#
# DEPS: data.table


## LOW-LEVEL CONSTRUCTOR (internal) ----
#
# Creates a validated NICHESObject from its 14 constituent slots. Called by
# create_NICHESObject() -- end users should never call this directly.
#
# All slot-level validation lives here so that every NICHESObject, regardless
# of how it was built, satisfies the invariants.

#' @noRd
.new_NICHESObject <- function(
    cell.list,
    cell.data,
    cell.meta,
    edge.list,
    edge.data,
    edge.meta,
    neighborhood.cell.list = NULL,
    neighborhood.edge.list = NULL,
    neighborhood.meta      = NULL,
    gene.list,
    gene.meta,
    LRM.list,
    LRM.meta,
    segmentation           = NULL
) {

  # -- Slot-type invariants ---------------------------------------------------
  stopifnot(
    is.character(cell.list),
    is.data.table(cell.data),
    is.data.table(cell.meta),
    is.data.table(edge.list),
    is.data.table(edge.data),
    is.data.table(edge.meta),
    is.character(gene.list),
    is.data.table(gene.meta),
    is.character(LRM.list),
    is.data.table(LRM.meta)
  )

  # -- Required column invariants ---------------------------------------------
  .require_cols <- function(dt, cols, slot.name) {
    missing <- setdiff(cols, colnames(dt))
    if (length(missing) > 0L)
      stop(sprintf("[NICHESObject] Slot '%s' is missing required columns: %s",
                   slot.name, paste(missing, collapse = ", ")))
  }

  .require_cols(cell.data, c("cell", "gene", "score", "score.norm"),    "cell.data")
  .require_cols(cell.meta, c("cell"),                                    "cell.meta")
  .require_cols(edge.list, c("edge", "SendingCell", "ReceivingCell"),    "edge.list")
  .require_cols(edge.data, c("edge", "LRM", "score", "score.norm"),     "edge.data")
  .require_cols(edge.meta, c("edge", "SendingCell", "ReceivingCell"),    "edge.meta")
  .require_cols(gene.meta, c("gene"),                                    "gene.meta")
  .require_cols(LRM.meta,  c("LRM", "ligand", "receptor"),              "LRM.meta")

  structure(
    list(
      cell.list              = cell.list,
      cell.data              = cell.data,
      cell.meta              = cell.meta,
      edge.list              = edge.list,
      edge.data              = edge.data,
      edge.meta              = edge.meta,
      neighborhood.cell.list = neighborhood.cell.list,
      neighborhood.edge.list = neighborhood.edge.list,
      neighborhood.meta      = neighborhood.meta,
      gene.list              = gene.list,
      gene.meta              = gene.meta,
      LRM.list               = LRM.list,
      LRM.meta               = LRM.meta,
      segmentation           = segmentation
    ),
    class = "NICHESObject"
  )
}


## PUBLIC IDENTITY CHECK ----

#' Test Whether an Object Is a NICHESObject
#'
#' Convenience wrapper around \code{inherits(x, "NICHESObject")}. Useful in
#' downstream pipelines that need to validate inputs before calling slot-aware
#' helpers.
#'
#' @param x Any R object.
#'
#' @return Logical scalar: \code{TRUE} if \code{x} inherits from class
#'   \code{"NICHESObject"}, otherwise \code{FALSE}.
#'
#' @examples
#' is_NICHESObject(list())   # FALSE
#' is_NICHESObject(1:10)     # FALSE
#'
#' @export
is_NICHESObject <- function(x) inherits(x, "NICHESObject")


## PRINT METHOD ----

#' Print a NICHESObject
#'
#' Compact one-screen summary of a \code{NICHESObject}: cell, gene, edge, and
#' LRM-pair counts (with non-empty / valid counts where applicable),
#' segmentation status, and a slot legend. Use \code{\link{summary.NICHESObject}}
#' for a longer breakdown including cell-type composition and signal density.
#'
#' @param x   A \code{NICHESObject}.
#' @param ... Unused. Present for S3 method signature conformance.
#'
#' @return The input \code{x}, returned invisibly. Called for the side effect
#'   of printing to the console.
#'
#' @seealso \code{\link{summary.NICHESObject}},
#'   \code{\link{dim.NICHESObject}}, \code{\link{is_NICHESObject}}.
#'
#' @examples
#' \donttest{
#' # Build a tiny NICHESObject for demonstration
#' set.seed(1L)
#' count.mtx <- Matrix::sparseMatrix(
#'   i = sample(6L, 30L, replace = TRUE),
#'   j = sample(8L, 30L, replace = TRUE),
#'   x = as.numeric(rpois(30L, lambda = 3) + 1L),
#'   dims = c(6L, 8L),
#'   dimnames = list(paste0("Gene", 1:6), paste0("Cell", 1:8))
#' )
#' meta.data <- data.frame(
#'   x = runif(8L), y = runif(8L),
#'   row.names = paste0("Cell", 1:8)
#' )
#' LRM.db <- data.frame(
#'   ligand   = c("Gene1", "Gene3"),
#'   receptor = c("Gene2", "Gene4")
#' )
#' obj <- create_NICHESObject(
#'   count.mtx = count.mtx, meta.data = meta.data, LRM.db = LRM.db,
#'   mode = "spatial", ligand.col = "ligand", receptor.col = "receptor",
#'   k = 3L, n.cores = 1L, verbose = FALSE
#' )
#'
#' print(obj)
#' obj            # implicit print
#' }
#'
#' @export
print.NICHESObject <- function(x, ...) {
  n.cells.nonempty <- uniqueN(x$cell.data$cell)
  n.genes.valid    <- uniqueN(x$cell.data$gene)
  n.edges.nonempty <- uniqueN(x$edge.data$edge)
  n.lrms.valid     <- uniqueN(x$edge.data$LRM)

  seg.status <- if (!is.null(x$segmentation)) {
    sprintf("attached (%d cells)", nrow(x$segmentation))
  } else {
    "none"
  }

  cat("-- NICHESObject --------------------------------------------\n")
  cat(sprintf("  Cells        : %d  (%d non-empty)\n",
              length(x$cell.list), n.cells.nonempty))
  cat(sprintf("  Genes        : %d  (%d valid)\n",
              length(x$gene.list), n.genes.valid))
  cat(sprintf("  Edges        : %d  (%d non-empty)\n",
              nrow(x$edge.list), n.edges.nonempty))
  cat(sprintf("  LRM pairs    : %d  (%d valid)\n",
              length(x$LRM.list), n.lrms.valid))
  cat(sprintf("  Segmentation : %s\n", seg.status))
  cat("\n  Slot summary:\n")
  cat("    $cell.list     -- character vector of cell barcodes\n")
  cat("    $cell.data     -- (cell, gene, score, score.norm) [sparse]\n")
  cat("    $cell.meta     -- per-cell metadata\n")
  cat("    $edge.list     -- (edge, SendingCell, ReceivingCell)\n")
  cat("    $edge.data     -- (edge, LRM, score, score.norm) [sparse]\n")
  cat("    $edge.meta     -- per-edge metadata with Sending/Receiving/Joint cols\n")
  cat("    $neighborhood.* -- reserved (NULL)\n")
  cat("    $gene.list     -- character vector of gene names\n")
  cat("    $gene.meta     -- per-gene metadata\n")
  cat("    $LRM.list      -- character vector of 'ligand|receptor' pairs\n")
  cat("    $LRM.meta      -- per-LRM metadata\n")
  cat("    $segmentation  -- sf object of cell polygons, or NULL\n")
  cat("------------------------------------------------------------\n")
  invisible(x)
}


## SUMMARY METHOD ----

#' Summarise a NICHESObject
#'
#' Long-form summary of a \code{NICHESObject}. Reports cell counts,
#' cell-type composition (when a \code{celltype} column is present in
#' \code{$cell.meta}), gene count, edge totals split into autocrine
#' (self-loop) versus inter-cell, LRM-pair count, edge-signal density
#' (filled / possible cells of the edge x LRM grid), how many of the
#' three neighborhood slots are populated, and segmentation status.
#'
#' Use \code{\link{print.NICHESObject}} for a shorter one-screen overview.
#'
#' @param object A \code{NICHESObject}.
#' @param ...    Unused. Present for S3 method signature conformance.
#'
#' @return The input \code{object}, returned invisibly. Called for the side
#'   effect of printing to the console.
#'
#' @seealso \code{\link{print.NICHESObject}}, \code{\link{dim.NICHESObject}}.
#'
#' @examples
#' \donttest{
#' set.seed(1L)
#' count.mtx <- Matrix::sparseMatrix(
#'   i = sample(6L, 30L, replace = TRUE),
#'   j = sample(8L, 30L, replace = TRUE),
#'   x = as.numeric(rpois(30L, lambda = 3) + 1L),
#'   dims = c(6L, 8L),
#'   dimnames = list(paste0("Gene", 1:6), paste0("Cell", 1:8))
#' )
#' meta.data <- data.frame(
#'   x = runif(8L), y = runif(8L),
#'   celltype = rep(c("A", "B"), each = 4L),
#'   row.names = paste0("Cell", 1:8)
#' )
#' LRM.db <- data.frame(
#'   ligand   = c("Gene1", "Gene3"),
#'   receptor = c("Gene2", "Gene4")
#' )
#' obj <- create_NICHESObject(
#'   count.mtx = count.mtx, meta.data = meta.data, LRM.db = LRM.db,
#'   mode = "spatial", ligand.col = "ligand", receptor.col = "receptor",
#'   k = 3L, n.cores = 1L, verbose = FALSE
#' )
#'
#' summary(obj)
#' }
#'
#' @export
summary.NICHESObject <- function(object, ...) {
  cat("-- NICHESObject Summary ------------------------------------\n")

  # Cell overview
  n.cells <- length(object$cell.list)
  cat(sprintf("  Cells: %d\n", n.cells))
  if ("celltype" %in% colnames(object$cell.meta)) {
    ct.tbl <- table(object$cell.meta$celltype)
    cat("  Cell types:\n")
    for (nm in names(ct.tbl))
      cat(sprintf("    %-20s : %d\n", nm, ct.tbl[[nm]]))
  }

  # Gene overview
  cat(sprintf("\n  Genes: %d\n", length(object$gene.list)))

  # Edge overview
  n.edges <- nrow(object$edge.list)
  cat(sprintf("\n  Edges: %d\n", n.edges))
  n.autocrine <- sum(object$edge.list$SendingCell == object$edge.list$ReceivingCell)
  cat(sprintf("    Autocrine (self-loop) : %d\n",  n.autocrine))
  cat(sprintf("    Inter-cell            : %d\n",  n.edges - n.autocrine))

  # LRM overview
  cat(sprintf("\n  LRM pairs: %d\n", length(object$LRM.list)))

  # Signal sparsity
  possible.ed <- as.numeric(n.edges) * length(object$LRM.list)
  actual.ed   <- nrow(object$edge.data)
  cat(sprintf("\n  Edge signal density: %d / %d (%.2f%%)\n",
              actual.ed, possible.ed,
              if (possible.ed > 0L) 100 * actual.ed / possible.ed else 0))

  # Neighborhood slots
  nb.filled <- sum(c(
    !is.null(object$neighborhood.cell.list),
    !is.null(object$neighborhood.edge.list),
    !is.null(object$neighborhood.meta)
  ))
  cat(sprintf("\n  Neighborhood slots populated: %d / 3\n", nb.filled))

  # Segmentation
  if (!is.null(object$segmentation)) {
    cat(sprintf("\n  Segmentation: attached (%d cells)\n", nrow(object$segmentation)))
  } else {
    cat("\n  Segmentation: none\n")
  }

  cat("------------------------------------------------------------\n")
  invisible(object)
}


## DIM METHOD ----

#' Dimensions of a NICHESObject
#'
#' Returns the four headline counts that characterise a \code{NICHESObject}:
#' the number of cells, genes, edges, and ligand--receptor (LRM) pairs.
#' Provided so that \code{dim(obj)} works on a \code{NICHESObject} the same
#' way it does on a matrix.
#'
#' @param x A \code{NICHESObject}.
#'
#' @return Named integer vector of length four with names \code{cells},
#'   \code{genes}, \code{edges}, \code{LRMs}.
#'
#' @seealso \code{\link{length.NICHESObject}} for cell count alone;
#'   \code{\link{summary.NICHESObject}} for a fuller breakdown.
#'
#' @examples
#' \donttest{
#' set.seed(1L)
#' count.mtx <- Matrix::sparseMatrix(
#'   i = sample(6L, 30L, replace = TRUE),
#'   j = sample(8L, 30L, replace = TRUE),
#'   x = as.numeric(rpois(30L, lambda = 3) + 1L),
#'   dims = c(6L, 8L),
#'   dimnames = list(paste0("Gene", 1:6), paste0("Cell", 1:8))
#' )
#' meta.data <- data.frame(
#'   x = runif(8L), y = runif(8L),
#'   row.names = paste0("Cell", 1:8)
#' )
#' LRM.db <- data.frame(
#'   ligand   = c("Gene1", "Gene3"),
#'   receptor = c("Gene2", "Gene4")
#' )
#' obj <- create_NICHESObject(
#'   count.mtx = count.mtx, meta.data = meta.data, LRM.db = LRM.db,
#'   mode = "spatial", ligand.col = "ligand", receptor.col = "receptor",
#'   k = 3L, n.cores = 1L, verbose = FALSE
#' )
#'
#' dim(obj)
#' }
#'
#' @export
dim.NICHESObject <- function(x) {
  c(
    cells = length(x$cell.list),
    genes = length(x$gene.list),
    edges = nrow(x$edge.list),
    LRMs  = length(x$LRM.list)
  )
}


## LENGTH METHOD ----

#' Number of Cells in a NICHESObject
#'
#' Returns the number of cells in a \code{NICHESObject}, identical to
#' \code{length(x$cell.list)}. Provided so that \code{length(obj)} works
#' the way most users expect for a per-cell collection.
#'
#' @param x A \code{NICHESObject}.
#'
#' @return Integer scalar: number of cells.
#'
#' @seealso \code{\link{dim.NICHESObject}} for cell, gene, edge, and LRM
#'   counts together.
#'
#' @examples
#' \donttest{
#' set.seed(1L)
#' count.mtx <- Matrix::sparseMatrix(
#'   i = sample(6L, 30L, replace = TRUE),
#'   j = sample(8L, 30L, replace = TRUE),
#'   x = as.numeric(rpois(30L, lambda = 3) + 1L),
#'   dims = c(6L, 8L),
#'   dimnames = list(paste0("Gene", 1:6), paste0("Cell", 1:8))
#' )
#' meta.data <- data.frame(
#'   x = runif(8L), y = runif(8L),
#'   row.names = paste0("Cell", 1:8)
#' )
#' LRM.db <- data.frame(
#'   ligand   = c("Gene1", "Gene3"),
#'   receptor = c("Gene2", "Gene4")
#' )
#' obj <- create_NICHESObject(
#'   count.mtx = count.mtx, meta.data = meta.data, LRM.db = LRM.db,
#'   mode = "spatial", ligand.col = "ligand", receptor.col = "receptor",
#'   k = 3L, n.cores = 1L, verbose = FALSE
#' )
#'
#' length(obj)   # 8
#' }
#'
#' @export
length.NICHESObject <- function(x) length(x$cell.list)


## SUBSET METHOD ----

#' Subset a NICHESObject to a Set of Cells
#'
#' Filters all cell-indexed slots (\code{$cell.list}, \code{$cell.data},
#' \code{$cell.meta}) and all edge-indexed slots
#' (\code{$edge.list}, \code{$edge.data}, \code{$edge.meta}) to keep only
#' edges where \strong{both} sending and receiving cells are in \code{i}.
#' If \code{$segmentation} is populated, it is filtered to retain only rows
#' whose \code{cell} is in the retained set. Slot-level invariants are
#' re-validated by the internal constructor before the subset is returned.
#'
#' Cell barcodes in \code{i} that are not present in \code{x$cell.list} are
#' silently dropped. If the intersection is empty the function errors rather
#' than returning an unusable empty object.
#'
#' @param x   A \code{NICHESObject}.
#' @param i   Character vector of cell barcodes to retain.
#' @param ... Unused. Present for S3 method signature conformance.
#'
#' @return A new \code{NICHESObject} restricted to the retained cells and
#'   their fully-internal edges. Gene and LRM slots are unchanged.
#'
#' @seealso \code{\link{is_NICHESObject}} to validate the result.
#'
#' @examples
#' \donttest{
#' set.seed(1L)
#' count.mtx <- Matrix::sparseMatrix(
#'   i = sample(6L, 30L, replace = TRUE),
#'   j = sample(8L, 30L, replace = TRUE),
#'   x = as.numeric(rpois(30L, lambda = 3) + 1L),
#'   dims = c(6L, 8L),
#'   dimnames = list(paste0("Gene", 1:6), paste0("Cell", 1:8))
#' )
#' meta.data <- data.frame(
#'   x = runif(8L), y = runif(8L),
#'   row.names = paste0("Cell", 1:8)
#' )
#' LRM.db <- data.frame(
#'   ligand   = c("Gene1", "Gene3"),
#'   receptor = c("Gene2", "Gene4")
#' )
#' obj <- create_NICHESObject(
#'   count.mtx = count.mtx, meta.data = meta.data, LRM.db = LRM.db,
#'   mode = "spatial", ligand.col = "ligand", receptor.col = "receptor",
#'   k = 3L, n.cores = 1L, verbose = FALSE
#' )
#'
#' # Subset to the first four cells; only fully-internal edges are kept
#' sub <- obj[paste0("Cell", 1:4)]
#' dim(sub)
#' }
#'
#' @export
`[.NICHESObject` <- function(x, i, ...) {

  if (!is.character(i))
    stop("[NICHESObject] Subset index 'i' must be a character vector of cell barcodes.")

  keep.cells <- intersect(i, x$cell.list)
  if (length(keep.cells) == 0L)
    stop("[NICHESObject] No cells from 'i' found in this object.")

  # Filter cell slots
  new.cell.list <- keep.cells
  new.cell.data <- x$cell.data[cell %in% keep.cells]
  new.cell.meta <- x$cell.meta[cell %in% keep.cells]

  # Filter edge slots -- keep only edges where BOTH endpoints are in keep.cells
  keep.edges <- x$edge.list[SendingCell %in% keep.cells & ReceivingCell %in% keep.cells, edge]
  new.edge.list <- x$edge.list[edge %in% keep.edges]
  new.edge.data <- x$edge.data[edge %in% keep.edges]
  new.edge.meta <- x$edge.meta[edge %in% keep.edges]

  # Filter segmentation if present -- keep only polygons for retained cells
  new.segmentation <- if (!is.null(x$segmentation)) {
    x$segmentation[x$segmentation$cell %in% keep.cells, ]
  } else {
    NULL
  }

  result <- .new_NICHESObject(
    cell.list              = new.cell.list,
    cell.data              = new.cell.data,
    cell.meta              = new.cell.meta,
    edge.list              = new.edge.list,
    edge.data              = new.edge.data,
    edge.meta              = new.edge.meta,
    neighborhood.cell.list = x$neighborhood.cell.list,
    neighborhood.edge.list = x$neighborhood.edge.list,
    neighborhood.meta      = x$neighborhood.meta,
    gene.list              = x$gene.list,
    gene.meta              = x$gene.meta,
    LRM.list               = x$LRM.list,
    LRM.meta               = x$LRM.meta,
    segmentation           = new.segmentation
  )
  attr(result, "mode")       <- attr(x, "mode")
  attr(result, "sample.col") <- attr(x, "sample.col")
  result
}
