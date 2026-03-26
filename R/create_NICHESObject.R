# ==============================================================================
# FILE:    create_NICHESObject.R
# PURPOSE: Top-level constructor for the NICHESObject S3 class. Orchestrates
#          the full pipeline: coordinate extraction → edge list → LR scoring →
#          slot assembly. Accepts Seurat-like inputs (count matrix + meta.data
#          data.frame) or raw R objects.
#
# FUNCTIONS:
#   create_NICHESObject()  — main function
#   add_cell_meta()        — add / update cell.meta columns post-construction
#
# DEPS: data.table (>= 1.14), Matrix (>= 1.5)
# ==============================================================================


# ==============================================================================
#' Create a NICHESObject from Spatial Single-Cell Data
#'
#' @description
#' The primary entry point for NICHESv2. Takes a gene × cell expression matrix
#' and a cell metadata table (which must contain spatial coordinates \code{x}
#' and \code{y}), builds a spatial neighbourhood graph, computes all valid
#' ligand-receptor scores, and assembles the result into a fully populated
#' \code{NICHESObject} with 13 slots.
#'
#' @section Slot construction summary:
#' \describe{
#'   \item{\code{$cell.list}}{Column names of \code{count.mtx}.}
#'   \item{\code{$cell.data}}{Long-format expression table (cell, gene, score,
#'     score.norm). Sparse: only rows where \code{score > 0} or
#'     \code{score.norm > 0} are retained. \code{score.norm} is per-cell
#'     proportional (gene count / total cell count).}
#'   \item{\code{$cell.meta}}{\code{meta.data} coerced to a data.table with a
#'     \code{cell} column added (= row names of \code{meta.data}).}
#'   \item{\code{$edge.list}}{Derived from \code{create_Edgelist_Spatial()};
#'     columns: \code{edge} (\code{"SendingCell|ReceivingCell"}),
#'     \code{SendingCell}, \code{ReceivingCell}.}
#'   \item{\code{$edge.data}}{Per-edge LRM scores; \code{score.norm} normalised
#'     \emph{per edge} (each LRM's score divided by that edge's total score).}
#'   \item{\code{$edge.meta}}{Sending / receiving cell metadata joined per edge.
#'     Categorical columns get a \code{Joint.*} column
#'     (\code{"SendingValue|ReceivingValue"}). Numeric metadata columns are
#'     added with \code{Sending.*} and \code{Receiving.*} prefixes only.}
#'   \item{\code{$neighborhood.*}}{Reserved as \code{NULL}.}
#'   \item{\code{$gene.list}}{Row names of \code{count.mtx}.}
#'   \item{\code{$gene.meta}}{Minimal table (\code{gene} column only) unless
#'     \code{gene.meta} is supplied.}
#'   \item{\code{$LRM.list}}{Character vector of \code{"ligand|receptor"}
#'     strings for all LRM pairs that pass gene-availability filtering.}
#'   \item{\code{$LRM.meta}}{Standardised LRM annotation table: \code{LRM},
#'     \code{ligand}, \code{receptor}, plus any extra columns from
#'     \code{LRM.db}.}
#' }
#'
#' @param count.mtx   Gene × cell expression matrix (sparse \code{dgCMatrix} or
#'   dense \code{matrix}). Row names = gene symbols; column names = cell
#'   barcodes.
#' @param meta.data   \code{data.frame} with one row per cell. Row names must
#'   match column names of \code{count.mtx}. \strong{Must contain columns
#'   \code{"x"} and \code{"y"}} (spatial coordinates). All other columns are
#'   propagated to \code{$cell.meta} and \code{$edge.meta}.
#' @param LRM.db      LR database: \code{"fantom5"}, \code{"omnipath"}, or a
#'   \code{data.frame} with at least \code{ligand.col} and \code{receptor.col}.
#' @param ligand.col  Column name for ligands in \code{LRM.db}. Default
#'   \code{"ligand"}.
#' @param receptor.col Column name for receptors in \code{LRM.db}. Default
#'   \code{"receptor"}.
#' @param k           \code{integer} or \code{NULL}. kNN neighbourhood size.
#'   Passed to \code{create_Edgelist_Spatial()}.
#' @param rad         \code{numeric} or \code{NULL}. Radius neighbourhood size.
#'   Passed to \code{create_Edgelist_Spatial()}.
#' @param method      Scoring method: \code{"product"} (default) or
#'   \code{"minimum"}.
#' @param normalize.method Normalisation for \code{score.norm}:
#'   \code{"prop"} (default, proportional within group) or \code{"none"}.
#' @param gene.meta   Optional \code{data.frame} with a \code{"gene"} column
#'   and any additional gene-level annotation. If \code{NULL}, a minimal
#'   table with only the \code{gene} column is created.
#' @param n.cores     CPU cores. Default 1. Fork-based on Unix; serial on
#'   Windows.
#' @param verbose     Print progress messages. Default \code{TRUE}.
#'
#' @return A \code{NICHESObject}.
#'
#' @seealso \code{\link{compute_CellToCell_Spatial}},
#'   \code{\link{create_Edgelist_Spatial}}, \code{\link{load_LRM_database}},
#'   \code{\link{Aggregate_CellToCell_Edges}}, \code{\link{CellToCell_InSitu}}
#'
#' @examples
#' \dontrun{
#' # From a Seurat object
#' count.mtx  <- Seurat::GetAssayData(seu, slot = "counts")
#' meta.data  <- seu@meta.data   # must contain x, y columns
#'
#' niches.obj <- create_NICHESObject(
#'   count.mtx = count.mtx,
#'   meta.data = meta.data,
#'   LRM.db    = "omnipath",
#'   k         = 6,
#'   method    = "product"
#' )
#'
#' # From raw matrices
#' niches.obj <- create_NICHESObject(
#'   count.mtx = my.count.matrix,
#'   meta.data = my.meta.data,
#'   LRM.db    = my.lr.dataframe,
#'   rad       = 200
#' )
#' }
#'
#' @export
create_NICHESObject <- function(count.mtx,
                                 meta.data,
                                 LRM.db,
                                 ligand.col        = "ligand",
                                 receptor.col      = "receptor",
                                 k                 = NULL,
                                 rad               = NULL,
                                 method            = c("product", "minimum"),
                                 normalize.method  = c("prop", "none"),
                                 gene.meta         = NULL,
                                 n.cores           = 1L,
                                 verbose           = TRUE) {

  method           <- match.arg(method)
  normalize.method <- match.arg(normalize.method)
  n.cores          <- as.integer(n.cores)

  # ════════════════════════════════════════════════════════════════════════════
  # SECTION 1 ── Validate core inputs
  # ════════════════════════════════════════════════════════════════════════════

  if (!is.matrix(count.mtx) && !inherits(count.mtx, "sparseMatrix"))
    stop("[create_NICHESObject] 'count.mtx' must be a matrix or sparseMatrix.")
  if (is.null(rownames(count.mtx)) || is.null(colnames(count.mtx)))
    stop("[create_NICHESObject] 'count.mtx' must have row names (genes) and column names (cells).")

  if (!is.data.frame(meta.data))
    stop("[create_NICHESObject] 'meta.data' must be a data.frame.")
  if (is.null(rownames(meta.data)))
    stop("[create_NICHESObject] 'meta.data' must have row names matching column names of 'count.mtx'.")

  # Spatial coordinate check
  if (!all(c("x", "y") %in% colnames(meta.data)))
    stop(
      "[create_NICHESObject] 'meta.data' must contain columns 'x' and 'y' ",
      "(spatial coordinates).\n",
      "  If your data uses different column names, please rename them to 'x' ",
      "and 'y' before calling create_NICHESObject().\n",
      "  Non-spatial analysis will be supported in a future version."
    )

  # Cell alignment: keep intersection of meta.data rownames and count.mtx colnames
  shared.cells <- intersect(rownames(meta.data), colnames(count.mtx))
  if (length(shared.cells) == 0L)
    stop("[create_NICHESObject] No shared cells between 'meta.data' row names ",
         "and 'count.mtx' column names.")

  n.dropped.cells <- (ncol(count.mtx) + nrow(meta.data)) -
    (length(colnames(count.mtx)) + length(rownames(meta.data))) +
    length(shared.cells)  # simplified: just warn if any mismatch

  if (length(shared.cells) < ncol(count.mtx) && verbose)
    message(sprintf(
      "[create_NICHESObject] %d cell(s) in count.mtx not found in meta.data — excluded.",
      ncol(count.mtx) - length(shared.cells)
    ))
  if (length(shared.cells) < nrow(meta.data) && verbose)
    message(sprintf(
      "[create_NICHESObject] %d cell(s) in meta.data not found in count.mtx — excluded.",
      nrow(meta.data) - length(shared.cells)
    ))

  count.mtx <- count.mtx[, shared.cells, drop = FALSE]
  meta.data <- meta.data[shared.cells, , drop = FALSE]

  # ════════════════════════════════════════════════════════════════════════════
  # SECTION 2 ── Build spatial edge list
  # ════════════════════════════════════════════════════════════════════════════

  if (verbose) message("[create_NICHESObject] Building spatial edge list ...")

  coord.df <- data.frame(
    x = meta.data[["x"]],
    y = meta.data[["y"]],
    row.names = shared.cells
  )

  el.raw <- create_Edgelist_Spatial(coord = coord.df, k = k, rad = rad, n.cores = n.cores)

  # ════════════════════════════════════════════════════════════════════════════
  # SECTION 3 ── Compute LR scores
  # ════════════════════════════════════════════════════════════════════════════

  if (verbose) message("[create_NICHESObject] Computing LR scores ...")

  elr <- compute_CellToCell_Spatial(
    count.mtx    = count.mtx,
    edge.list    = el.raw,
    LRM.db       = LRM.db,
    ligand.col   = ligand.col,
    receptor.col = receptor.col,
    method       = method,
    n.cores      = n.cores
  )

  # ════════════════════════════════════════════════════════════════════════════
  # SECTION 4 ── Build $cell.list and $gene.list
  # ════════════════════════════════════════════════════════════════════════════

  cell.list <- shared.cells
  gene.list <- rownames(count.mtx)

  # ════════════════════════════════════════════════════════════════════════════
  # SECTION 5 ── Build $cell.data (sparse, per-cell normalised)
  #
  # Convert count matrix to long format, filter zeros, add score.norm.
  # For large matrices we iterate over cells in batches to avoid
  # a single dense materialisation.
  # ════════════════════════════════════════════════════════════════════════════

  if (verbose) message("[create_NICHESObject] Building cell.data ...")

  cell.data <- .build_cell_data(count.mtx, normalize.method)

  # ════════════════════════════════════════════════════════════════════════════
  # SECTION 6 ── Build $cell.meta
  # ════════════════════════════════════════════════════════════════════════════

  cell.meta <- as.data.table(meta.data, keep.rownames = "cell")
  setkeyv(cell.meta, "cell")

  # ════════════════════════════════════════════════════════════════════════════
  # SECTION 7 ── Build $edge.list
  # ════════════════════════════════════════════════════════════════════════════

  edge.list.dt <- data.table(
    edge          = paste(el.raw$from, el.raw$to, sep = "|"),
    SendingCell   = el.raw$from,
    ReceivingCell = el.raw$to
  )
  setkeyv(edge.list.dt, "edge")

  # ════════════════════════════════════════════════════════════════════════════
  # SECTION 8 ── Build $edge.data (sparse, per-edge normalised)
  # ════════════════════════════════════════════════════════════════════════════

  if (verbose) message("[create_NICHESObject] Building edge.data ...")

  edge.data <- data.table(
    edge     = paste(elr$edges$from, elr$edges$to, sep = "|"),
    LRM      = elr$edges$LRM,
    score    = elr$edges$score
  )

  # Per-edge normalisation of score.norm
  if (normalize.method == "prop") {
    edge.data[, score.norm := .prop_norm(score), by = edge]
  } else {
    edge.data[, score.norm := score]
  }

  setkeyv(edge.data, c("edge", "LRM"))

  # ════════════════════════════════════════════════════════════════════════════
  # SECTION 9 ── Build $edge.meta
  # ════════════════════════════════════════════════════════════════════════════

  if (verbose) message("[create_NICHESObject] Building edge.meta ...")

  edge.meta <- .build_edge_meta(edge.list.dt, cell.meta)

  # ════════════════════════════════════════════════════════════════════════════
  # SECTION 10 ── Build $gene.meta
  # ════════════════════════════════════════════════════════════════════════════

  if (!is.null(gene.meta)) {
    gene.meta.dt <- as.data.table(gene.meta)
    if (!"gene" %in% colnames(gene.meta.dt))
      stop("[create_NICHESObject] 'gene.meta' must have a 'gene' column.")
    # Fill in any genes missing from the supplied table
    missing.genes <- setdiff(gene.list, gene.meta.dt$gene)
    if (length(missing.genes) > 0L)
      gene.meta.dt <- rbindlist(
        list(gene.meta.dt, data.table(gene = missing.genes)),
        fill = TRUE
      )
    setkeyv(gene.meta.dt, "gene")
  } else {
    gene.meta.dt <- data.table(gene = gene.list)
    setkeyv(gene.meta.dt, "gene")
  }

  # ════════════════════════════════════════════════════════════════════════════
  # SECTION 11 ── Build $LRM.list and $LRM.meta
  # ════════════════════════════════════════════════════════════════════════════

  LRM.meta.dt <- copy(elr$lr.meta)
  LRM.list    <- LRM.meta.dt$LRM

  # Ensure LRM column is first, followed by ligand, receptor, then extras
  other.cols <- setdiff(colnames(LRM.meta.dt), c("LRM", "ligand", "receptor"))
  setcolorder(LRM.meta.dt, c("LRM", "ligand", "receptor", other.cols))
  setkeyv(LRM.meta.dt, "LRM")

  # ════════════════════════════════════════════════════════════════════════════
  # SECTION 12 ── Assemble NICHESObject
  # ════════════════════════════════════════════════════════════════════════════

  obj <- .new_NICHESObject(
    cell.list              = cell.list,
    cell.data              = cell.data,
    cell.meta              = cell.meta,
    edge.list              = edge.list.dt,
    edge.data              = edge.data,
    edge.meta              = edge.meta,
    neighborhood.cell.list = NULL,
    neighborhood.edge.list = NULL,
    neighborhood.meta      = NULL,
    gene.list              = gene.list,
    gene.meta              = gene.meta.dt,
    LRM.list               = LRM.list,
    LRM.meta               = LRM.meta.dt
  )

  if (verbose) {
    message(sprintf(
      "[create_NICHESObject] Done. %d cells | %d edges | %d LRM pairs | %d non-zero edge signals.",
      length(cell.list), nrow(edge.list.dt), length(LRM.list), nrow(edge.data)
    ))
  }

  obj
}


# ==============================================================================
#' Add or Update Columns in a NICHESObject's cell.meta Slot
#'
#' @description
#' Joins additional metadata columns onto \code{$cell.meta} and rebuilds
#' \code{$edge.meta} to reflect the change. Useful for adding cell type
#' annotations, clustering results, or QC metrics after initial construction.
#'
#' @param niches.obj A \code{NICHESObject}.
#' @param new.meta   \code{data.frame} or \code{data.table} with a \code{"cell"}
#'   column and one or more additional columns to add or overwrite.
#' @param rebuild.edge.meta Logical. Whether to rebuild \code{$edge.meta} after
#'   the join. Default \code{TRUE}. Set \code{FALSE} to skip if you plan to
#'   call \code{add_cell_meta()} multiple times before needing \code{$edge.meta}.
#'
#' @return An updated \code{NICHESObject}.
#' @export
add_cell_meta <- function(niches.obj, new.meta, rebuild.edge.meta = TRUE) {

  if (!is_NICHESObject(niches.obj))
    stop("[add_cell_meta] 'niches.obj' must be a NICHESObject.")

  new.meta <- as.data.table(new.meta)
  if (!"cell" %in% colnames(new.meta))
    stop("[add_cell_meta] 'new.meta' must have a 'cell' column.")

  # Overwrite / add columns via update join
  new.cols <- setdiff(colnames(new.meta), "cell")
  niches.obj$cell.meta <- new.meta[niches.obj$cell.meta, on = "cell"]

  if (rebuild.edge.meta)
    niches.obj$edge.meta <- .build_edge_meta(niches.obj$edge.list, niches.obj$cell.meta)

  niches.obj
}


# ==============================================================================
# INTERNAL CONSTRUCTION HELPERS
# ==============================================================================

#' Build $cell.data: sparse long-format expression table with normalised scores.
#' @noRd
.build_cell_data <- function(count.mtx, normalize.method) {

  # Convert to TsparseMatrix for fast triplet extraction
  tsm  <- as(count.mtx, "TsparseMatrix")
  rows <- rownames(count.mtx)
  cols <- colnames(count.mtx)

  dt <- data.table(
    cell  = cols[tsm@j + 1L],
    gene  = rows[tsm@i + 1L],
    score = tsm@x
  )

  # Per-cell proportional normalisation
  if (normalize.method == "prop") {
    dt[, score.norm := .prop_norm(score), by = cell]
  } else {
    dt[, score.norm := score]
  }

  # Keep only rows where at least one value is non-zero (already guaranteed
  # for score by sparse format, but defensive for score.norm == 0)
  dt <- dt[score > 0 | score.norm > 0]

  setkeyv(dt, c("cell", "gene"))
  dt
}


#' Build $edge.meta: per-edge metadata with Sending/Receiving/Joint columns.
#'
#' Categorical meta columns get both Sending.*, Receiving.* and Joint.*
#' (value concatenated with "|"). Numeric columns get Sending.* and
#' Receiving.* only (Joint.* deferred pending future discussion).
#' @noRd
.build_edge_meta <- function(edge.list.dt, cell.meta) {

  # Identify metadata columns (exclude the key "cell" column)
  meta.cols <- setdiff(colnames(cell.meta), "cell")

  # Start with the edge identity columns
  edge.meta <- copy(edge.list.dt)   # edge, SendingCell, ReceivingCell

  if (length(meta.cols) == 0L) return(edge.meta)

  # Join sending-cell metadata
  send.meta <- cell.meta[edge.meta[, .(cell = SendingCell)], on = "cell"]
  for (col in meta.cols) {
    edge.meta[[paste0("Sending.", col)]] <- send.meta[[col]]
  }

  # Join receiving-cell metadata
  recv.meta <- cell.meta[edge.meta[, .(cell = ReceivingCell)], on = "cell"]
  for (col in meta.cols) {
    edge.meta[[paste0("Receiving.", col)]] <- recv.meta[[col]]
  }

  # Joint columns for categorical variables only
  cat.cols <- .categorical_cols(cell.meta, exclude = "cell")
  for (col in cat.cols) {
    s.col <- paste0("Sending.",   col)
    r.col <- paste0("Receiving.", col)
    edge.meta[[paste0("Joint.", col)]] <- paste(
      edge.meta[[s.col]], edge.meta[[r.col]], sep = "|"
    )
  }

  setkeyv(edge.meta, "edge")
  edge.meta
}
