# ==============================================================================
# FILE:    create_NICHESObject.R
# PURPOSE: Top-level constructor for the NICHESObject S3 class. Orchestrates
#          the full pipeline: edge list construction → LR scoring → slot
#          assembly. Three analysis modes are supported:
#
#   "spatial"    — physical kNN or radius neighbourhood from x/y coordinates
#   "pseudobulk" — one directed edge per ordered cell type pair; scoring uses
#                  per-type mean expression (genes × types matrix)
#   "sampled"    — NICHES v1-style stratified random sampling of cell barcodes
#                  within each directed type pair crossing
#
# FUNCTIONS:
#   create_NICHESObject()  — main function
#   add_cell_meta()        — add / update cell.meta columns post-construction
#
# DEPS: data.table (>= 1.14), Matrix (>= 1.5)
# ==============================================================================


# ==============================================================================
#' Create a NICHESObject from Single-Cell Data
#'
#' @description
#' The primary entry point for NICHESv2. Supports three analysis modes that
#' differ in how the cell-to-cell edge list is defined:
#'
#' \describe{
#'   \item{\code{"spatial"}}{Edges are defined by physical proximity using
#'     k-nearest neighbours or a fixed radius in x/y coordinate space.
#'     Requires \code{x} and \code{y} columns in \code{meta.data}.}
#'   \item{\code{"pseudobulk"}}{One directed edge per ordered cell type pair
#'     (T² edges total, including autocrine). Scoring uses the per-type mean
#'     expression profile. \code{$cell.list} and \code{$cell.data} retain
#'     original single-cell resolution; only edge-level data are type-level.
#'     Requires \code{cell.type.col}.}
#'   \item{\code{"sampled"}}{NICHES v1-style stratified random sampling. For
#'     each directed type pair (A → B), up to \code{min(n_A, n_B)} cell
#'     barcodes are independently drawn from each side and paired positionally.
#'     All T² crossings are represented. Requires \code{cell.type.col}.}
#' }
#'
#' @section Slot construction summary:
#' \describe{
#'   \item{\code{$cell.list}}{Column names of \code{count.mtx} (original cell
#'     barcodes in all modes).}
#'   \item{\code{$cell.data}}{Long-format expression table (cell, gene, score,
#'     score.norm). Sparse: only rows where \code{score > 0} retained.
#'     Always built from original single-cell \code{count.mtx}.}
#'   \item{\code{$cell.meta}}{\code{meta.data} coerced to a data.table with a
#'     \code{cell} column. Always contains original cell barcodes.}
#'   \item{\code{$edge.list}}{Columns: \code{edge}
#'     (\code{"Sender|Receiver"}), \code{SendingCell}, \code{ReceivingCell}.
#'     In \code{"pseudobulk"} mode, SendingCell and ReceivingCell hold type
#'     labels rather than barcodes.}
#'   \item{\code{$edge.data}}{Per-edge LRM scores with \code{score.norm}
#'     normalised per edge (proportional).}
#'   \item{\code{$edge.meta}}{Sending / receiving metadata joined per edge.
#'     In \code{"pseudobulk"} mode, contains only edge identity columns since
#'     the edge nodes are type labels, not individual cells.}
#'   \item{\code{$neighborhood.*}}{Reserved as \code{NULL}.}
#'   \item{\code{$gene.list}}{Row names of \code{count.mtx}.}
#'   \item{\code{$gene.meta}}{Minimal table unless \code{gene.meta} supplied.}
#'   \item{\code{$LRM.list}}{Valid \code{"ligand|receptor"} strings after gene
#'     filtering.}
#'   \item{\code{$LRM.meta}}{LRM annotation table: LRM, ligand, receptor,
#'     plus any extra columns from \code{LRM.db}.}
#' }
#' The assembled object carries a \code{mode} attribute recording which of the
#' three modes was used.
#'
#' @param count.mtx   Gene × cell expression matrix (sparse \code{dgCMatrix} or
#'   dense \code{matrix}). Row names = gene symbols; column names = cell
#'   barcodes.
#' @param meta.data   \code{data.frame} with one row per cell. Row names must
#'   match column names of \code{count.mtx}. For \code{mode = "spatial"},
#'   must contain columns \code{"x"} and \code{"y"} (spatial coordinates).
#'   For non-spatial modes, must contain \code{cell.type.col}.
#' @param LRM.db      LR database: \code{"fantom5"}, \code{"omnipath"}, or a
#'   \code{data.frame} with at least \code{ligand.col} and \code{receptor.col}.
#' @param mode        Analysis mode. One of \code{"spatial"} (default),
#'   \code{"pseudobulk"}, or \code{"sampled"}.
#' @param cell.type.col \code{character(1)} or \code{NULL}. Column name in
#'   \code{meta.data} holding cell type labels. Required for
#'   \code{mode = "pseudobulk"} and \code{mode = "sampled"}; ignored for
#'   \code{mode = "spatial"}.
#' @param ligand.col  Column name for ligands in \code{LRM.db}. Default
#'   \code{"ligand"}.
#' @param receptor.col Column name for receptors in \code{LRM.db}. Default
#'   \code{"receptor"}.
#' @param k           \code{integer} or \code{NULL}. kNN neighbourhood size.
#'   Spatial mode only. Passed to \code{create_Edgelist_Spatial()}.
#' @param rad         \code{numeric} or \code{NULL}. Radius neighbourhood size.
#'   Spatial mode only. Passed to \code{create_Edgelist_Spatial()}.
#' @param seed        \code{integer(1)}. Random seed for reproducibility.
#'   Sampled mode only. Default \code{42L}.
#' @param max.cells   \code{integer(1)} or \code{NULL}. Per-crossing edge cap.
#'   Sampled mode only. Default \code{NULL} (cap = \code{min(n_i, n_j)}).
#' @param method      Scoring method: \code{"product"} (default) or
#'   \code{"minimum"}.
#' @param normalize.method Normalisation for \code{score.norm}:
#'   \code{"prop"} (default) or \code{"none"}.
#' @param gene.meta   Optional \code{data.frame} with a \code{"gene"} column
#'   and any additional gene-level annotation.
#' @param n.cores     CPU cores. Default 1. Uses socket clusters
#'   (\code{parallel::makeCluster}) — works on all platforms.
#' @param verbose     Print progress messages. Default \code{TRUE}.
#'
#' @return A \code{NICHESObject} with a \code{mode} attribute.
#'
#' @seealso \code{\link{compute_CellToCell}},
#'   \code{\link{create_Edgelist_Spatial}},
#'   \code{\link{create_Edgelist_PseudoBulk}},
#'   \code{\link{create_Edgelist_SamplingByType}},
#'   \code{\link{load_LRM_database}},
#'   \code{\link{Aggregate_CellToCell_Edges}},
#'   \code{\link{CellToCell_InSitu}}
#'
#' @examples
#' \dontrun{
#' # Spatial mode (from a Seurat object)
#' count.mtx <- Seurat::GetAssayData(seu, slot = "counts")
#' meta.data  <- seu@meta.data   # must contain x, y columns
#' niches.obj <- create_NICHESObject(count.mtx, meta.data,
#'                                   LRM.db = "omnipath", k = 6)
#'
#' # Pseudobulk mode
#' niches.obj <- create_NICHESObject(count.mtx, meta.data,
#'                                   LRM.db = "omnipath",
#'                                   mode = "pseudobulk",
#'                                   cell.type.col = "cell.type")
#'
#' # Sampled mode
#' niches.obj <- create_NICHESObject(count.mtx, meta.data,
#'                                   LRM.db = "omnipath",
#'                                   mode = "sampled",
#'                                   cell.type.col = "cell.type",
#'                                   seed = 42L, max.cells = 100L)
#' }
#'
#' @export
create_NICHESObject <- function(count.mtx,
                                 meta.data,
                                 LRM.db,
                                 mode             = c("spatial", "pseudobulk", "sampled"),
                                 cell.type.col    = NULL,
                                 ligand.col       = "ligand",
                                 receptor.col     = "receptor",
                                 k                = NULL,
                                 rad              = NULL,
                                 seed             = 42L,
                                 max.cells        = NULL,
                                 method           = c("product", "minimum"),
                                 normalize.method = c("prop", "none"),
                                 gene.meta        = NULL,
                                 n.cores          = 1L,
                                 verbose          = TRUE) {

  mode             <- match.arg(mode)
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

  # Mode-specific input requirements
  if (mode == "spatial") {
    if (!all(c("x", "y") %in% colnames(meta.data)))
      stop(
        "[create_NICHESObject] mode = 'spatial' requires columns 'x' and 'y' in 'meta.data'.\n",
        "  If your data uses different column names, please rename them to 'x' and 'y'."
      )
  } else {
    if (is.null(cell.type.col))
      stop(sprintf(
        "[create_NICHESObject] mode = '%s' requires 'cell.type.col' to be specified.", mode
      ))
    if (!is.character(cell.type.col) || length(cell.type.col) != 1L)
      stop("[create_NICHESObject] 'cell.type.col' must be a single character string.")
    if (!cell.type.col %in% colnames(meta.data))
      stop(sprintf(
        "[create_NICHESObject] 'cell.type.col' ('%s') not found in 'meta.data'.",
        cell.type.col
      ))
  }

  # Cell alignment: keep intersection of meta.data rownames and count.mtx colnames
  shared.cells <- intersect(rownames(meta.data), colnames(count.mtx))
  if (length(shared.cells) == 0L)
    stop("[create_NICHESObject] No shared cells between 'meta.data' row names ",
         "and 'count.mtx' column names.")

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
  meta.data  <- meta.data[shared.cells, , drop = FALSE]

  # ════════════════════════════════════════════════════════════════════════════
  # SECTION 2 ── Build edge list (mode-dependent)
  # ════════════════════════════════════════════════════════════════════════════

  if (verbose) message(sprintf(
    "[create_NICHESObject] Building edge list (mode = '%s') ...", mode
  ))

  if (mode == "spatial") {

    coord.df <- data.frame(
      x = meta.data[["x"]],
      y = meta.data[["y"]],
      row.names = shared.cells
    )
    el.raw <- create_Edgelist_Spatial(coord = coord.df, k = k, rad = rad,
                                      n.cores = n.cores)

  } else if (mode == "pseudobulk") {

    el.raw <- create_Edgelist_PseudoBulk(meta.data     = meta.data,
                                          cell.type.col = cell.type.col)

  } else {  # mode == "sampled"

    el.raw <- create_Edgelist_SamplingByType(meta.data     = meta.data,
                                              cell.type.col = cell.type.col,
                                              seed          = seed,
                                              max.cells     = max.cells)
  }

  # ════════════════════════════════════════════════════════════════════════════
  # SECTION 3 ── Compute LR scores
  #
  # Pseudobulk mode: pass a genes × cell_types mean-expression matrix so that
  # compute_CellToCell() can look up type-label "cells" in colnames().
  # All other modes: pass the original single-cell count matrix.
  # ════════════════════════════════════════════════════════════════════════════

  if (verbose) message("[create_NICHESObject] Computing LR scores ...")

  if (mode == "pseudobulk") {
    # Compute per-type mean expression: genes × cell_types dense matrix.
    # Types are derived from the filtered meta.data (after cell alignment).
    cell.types.ordered <- attr(el.raw, "cell.names")   # sorted type labels
    pb.mtx <- do.call(cbind, lapply(cell.types.ordered, function(ct) {
      ct.barcodes <- rownames(meta.data)[meta.data[[cell.type.col]] == ct]
      Matrix::rowMeans(count.mtx[, ct.barcodes, drop = FALSE])
    }))
    colnames(pb.mtx) <- cell.types.ordered
    score.mtx <- pb.mtx
  } else {
    score.mtx <- count.mtx
  }

  elr <- compute_CellToCell(
    count.mtx    = score.mtx,
    edge.list    = el.raw,
    LRM.db       = LRM.db,
    ligand.col   = ligand.col,
    receptor.col = receptor.col,
    method       = method,
    n.cores      = n.cores
  )

  # ════════════════════════════════════════════════════════════════════════════
  # SECTION 4 ── Build $cell.list and $gene.list
  #
  # Always derived from the original count.mtx (cell barcodes), not the
  # pseudobulk matrix, so single-cell resolution is preserved in all modes.
  # ════════════════════════════════════════════════════════════════════════════

  cell.list <- shared.cells
  gene.list <- rownames(count.mtx)

  # ════════════════════════════════════════════════════════════════════════════
  # SECTION 5 ── Build $cell.data (always from original single-cell matrix)
  # ════════════════════════════════════════════════════════════════════════════

  if (verbose) message("[create_NICHESObject] Building cell.data ...")

  cell.data <- .build_cell_data(count.mtx, normalize.method)

  # ════════════════════════════════════════════════════════════════════════════
  # SECTION 6 ── Build $cell.meta (always original cell barcodes)
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
    edge  = paste(elr$edges$from, elr$edges$to, sep = "|"),
    LRM   = elr$edges$LRM,
    score = elr$edges$score
  )

  if (normalize.method == "prop") {
    edge.data[, score.norm := .prop_norm(score), by = edge]
  } else {
    edge.data[, score.norm := score]
  }

  setkeyv(edge.data, c("edge", "LRM"))

  # ════════════════════════════════════════════════════════════════════════════
  # SECTION 9 ── Build $edge.meta
  #
  # Pseudobulk: edges connect type labels, not barcodes. We pass a minimal
  # type-label table (no extra columns) to .build_edge_meta() so the result
  # contains edge, SendingCell, and ReceivingCell columns only — which is all
  # that is meaningful when nodes are type names rather than individual cells.
  #
  # Spatial and sampled: edges connect cell barcodes, so we join the full
  # cell.meta as usual to get Sending.*, Receiving.*, and Joint.* columns.
  # ════════════════════════════════════════════════════════════════════════════

  if (verbose) message("[create_NICHESObject] Building edge.meta ...")

  if (mode == "pseudobulk") {
    type.meta.dt <- data.table(cell = attr(el.raw, "cell.names"))
    setkeyv(type.meta.dt, "cell")
    edge.meta <- .build_edge_meta(edge.list.dt, type.meta.dt)
  } else {
    edge.meta <- .build_edge_meta(edge.list.dt, cell.meta)
  }

  # ════════════════════════════════════════════════════════════════════════════
  # SECTION 10 ── Build $gene.meta
  # ════════════════════════════════════════════════════════════════════════════

  if (!is.null(gene.meta)) {
    gene.meta.dt <- as.data.table(gene.meta)
    if (!"gene" %in% colnames(gene.meta.dt))
      stop("[create_NICHESObject] 'gene.meta' must have a 'gene' column.")
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

  # Stamp the mode as an object-level attribute for downstream inspection
  attr(obj, "mode") <- mode

  if (verbose) {
    message(sprintf(
      "[create_NICHESObject] Done. mode = '%s' | %d cells | %d edges | %d LRM pairs | %d non-zero edge signals.",
      mode, length(cell.list), nrow(edge.list.dt), length(LRM.list), nrow(edge.data)
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

  tsm  <- as(count.mtx, "TsparseMatrix")
  rows <- rownames(count.mtx)
  cols <- colnames(count.mtx)

  dt <- data.table(
    cell  = cols[tsm@j + 1L],
    gene  = rows[tsm@i + 1L],
    score = tsm@x
  )

  if (normalize.method == "prop") {
    dt[, score.norm := .prop_norm(score), by = cell]
  } else {
    dt[, score.norm := score]
  }

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

  meta.cols <- setdiff(colnames(cell.meta), "cell")
  edge.meta <- copy(edge.list.dt)

  if (length(meta.cols) == 0L) return(edge.meta)

  send.meta <- cell.meta[edge.meta[, .(cell = SendingCell)], on = "cell"]
  for (col in meta.cols) {
    edge.meta[[paste0("Sending.", col)]] <- send.meta[[col]]
  }

  recv.meta <- cell.meta[edge.meta[, .(cell = ReceivingCell)], on = "cell"]
  for (col in meta.cols) {
    edge.meta[[paste0("Receiving.", col)]] <- recv.meta[[col]]
  }

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
