# create_NICHESObject_single.R --------------------------------------------
# Single-partition implementation of create_NICHESObject(). The public
# dispatcher in create_NICHESObject.R routes to this function (or to
# .create_NICHESObject_by_sample()) depending on whether sample.col
# is supplied.
#
#          Orchestrates the full pipeline: edge list construction ->
#          LR scoring -> slot assembly. Two analysis modes are supported:
#
#   "spatial" - physical kNN or radius neighbourhood from x/y coordinates
#   "nonspatial" - NICHES v1-style stratified random sampling of cell barcodes
#               within each directed type pair crossing
#
#   Pseudobulk analysis is handled by the separate lightweight function
#   compute_PseudoBulk_Connectivity(), which returns a plain data.table
#   rather than a full NICHESObject.
#
# FUNCTIONS:
#   .create_NICHESObject_single()  - internal single-partition implementation
#
# DEPS: data.table (>= 1.14), Matrix (>= 1.5)


#' Single-partition implementation of create_NICHESObject
#'
#' @description
#' Internal. Implements the full NICHESObject construction pipeline for one
#' partition of cells. Called by the public dispatcher
#' \code{create_NICHESObject()} (or, for sample-aware construction, by
#' \code{.create_NICHESObject_by_sample()}).
#'
#' @param count.mtx   Gene x cell matrix. Row names = genes; col names = cells.
#' @param meta.data   \code{data.frame}, row names = cell barcodes.
#' @param LRM.db      \code{"fantom5"}, \code{"omnipath"}, or a data.frame.
#' @param mode        \code{"spatial"} (default) or \code{"nonspatial"}.
#' @param species     Species code passed through to \code{load_LRM_database()}.
#' @param cell.type.col Column in \code{meta.data} with cell type labels.
#'   Required for \code{mode = "nonspatial"}.
#' @param ligand.col  Default \code{"ligand"}.
#' @param receptor.col Default \code{"receptor"}.
#' @param k           kNN size. Spatial mode only.
#' @param rad         Radius. Spatial mode only.
#' @param seed        Random seed. Nonspatial mode only. Default \code{42L}.
#' @param max.cells   Per-crossing edge cap. Nonspatial mode only. Default
#'   \code{NULL}.
#' @param method      \code{"product"} (default) or \code{"minimum"}.
#' @param normalize.method \code{"prop"} (default) or \code{"none"}.
#' @param gene.meta   Optional data.frame with \code{"gene"} column.
#' @param n.cores     CPU cores. Default 1. Socket cluster -- all platforms.
#' @param verbose     Default \code{TRUE}.
#'
#' @return \code{NICHESObject} with \code{attr(obj, "mode")} set.
#'
#' @noRd
.create_NICHESObject_single <- function(count.mtx,
                                 meta.data,
                                 LRM.db,
                                 mode             = c("spatial", "nonspatial"),
                                 species          = "human",
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

  ## SECTION 1: Validate core inputs ----

  if (!is.matrix(count.mtx) && !inherits(count.mtx, "sparseMatrix"))
    stop("[create_NICHESObject] 'count.mtx' must be a matrix or sparseMatrix.")
  if (is.null(rownames(count.mtx)) || is.null(colnames(count.mtx)))
    stop("[create_NICHESObject] 'count.mtx' must have row names (genes) and column names (cells).")
  if (!is.data.frame(meta.data))
    stop("[create_NICHESObject] 'meta.data' must be a data.frame.")
  if (is.null(rownames(meta.data)))
    stop("[create_NICHESObject] 'meta.data' must have row names matching colnames of 'count.mtx'.")

  if (mode == "spatial") {
    if (!all(c("x", "y") %in% colnames(meta.data)))
      stop(
        "[create_NICHESObject] mode = 'spatial' requires columns 'x' and 'y' in 'meta.data'.\n",
        "  Rename your coordinate columns to 'x' and 'y' before calling this function."
      )
    if (!is.null(k) && !is.null(rad))
      stop("[create_NICHESObject] Provide either 'k' or 'rad', not both.")
    if (is.null(k) && is.null(rad))
      stop("[create_NICHESObject] mode = 'spatial' requires either 'k' or 'rad'.")
  } else {
    if (is.null(cell.type.col))
      stop("[create_NICHESObject] mode = 'nonspatial' requires 'cell.type.col' to be specified.")
    if (!is.character(cell.type.col) || length(cell.type.col) != 1L)
      stop("[create_NICHESObject] 'cell.type.col' must be a single character string.")
    if (!cell.type.col %in% colnames(meta.data))
      stop(sprintf(
        "[create_NICHESObject] 'cell.type.col' ('%s') not found in 'meta.data'.",
        cell.type.col
      ))
  }

  shared.cells <- intersect(rownames(meta.data), colnames(count.mtx))
  if (length(shared.cells) == 0L)
    stop("[create_NICHESObject] No shared cells between 'meta.data' and 'count.mtx'.")

  if (length(shared.cells) < ncol(count.mtx) && verbose)
    message(sprintf("[create_NICHESObject] %d cell(s) in count.mtx not in meta.data -- excluded.",
                    ncol(count.mtx) - length(shared.cells)))
  if (length(shared.cells) < nrow(meta.data) && verbose)
    message(sprintf("[create_NICHESObject] %d cell(s) in meta.data not in count.mtx -- excluded.",
                    nrow(meta.data) - length(shared.cells)))

  count.mtx <- count.mtx[, shared.cells, drop = FALSE]
  meta.data  <- meta.data[shared.cells, , drop = FALSE]

  ## SECTION 2: Build edge list (mode-dependent) ----

  if (verbose) message(sprintf(
    "[create_NICHESObject] Building edge list (mode = '%s') ...", mode
  ))

  if (mode == "spatial") {
    coord.df <- data.frame(x = meta.data[["x"]], y = meta.data[["y"]],
                            row.names = shared.cells)
    el.raw <- create_Edgelist_Spatial(coord = coord.df, k = k, rad = rad,
                                      n.cores = n.cores)
  } else {
    el.raw <- create_Edgelist_SamplingByType(meta.data     = meta.data,
                                              cell.type.col = cell.type.col,
                                              seed          = seed,
                                              max.cells     = max.cells)
  }

  ## SECTION 3: Prepare scoring matrix ----

  score.mtx <- count.mtx

  ## SECTION 4: Compute LR scores ----

  if (verbose) message("[create_NICHESObject] Computing LR scores ...")

  # Resolve string shorthand to a data.frame before scoring (handles species)
  if(is.character(LRM.db)){
    LRM.db <- load_LRM_database(
      db = LRM.db,
      species = species,
      ligand.col = ligand.col,
      receptor.col = receptor.col,
      verbose = verbose
    )
  }

  # `elr` = edge-LR result. Plain 2-element list from compute_CellToCell():
  #   $edges   - non-zero (edge x LRM) scores (from, to, ligand, receptor, LRM, score [, coords])
  #   $lr.meta - one row per valid LR pair (ligand, receptor, LRM, + extras from LRM.db)
  # Unpacked into $edge.data (Section 9) and $LRM.meta (Section 12) below.
  elr <- compute_CellToCell(
    count.mtx    = score.mtx,
    edge.list    = el.raw,
    LRM.db       = LRM.db,
    ligand.col   = ligand.col,
    receptor.col = receptor.col,
    method       = method,
    n.cores      = n.cores
  )

  ## SECTION 5: Build $cell.list and $gene.list (always from original matrix) ----

  cell.list <- shared.cells
  gene.list <- rownames(count.mtx)

  ## SECTION 6: Build $cell.data (always from original count.mtx) ----

  if (verbose) message("[create_NICHESObject] Building cell.data ...")
  cell.data <- .build_cell_data(count.mtx, normalize.method)

  ## SECTION 7: Build $cell.meta ----

  cell.meta <- as.data.table(meta.data, keep.rownames = "cell")
  setkeyv(cell.meta, "cell")

  ## SECTION 8: Build $edge.list ----

  edge.list.dt <- data.table(
    edge          = paste(el.raw$from, el.raw$to, sep = "|"),
    SendingCell   = el.raw$from,
    ReceivingCell = el.raw$to
  )
  setkeyv(edge.list.dt, "edge")

  ## SECTION 9: Build $edge.data ----

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

  ## SECTION 10: Build $edge.meta ----

  if (verbose) message("[create_NICHESObject] Building edge.meta ...")
  edge.meta <- .build_edge_meta(edge.list.dt, cell.meta)

  ## SECTION 11: Build $gene.meta ----

  if (!is.null(gene.meta)) {
    gene.meta.dt <- as.data.table(gene.meta)
    if (!"gene" %in% colnames(gene.meta.dt))
      stop("[create_NICHESObject] 'gene.meta' must have a 'gene' column.")
    missing.genes <- setdiff(gene.list, gene.meta.dt$gene)
    if (length(missing.genes) > 0L)
      gene.meta.dt <- rbindlist(list(gene.meta.dt,
                                      data.table(gene = missing.genes)),
                                 fill = TRUE)
    setkeyv(gene.meta.dt, "gene")
  } else {
    gene.meta.dt <- data.table(gene = gene.list)
    setkeyv(gene.meta.dt, "gene")
  }

  ## SECTION 12: Build $LRM.list and $LRM.meta ----

  LRM.meta.dt <- copy(elr$lr.meta)
  LRM.list    <- LRM.meta.dt$LRM
  other.cols  <- setdiff(colnames(LRM.meta.dt), c("LRM", "ligand", "receptor"))
  setcolorder(LRM.meta.dt, c("LRM", "ligand", "receptor", other.cols))
  setkeyv(LRM.meta.dt, "LRM")

  ## SECTION 13: Assemble and return NICHESObject ----

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

  attr(obj, "mode") <- mode

  if (verbose) {
    message(sprintf(
      "[create_NICHESObject] Done. mode='%s' | %d cells | %d edges | %d LRMs | %d non-zero signals.",
      mode, length(cell.list), nrow(edge.list.dt),
      length(LRM.list), nrow(edge.data)
    ))
  }

  obj
}
