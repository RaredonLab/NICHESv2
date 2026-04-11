# ==============================================================================
# FILE:    create_NICHESObject.R
# PURPOSE: Public wrapper for NICHESObject construction.
#          Routes to .create_NICHESObject_single() (single-pass, original logic)
#          or .create_NICHESObject_by_sample() (sample-aware split/merge path).
#
# FUNCTIONS (exported):
#   create_NICHESObject()   — main entry point; gains sample.col argument
#
# INTERNAL HELPERS (unexported, in this file):
#   create_NICHESObject_by_sample()  — splits meta.data, loops, merges
#
# NOTE: create_NICHESObject_single() lives in create_NICHESObject_single.R
#       alongside add_cell_meta() and the slot-building helpers. It is the
#       original create_NICHESObject() body with only the function name changed.
# ==============================================================================

#' Create a NICHESObject from Single-Cell Data
#'
#' @description
#' The primary entry point for NICHESv2. Supports two analysis modes:
#'
#' \describe{
#'   \item{\code{"spatial"}}{Edges defined by physical proximity using kNN or
#'     radius in x/y space. Requires \code{x} and \code{y} in \code{meta.data}.}
#'   \item{\code{"sampled"}}{NICHES v1-style type-stratified sampling. Requires
#'     \code{cell.type.col}. Optional PCA-based neighbourhood smoothing via
#'     \code{smooth = TRUE}.}
#' }
#'
#' For pseudobulk (type-averaged) connectivity use
#' \code{\link{compute_PseudoBulk_Connectivity}} — it returns a plain
#' \code{data.table}, not a \code{NICHESObject}.
#'
#' @param count.mtx   Gene x cell matrix. Row names = genes; col names = cells.
#' @param meta.data   \code{data.frame}, row names = cell barcodes.
#' @param LRM.db      \code{"fantom5"}, \code{"omnipath"}, or a data.frame.
#' @param mode        \code{"spatial"} (default) or \code{"sampled"}.
#' @param species     Species for built-in LRM databases. One of \code{"human"}
#'   (default), \code{"mouse"}, \code{"rat"}, or \code{"pig"} (FANTOM5 only).
#'   Ignored when \code{LRM.db} is a \code{data.frame}.
#' @param sample.col  Character or \code{NULL} (default). Name of a column in
#'   \code{meta.data} identifying the biological sample of origin for each cell.
#'   When supplied, cell-cell connectivity is computed \emph{independently within
#'   each sample}, preventing spurious cross-sample edges. Results are then
#'   merged into a single \code{NICHESObject}.
#'   \strong{Barcodes must be globally unique across samples.}
#'   Prepend a sample-specific prefix if needed
#'   (e.g. \code{paste0(sample.id, "_", barcode)}).
#' @param cell.type.col Column in \code{meta.data} with cell type labels.
#'   Required for \code{mode = "sampled"}.
#' @param ligand.col  Default \code{"ligand"}.
#' @param receptor.col Default \code{"receptor"}.
#' @param k           kNN size. Spatial mode only.
#' @param rad         Radius. Spatial mode only.
#' @param seed        Random seed. Sampled mode only. Default \code{42L}.
#' @param max.cells   Per-crossing edge cap. Sampled mode only. Default
#'   \code{NULL}.
#' @param smooth      Logical. Enable PCA neighbourhood smoothing before
#'   scoring. Sampled mode only. Default \code{FALSE}.
#' @param smooth.k    Neighbours for smoothing (k+1 cells averaged). Default
#'   \code{9L}. Sampled + smooth only.
#' @param n.pcs       PCA components for neighbour search. Default \code{30L}.
#'   Capped at \code{min(n.pcs, ncells-1, ngenes-1)}. Sampled + smooth only.
#' @param method      \code{"product"} (default) or \code{"minimum"}.
#' @param normalize.method \code{"prop"} (default) or \code{"none"}.
#' @param gene.meta   Optional data.frame with \code{"gene"} column.
#' @param n.cores     CPU cores. Default 1. Socket cluster — all platforms.
#' @param verbose     Default \code{TRUE}.
#'
#' @return \code{NICHESObject} with \code{attr(obj, "mode")} set. If
#'   \code{sample.col} is provided, \code{attr(obj, "sample.col")} is also set.
#'
#' @seealso \code{\link{merge_NICHESObjects}},
#'   \code{\link{compute_PseudoBulk_Connectivity}},
#'   \code{\link{compute_CellToCell}},
#'   \code{\link{create_Edgelist_Spatial}},
#'   \code{\link{create_Edgelist_SamplingByType}},
#'   \code{\link{load_LRM_database}},
#'   \code{\link{Aggregate_CellToCell_Edges}},
#'   \code{\link{CellToCell_InSitu}}
#'
#' @export
create_NICHESObject <- function(count.mtx,
                                meta.data,
                                LRM.db,
                                mode             = c("spatial", "sampled"),
                                species          = "human",
                                sample.col       = NULL,
                                cell.type.col    = NULL,
                                ligand.col       = "ligand",
                                receptor.col     = "receptor",
                                k                = NULL,
                                rad              = NULL,
                                seed             = 42L,
                                max.cells        = NULL,
                                smooth           = FALSE,
                                smooth.k         = 9L,
                                n.pcs            = 30L,
                                method           = c("product", "minimum"),
                                normalize.method = c("prop", "none"),
                                gene.meta        = NULL,
                                n.cores          = 1L,
                                verbose          = TRUE) {

  # ----------------------------------------------------------
  # Validate sample.col before anything else
  # ----------------------------------------------------------
  if (!is.null(sample.col)) {
    if (!is.character(sample.col) || length(sample.col) != 1L)
      stop("[create_NICHESObject] 'sample.col' must be a single character string or NULL.",
           call. = FALSE)
    if (!is.data.frame(meta.data) || !sample.col %in% colnames(meta.data))
      stop(sprintf(
        "[create_NICHESObject] 'sample.col' \"%s\" not found in colnames(meta.data).\nAvailable columns: %s",
        sample.col, paste(colnames(meta.data), collapse = ", ")
      ), call. = FALSE)
  }

  # ----------------------------------------------------------
  # BRANCH: sample.col provided → split, compute per sample, merge
  # ----------------------------------------------------------
  if (!is.null(sample.col)) {
    return(.create_NICHESObject_by_sample(
      count.mtx        = count.mtx,
      meta.data        = meta.data,
      LRM.db           = LRM.db,
      mode             = mode,
      species          = species,
      sample.col       = sample.col,
      cell.type.col    = cell.type.col,
      ligand.col       = ligand.col,
      receptor.col     = receptor.col,
      k                = k,
      rad              = rad,
      seed             = seed,
      max.cells        = max.cells,
      smooth           = smooth,
      smooth.k         = smooth.k,
      n.pcs            = n.pcs,
      method           = method,
      normalize.method = normalize.method,
      gene.meta        = gene.meta,
      n.cores          = n.cores,
      verbose          = verbose
    ))
  }

  # ----------------------------------------------------------
  # BRANCH: no sample.col → original single-pass logic
  # ----------------------------------------------------------
  .create_NICHESObject_single(
    count.mtx        = count.mtx,
    meta.data        = meta.data,
    LRM.db           = LRM.db,
    mode             = mode,
    species          = species,
    cell.type.col    = cell.type.col,
    ligand.col       = ligand.col,
    receptor.col     = receptor.col,
    k                = k,
    rad              = rad,
    seed             = seed,
    max.cells        = max.cells,
    smooth           = smooth,
    smooth.k         = smooth.k,
    n.pcs            = n.pcs,
    method           = method,
    normalize.method = normalize.method,
    gene.meta        = gene.meta,
    n.cores          = n.cores,
    verbose          = verbose
  )
}
