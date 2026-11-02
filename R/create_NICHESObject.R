# create_NICHESObject.R ---------------------------------------------------
# Public wrapper for NICHESObject construction.
# Routes to .create_NICHESObject_single() (single-pass, original logic)
# or .create_NICHESObject_by_sample() (sample-aware split/merge path).
#
# FUNCTIONS (exported):
#   create_NICHESObject()   -- main entry point; gains sample.col argument
#
# INTERNAL HELPERS (unexported, in this file):
#   create_NICHESObject_by_sample()  -- splits meta.data, loops, merges
#
# NOTE: create_NICHESObject_single() lives in create_NICHESObject_single.R
#       alongside add_cell_meta() and the slot-building helpers. It is the
#       original create_NICHESObject() body with only the function name changed.

#' Create a NICHESObject from Single-Cell Data
#'
#' @description
#' The primary entry point for NICHESv2. Supports two analysis modes:
#'
#' \describe{
#'   \item{\code{"spatial"}}{Edges defined by physical proximity using kNN or
#'     radius in x/y space. Requires \code{x} and \code{y} in \code{meta.data}.}
#'   \item{\code{"nonspatial"}}{NICHES v1-style type-stratified sampling. Requires
#'     \code{cell.type.col}.}
#' }
#'
#' @section Pipeline overview:
#' Internally, \code{create_NICHESObject()} chains the following steps:
#' \enumerate{
#'   \item Load the LR database (\code{\link{load_LRM_database}}) when
#'         \code{LRM.db} is a character string; otherwise the supplied
#'         \code{data.frame} is used directly.
#'   \item Build the edge list -- \code{\link{create_Edgelist_Spatial}} for
#'         \code{mode = "spatial"} or \code{\link{create_Edgelist_SamplingByType}}
#'         for \code{mode = "nonspatial"}.
#'   \item Score every edge x LRM combination via
#'         \code{\link{compute_CellToCell}}, sparsifying as it goes.
#'   \item Assemble the 15-slot \code{NICHESObject} (cell, edge, gene, LRM,
#'         segmentation, and aggregation slots).
#'   \item Extract neighborhoods via \code{\link{add_neighborhoods}} using
#'         \code{edge.filter.mode = "both_in_neighborhood"} and
#'         \code{remove.autocrine = FALSE} (fixed by design).
#' }
#'
#' Each step is also exported individually for users who need finer control.
#'
#' @section Sample-aware processing:
#' When \code{sample.col} is supplied, \code{meta.data} is split by that
#' column and the full pipeline runs \emph{independently within each sample}.
#' Per-sample objects are then combined via \code{\link{merge_NICHESObjects}}.
#' This prevents physically meaningless cross-sample edges (a kNN search on
#' coordinates pooled across slides will happily connect cells that were
#' never on the same tissue).
#'
#' \strong{Barcodes must be globally unique across samples.} A duplicate
#' barcode is a hard error during the merge step. Prepend a sample-specific
#' prefix before construction (e.g.\ \code{paste0(sample.id, "_", barcode)},
#' which is what \code{Seurat::RenameCells()} does automatically).
#'
#' Outer iteration over samples is always sequential; the inner pipeline
#' uses \code{n.cores} for parallelism, so socket clusters are never nested.
#'
#' @param count.mtx   Gene x cell matrix (sparse \code{dgCMatrix} or dense).
#'   Row names = gene symbols; column names = cell barcodes. Both name sets
#'   must match \code{rownames(meta.data)} and the columns of
#'   \code{LRM.db} respectively.
#' @param meta.data   \code{data.frame} with row names = cell barcodes. Must
#'   contain \code{x} and \code{y} columns when \code{mode = "spatial"}, and
#'   the column named by \code{cell.type.col} when \code{mode = "nonspatial"}.
#' @param LRM.db      Either a character string naming a built-in database
#'   (\code{"connectomedb2025"}, \code{"fantom5"} or \code{"omnipath"}) or a
#'   \code{data.frame} with at
#'   least the columns named by \code{ligand.col} and \code{receptor.col}.
#'   See \code{\link{load_LRM_database}} for details and supported species.
#' @param mode        Character. \code{"spatial"} (default, kNN/radius edges)
#'   or \code{"nonspatial"} (type-stratified random sampling).
#' @param species     Species for built-in LRM databases. One of
#'   \code{"human"} (default), \code{"mouse"}, \code{"rat"}, or \code{"pig"}
#'   (\code{"pig"} is not available from OmniPath). Silently ignored when
#'   \code{LRM.db} is a \code{data.frame}.
#' @param sample.col  Character or \code{NULL} (default). Name of a column
#'   in \code{meta.data} identifying the biological sample of origin for
#'   each cell. When supplied, the pipeline runs independently within each
#'   sample and the per-sample objects are merged. See the Sample-aware
#'   processing section above for the barcode-uniqueness requirement.
#' @param cell.type.col Character or \code{NULL}. Column in \code{meta.data}
#'   holding cell type labels. \strong{Required for \code{mode = "nonspatial"}};
#'   ignored otherwise.
#' @param ligand.col  Character. Column name for ligands in \code{LRM.db}.
#'   Default \code{"ligand"}.
#' @param receptor.col Character. Column name for receptors in \code{LRM.db}.
#'   Default \code{"receptor"}.
#' @param k           Integer or \code{NULL}. Number of nearest neighbors
#'   for kNN edges. \strong{Spatial mode only; supply exactly one of
#'   \code{k} or \code{rad}.}
#' @param rad         Numeric or \code{NULL}. Neighborhood radius for
#'   fixed-radius edges (same units as the \code{x}/\code{y} columns).
#'   \strong{Spatial mode only; supply exactly one of \code{k} or \code{rad}.}
#'   Use \code{\link{Explore_Radius_Neighborhood}} to choose a sensible value.
#' @param seed        Integer. Random seed for reproducible type-stratified
#'   sampling. Default \code{42L}. Nonspatial mode only.
#' @param max.cells   Integer or \code{NULL}. Per cell-type-pair edge cap.
#'   Default \code{NULL} (no cap). Nonspatial mode only.
#' @param method      Per-edge LR scoring method. \code{"product"} (default)
#'   multiplies sender ligand expression by receiver receptor expression;
#'   \code{"minimum"} takes the elementwise minimum. Both handle gene
#'   complexes (underscore-joined symbols).
#' @param normalize.method Edge-level normalization for \code{score.norm}
#'   in \code{$edge.data}. \code{"prop"} (default) scales each edge's
#'   per-LRM scores to sum to 1 across LRMs; \code{"none"} leaves
#'   \code{score.norm} equal to \code{score}.
#' @param gene.meta   Optional \code{data.frame} with a \code{"gene"} column
#'   (and any number of additional annotation columns) joined onto the
#'   \code{$gene.meta} slot. Default \code{NULL}.
#' @param n.cores     Integer. Number of CPU cores for the inner scoring
#'   pipeline. Default \code{1L}. Always uses socket clusters
#'   (\code{parallel::makeCluster}) for cross-platform compatibility.
#' @param verbose     Logical. Print progress messages. Default \code{TRUE}.
#'
#' @return A \code{NICHESObject} with attribute \code{"mode"} set to either
#'   \code{"spatial"} or \code{"nonspatial"}. When \code{sample.col} is
#'   provided, attribute \code{"sample.col"} is also set to the column name.
#'   Slots \code{$cell.list}, \code{$cell.data}, \code{$cell.meta},
#'   \code{$edge.list}, \code{$edge.data}, \code{$edge.meta},
#'   \code{$gene.list}, \code{$gene.meta}, \code{$LRM.list}, and
#'   \code{$LRM.meta} are always populated. The three neighborhood slots are
#'   always populated and \code{attr(obj, "edge.filter.mode")} and
#'   \code{attr(obj, "remove.autocrine")} are set on the object.
#'   \code{$segmentation} always remains
#'   \code{NULL} until \code{\link{add_segmentation}} is called.
#'
#' @seealso
#'   \code{\link{NICHESv2_example}} for an example object produced by this
#'     function;
#'   \code{\link{merge_NICHESObjects}} for combining objects produced
#'     independently;
#'   \code{\link{compute_CellToCell}} for the underlying scoring engine;
#'   \code{\link{create_Edgelist_Spatial}},
#'   \code{\link{create_Edgelist_SamplingByType}} for the edge-list
#'     constructors invoked internally;
#'   \code{\link{load_LRM_database}} for the LR-database loader;
#'   \code{\link{aggregate_NICHESObject}} for downstream summarization.
#'
#' @examples
#' \donttest{
#' data(NICHESv2_inputs)
#'
#' # ---- Basic spatial mode (kNN edges) ----
#' obj <- create_NICHESObject(
#'   count.mtx    = NICHESv2_inputs$count.mtx,
#'   meta.data    = NICHESv2_inputs$meta.data,
#'   LRM.db       = NICHESv2_inputs$LRM.db,
#'   mode         = "spatial",
#'   ligand.col   = "ligand",
#'   receptor.col = "receptor",
#'   k            = 5L,
#'   n.cores      = 1L,
#'   verbose      = FALSE
#' )
#' print(obj)
#' attr(obj, "mode")           # "spatial"
#'
#' # ---- Nonspatial mode (NICHES v1-style type-stratified sampling) ----
#' obj.sampled <- create_NICHESObject(
#'   count.mtx     = NICHESv2_inputs$count.mtx,
#'   meta.data     = NICHESv2_inputs$meta.data,
#'   LRM.db        = NICHESv2_inputs$LRM.db,
#'   mode          = "nonspatial",
#'   cell.type.col = "celltype",
#'   ligand.col    = "ligand",
#'   receptor.col  = "receptor",
#'   n.cores       = 1L,
#'   verbose       = FALSE
#' )
#' attr(obj.sampled, "mode")   # "nonspatial"
#'
#' # ---- Sample-aware processing (sample.col) ----
#' # Synthesise a 'sample' column to demonstrate the split-and-merge path.
#' # Real workflows: ensure barcodes are globally unique across samples.
#' md <- NICHESv2_inputs$meta.data
#' md$sample <- rep(c("S1", "S2"), each = nrow(md) / 2L)
#'
#' obj.split <- create_NICHESObject(
#'   count.mtx    = NICHESv2_inputs$count.mtx,
#'   meta.data    = md,
#'   LRM.db       = NICHESv2_inputs$LRM.db,
#'   mode         = "spatial",
#'   sample.col   = "sample",
#'   ligand.col   = "ligand",
#'   receptor.col = "receptor",
#'   k            = 5L,
#'   n.cores      = 1L,
#'   verbose      = FALSE
#' )
#' attr(obj.split, "sample.col")   # "sample"
#' }
#'
#' @export
create_NICHESObject <- function(count.mtx,
                                meta.data,
                                LRM.db,
                                mode             = c("spatial", "nonspatial"),
                                species          = "human",
                                sample.col       = NULL,
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
  ## Validate sample.col before anything else ----
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

  ## BRANCH: sample.col provided -> split, compute per sample, merge ----
  if (!is.null(sample.col)) {
    obj <- .create_NICHESObject_by_sample(
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
      method           = method,
      normalize.method = normalize.method,
      gene.meta        = gene.meta,
      n.cores          = n.cores,
      verbose          = verbose
    )
  } else {
    ## BRANCH: no sample.col -> original single-pass logic ----
    obj <- .create_NICHESObject_single(
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
      method           = method,
      normalize.method = normalize.method,
      gene.meta        = gene.meta,
      n.cores          = n.cores,
      verbose          = verbose
    )
  }
  
  ## Neighborhood extraction ----
  # Always uses "both_in_neighborhood" with remove.autocrine = FALSE.
  # These are fixed by design -- use add_neighborhoods() directly to
  # experiment with other parameter combinations.
  obj <- add_neighborhoods(
    niches.obj       = obj,
    edge.filter.mode = "both_in_neighborhood",
    remove.autocrine = FALSE,
    overwrite        = FALSE,
    n.cores          = n.cores,
    verbose          = verbose
  )
  
  obj
}