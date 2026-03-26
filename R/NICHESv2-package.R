# ==============================================================================
# FILE:    NICHESv2-package.R
# PURPOSE: Package-level documentation and .onLoad hook.
#          This file produces the ?NICHESv2 help page and registers the
#          data.table aware-of-package declaration so that := and friends
#          work without triggering R CMD check NOTEs.
# ==============================================================================


# ==============================================================================
#' NICHESv2: Spatial Cell-to-Cell Ligand-Receptor Communication Analysis
#'
#' @description
#' NICHESv2 computes ligand-receptor (LR) signalling scores between spatially
#' neighbouring cells in single-cell spatial transcriptomics data. It improves
#' on the original NICHES package with:
#'
#' \itemize{
#'   \item A memory-efficient long-format output object (\code{NICHESObject},
#'     13 slots).
#'   \item Support for multi-gene \strong{ligand and receptor complexes}
#'     (underscore-delimited, e.g. \code{"COL1A1_COL1A2"}).
#'   \item Built-in loaders for the \strong{FANTOM5} and \strong{OmniPath}
#'     LR databases.
#'   \item Per-edge and per-cell \strong{normalised scores} alongside raw
#'     counts.
#'   \item Parallel LR computation via \code{mclapply} on Unix/macOS.
#'   \item Flexible aggregation (\code{Aggregate_CellToCell_Edges}) and
#'     batch spatial visualisation (\code{CellToCell_InSitu}).
#' }
#'
#' @section Typical workflow:
#' \preformatted{
#' library(NICHESv2)
#'
#' # 1. Load / prepare inputs
#' #    count.mtx : gene x cell sparse matrix  (rownames = genes, colnames = cells)
#' #    meta.data : data.frame                  (rownames = cells, must have "x", "y")
#' #    LRM.db    : "omnipath" | "fantom5" | user data.frame
#'
#' # 2. Build the NICHESObject (one call does everything)
#' niches.obj <- create_NICHESObject(
#'   count.mtx = count.mtx,
#'   meta.data = meta.data,
#'   LRM.db    = "omnipath",
#'   k         = 6L,           # 6-nearest-neighbour graph
#'   method    = "product"
#' )
#'
#' # 3. Inspect
#' print(niches.obj)
#' summary(niches.obj)
#' dim(niches.obj)            # cells / genes / edges / LRMs
#'
#' # 4. Aggregate per-cell sending / receiving scores
#' agg <- Aggregate_CellToCell_Edges(niches.obj, stat = "mean")
#' head(agg$sending)
#'
#' # 5. Visualise a single LRM pair in situ
#' p <- CellToCell_InSitu_Single_LRM(niches.obj, LRM = "MIF|CD74")
#' print(p)
#'
#' # 6. Batch-save PNGs for all LRM pairs
#' CellToCell_InSitu(niches.obj, folder = "niches_plots")
#'
#' # 7. Subset to a cell population
#' b.cells <- niches.obj$cell.meta[celltype == "B_cell", cell]
#' sub.obj  <- niches.obj[b.cells]
#' }
#'
#' @section NICHESObject slot reference:
#' \tabular{lll}{
#'   \strong{Slot} \tab \strong{Type} \tab \strong{Contents} \cr
#'   \code{$cell.list} \tab character \tab Cell barcodes \cr
#'   \code{$cell.data} \tab data.table \tab (cell, gene, score, score.norm) — sparse \cr
#'   \code{$cell.meta} \tab data.table \tab Per-cell metadata (x, y, celltype, …) \cr
#'   \code{$edge.list} \tab data.table \tab (edge, SendingCell, ReceivingCell) \cr
#'   \code{$edge.data} \tab data.table \tab (edge, LRM, score, score.norm) — sparse \cr
#'   \code{$edge.meta} \tab data.table \tab Per-edge metadata with Sending/Receiving/Joint cols \cr
#'   \code{$neighborhood.*} \tab reserved \tab NULL — populated by collaborator functions \cr
#'   \code{$gene.list} \tab character \tab Gene names \cr
#'   \code{$gene.meta} \tab data.table \tab Per-gene annotation \cr
#'   \code{$LRM.list} \tab character \tab "ligand|receptor" strings \cr
#'   \code{$LRM.meta} \tab data.table \tab Per-LRM annotation (pathway, class, …) \cr
#' }
#'
#' @section Key separators:
#' \itemize{
#'   \item Edge IDs: \code{"SendingCell|ReceivingCell"}
#'   \item LRM IDs:  \code{"ligand|receptor"}
#'   \item Gene complexes in database: \code{"gene1_gene2"} (underscore)
#' }
#'
#' @section Key naming conventions:
#' \itemize{
#'   \item Function names use underscores: \code{create_NICHESObject()}
#'   \item Variable / argument names use dots: \code{n.cores}, \code{bg.color}
#'   \item Score columns: \code{score} (raw), \code{score.norm} (normalised)
#' }
#'
#' @author NICHESv2 Authors
#'
#' @references
#' Ramilowski, J.A. et al. (2015). A draft network of ligand-receptor-mediated
#' multicellular signalling in human. \emph{Nature Communications}, 6, 7866.
#'
#' Türei, D. et al. (2021). Integrated intra- and intercellular signaling
#' knowledge for multicellular omics analysis. \emph{Molecular Systems Biology},
#' 17(3), e9923.
#'
#' @seealso
#' \code{\link{create_NICHESObject}}, \code{\link{Aggregate_CellToCell_Edges}},
#' \code{\link{CellToCell_InSitu}}, \code{\link{load_LRM_database}}
#'
#' @keywords internal
"_PACKAGE"


# ==============================================================================
# .onLoad — package startup tasks
# ==============================================================================

.onLoad <- function(libname, pkgname) {
  # Suppress the data.table "no visible binding" NOTE for := and .SD etc.
  # This is the recommended pattern from data.table's own vignette.
  if (requireNamespace("data.table", quietly = TRUE)) {
    data.table::setDTthreads(0L)  # use all available threads by default
  }
  invisible(NULL)
}

.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    "NICHESv2 ", utils::packageVersion("NICHESv2"),
    " — Spatial LR communication analysis. ",
    "See ?NICHESv2 or vignette('NICHESv2') to get started."
  )
}
