# NICHESv2-package.R ------------------------------------------------------
# Package-level documentation and .onAttach hook.
# This file produces the ?NICHESv2 help page and the startup message.


# ==============================================================================
#' NICHESv2: Spatial Cell-to-Cell Ligand-Receptor Communication Analysis
#'
#' @description
#' NICHESv2 computes ligand-receptor (LR) signaling scores between pairs of
#' cells in spatial and non-spatial single-cell transcriptomics data. It
#' improves on the original NICHES package with:
#'
#' \itemize{
#'   \item A memory-efficient long-format output object (\code{NICHESObject},
#'     15 slots).
#'   \item Support for multi-gene \strong{ligand and receptor complexes}
#'     (underscore-delimited, e.g. \code{"COL1A1_COL1A2"}).
#'   \item Built-in loaders for the \strong{ConnectomeDB2025} (default),
#'     \strong{FANTOM5}, and \strong{OmniPath} LR databases.
#'   \item Per-edge and per-cell \strong{normalized scores} alongside raw
#'     counts.
#'   \item Parallel LR computation via socket clusters on all platforms.
#'   \item Flexible aggregation (\code{aggregate_NICHESObject}).
#' }
#'
#' @section Typical workflow:
#' \preformatted{
#' library(NICHESv2)
#'
#' # 1. Load / prepare inputs
#' #    count.mtx : gene x cell sparse matrix  (rownames = genes, colnames = cells)
#' #    meta.data : data.frame                  (rownames = cells; spatial mode
#' #                                             needs "x" and "y" columns)
#' #    LRM.db    : "connectomedb2025" | "fantom5" | "omnipath" | user data.frame
#'
#' # 2. Build the NICHESObject (one call does everything)
#' niches.obj <- create_NICHESObject(
#'   count.mtx = count.mtx,
#'   meta.data = meta.data,
#'   LRM.db    = "connectomedb2025",
#'   k         = 6L,           # 6-nearest-neighbor graph
#'   method    = "product"
#' )
#'
#' # 3. Inspect
#' print(niches.obj)
#' summary(niches.obj)
#' dim(niches.obj)            # cells / genes / edges / LRMs
#'
#' # 4. Aggregate neighborhood edge signals
#' niches.obj <- aggregate_NICHESObject(niches.obj, cell.type.col = "celltype")
#'
#' # 5. Subset to a cell population (aggregation rows carry over unchanged)
#' b.cells <- niches.obj$cell.meta[celltype == "B_cell", cell]
#' sub.obj  <- niches.obj[b.cells]
#' }
#'
#' @section NICHESObject slot reference:
#' \tabular{lll}{
#'   \strong{Slot} \tab \strong{Type} \tab \strong{Contents} \cr
#'   \code{$cell.list} \tab character \tab Cell barcodes \cr
#'   \code{$cell.data} \tab data.table \tab (cell, gene, score, score.norm) -- sparse \cr
#'   \code{$cell.meta} \tab data.table \tab Per-cell metadata (x, y, celltype, ...) \cr
#'   \code{$edge.list} \tab data.table \tab (edge, SendingCell, ReceivingCell) \cr
#'   \code{$edge.data} \tab data.table \tab (edge, LRM, score, score.norm) -- sparse \cr
#'   \code{$edge.meta} \tab data.table \tab Per-edge metadata with Sending/Receiving/Joint cols \cr
#'   \code{$neighborhood.cell.list} \tab list \tab Per-cell neighbor barcodes \cr
#'   \code{$neighborhood.edge.list} \tab list \tab Per-cell neighborhood edge IDs \cr
#'   \code{$neighborhood.meta} \tab data.table \tab (cell, n.neighbors, n.edges) \cr
#'   \code{$gene.list} \tab character \tab Gene names \cr
#'   \code{$gene.meta} \tab data.table \tab Per-gene annotation \cr
#'   \code{$LRM.list} \tab character \tab "ligand|receptor" strings \cr
#'   \code{$LRM.meta} \tab data.table \tab Per-LRM annotation (pathway, class, ...) \cr
#'   \code{$segmentation} \tab sf or NULL \tab Optional cell polygons \cr
#'   \code{$aggregations} \tab list \tab Aggregation tables computed on demand \cr
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
#'   \item Score columns: \code{score} (raw), \code{score.norm} (normalized)
#' }
#'
#' @author NICHESv2 Authors
#'
#' @references
#' Ramilowski, J.A. et al. (2015). A draft network of ligand-receptor-mediated
#' multicellular signalling in human. \emph{Nature Communications}, 6, 7866.
#'
#' Turei, D. et al. (2021). Integrated intra- and intercellular signaling
#' knowledge for multicellular omics analysis. \emph{Molecular Systems Biology},
#' 17(3), e9923.
#'
#' @seealso
#' \code{\link{create_NICHESObject}}, \code{\link{aggregate_NICHESObject}},
#' \code{\link{load_LRM_database}}
#'
#' @keywords internal
#' @import data.table
#' @import ggplot2
#' @importFrom methods as
#' @importFrom stats median quantile sd setNames
#' @importFrom utils head
"_PACKAGE"


## .onAttach -- startup message ----

.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    "NICHESv2 ", utils::packageVersion("NICHESv2"),
    " -- Spatial LR communication analysis. ",
    "See ?NICHESv2 or vignette('NICHESv2') to get started."
  )
}
