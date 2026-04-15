# compute_PseudoBulk_Connectivity.R ---------------------------------------
# Compute pseudobulk cell-type-to-cell-type ligand-receptor
# connectivity. Returns a lightweight data.table (not a NICHESObject)
# summarising average signal from each cell type to each other type.
#
# FUNCTIONS:
#   compute_PseudoBulk_Connectivity()  -- main function
#
# DEPS: data.table (>= 1.14), Matrix (>= 1.5), parallel


# ==============================================================================
#' Compute Pseudobulk Cell-Type Connectivity
#'
#' @description
#' Computes ligand-receptor signal scores between all ordered cell type pairs
#' using per-type mean expression profiles (pseudobulk). Covers all T^2 directed
#' type-to-type crossings including autocrine (A -> A).
#'
#' Unlike \code{create_NICHESObject()}, this function does \strong{not} build
#' a \code{NICHESObject}. It returns a compact \code{data.table} directly
#' suitable for heatmap or network visualisation of type-level communication.
#' This is intentional -- the pseudobulk output is type-level rather than
#' cell-level and does not benefit from the full NICHESObject slot structure.
#'
#' @section Output schema:
#' One row per (type-pair x non-zero LRM):
#' \describe{
#'   \item{\code{SendingType}}{Sending cell type label.}
#'   \item{\code{ReceivingType}}{Receiving cell type label.}
#'   \item{\code{LRM}}{\code{"ligand|receptor"} string.}
#'   \item{\code{ligand}}{Ligand gene or complex.}
#'   \item{\code{receptor}}{Receptor gene or complex.}
#'   \item{\code{score}}{Raw LR signal (product or minimum of mean expressions).}
#'   \item{\code{score.norm}}{Score normalised proportionally per type-pair edge
#'     (sums to 1.0 across all LRMs for a given sending-receiving pair), or equal
#'     to \code{score} when \code{normalize.method = "none"}.}
#' }
#'
#' @param count.mtx      Gene x cell expression matrix (sparse or dense).
#'   Row names = gene symbols; column names = cell barcodes.
#' @param meta.data      \code{data.frame} with row names = cell barcodes.
#'   Must contain the column named by \code{cell.type.col}.
#' @param LRM.db         LR database: \code{"fantom5"}, \code{"omnipath"}, or
#'   a \code{data.frame} with at least \code{ligand.col} and
#'   \code{receptor.col}.
#' @param cell.type.col  \code{character(1)}. Column in \code{meta.data}
#'   holding cell type labels.
#' @param ligand.col     Column name for ligands in \code{LRM.db}.
#'   Default \code{"ligand"}.
#' @param receptor.col   Column name for receptors in \code{LRM.db}.
#'   Default \code{"receptor"}.
#' @param method         Scoring method: \code{"product"} (default) or
#'   \code{"minimum"}.
#' @param normalize.method Normalisation for \code{score.norm}:
#'   \code{"prop"} (default, proportional per type-pair) or \code{"none"}.
#' @param n.cores        CPU cores. Default 1. Socket cluster -- works on all
#'   platforms.
#' @param verbose        Print progress messages. Default \code{TRUE}.
#'
#' @return A \code{data.table} keyed on \code{(SendingType, ReceivingType, LRM)}
#'   with columns as described above. Only non-zero score rows are included
#'   (implicit sparsity).
#'
#' @seealso \code{\link{create_NICHESObject}},
#'   \code{\link{create_Edgelist_PseudoBulk}},
#'   \code{\link{compute_CellToCell}}
#'
#' @examples
#' \dontrun{
#' pb <- compute_PseudoBulk_Connectivity(
#'   count.mtx     = count.mtx,
#'   meta.data     = meta.data,
#'   LRM.db        = "omnipath",
#'   cell.type.col = "cell.type"
#' )
#' # Quick heatmap of total signal per type pair
#' agg <- pb[, .(total.score = sum(score)), by = .(SendingType, ReceivingType)]
#' }
#'
#' @export
compute_PseudoBulk_Connectivity <- function(count.mtx,
                                             meta.data,
                                             LRM.db,
                                             cell.type.col,
                                             ligand.col       = "ligand",
                                             receptor.col     = "receptor",
                                             method           = c("product", "minimum"),
                                             normalize.method = c("prop", "none"),
                                             n.cores          = 1L,
                                             verbose          = TRUE) {

  method           <- match.arg(method)
  normalize.method <- match.arg(normalize.method)
  n.cores          <- as.integer(n.cores)

  ## SECTION 1: Input validation ----

  if (!is.matrix(count.mtx) && !inherits(count.mtx, "sparseMatrix"))
    stop("[compute_PseudoBulk_Connectivity] 'count.mtx' must be a matrix or sparseMatrix.")
  if (is.null(rownames(count.mtx)) || is.null(colnames(count.mtx)))
    stop("[compute_PseudoBulk_Connectivity] 'count.mtx' must have row names (genes) and column names (cells).")
  if (!is.data.frame(meta.data))
    stop("[compute_PseudoBulk_Connectivity] 'meta.data' must be a data.frame.")
  if (is.null(rownames(meta.data)))
    stop("[compute_PseudoBulk_Connectivity] 'meta.data' must have row names (cell barcodes).")
  if (!is.character(cell.type.col) || length(cell.type.col) != 1L)
    stop("[compute_PseudoBulk_Connectivity] 'cell.type.col' must be a single character string.")
  if (!cell.type.col %in% colnames(meta.data))
    stop(sprintf(
      "[compute_PseudoBulk_Connectivity] 'cell.type.col' ('%s') not found in 'meta.data'.",
      cell.type.col
    ))

  ## SECTION 2: Cell alignment ----

  shared.cells <- intersect(rownames(meta.data), colnames(count.mtx))
  if (length(shared.cells) == 0L)
    stop("[compute_PseudoBulk_Connectivity] No shared cells between 'meta.data' and 'count.mtx'.")

  if (length(shared.cells) < ncol(count.mtx) && verbose)
    message(sprintf(
      "[compute_PseudoBulk_Connectivity] %d cell(s) in count.mtx not in meta.data -- excluded.",
      ncol(count.mtx) - length(shared.cells)
    ))
  if (length(shared.cells) < nrow(meta.data) && verbose)
    message(sprintf(
      "[compute_PseudoBulk_Connectivity] %d cell(s) in meta.data not in count.mtx -- excluded.",
      nrow(meta.data) - length(shared.cells)
    ))

  count.mtx <- count.mtx[, shared.cells, drop = FALSE]
  meta.data  <- meta.data[shared.cells, , drop = FALSE]

  ## SECTION 3: Build pseudobulk edge list and mean-expression matrix ----

  if (verbose) message("[compute_PseudoBulk_Connectivity] Building type-pair edge list ...")

  el.raw <- create_Edgelist_PseudoBulk(meta.data     = meta.data,
                                        cell.type.col = cell.type.col)

  # Per-type mean expression: genes x cell_types dense matrix.
  # Column order follows the sorted type labels in attr(el.raw, "cell.names")
  # so colnames() matches the from/to values in el.raw exactly.
  cell.types.ordered <- attr(el.raw, "cell.names")

  if (verbose) message("[compute_PseudoBulk_Connectivity] Computing per-type mean expression ...")

  pb.mtx <- do.call(cbind, lapply(cell.types.ordered, function(ct) {
    ct.barcodes <- rownames(meta.data)[meta.data[[cell.type.col]] == ct]
    Matrix::rowMeans(count.mtx[, ct.barcodes, drop = FALSE])
  }))
  colnames(pb.mtx) <- cell.types.ordered

  ## SECTION 4: Score LRM pairs across all type-pair edges ----

  if (verbose) message("[compute_PseudoBulk_Connectivity] Scoring LRM pairs ...")

  elr <- compute_CellToCell(
    count.mtx    = pb.mtx,
    edge.list    = el.raw,
    LRM.db       = LRM.db,
    ligand.col   = ligand.col,
    receptor.col = receptor.col,
    method       = method,
    n.cores      = n.cores
  )

  ## SECTION 5: Assemble and normalise output data.table ----

  result <- elr$edges[, .(
    SendingType   = from,
    ReceivingType = to,
    LRM, ligand, receptor, score
  )]

  if (normalize.method == "prop") {
    result[, score.norm := .prop_norm(score), by = .(SendingType, ReceivingType)]
  } else {
    result[, score.norm := score]
  }

  setcolorder(result, c("SendingType", "ReceivingType", "LRM",
                         "ligand", "receptor", "score", "score.norm"))
  setkeyv(result, c("SendingType", "ReceivingType", "LRM"))

  if (verbose) {
    n.types <- length(cell.types.ordered)
    message(sprintf(
      "[compute_PseudoBulk_Connectivity] Done. %d type pairs | %d valid LRMs | %d non-zero signals.",
      n.types^2, length(elr$lr.meta$LRM), nrow(result)
    ))
  }

  result
}
