# compute_LRM_Score.R -----------------------------------------------------
# Compute the per-edge ligand-receptor score vector for a single LRM.
# Internal helper used by compute_CellToCell(); also exported to
# socket-cluster workers via clusterExport(asNamespace("NICHESv2")).

#' Compute the per-edge score vector for one ligand-receptor pair
#'
#' Internal. Computes a score for one ligand-receptor pair across a vector of
#' sending/receiving cell pairs. Supports gene complexes via underscore
#' notation (e.g. \code{"CXCL12_PF4"}).
#'
#' For each edge, the ligand expression is read from the \emph{sending} cell
#' and the receptor expression from the \emph{receiving} cell. When either
#' the ligand or receptor is a complex, every component is extracted, and all
#' components are combined via the chosen \code{method}.
#'
#' @param lig.str       Character. Ligand gene or complex (e.g. \code{"VEGFA"}
#'   or \code{"CXCL12_PF4"}).
#' @param rec.str       Character. Receptor gene or complex.
#' @param count.mtx     Gene x cell matrix (sparse or dense). Rownames are
#'   gene symbols; colnames are cell barcodes.
#' @param sending.cells Character vector of length \code{n.edges}.
#' @param receiv.cells  Character vector of length \code{n.edges}.
#' @param method        \code{"product"} or \code{"minimum"} -- how to combine
#'   ligand and receptor components.
#' @return Numeric vector of length \code{n.edges} (may contain zeros).
#' @noRd
.compute_LRM_Score <- function(lig.str, rec.str,
                               count.mtx,
                               sending.cells, receiv.cells,
                               method = c("product", "minimum")) {
  method <- match.arg(method)

  lig.parts <- .parse_complex(lig.str)
  rec.parts <- .parse_complex(rec.str)

  lig.vecs <- lapply(lig.parts, function(g)
    as.numeric(count.mtx[g, sending.cells, drop = TRUE]))

  rec.vecs <- lapply(rec.parts, function(g)
    as.numeric(count.mtx[g, receiv.cells, drop = TRUE]))

  all.vecs <- c(lig.vecs, rec.vecs)

  if (method == "product") {
    Reduce("*", all.vecs)
  } else {
    Reduce(pmin, all.vecs)
  }
}
