# build_cell_data.R -------------------------------------------------------
# Convert a sparse count matrix to long-format cell.data data.table.
# Internal helper used by .create_NICHESObject_single().

#' Build the long-format cell.data table from a count matrix
#'
#' Internal. Converts a gene x cell sparse count matrix into the long-format
#' \code{data.table} stored in \code{NICHESObject$cell.data}, with columns
#' \code{cell}, \code{gene}, \code{score}, and \code{score.norm}.
#'
#' When \code{normalize.method = "prop"}, \code{score.norm} is the per-cell
#' proportional normalization (each gene divided by the total counts in its
#' cell). Otherwise \code{score.norm} equals \code{score}.
#'
#' Zero-only rows are dropped. The returned table is keyed on
#' \code{c("cell", "gene")}.
#'
#' @param count.mtx        Gene x cell sparse or dense matrix. Rownames are
#'   gene symbols; colnames are cell barcodes.
#' @param normalize.method \code{"prop"} for per-cell proportional
#'   normalization, or any other string to leave \code{score.norm} equal to
#'   \code{score}.
#' @return A keyed \code{data.table} with columns \code{cell}, \code{gene},
#'   \code{score}, \code{score.norm}.
#' @noRd
#' @import data.table
.build_cell_data <- function(count.mtx, normalize.method = "prop") {
  tsm  <- methods::as(count.mtx, "TsparseMatrix")
  rows <- rownames(count.mtx)
  cols <- colnames(count.mtx)
  dt   <- data.table(cell  = cols[tsm@j + 1L],
                     gene  = rows[tsm@i + 1L],
                     score = tsm@x)
  if (normalize.method == "prop") {
    dt[, score.norm := .prop_norm(score), by = cell]
  } else {
    dt[, score.norm := score]
  }
  dt <- dt[score > 0 | score.norm > 0]
  setkeyv(dt, c("cell", "gene"))
  dt
}
