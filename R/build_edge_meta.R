# build_edge_meta.R -------------------------------------------------------
# Join sending/receiving cell metadata onto an edge list.
# Internal helper used by .create_NICHESObject_single() and
# add_cell_meta().

#' Join cell metadata onto an edge list to produce edge.meta
#'
#' Internal. Takes an edge list (with columns \code{edge}, \code{SendingCell},
#' \code{ReceivingCell}) and a keyed cell metadata table, and produces the
#' \code{edge.meta} table by joining sending and receiving cell metadata onto
#' each edge.
#'
#' For every metadata column in \code{cell.meta} (other than \code{cell}),
#' two columns are added: \code{Sending.<col>} and \code{Receiving.<col>}.
#' For categorical columns (character or factor), an additional
#' \code{Joint.<col>} column is added with the format
#' \code{"<sending>|<receiving>"}.
#'
#' The returned table is keyed on \code{edge}.
#'
#' @param edge.list.dt A \code{data.table} with at least \code{edge},
#'   \code{SendingCell}, and \code{ReceivingCell}.
#' @param cell.meta    A \code{data.table} keyed on \code{cell}.
#' @return A keyed \code{data.table} containing the edge list plus
#'   \code{Sending.*}, \code{Receiving.*}, and \code{Joint.*} columns.
#' @noRd
#' @import data.table
.build_edge_meta <- function(edge.list.dt, cell.meta) {
  meta.cols <- setdiff(colnames(cell.meta), "cell")
  edge.meta <- copy(edge.list.dt)
  if (length(meta.cols) == 0L) return(edge.meta)
  send.meta <- cell.meta[edge.meta[, .(cell = SendingCell)], on = "cell"]
  for (col in meta.cols)
    edge.meta[[paste0("Sending.", col)]] <- send.meta[[col]]
  recv.meta <- cell.meta[edge.meta[, .(cell = ReceivingCell)], on = "cell"]
  for (col in meta.cols)
    edge.meta[[paste0("Receiving.", col)]] <- recv.meta[[col]]
  cat.cols <- .categorical_cols(cell.meta, exclude = "cell")
  for (col in cat.cols) {
    edge.meta[[paste0("Joint.", col)]] <- paste(
      edge.meta[[paste0("Sending.",   col)]],
      edge.meta[[paste0("Receiving.", col)]], sep = "|"
    )
  }
  setkeyv(edge.meta, "edge")
  edge.meta
}
