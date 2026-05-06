#' Export a NICHESObject to TissuePlex edges.parquet
#'
#' Writes a Parquet file that TissuePlex can load directly as its edge layer.
#' Each row represents one (directed edge, LRM) pair.  Only edges present in
#' \code{$edge.data} (i.e. edges with at least one non-zero LR score) are
#' exported.
#'
#' @param object      A \code{NICHESObject}.
#' @param output.path Full path to the output \code{.parquet} file.
#'   Parent directories are created automatically.
#' @param x.col       Name of the x-coordinate column in \code{$cell.meta}
#'   (default \code{"x"}).
#' @param y.col       Name of the y-coordinate column in \code{$cell.meta}
#'   (default \code{"y"}).
#' @param celltype.col Name of the cell-type column in \code{$cell.meta}, or
#'   \code{NULL} to omit cell-type columns (default \code{"celltype"}).
#'
#' @return Invisibly returns \code{output.path}.
#' @export
#' @import data.table
export_to_TissuePlex <- function(object,
                                  output.path,
                                  x.col        = "x",
                                  y.col        = "y",
                                  celltype.col = "celltype") {

  if (!inherits(object, "NICHESObject"))
    stop("'object' must be a NICHESObject")

  if (!requireNamespace("arrow", quietly = TRUE))
    stop("Package 'arrow' is required. Install with: install.packages('arrow')")

  edge.data <- object$edge.data
  edge.meta <- object$edge.meta

  if (is.null(edge.data) || nrow(edge.data) == 0L)
    stop("NICHESObject$edge.data is empty — nothing to export")
  if (is.null(edge.meta) || nrow(edge.meta) == 0L)
    stop("NICHESObject$edge.meta is empty — nothing to export")

  # Derived column names in edge.meta (Sending.<col>, Receiving.<col>)
  send.x <- paste0("Sending.",   x.col)
  recv.x <- paste0("Receiving.", x.col)
  send.y <- paste0("Sending.",   y.col)
  recv.y <- paste0("Receiving.", y.col)

  missing.coord <- setdiff(c(send.x, recv.x, send.y, recv.y), names(edge.meta))
  if (length(missing.coord) > 0L)
    stop(sprintf(
      "Coordinate columns not found in $edge.meta: %s\n  Available: %s",
      paste(missing.coord, collapse = ", "),
      paste(names(edge.meta), collapse = ", ")
    ))

  # Optional cell-type columns
  has.celltype <- FALSE
  send.type    <- NULL
  recv.type    <- NULL
  if (!is.null(celltype.col)) {
    send.type <- paste0("Sending.",   celltype.col)
    recv.type <- paste0("Receiving.", celltype.col)
    if (all(c(send.type, recv.type) %in% names(edge.meta))) {
      has.celltype <- TRUE
    } else {
      message(sprintf(
        "Cell-type column '%s' not found in $edge.meta; sending_type/receiving_type will be NA",
        celltype.col
      ))
    }
  }

  # ── Build LRM catalogue (integer index, stable alphabetical sort) ──────────
  unique.lrms <- sort(unique(edge.data[["LRM"]]))
  lrm.catalogue <- data.table::data.table(
    LRM    = unique.lrms,
    lrm_id = seq_along(unique.lrms)
  )
  data.table::setkeyv(lrm.catalogue, "LRM")

  # ── Pull the subset of edge.meta columns we need ──────────────────────────
  meta.cols <- c("edge", "SendingCell", "ReceivingCell", send.x, send.y, recv.x, recv.y)
  if (has.celltype) meta.cols <- c(meta.cols, send.type, recv.type)
  meta.sub  <- edge.meta[, .SD, .SDcols = meta.cols]

  # ── Join: edge.data + lrm_id + edge.meta ──────────────────────────────────
  # lrm.catalogue[edge.data] = right join → all edge.data rows, adds lrm_id
  out <- lrm.catalogue[edge.data, on = "LRM"]
  # meta.sub[out] = right join → all out rows, adds coordinate/type columns
  out <- meta.sub[out, on = "edge"]

  # ── Derived columns ────────────────────────────────────────────────────────
  out[, is_autocrine := (SendingCell == ReceivingCell)]

  lrm.parts <- data.table::tstrsplit(out[["LRM"]], "|", fixed = TRUE)
  out[, ligand   := lrm.parts[[1L]]]
  out[, receptor := lrm.parts[[2L]]]

  # ── Rename to TissuePlex schema ────────────────────────────────────────────
  data.table::setnames(out,
    old = c("SendingCell",  "ReceivingCell",  "LRM", "score.norm",
            send.x, send.y, recv.x, recv.y),
    new = c("sending_cell", "receiving_cell", "lrm", "score_norm",
            "x1",   "y1",   "x2",   "y2")
  )

  if (has.celltype) {
    data.table::setnames(out, c(send.type, recv.type), c("sending_type", "receiving_type"))
  } else {
    out[, sending_type   := NA_character_]
    out[, receiving_type := NA_character_]
  }

  # ── Select and order final columns ────────────────────────────────────────
  final.cols <- c(
    "edge", "sending_cell", "receiving_cell", "is_autocrine",
    "lrm", "lrm_id", "ligand", "receptor",
    "score", "score_norm",
    "x1", "y1", "x2", "y2",
    "sending_type", "receiving_type"
  )
  out <- out[, .SD, .SDcols = final.cols]

  # ── Write ─────────────────────────────────────────────────────────────────
  dir.create(dirname(normalizePath(output.path, mustWork = FALSE)),
             showWarnings = FALSE, recursive = TRUE)

  arrow::write_parquet(out, output.path)

  message(sprintf(
    "Wrote %d rows | %d directed edges | %d LRMs → %s",
    nrow(out), data.table::uniqueN(out[["edge"]]),
    length(unique.lrms), output.path
  ))

  invisible(output.path)
}
