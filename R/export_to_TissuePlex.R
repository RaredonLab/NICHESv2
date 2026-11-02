#' Export a NICHESObject to TissuePlex edges.parquet
#'
#' Writes a Parquet file that TissuePlex can load directly as its edge layer.
#' Every edge in \code{$edge.list} is exported.  Edges present in
#' \code{$edge.data} produce one row per (edge, LRM) pair with real scores;
#' edges absent from \code{$edge.data} (no active LR signal) produce a single
#' placeholder row with \code{NA} for all LRM fields.  This ensures TissuePlex
#' can render the full tissue graph even for edges with no scored LR pairs.
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
#' @param n.threads   Number of threads for \code{data.table} operations
#'   (\code{NULL} uses the current \code{getDTthreads()} setting).  On a
#'   MacBook Pro, setting this to the number of performance cores (e.g. 8--12)
#'   significantly speeds up joins and string operations on large datasets.
#'
#' @return Invisibly returns \code{output.path}.
#'
#' @examples
#' if (requireNamespace("arrow", quietly = TRUE)) {
#'   data(NICHESv2_example)
#'   out.file <- tempfile(fileext = ".parquet")
#'   export_to_TissuePlex(NICHESv2_example, output.path = out.file)
#'   unlink(out.file)
#' }
#'
#' @export
export_to_TissuePlex <- function(object,
                                  output.path,
                                  x.col        = "x",
                                  y.col        = "y",
                                  celltype.col = "celltype",
                                  n.threads    = NULL) {

  if (!inherits(object, "NICHESObject"))
    stop("'object' must be a NICHESObject")

  if (!requireNamespace("arrow", quietly = TRUE))
    stop("Package 'arrow' is required. Install with: install.packages('arrow')")

  # Temporarily increase data.table thread count if requested
  if (!is.null(n.threads)) {
    old.threads <- data.table::getDTthreads()
    data.table::setDTthreads(n.threads)
    on.exit(data.table::setDTthreads(old.threads), add = TRUE)
  }

  edge.list <- object$edge.list
  edge.data <- object$edge.data
  edge.meta <- object$edge.meta

  if (is.null(edge.list) || nrow(edge.list) == 0L)
    stop("NICHESObject$edge.list is empty -- nothing to export")
  if (is.null(edge.meta) || nrow(edge.meta) == 0L)
    stop("NICHESObject$edge.meta is empty -- nothing to export")

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

  ## LRM catalogue: unique LRMs with lrm_id + pre-split ligand/receptor ----
  has.edge.data <- !is.null(edge.data) && nrow(edge.data) > 0L

  if (has.edge.data) {
    unique.lrms <- sort(unique(edge.data[["LRM"]]))
    lrm.parts   <- data.table::tstrsplit(unique.lrms, "|", fixed = TRUE)
    lrm.catalogue <- data.table::data.table(
      LRM      = unique.lrms,
      lrm_id   = seq_along(unique.lrms),
      ligand   = lrm.parts[[1L]],
      receptor = lrm.parts[[2L]]
    )
    data.table::setkeyv(lrm.catalogue, "LRM")
  } else {
    unique.lrms   <- character(0L)
    lrm.catalogue <- data.table::data.table(
      LRM      = character(0L),
      lrm_id   = integer(0L),
      ligand   = character(0L),
      receptor = character(0L)
    )
  }

  ## edge.meta: pull only the columns we need ----
  meta.cols <- c("edge", "SendingCell", "ReceivingCell",
                 send.x, send.y, recv.x, recv.y)
  if (has.celltype) meta.cols <- c(meta.cols, send.type, recv.type)
  meta.sub  <- edge.meta[, .SD, .SDcols = meta.cols]
  data.table::setkeyv(meta.sub, "edge")

  ## Combine scored rows (edge.data) with placeholders for unscored edges ----
  if (has.edge.data) {
    scored.edge.dt   <- unique(edge.data[, .(edge)])
    all.edge.dt      <- data.table::data.table(edge = edge.list[["edge"]])
    unscored         <- all.edge.dt[!scored.edge.dt, on = "edge"]
    n.unscored       <- nrow(unscored)

    if (n.unscored > 0L) {
      placeholder <- data.table::data.table(
        edge       = unscored[["edge"]],
        LRM        = NA_character_,
        score      = NA_real_,
        score.norm = NA_real_
      )
      combined <- data.table::rbindlist(
        list(edge.data[, .(edge, LRM, score, score.norm)], placeholder)
      )
    } else {
      combined   <- edge.data[, .(edge, LRM, score, score.norm)]
      n.unscored <- 0L
    }
  } else {
    # No edge.data at all -- every edge gets a placeholder row
    combined <- data.table::data.table(
      edge       = edge.list[["edge"]],
      LRM        = NA_character_,
      score      = NA_real_,
      score.norm = NA_real_
    )
    n.unscored <- nrow(edge.list)
  }

  ## Join lrm_id + pre-split ligand/receptor from catalogue ----
  out <- lrm.catalogue[combined, on = "LRM"]

  ## Join coordinates + SendingCell/ReceivingCell from edge.meta ----
  out <- meta.sub[out, on = "edge"]

  ## Derived column ----
  out[, is_autocrine := (SendingCell == ReceivingCell)]

  ## Rename to TissuePlex schema ----
  data.table::setnames(out,
    old = c("SendingCell",  "ReceivingCell",  "LRM", "score.norm",
            send.x, send.y, recv.x, recv.y),
    new = c("sending_cell", "receiving_cell", "lrm", "score_norm",
            "x1",   "y1",   "x2",   "y2")
  )

  if (has.celltype) {
    data.table::setnames(out, c(send.type, recv.type),
                              c("sending_type", "receiving_type"))
  } else {
    out[, sending_type   := NA_character_]
    out[, receiving_type := NA_character_]
  }

  ## Select and order final columns ----
  final.cols <- c(
    "edge", "sending_cell", "receiving_cell", "is_autocrine",
    "lrm", "lrm_id", "ligand", "receptor",
    "score", "score_norm",
    "x1", "y1", "x2", "y2",
    "sending_type", "receiving_type"
  )
  out <- out[, .SD, .SDcols = final.cols]

  ## Write ----
  dir.create(dirname(normalizePath(output.path, mustWork = FALSE)),
             showWarnings = FALSE, recursive = TRUE)

  arrow::write_parquet(out, output.path)

  n.scored <- nrow(edge.list) - n.unscored
  message(sprintf(
    "Wrote %d rows | %d scored edges (%d LRMs) | %d unscored edges -> %s",
    nrow(out), n.scored, length(unique.lrms), n.unscored, output.path
  ))

  invisible(output.path)
}
