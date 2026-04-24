# add_cell_centric_edge_agg.R ------------------------------------------------
# Compute per-(cell, LRM) aggregation over edges directly linked to each cell.
# For each cell, edges are split into out-edges (cell is SendingCell) and
# in-edges (cell is ReceivingCell). Three aggregation methods are supported:
#   sum   -- total score across all in- or out-edges for this LRM
#   mean  -- sum / total edges in that direction (LRM-agnostic denominator)
#   cmean -- sum / count of non-zero edges only (conditional mean)
#
# Result stored at $aggregations$cell.centric.edge.agg as a sparse data.table
# keyed on (cell, LRM). Only (cell, LRM) pairs with at least one non-zero
# edge in either direction are included.
#
# Sources: $edge.list, $edge.data  -- no neighborhood slot dependency.
#
# DEPS: data.table, parallel


## CELL-CENTRIC EDGE AGGREGATION ----

#' Add Cell-Centric Edge Aggregation to a NICHESObject
#'
#' For each cell and each LRM pair, aggregates scores across all edges where
#' this cell directly participates as a sender (\emph{out-edges}) or receiver
#' (\emph{in-edges}). Results are stored as a sparse data.table keyed on
#' \code{(cell, LRM)}: only pairs where at least one non-zero edge exists in
#' either direction are included.
#'
#' @section Column definitions:
#' \describe{
#'   \item{\code{count.out}, \code{count.in}}{Number of non-zero-score edges
#'     in each direction for this LRM. Always computed regardless of
#'     \code{methods}. Because \code{$edge.data} is stored sparse (only
#'     non-zero scores), this is simply the row count per group.}
#'   \item{\code{score.out.sum}, \code{score.in.sum}}{Total score across all
#'     non-zero edges in each direction. Requested via \code{"sum"}.}
#'   \item{\code{score.out.mean}, \code{score.in.mean}}{Sum divided by the
#'     \emph{total} number of edges in that direction for this cell, including
#'     edges with zero score for this LRM. The denominator is LRM-agnostic
#'     (same for every LRM). This penalises cells with many silent connections.
#'     Requested via \code{"mean"}.}
#'   \item{\code{score.out.cmean}, \code{score.in.cmean}}{Conditional mean:
#'     sum divided by \code{count.out} or \code{count.in} respectively.
#'     \code{NA} when the corresponding count is zero (no non-zero edge exists).
#'     This reports signal strength \emph{given} that signalling occurs.
#'     Requested via \code{"cmean"}.}
#' }
#'
#' @section Incremental updates:
#' Re-calling this function with additional \code{methods} and
#' \code{overwrite = FALSE} appends only the newly requested score columns
#' to the existing table; already-computed columns are not touched.
#'
#' @param niches.obj A \code{NICHESObject}.
#' @param methods    Character vector selecting which aggregation statistics to
#'   compute. Any non-empty subset of \code{c("sum", "mean", "cmean")}.
#'   Default: all three.
#' @param overwrite  Logical. When \code{FALSE} (default), methods whose
#'   columns are already present in the slot are skipped; only genuinely new
#'   columns are computed and joined in. When \code{TRUE}, the slot is fully
#'   recomputed from scratch regardless of existing content.
#' @param n.cores    Positive integer. Number of parallel worker processes.
#'   Parallelises over LRM pairs using socket clusters (\code{parLapply}).
#'   Default \code{1L} (serial).
#' @param verbose    Logical. Emit progress messages. Default \code{TRUE}.
#'
#' @return The input \code{niches.obj} with
#'   \code{$aggregations$cell.centric.edge.agg} populated or updated.
#'
#' @seealso \code{\link{create_NICHESObject}}, \code{\link{add_neighborhoods}}
#'
#' @examples
#' \donttest{
#' data(NICHESv2_example)
#' obj <- add_cell_centric_edge_agg(
#'   NICHESv2_example,
#'   methods = c("sum", "mean", "cmean"),
#'   n.cores = 1L
#' )
#' head(obj$aggregations$cell.centric.edge.agg)
#'
#' # Incremental: add cmean without recomputing sum/mean
#' obj <- add_cell_centric_edge_agg(obj, methods = "cmean", overwrite = FALSE)
#' }
#'
#' @export
add_cell_centric_edge_agg <- function(
    niches.obj,
    methods   = c("sum", "mean", "cmean"),
    overwrite = FALSE,
    n.cores   = 1L,
    verbose   = TRUE
) {

  # -- Input validation -------------------------------------------------------
  if (!inherits(niches.obj, "NICHESObject"))
    stop("[add_cell_centric_edge_agg] 'niches.obj' must be a NICHESObject.")

  methods <- match.arg(methods, c("sum", "mean", "cmean"), several.ok = TRUE)

  if (!is.logical(overwrite) || length(overwrite) != 1L)
    stop("[add_cell_centric_edge_agg] 'overwrite' must be a single logical.")

  n.cores <- as.integer(n.cores)
  if (is.na(n.cores) || n.cores < 1L)
    stop("[add_cell_centric_edge_agg] 'n.cores' must be a positive integer.")

  # -- Resolve which methods genuinely need computing ------------------------
  # Presence of score.in.sum is a reliable proxy: if "sum" was computed before,
  # this column exists. Same logic applies to "mean" and "cmean".
  existing <- niches.obj$aggregations$cell.centric.edge.agg

  if (!overwrite && !is.null(existing)) {
    existing.cols <- colnames(existing)
    methods.todo <- methods[c(
      "sum"   %in% methods & !("score.in.sum"   %in% existing.cols),
      "mean"  %in% methods & !("score.in.mean"  %in% existing.cols),
      "cmean" %in% methods & !("score.in.cmean" %in% existing.cols)
    )]
    if (length(methods.todo) == 0L) {
      if (verbose)
        message("[add_cell_centric_edge_agg] All requested methods already ",
                "present. Use overwrite = TRUE to recompute.")
      return(niches.obj)
    }
    if (verbose && length(methods.todo) < length(methods))
      message("[add_cell_centric_edge_agg] Skipping already-computed methods; ",
              "adding: ", paste(methods.todo, collapse = ", "), ".")
    methods      <- methods.todo
    incremental  <- TRUE
  } else {
    incremental  <- FALSE
  }

  # -- Build flat table: one row per (edge, LRM) with cell endpoints ----------
  # edge.data is SPARSE by construction -- every stored row has score > 0.
  # No additional score > 0 filter is required.
  edge.list <- niches.obj$edge.list                          # edge|SendingCell|ReceivingCell
  edge.data <- niches.obj$edge.data[, .(edge, LRM, score)]  # drop score.norm; not needed

  # Left join: attach SendingCell / ReceivingCell to each scored (edge, LRM)
  flat <- edge.list[edge.data, on = "edge", nomatch = NULL]
  # flat: edge | SendingCell | ReceivingCell | LRM | score

  # -- Precompute total edge counts per cell (LRM-agnostic) ------------------
  # These are the denominators for score.*.mean. Computed once and passed to
  # every worker, avoiding redundant recomputation inside the parallel loop.
  total.out <- edge.list[, .(total.out = .N), by = .(cell = SendingCell)]
  total.in  <- edge.list[, .(total.in  = .N), by = .(cell = ReceivingCell)]

  if (verbose)
    message(sprintf(
      "[add_cell_centric_edge_agg] Aggregating %d LRM pairs across %d edges ...",
      length(niches.obj$LRM.list), nrow(edge.list)
    ))

  # -- Worker function: aggregates one LRM's subset of flat ------------------
  # Defined in local scope; exported to socket workers via clusterExport.
  .cell_centric_worker <- function(lrm.dt, total.out, total.in, methods) {

    lrm.id <- lrm.dt$LRM[1L]

    # Aggregate out-edges (this cell is the sender)
    out.agg <- lrm.dt[,
      .(count.out = .N, score.out.sum = sum(score)),
      by = .(cell = SendingCell)
    ]

    # Aggregate in-edges (this cell is the receiver)
    in.agg <- lrm.dt[,
      .(count.in = .N, score.in.sum = sum(score)),
      by = .(cell = ReceivingCell)
    ]

    # Full outer merge: a cell may send but not receive for this LRM, or vice versa
    result <- merge(out.agg, in.agg, by = "cell", all = TRUE)

    # Fill NAs introduced by the full outer join
    result[is.na(count.out),     count.out     := 0L]
    result[is.na(count.in),      count.in      := 0L]
    result[is.na(score.out.sum), score.out.sum := 0]
    result[is.na(score.in.sum),  score.in.sum  := 0]
    result[, count.out := as.integer(count.out)]
    result[, count.in  := as.integer(count.in)]

    # mean: denominator is total edges per cell in each direction, all LRMs
    if ("mean" %in% methods) {
      idx.out <- match(result$cell, total.out$cell)
      idx.in  <- match(result$cell, total.in$cell)
      result[, score.out.mean := score.out.sum / total.out$total.out[idx.out]]
      result[, score.in.mean  := score.in.sum  / total.in$total.in[idx.in]]
      # Guard against NA (cell appears in result but not in total.* -- should
      # not happen in a well-formed object, but safe to zero-fill)
      result[is.na(score.out.mean), score.out.mean := 0]
      result[is.na(score.in.mean),  score.in.mean  := 0]
    }

    # cmean: denominator is count of non-zero edges for this specific LRM
    if ("cmean" %in% methods) {
      result[, score.out.cmean := fifelse(
        count.out > 0L, score.out.sum / count.out, NA_real_
      )]
      result[, score.in.cmean := fifelse(
        count.in  > 0L, score.in.sum  / count.in,  NA_real_
      )]
    }

    # Drop internal sum columns if "sum" was not requested by the user
    if (!"sum" %in% methods)
      result[, c("score.out.sum", "score.in.sum") := NULL]

    result[, LRM := lrm.id]
    result
  }

  # -- Split by LRM and dispatch to workers ----------------------------------
  lrm.split <- split(flat, by = "LRM", keep.by = TRUE)

  if (n.cores > 1L) {
    cl <- parallel::makeCluster(n.cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterEvalQ(cl, {
      library(data.table)
      data.table::setDTthreads(1L)
    })
    parallel::clusterExport(
      cl,
      varlist = c(".cell_centric_worker", "total.out", "total.in", "methods"),
      envir   = environment()
    )
    results <- parallel::parLapply(
      cl, lrm.split,
      function(lrm.dt) .cell_centric_worker(lrm.dt, total.out, total.in, methods)
    )
  } else {
    results <- lapply(
      lrm.split,
      function(lrm.dt) .cell_centric_worker(lrm.dt, total.out, total.in, methods)
    )
  }

  # -- Assemble final table --------------------------------------------------
  agg.dt <- rbindlist(results)
  setkeyv(agg.dt, c("cell", "LRM"))

  # Enforce column order: identity cols, counts, then score cols in method order
  score.col.order <- c(
    "score.out.sum",  "score.in.sum",
    "score.out.mean", "score.in.mean",
    "score.out.cmean","score.in.cmean"
  )
  present.score.cols <- intersect(score.col.order, colnames(agg.dt))
  setcolorder(agg.dt, c("cell", "LRM", "count.out", "count.in", present.score.cols))

  # -- Store result: full replace or incremental column join -----------------
  if (incremental) {
    # Join only the newly computed score columns onto the existing keyed table.
    # count.out / count.in are recomputed but not re-stored (existing values
    # are authoritative when overwrite = FALSE).
    new.cols <- setdiff(colnames(agg.dt), c("cell", "LRM", "count.out", "count.in"))
    join.dt  <- agg.dt[, c("cell", "LRM", new.cols), with = FALSE]
    existing <- existing[join.dt, on = c("cell", "LRM")]
    niches.obj$aggregations$cell.centric.edge.agg <- existing
  } else {
    niches.obj$aggregations$cell.centric.edge.agg <- agg.dt
  }

  if (verbose)
    message(sprintf(
      "[add_cell_centric_edge_agg] Done. %d rows | %d cells | %d LRM pairs.",
      nrow(niches.obj$aggregations$cell.centric.edge.agg),
      uniqueN(niches.obj$aggregations$cell.centric.edge.agg$cell),
      uniqueN(niches.obj$aggregations$cell.centric.edge.agg$LRM)
    ))

  niches.obj
}
