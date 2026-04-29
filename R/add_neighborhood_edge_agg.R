# add_neighborhood_edge_agg.R ------------------------------------------------
# Unified neighborhood edge aggregation for NICHESv2.
#
# For each focal cell and each LRM pair, every edge in the cell's neighborhood
# is classified into exactly one of four mutually exclusive categories relative
# to the focal cell:
#
#   out       -- focal cell is the SendingCell   (non-autocrine)
#   in        -- focal cell is the ReceivingCell (non-autocrine)
#   cross     -- focal cell is neither endpoint  (non-autocrine)
#   autocrine -- SendingCell == ReceivingCell     (any cell in neighborhood)
#
# Note: the autocrine category captures ALL self-loop edges in the
# neighborhood, not only those where the focal cell is the self-looping cell.
# A neighbor cell B with a self-loop B->B contributes to the focal cell's
# autocrine counts.
#
# All three score methods (sum, mean, cmean) and score.inout.ratio are always
# computed. The slot is always fully replaced on each call (no incremental
# updates). This mirrors the Seurat convention: call the function when you
# want the result recomputed.
#
# Neighborhood definition is fixed: "both_in_neighborhood", remove.autocrine =
# FALSE. These are determined at create_NICHESObject() time and stored as
# object attributes. Use add_neighborhoods(overwrite = TRUE) to change them.
#
# Sources: $neighborhood.edge.list, $neighborhood.meta, $edge.list, $edge.data
#
# DEPS: data.table, parallel


## NEIGHBORHOOD EDGE AGGREGATION ----

#' Add Unified Neighborhood Edge Aggregation to a NICHESObject
#'
#' For each focal cell and each LRM pair, aggregates LR communication scores
#' across all edges in the cell's neighborhood. Each neighborhood edge is
#' classified into one of four mutually exclusive categories relative to the
#' focal cell: \code{out} (focal cell sends), \code{in} (focal cell receives),
#' \code{cross} (neither endpoint is the focal cell), or \code{autocrine}
#' (SendingCell == ReceivingCell, any cell in the neighborhood).
#'
#' All three score methods (\code{sum}, \code{mean}, \code{cmean}) and the
#' in/out ratio are always computed. The slot is always fully replaced on each
#' call -- there is no incremental update mode. This follows the same
#' convention as Seurat slot-populating functions: call the function when you
#' want a fresh result.
#'
#' Results are stored as a sparse data.table keyed on \code{(cell, LRM)} at
#' \code{$aggregations$neighborhood.edge.agg}. Only \code{(cell, LRM)} pairs
#' with at least one non-zero-score neighborhood edge are included.
#'
#' @section Four-way edge classification:
#' For a focal cell F and an edge (S -> R) in F's neighborhood:
#' \describe{
#'   \item{\code{autocrine}}{S == R. Self-loop edge involving any cell in the
#'     neighborhood. Checked first, before directional classification.}
#'   \item{\code{out}}{S == F and S != R. Focal cell is the sender.}
#'   \item{\code{in}}{R == F and S != R. Focal cell is the receiver.}
#'   \item{\code{cross}}{S != F and R != F and S != R. Edge between two
#'     neighbors; focal cell is not involved.}
#' }
#'
#' @section Columns produced:
#' Count columns (always produced):
#' \describe{
#'   \item{\code{count.in}, \code{count.out}, \code{count.cross},
#'     \code{count.autocrine}}{Number of non-zero-score edges in each category
#'     for this LRM. Because \code{$edge.data} is sparse, every stored row
#'     represents an active edge.}
#'   \item{\code{count.total}}{Sum of all four category counts.}
#' }
#' Score columns (all three methods always produced, for each category):
#' \describe{
#'   \item{\code{score.\{cat\}.sum}}{Total score across category edges.}
#'   \item{\code{score.\{cat\}.mean}}{Sum divided by total neighborhood edges
#'     (\code{$neighborhood.meta$n.edges}), including zero-score edges.
#'     LRM-agnostic denominator.}
#'   \item{\code{score.\{cat\}.cmean}}{Conditional mean: sum divided by the
#'     category count. \code{NA} when count is zero.}
#'   \item{\code{score.inout.ratio}}{\code{score.in.sum / score.out.sum}.
#'     \code{NA} when \code{score.out.sum == 0}.}
#' }
#'
#' @param niches.obj A \code{NICHESObject} with neighborhood slots populated.
#' @param n.cores    Positive integer. Socket-cluster workers for parallel
#'   computation over LRM pairs. Default \code{1L}.
#' @param verbose    Logical. Emit progress messages. Default \code{TRUE}.
#'
#' @return The input \code{niches.obj} with
#'   \code{$aggregations$neighborhood.edge.agg} populated (always replaced).
#'
#' @seealso \code{\link{add_neighborhoods}},
#'   \code{\link{add_neighborhood_composition}},
#'   \code{\link{aggregate_NICHESObject}}
#'
#' @examples
#' \donttest{
#' data(NICHESv2_example)
#'
#' obj <- add_neighborhood_edge_agg(NICHESv2_example, n.cores = 1L)
#' head(obj$aggregations$neighborhood.edge.agg)
#' }
#'
#' @export
add_neighborhood_edge_agg <- function(
    niches.obj,
    n.cores = 1L,
    verbose = TRUE
) {

  ## Input validation ----------------------------------------------------------
  if (!inherits(niches.obj, "NICHESObject"))
    stop("[add_neighborhood_edge_agg] 'niches.obj' must be a NICHESObject.",
         call. = FALSE)

  n.cores <- as.integer(n.cores)
  if (is.na(n.cores) || n.cores < 1L)
    stop("[add_neighborhood_edge_agg] 'n.cores' must be a positive integer.",
         call. = FALSE)

  ## Neighborhood slot guard ---------------------------------------------------
  if (is.null(niches.obj$neighborhood.edge.list))
    stop(
      "[add_neighborhood_edge_agg] Neighborhood slots are NULL. ",
      "Run add_neighborhoods() first.",
      call. = FALSE
    )

  ## Build (cell, edge, edge.type) flat table ----------------------------------
  # Unroll neighborhood.edge.list: one row per (focal cell, neighborhood edge).
  nbr.lens     <- lengths(niches.obj$neighborhood.edge.list)
  cell.edge.dt <- data.table(
    cell = rep(names(niches.obj$neighborhood.edge.list), nbr.lens),
    edge = unlist(niches.obj$neighborhood.edge.list, use.names = FALSE)
  )

  if (nrow(cell.edge.dt) == 0L)
    stop("[add_neighborhood_edge_agg] No edges found in any neighborhood.",
         call. = FALSE)

  # Attach SendingCell and ReceivingCell via match-based lookup (avoids
  # composite-key join gotcha).
  el     <- niches.obj$edge.list
  el.idx <- match(cell.edge.dt$edge, el$edge)
  cell.edge.dt[, SendingCell   := el$SendingCell[el.idx]]
  cell.edge.dt[, ReceivingCell := el$ReceivingCell[el.idx]]

  # Classify each (focal cell, edge) into one of four categories.
  # Autocrine is checked first: SendingCell == ReceivingCell defines a
  # self-loop regardless of which cell is focal.
  cell.edge.dt[, edge.type := data.table::fcase(
    SendingCell == ReceivingCell,                         "autocrine",
    SendingCell == cell & SendingCell != ReceivingCell,   "out",
    ReceivingCell == cell & SendingCell != ReceivingCell, "in",
    default = "cross"
  )]
  cell.edge.dt[, c("SendingCell", "ReceivingCell") := NULL]

  ## Join edge.data to get per-(cell, edge, LRM) scores -----------------------
  # edge.data is SPARSE: every stored row has score > 0. No filter needed.
  edge.data.sub <- niches.obj$edge.data[, .(edge, LRM, score)]
  flat <- merge(cell.edge.dt, edge.data.sub,
                by = "edge", allow.cartesian = TRUE)
  # flat: edge | cell | edge.type | LRM | score

  ## Precompute total neighborhood edges per cell (LRM-agnostic denominator) --
  # Used as denominator for score.*.mean. Precomputed once, passed to workers.
  total.nbr.edges <- niches.obj$neighborhood.meta[, .(cell, n.edges)]

  if (verbose)
    message(sprintf(
      "[add_neighborhood_edge_agg] Aggregating %d LRM pairs | %d neighborhood edges ...",
      length(niches.obj$LRM.list),
      nrow(cell.edge.dt)
    ))

  ## Worker function -----------------------------------------------------------
  # Processes one LRM's slice. Always computes all counts, all three score
  # methods (sum/mean/cmean) for all four categories, and score.inout.ratio.
  .nbr_edge_worker <- function(lrm.dt, total.nbr.edges) {

    lrm.id   <- lrm.dt$LRM[1L]
    all.cats <- c("in", "out", "cross", "autocrine")

    ## Aggregate counts and sums per (cell, edge.type) -----------------------
    agg <- lrm.dt[,
      .(n.cat = .N, score.cat.sum = sum(score)),
      by = .(cell, edge.type)
    ]

    ## Pivot counts to wide --------------------------------------------------
    cnt.wide <- data.table::dcast(
      agg, cell ~ edge.type, value.var = "n.cat", fill = 0L
    )
    for (cat in setdiff(all.cats, colnames(cnt.wide)))
      cnt.wide[, (cat) := 0L]
    data.table::setnames(cnt.wide,
      old = all.cats, new = paste0("count.", all.cats)
    )
    cnt.wide[, count.total :=
               count.in + count.out + count.cross + count.autocrine]

    ## Pivot score sums to wide ----------------------------------------------
    sum.wide <- data.table::dcast(
      agg, cell ~ edge.type, value.var = "score.cat.sum", fill = 0
    )
    for (cat in setdiff(all.cats, colnames(sum.wide)))
      sum.wide[, (cat) := 0]
    data.table::setnames(sum.wide,
      old = all.cats, new = paste0("score.", all.cats, ".sum")
    )
    sum.wide[, score.total.sum :=
               score.in.sum + score.out.sum +
               score.cross.sum + score.autocrine.sum]

    ## Merge counts and sums -------------------------------------------------
    result <- merge(cnt.wide, sum.wide, by = "cell", all = TRUE)

    ## Mean: LRM-agnostic denominator (total neighborhood edges) -------------
    idx     <- match(result$cell, total.nbr.edges$cell)
    n.edges <- total.nbr.edges$n.edges[idx]
    for (cat in all.cats)
      result[, (paste0("score.", cat, ".mean")) :=
               get(paste0("score.", cat, ".sum")) / n.edges]
    result[, score.total.mean := score.total.sum / n.edges]
    # Guard: cells absent from neighborhood.meta (should not occur)
    for (col in grep("\\.mean$", colnames(result), value = TRUE))
      result[is.na(get(col)), (col) := 0]

    ## Cmean: denominator is the category count (conditional mean) -----------
    for (cat in all.cats)
      result[, (paste0("score.", cat, ".cmean")) :=
               data.table::fifelse(
                 get(paste0("count.", cat)) > 0L,
                 get(paste0("score.", cat, ".sum")) /
                   get(paste0("count.", cat)),
                 NA_real_
               )]
    result[, score.total.cmean := data.table::fifelse(
      count.total > 0L, score.total.sum / count.total, NA_real_
    )]

    ## In/out ratio ----------------------------------------------------------
    result[, score.inout.ratio := data.table::fifelse(
      score.out.sum > 0, score.in.sum / score.out.sum, NA_real_
    )]

    result[, LRM := lrm.id]
    result
  }

  ## Split by LRM and dispatch to workers -------------------------------------
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
      varlist = c(".nbr_edge_worker", "total.nbr.edges"),
      envir   = environment()
    )
    results <- parallel::parLapply(
      cl, lrm.split,
      function(lrm.dt) .nbr_edge_worker(lrm.dt, total.nbr.edges)
    )
  } else {
    results <- lapply(
      lrm.split,
      function(lrm.dt) .nbr_edge_worker(lrm.dt, total.nbr.edges)
    )
  }

  ## Assemble final table -----------------------------------------------------
  agg.dt <- rbindlist(results, fill = TRUE)
  setkeyv(agg.dt, c("cell", "LRM"))

  # Enforce column order: identity | counts | sums | means | cmeans | ratio
  count.order <- c("count.in", "count.out", "count.cross",
                   "count.autocrine", "count.total")
  sum.order   <- c("score.in.sum",   "score.out.sum",   "score.cross.sum",
                   "score.autocrine.sum",   "score.total.sum")
  mean.order  <- c("score.in.mean",  "score.out.mean",  "score.cross.mean",
                   "score.autocrine.mean",  "score.total.mean")
  cmean.order <- c("score.in.cmean", "score.out.cmean", "score.cross.cmean",
                   "score.autocrine.cmean", "score.total.cmean")
  col.order <- c(
    "cell", "LRM",
    intersect(count.order,         colnames(agg.dt)),
    intersect(sum.order,           colnames(agg.dt)),
    intersect(mean.order,          colnames(agg.dt)),
    intersect(cmean.order,         colnames(agg.dt)),
    intersect("score.inout.ratio", colnames(agg.dt))
  )
  setcolorder(agg.dt, col.order)

  ## Always replace the slot --------------------------------------------------
  niches.obj$aggregations$neighborhood.edge.agg <- agg.dt

  if (verbose)
    message(sprintf(
      "[add_neighborhood_edge_agg] Done. %d rows | %d cells | %d LRM pairs.",
      nrow(agg.dt),
      uniqueN(agg.dt$cell),
      uniqueN(agg.dt$LRM)
    ))

  niches.obj
}
