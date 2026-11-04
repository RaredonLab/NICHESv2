# add_neighborhood_edge_agg.R ------------------------------------------------
# Unified neighborhood edge aggregation for NICHESv2.
#
# For each focal cell and each LRM pair, every edge in the cell's neighborhood
# is classified into exactly one of five mutually exclusive categories relative
# to the focal cell:
#
#   out             -- focal cell is the SendingCell   (non-autocrine)
#   in              -- focal cell is the ReceivingCell (non-autocrine)
#   cross           -- focal cell is neither endpoint  (non-autocrine)
#   self.autocrine  -- SendingCell == ReceivingCell == focal cell (own self-loop)
#   other.autocrine -- SendingCell == ReceivingCell != focal cell (neighbor self-loop)
#
# self.autocrine and other.autocrine together replace the former single
# "autocrine" category. count.self.autocrine is always 0 or 1 per (cell, LRM)
# because a cell can have at most one self-loop edge.
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
# DEPS: data.table, Matrix


## NEIGHBORHOOD EDGE AGGREGATION ----

#' Add Unified Neighborhood Edge Aggregation to a NICHESObject
#'
#' For each focal cell and each LRM pair, aggregates LR communication scores
#' across all edges in the cell's neighborhood. Each neighborhood edge is
#' classified into one of five mutually exclusive categories relative to the
#' focal cell: \code{out} (focal cell sends), \code{in} (focal cell receives),
#' \code{cross} (neither endpoint is the focal cell), \code{self.autocrine}
#' (focal cell's own self-loop: SendingCell == ReceivingCell == focal cell), or
#' \code{other.autocrine} (a neighbor cell's self-loop: SendingCell ==
#' ReceivingCell != focal cell).
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
#' @section Five-way edge classification:
#' For a focal cell F and an edge (S -> R) in F's neighborhood:
#' \describe{
#'   \item{\code{self.autocrine}}{S == R == F. F's own self-loop edge. Checked
#'     first. \code{count.self.autocrine} is always 0 or 1 per (cell, LRM).}
#'   \item{\code{other.autocrine}}{S == R != F. A neighbor cell's self-loop
#'     edge; F is not the self-looping cell.}
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
#'     \code{count.self.autocrine}, \code{count.other.autocrine}}{Number of
#'     non-zero-score edges in each category for this LRM. Because
#'     \code{$edge.data} is sparse, every stored row represents an active edge.
#'     \code{count.self.autocrine} is always 0 or 1 per (cell, LRM).}
#'   \item{\code{count.total}}{Sum of all five category counts.}
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
#' @param n.cores    Positive integer. Retained for API compatibility; not used
#'   in the current sparse matrix implementation. Default \code{1L}.
#' @param chunk.size Positive integer. Number of LRM pairs processed per sparse
#'   matrix batch. Smaller values reduce peak memory at the cost of slightly
#'   more overhead; larger values are faster but require more RAM. Default
#'   \code{50L} (about 12 batches for 558 LRMs). Reduce to between \code{20L}
#'   and \code{30L} on memory-constrained machines.
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
#' data(NICHESv2_example)
#'
#' obj <- add_neighborhood_edge_agg(NICHESv2_example, n.cores = 1L)
#' head(obj$aggregations$neighborhood.edge.agg)
#'
#' @export
add_neighborhood_edge_agg <- function(
    niches.obj,
    n.cores  = 1L,
    chunk.size = 50L,
    verbose  = TRUE
) {
  
  ## Input validation ----------------------------------------------------------
  if (!inherits(niches.obj, "NICHESObject"))
    stop("[add_neighborhood_edge_agg] 'niches.obj' must be a NICHESObject.",
         call. = FALSE)
  
  n.cores    <- as.integer(n.cores)
  chunk.size <- as.integer(chunk.size)
  if (is.na(n.cores) || n.cores < 1L)
    stop("[add_neighborhood_edge_agg] 'n.cores' must be a positive integer.",
         call. = FALSE)
  if (is.na(chunk.size) || chunk.size < 1L)
    stop("[add_neighborhood_edge_agg] 'chunk.size' must be a positive integer.",
         call. = FALSE)
  
  ## Neighborhood slot guard ---------------------------------------------------
  if (is.null(niches.obj$neighborhood.edge.list))
    stop(
      "[add_neighborhood_edge_agg] Neighborhood slots are NULL. ",
      "Run add_neighborhoods() first.",
      call. = FALSE
    )
  .stop_if_neighborhood_edges_missing(niches.obj, "add_neighborhood_edge_agg")

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
  
  # Classify each (focal cell, edge) into one of five categories.
  # Self-loop detection (SendingCell == ReceivingCell) is checked first and
  # immediately sub-divided: self.autocrine when the focal cell is the
  # self-looping cell; other.autocrine when a neighbor cell is.
  cell.edge.dt[, edge.type := data.table::fcase(
    SendingCell == ReceivingCell & SendingCell == cell,   "self.autocrine",
    SendingCell == ReceivingCell & SendingCell != cell,   "other.autocrine",
    SendingCell == cell & SendingCell != ReceivingCell,   "out",
    ReceivingCell == cell & SendingCell != ReceivingCell, "in",
    default = "cross"
  )]
  cell.edge.dt[, c("SendingCell", "ReceivingCell") := NULL]
  
  ## Extract edge scores (sparse: every stored row has score > 0) -------------
  edge.data.sub <- niches.obj$edge.data[, .(edge, LRM, score)]
  
  ## Precompute total neighborhood edges per cell (LRM-agnostic denominator) --
  # Used as denominator for score.*.mean. Precomputed once, passed to workers.
  total.nbr.edges <- niches.obj$neighborhood.meta[, .(cell, n.edges)]
  
  if (verbose)
    message(sprintf(
      "[add_neighborhood_edge_agg] Aggregating %d cells x %d LRMs x %d neighborhood edges (chunk.size = %d) ...",
      length(unique(cell.edge.dt$cell)),
      length(niches.obj$LRM.list),
      nrow(cell.edge.dt),
      chunk.size
    ))

  ## Chunked sparse matrix multiplication approach ----------------------------
  # Build 5 binary A matrices (cells x edges) once, then process LRMs in
  # batches of chunk.size. Each chunk multiplies A_t %*% B_chunk (score) and
  # A_t %*% B.ind_chunk (count), extracts the non-zero (cell, LRM) entries,
  # and discards the chunk matrices before processing the next batch.
  #
  # This bounds peak memory to one chunk at a time while still reducing the
  # number of iterations from n.lrms to ceil(n.lrms / chunk.size).
  #
  # n.cores is retained for API compatibility; not used (sparse matmul is
  # single-threaded; socket-cluster overhead exceeds gain for 5 tasks).

  all.cats  <- c("in", "out", "cross", "self.autocrine", "other.autocrine")

  # Index mappings
  all.cells <- unique(cell.edge.dt$cell)
  all.lrms  <- niches.obj$LRM.list
  n.cells   <- length(all.cells)
  n.lrms    <- length(all.lrms)

  # Restrict edge.data to edges that actually appear in neighbourhoods
  nbr.edge.set <- unique(cell.edge.dt$edge)
  ed.sub       <- edge.data.sub[edge %in% nbr.edge.set]
  common.edges <- unique(ed.sub$edge)
  n.e          <- length(common.edges)

  if (n.e == 0L) {
    if (verbose)
      message("[add_neighborhood_edge_agg] No scored edges in any neighborhood. Returning empty result.")
    niches.obj$aggregations$neighborhood.edge.agg <- data.table::data.table()
    return(niches.obj)
  }

  # Restrict cell.edge.dt to common edges and build index vectors
  ce.sub    <- cell.edge.dt[edge %in% common.edges]
  c.idx.vec <- match(ce.sub$cell, all.cells)
  e.idx.vec <- match(ce.sub$edge, common.edges)

  # Build A matrices once -- one per edge type (cells x edges), binary.
  # These are reused across all LRM chunks.
  A.list <- vector("list", 5L); names(A.list) <- all.cats
  for (cat in all.cats) {
    mask <- ce.sub$edge.type == cat
    if (any(mask))
      A.list[[cat]] <- Matrix::sparseMatrix(
        i = c.idx.vec[mask], j = e.idx.vec[mask], x = 1L,
        dims = c(n.cells, n.e)
      )
  }

  # Precompute edge-level index for fast chunk subsetting
  ed.sub.edge.idx <- match(ed.sub$edge, common.edges)
  ed.sub.lrm      <- ed.sub$LRM
  ed.sub.score    <- ed.sub$score

  # Split LRMs into chunks and process each independently
  lrm.chunks   <- split(seq_len(n.lrms),
                        ceiling(seq_len(n.lrms) / chunk.size))
  n.chunks     <- length(lrm.chunks)
  chunk.results <- vector("list", n.chunks)

  for (k in seq_len(n.chunks)) {

    chunk.l.idx  <- lrm.chunks[[k]]           # column indices into all.lrms
    chunk.lrms   <- all.lrms[chunk.l.idx]
    n.chunk      <- length(chunk.lrms)

    # Build B_chunk (n.e x n.chunk) and B.ind_chunk for this set of LRMs
    rows.k  <- ed.sub.lrm %in% chunk.lrms
    e.idx.k <- ed.sub.edge.idx[rows.k]
    l.idx.k <- match(ed.sub.lrm[rows.k], chunk.lrms)
    keep.k  <- !is.na(e.idx.k) & !is.na(l.idx.k)

    if (!any(keep.k)) {
      chunk.results[[k]] <- NULL
      next
    }

    B.k <- Matrix::sparseMatrix(
      i = e.idx.k[keep.k], j = l.idx.k[keep.k], x = ed.sub.score[rows.k][keep.k],
      dims = c(n.e, n.chunk)
    )
    B.ind.k <- Matrix::sparseMatrix(
      i = e.idx.k[keep.k], j = l.idx.k[keep.k], x = 1L,
      dims = c(n.e, n.chunk)
    )

    # 10 matrix products for this chunk
    mat.sum.k <- vector("list", 5L); names(mat.sum.k) <- all.cats
    mat.cnt.k <- vector("list", 5L); names(mat.cnt.k) <- all.cats
    zero.k    <- Matrix::sparseMatrix(i = integer(0L), j = integer(0L),
                                      dims = c(n.cells, n.chunk))
    for (cat in all.cats) {
      if (is.null(A.list[[cat]])) {
        mat.sum.k[[cat]] <- zero.k
        mat.cnt.k[[cat]] <- zero.k
      } else {
        mat.sum.k[[cat]] <- A.list[[cat]] %*% B.k
        mat.cnt.k[[cat]] <- A.list[[cat]] %*% B.ind.k
      }
    }

    # Identify active (cell, LRM) pairs in this chunk
    cnt.total.k <- Reduce("+", mat.cnt.k)
    nz.k        <- Matrix::summary(cnt.total.k)   # data.frame: i, j, x

    if (nrow(nz.k) == 0L) {
      chunk.results[[k]] <- NULL
      next
    }

    nz.mat.k <- cbind(nz.k$i, nz.k$j)

    # Build partial result table for this chunk
    cdt <- data.table::data.table(
      cell = all.cells[nz.mat.k[, 1L]],
      LRM  = chunk.lrms[nz.mat.k[, 2L]]
    )
    for (cat in all.cats) {
      cdt[, paste0("count.", cat) := as.integer(mat.cnt.k[[cat]][nz.mat.k])]
      cdt[, paste0("score.", cat, ".sum") := as.numeric(mat.sum.k[[cat]][nz.mat.k])]
    }

    chunk.results[[k]] <- cdt

    # Explicitly drop chunk matrices to free memory before next iteration
    rm(B.k, B.ind.k, mat.sum.k, mat.cnt.k, cnt.total.k, nz.k, nz.mat.k, cdt)
    gc(verbose = FALSE)

    if (verbose && (k %% max(1L, n.chunks %/% 5L) == 0L || k == n.chunks))
      message(sprintf("  [add_neighborhood_edge_agg] chunk %d / %d done.", k, n.chunks))
  }

  ## Assemble final table from chunk results -----------------------------------
  agg.dt <- data.table::rbindlist(
    Filter(Negate(is.null), chunk.results),
    fill = TRUE, use.names = TRUE
  )
  rm(chunk.results, A.list); gc(verbose = FALSE)

  if (nrow(agg.dt) == 0L) {
    if (verbose)
      message("[add_neighborhood_edge_agg] No active (cell, LRM) pairs found. Returning empty result.")
    niches.obj$aggregations$neighborhood.edge.agg <- agg.dt
    return(niches.obj)
  }

  data.table::setkeyv(agg.dt, c("cell", "LRM"))

  # Fill any category columns absent from all chunks with zero
  for (cat in all.cats) {
    cnt.col <- paste0("count.", cat);      if (is.null(agg.dt[[cnt.col]])) agg.dt[, (cnt.col) := 0L]
    sum.col <- paste0("score.", cat, ".sum"); if (is.null(agg.dt[[sum.col]])) agg.dt[, (sum.col) := 0]
  }

  agg.dt[, count.total :=
           count.in + count.out + count.cross +
           count.self.autocrine + count.other.autocrine]
  agg.dt[, score.total.sum :=
           score.in.sum + score.out.sum + score.cross.sum +
           score.self.autocrine.sum + score.other.autocrine.sum]

  ## Mean: LRM-agnostic denominator (total neighbourhood edges) ---------------
  idx     <- match(agg.dt$cell, total.nbr.edges$cell)
  n.edges <- total.nbr.edges$n.edges[idx]
  for (cat in all.cats)
    agg.dt[, paste0("score.", cat, ".mean") :=
             get(paste0("score.", cat, ".sum")) / n.edges]
  agg.dt[, score.total.mean := score.total.sum / n.edges]
  # Guard: cells absent from neighborhood.meta (should not occur)
  for (col in grep("\\.mean$", colnames(agg.dt), value = TRUE))
    agg.dt[is.na(get(col)), (col) := 0]

  ## Cmean: conditional mean (denominator = category count) -------------------
  for (cat in all.cats)
    agg.dt[, paste0("score.", cat, ".cmean") :=
             data.table::fifelse(
               get(paste0("count.", cat)) > 0L,
               get(paste0("score.", cat, ".sum")) / get(paste0("count.", cat)),
               NA_real_
             )]
  agg.dt[, score.total.cmean := data.table::fifelse(
    count.total > 0L, score.total.sum / count.total, NA_real_
  )]

  ## In/out ratio --------------------------------------------------------------
  agg.dt[, score.inout.ratio := data.table::fifelse(
    score.out.sum > 0, score.in.sum / score.out.sum, NA_real_
  )]
  
  # Enforce column order: identity | counts | sums | means | cmeans | ratio
  count.order <- c("count.in", "count.out", "count.cross",
                   "count.self.autocrine", "count.other.autocrine", "count.total")
  sum.order   <- c("score.in.sum",   "score.out.sum",   "score.cross.sum",
                   "score.self.autocrine.sum",   "score.other.autocrine.sum",
                   "score.total.sum")
  mean.order  <- c("score.in.mean",  "score.out.mean",  "score.cross.mean",
                   "score.self.autocrine.mean",  "score.other.autocrine.mean",
                   "score.total.mean")
  cmean.order <- c("score.in.cmean", "score.out.cmean", "score.cross.cmean",
                   "score.self.autocrine.cmean", "score.other.autocrine.cmean",
                   "score.total.cmean")
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