# add_neighborhood_centric_edge_agg.R ----------------------------------------
# For each focal cell and each LRM pair, aggregate scores across all edges
# within the cell's neighborhood (as defined by the neighborhood slots).
#
# Supports three aggregation methods:
#   sum   -- total score across neighborhood edges for this LRM
#   mean  -- sum / total neighborhood edges (LRM-agnostic denominator)
#   cmean -- sum / count of non-zero neighborhood edges (conditional mean)
#
# Optional directional decomposition (by.direction = TRUE) splits each
# cell's neighborhood edges into three categories relative to the focal cell:
#   focal.in  -- focal cell is the ReceivingCell
#   focal.out -- focal cell is the SendingCell
#   cross     -- focal cell is neither endpoint
#
# Results stored in $aggregations$neighborhood.centric.edge.agg as a named
# list of data.tables. Each element corresponds to one combination of
# neighborhood parameters. Key format:
#   "{edge.filter.mode}.ac={!remove.autocrine}"
#   e.g. "both_in_neighborhood.ac=TRUE"  (autocrine included)
#        "central_cell_only.ac=FALSE"    (autocrine removed)
#
# Sources: $neighborhood.edge.list, $neighborhood.meta, $edge.data, $edge.list
#
# DEPS: data.table, parallel


## NEIGHBORHOOD-CENTRIC EDGE AGGREGATION ----

#' Add Neighborhood-Centric Edge Aggregation to a NICHESObject
#'
#' For each focal cell and each LRM pair, aggregates LR communication scores
#' across all edges that fall within the cell's neighborhood (as defined by the
#' neighborhood slots). Results capture the overall signalling environment
#' around each cell, rather than only edges directly involving it.
#'
#' Results are stored in \code{$aggregations$neighborhood.centric.edge.agg} as
#' a named list of data.tables. Each element is keyed by the neighborhood
#' parameter combination used (see Key format below), so multiple versions
#' computed under different \code{edge.filter.mode} or
#' \code{remove.autocrine} settings coexist without overwriting each other.
#'
#' @section Key format:
#' The storage key encodes both neighborhood parameters:
#' \code{"{edge.filter.mode}.ac={!remove.autocrine}"}
#' where \code{ac=TRUE} means autocrine edges \emph{are included} (i.e.\
#' \code{remove.autocrine = FALSE}).
#' Example keys:
#' \itemize{
#'   \item \code{"both_in_neighborhood.ac=TRUE"}
#'   \item \code{"central_cell_only.ac=FALSE"}
#' }
#'
#' @section Column definitions:
#' The following columns are always computed, regardless of \code{methods}:
#' \describe{
#'   \item{\code{count.total}}{Number of non-zero-score neighborhood edges for
#'     this LRM. Because \code{$edge.data} is sparse, every stored row has
#'     score > 0, so this is simply the row count per focal-cell group.}
#' }
#' Score columns (produced according to \code{methods}):
#' \describe{
#'   \item{\code{score.total.sum}}{Sum of scores across all neighborhood edges
#'     for this LRM. Requested via \code{"sum"}.}
#'   \item{\code{score.total.mean}}{Sum divided by the \emph{total} number of
#'     edges in the neighborhood (from \code{$neighborhood.meta$n.edges}),
#'     including edges with zero score for this LRM. Penalises cells with
#'     large but silent neighbourhoods. Requested via \code{"mean"}.}
#'   \item{\code{score.total.cmean}}{Conditional mean: sum divided by
#'     \code{count.total}. \code{NA} when \code{count.total} is zero.
#'     Reports signal strength \emph{given} that signalling occurs.
#'     Requested via \code{"cmean"}.}
#' }
#' When \code{by.direction = TRUE}, the three directional variants of each
#' count and score column are also produced:
#' \describe{
#'   \item{\code{count.focal.in}, \code{score.focal.in.*}}{Edges where the
#'     focal cell is the \emph{receiver}.}
#'   \item{\code{count.focal.out}, \code{score.focal.out.*}}{Edges where the
#'     focal cell is the \emph{sender}.}
#'   \item{\code{count.cross}, \code{score.cross.*}}{Edges between two
#'     \emph{neighbours} of the focal cell (neither endpoint is the focal cell
#'     itself). Always zero under \code{edge.filter.mode =
#'     "central_cell_only"}.}
#' }
#'
#' @section Incremental updates:
#' Calling this function again with additional \code{methods} or with
#' \code{by.direction = TRUE} on a key that already exists appends only the
#' missing columns; existing columns are not recomputed unless
#' \code{overwrite = TRUE}.
#'
#' @section Non-spatial mode:
#' In \code{"sampled"} mode the neighbourhood edges are statistically paired
#' rather than spatially proximate. The computation is identical but the
#' biological interpretation shifts: results reflect signalling activity among
#' co-occurring cell pairs rather than a local tissue microenvironment. Use
#' and interpret with care.
#'
#' @param niches.obj   A \code{NICHESObject} with neighborhood slots populated.
#' @param neighborhoods Optional list returned by
#'   \code{\link{extract_neighborhoods}}, or \code{NULL} (default). When
#'   \code{NULL}, the object's own neighborhood slots and attributes are used.
#'   Supply an external list to aggregate a neighbourhood result that was not
#'   stored back into the object (e.g.\ produced with a custom edge list).
#' @param methods      Character vector. Any non-empty subset of
#'   \code{c("sum", "mean", "cmean")}. Default: all three.
#' @param by.direction Logical. When \code{TRUE}, additionally decompose each
#'   cell's neighborhood aggregation into three directional categories relative
#'   to the focal cell (\code{focal.in}, \code{focal.out}, \code{cross}).
#'   Default \code{FALSE}. A warning is emitted when \code{TRUE} and
#'   \code{edge.filter.mode = "central_cell_only"} because \code{cross} edges
#'   will always be empty in that mode.
#' @param overwrite    Logical. When \code{FALSE} (default), only missing
#'   columns are appended to the existing entry for this key; already-computed
#'   columns are preserved. When \code{TRUE}, the entry for this key is fully
#'   recomputed from scratch.
#' @param n.cores      Positive integer. Number of parallel worker processes.
#'   Parallelises over LRM pairs via socket clusters. Default \code{1L}.
#' @param verbose      Logical. Emit progress messages. Default \code{TRUE}.
#'
#' @return The input \code{niches.obj} with
#'   \code{$aggregations$neighborhood.centric.edge.agg} populated or updated.
#'
#' @seealso \code{\link{add_neighborhoods}}, \code{\link{extract_neighborhoods}},
#'   \code{\link{add_cell_centric_edge_agg}}
#'
#' @examples
#' \donttest{
#' data(NICHESv2_example)
#'
#' # Default: total aggregation only
#' obj <- add_neighborhood_centric_edge_agg(
#'   NICHESv2_example,
#'   methods  = c("sum", "mean", "cmean"),
#'   n.cores  = 1L
#' )
#' # Access via the stored key
#' key <- "both_in_neighborhood.ac=TRUE"
#' head(obj$aggregations$neighborhood.centric.edge.agg[[key]])
#'
#' # Add directional breakdown without recomputing existing columns
#' obj <- add_neighborhood_centric_edge_agg(
#'   obj,
#'   by.direction = TRUE,
#'   overwrite    = FALSE
#' )
#'
#' # Compute under a second neighbourhood mode -- stored alongside the first
#' obj2 <- add_neighborhoods(NICHESv2_example,
#'                           edge.filter.mode = "central_cell_only",
#'                           overwrite        = TRUE)
#' obj2 <- add_neighborhood_centric_edge_agg(obj2, methods = "sum")
#' names(obj2$aggregations$neighborhood.centric.edge.agg)  # two keys
#' }
#'
#' @export
add_neighborhood_centric_edge_agg <- function(
    niches.obj,
    neighborhoods  = NULL,
    methods        = c("sum", "mean", "cmean"),
    by.direction   = FALSE,
    overwrite      = FALSE,
    n.cores        = 1L,
    verbose        = TRUE
) {

  ## Input validation ----------------------------------------------------------
  if (!inherits(niches.obj, "NICHESObject"))
    stop("[add_neighborhood_centric_edge_agg] 'niches.obj' must be a NICHESObject.",
         call. = FALSE)

  methods <- match.arg(methods, c("sum", "mean", "cmean"), several.ok = TRUE)

  if (!is.logical(by.direction) || length(by.direction) != 1L)
    stop("[add_neighborhood_centric_edge_agg] 'by.direction' must be a single logical.",
         call. = FALSE)

  if (!is.logical(overwrite) || length(overwrite) != 1L)
    stop("[add_neighborhood_centric_edge_agg] 'overwrite' must be a single logical.",
         call. = FALSE)

  n.cores <- as.integer(n.cores)
  if (is.na(n.cores) || n.cores < 1L)
    stop("[add_neighborhood_centric_edge_agg] 'n.cores' must be a positive integer.",
         call. = FALSE)

  ## Resolve neighborhoods source ---------------------------------------------
  if (is.null(neighborhoods)) {
    if (is.null(niches.obj$neighborhood.edge.list))
      stop("[add_neighborhood_centric_edge_agg] Neighborhood slots are NULL. ",
           "Run add_neighborhoods() first.", call. = FALSE)
    nbr.edge.list <- niches.obj$neighborhood.edge.list
    nbr.meta      <- niches.obj$neighborhood.meta
    efm           <- attr(niches.obj, "edge.filter.mode")
    rem.auto      <- attr(niches.obj, "remove.autocrine")
    if (is.null(efm) || is.null(rem.auto))
      stop("[add_neighborhood_centric_edge_agg] Object is missing 'edge.filter.mode' ",
           "or 'remove.autocrine' attributes. Re-run add_neighborhoods() to update.",
           call. = FALSE)
  } else {
    if (!is.list(neighborhoods) || is.null(neighborhoods$neighborhood.edge.list))
      stop("[add_neighborhood_centric_edge_agg] 'neighborhoods' must be a list ",
           "returned by extract_neighborhoods().", call. = FALSE)
    nbr.edge.list <- neighborhoods$neighborhood.edge.list
    nbr.meta      <- neighborhoods$neighborhood.meta
    efm           <- attr(neighborhoods, "edge.filter.mode")
    rem.auto      <- attr(neighborhoods, "remove.autocrine")
    if (is.null(efm) || is.null(rem.auto))
      stop("[add_neighborhood_centric_edge_agg] 'neighborhoods' is missing ",
           "'edge.filter.mode' or 'remove.autocrine' attributes. Ensure it was ",
           "produced by extract_neighborhoods().", call. = FALSE)
  }

  ## Build storage key ---------------------------------------------------------
  # Key encodes both neighborhood parameters.
  # ac=TRUE means autocrine edges ARE included (!remove.autocrine).
  agg.key <- paste0(efm, ".ac=", !rem.auto)

  ## Warn for uninformative by.direction combination --------------------------
  if (by.direction && efm == "central_cell_only")
    warning(
      "[add_neighborhood_centric_edge_agg] 'by.direction = TRUE' with ",
      "edge.filter.mode = 'central_cell_only': 'cross' direction edges will ",
      "always be empty because every neighborhood edge involves the focal cell.",
      call. = FALSE
    )

  ## Determine which methods and direction components still need computing -----
  existing.slot  <- niches.obj$aggregations$neighborhood.centric.edge.agg
  existing.entry <- existing.slot[[agg.key]]

  if (!overwrite && !is.null(existing.entry)) {
    existing.cols <- colnames(existing.entry)

    methods.todo <- methods[c(
      "sum"   %in% methods & !("score.total.sum"   %in% existing.cols),
      "mean"  %in% methods & !("score.total.mean"  %in% existing.cols),
      "cmean" %in% methods & !("score.total.cmean" %in% existing.cols)
    )]
    direction.needed <- by.direction && !("count.focal.in" %in% existing.cols)

    if (length(methods.todo) == 0L && !direction.needed) {
      if (verbose)
        message("[add_neighborhood_centric_edge_agg] All requested methods ",
                "already present for key '", agg.key,
                "'. Use overwrite = TRUE to recompute.")
      return(niches.obj)
    }
    if (verbose && length(methods.todo) < length(methods))
      message("[add_neighborhood_centric_edge_agg] Skipping already-computed ",
              "methods; adding: ",
              paste(c(methods.todo,
                      if (direction.needed) "direction"), collapse = ", "), ".")

    # methods.total  : new total-score methods to compute
    # methods.dir    : direction score methods to compute (all user methods,
    #                  first time direction columns appear)
    methods.total  <- methods.todo
    methods.dir    <- if (direction.needed) methods else character(0)
    incremental    <- TRUE
    direction.todo <- direction.needed
  } else {
    methods.total  <- methods
    methods.dir    <- if (by.direction) methods else character(0)
    incremental    <- FALSE
    direction.todo <- by.direction
  }

  ## Build cell→edge long table ------------------------------------------------
  # Unroll neighborhood.edge.list: one row per (focal cell, neighborhood edge).
  nbr.lens     <- lengths(nbr.edge.list)
  cell.edge.dt <- data.table(
    cell = rep(names(nbr.edge.list), nbr.lens),
    edge = unlist(nbr.edge.list, use.names = FALSE)
  )

  if (nrow(cell.edge.dt) == 0L)
    stop("[add_neighborhood_centric_edge_agg] No edges found in any neighborhood.",
         call. = FALSE)

  ## Join edge.data to get per-(focal cell, edge, LRM) scores -----------------
  # edge.data is SPARSE: every stored row has score > 0. No filter needed.
  # merge() handles the many-to-many expansion (same edge in multiple hoods)
  # safely without the composite-key join gotcha.
  edge.data.sub <- niches.obj$edge.data[, .(edge, LRM, score)]
  flat <- merge(cell.edge.dt, edge.data.sub, by = "edge")
  # flat: edge | cell | LRM | score

  ## Attach direction classification (pre-loop; LRM-agnostic) -----------------
  if (direction.todo) {
    # Use match-based lookup to avoid composite-key join issue.
    el   <- niches.obj$edge.list
    idx  <- match(flat$edge, el$edge)
    flat[, SendingCell   := el$SendingCell[idx]]
    flat[, ReceivingCell := el$ReceivingCell[idx]]
    flat[, direction := data.table::fcase(
      SendingCell == cell,   "focal.out",
      ReceivingCell == cell, "focal.in",
      default = "cross"
    )]
    flat[, c("SendingCell", "ReceivingCell") := NULL]
  }

  ## Precompute total neighborhood edge counts per cell (LRM-agnostic) --------
  # nbr.meta$n.edges is the denominator for score.total.mean /
  # score.focal.*.mean. Extracted once and passed to every worker.
  total.nbr.edges <- nbr.meta[, .(cell, n.edges)]

  if (verbose)
    message(sprintf(
      "[add_neighborhood_centric_edge_agg] Aggregating %d LRM pairs | key: '%s' ...",
      length(niches.obj$LRM.list), agg.key
    ))

  ## Worker function -----------------------------------------------------------
  # Processes one LRM's slice of the flat table.
  # All named arguments are explicit to allow safe clusterExport.
  .nbr_centric_worker <- function(lrm.dt, total.nbr.edges,
                                   methods.total, methods.dir,
                                   direction.todo) {
    lrm.id <- lrm.dt$LRM[1L]

    ## Total aggregation ----
    total.agg <- lrm.dt[,
      .(count.total = .N, score.total.sum = sum(score)),
      by = cell
    ]

    if ("mean" %in% methods.total) {
      idx <- match(total.agg$cell, total.nbr.edges$cell)
      total.agg[, score.total.mean :=
                  score.total.sum / total.nbr.edges$n.edges[idx]]
      total.agg[is.na(score.total.mean), score.total.mean := 0]
    }
    if ("cmean" %in% methods.total) {
      total.agg[, score.total.cmean := data.table::fifelse(
        count.total > 0L, score.total.sum / count.total, NA_real_
      )]
    }
    if (!"sum" %in% methods.total)
      total.agg[, score.total.sum := NULL]

    ## Directional aggregation ----
    if (direction.todo) {

      dir.agg <- lrm.dt[,
        .(n.dir = .N, score.sum = sum(score)),
        by = .(cell, direction)
      ]

      # Pivot to wide: columns for each direction × value combination.
      # dcast names: "n.dir_focal.in", "n.dir_focal.out", "n.dir_cross",
      #              "score.sum_focal.in", etc.
      dir.wide <- data.table::dcast(
        dir.agg, cell ~ direction,
        value.var = c("n.dir", "score.sum"),
        fill      = 0
      )

      # Ensure all six direction columns exist even if a direction never
      # occurs for this LRM (dcast omits absent levels).
      expected.raw <- c("n.dir_focal.in",    "n.dir_focal.out",    "n.dir_cross",
                        "score.sum_focal.in", "score.sum_focal.out", "score.sum_cross")
      for (col in setdiff(expected.raw, colnames(dir.wide))) {
        dir.wide[, (col) := if (startsWith(col, "n.dir")) 0L else 0]
      }

      # Rename to final NICHESv2 column convention.
      data.table::setnames(dir.wide,
        old = c("n.dir_focal.in",    "n.dir_focal.out",    "n.dir_cross",
                "score.sum_focal.in", "score.sum_focal.out", "score.sum_cross"),
        new = c("count.focal.in",    "count.focal.out",    "count.cross",
                "score.focal.in.sum", "score.focal.out.sum", "score.cross.sum")
      )
      # Ensure count columns are integer
      dir.wide[, count.focal.in  := as.integer(count.focal.in)]
      dir.wide[, count.focal.out := as.integer(count.focal.out)]
      dir.wide[, count.cross     := as.integer(count.cross)]

      if ("mean" %in% methods.dir) {
        idx     <- match(dir.wide$cell, total.nbr.edges$cell)
        n.edges <- total.nbr.edges$n.edges[idx]
        dir.wide[, score.focal.in.mean  := score.focal.in.sum  / n.edges]
        dir.wide[, score.focal.out.mean := score.focal.out.sum / n.edges]
        dir.wide[, score.cross.mean     := score.cross.sum     / n.edges]
        dir.wide[is.na(score.focal.in.mean),  score.focal.in.mean  := 0]
        dir.wide[is.na(score.focal.out.mean), score.focal.out.mean := 0]
        dir.wide[is.na(score.cross.mean),     score.cross.mean     := 0]
      }
      if ("cmean" %in% methods.dir) {
        dir.wide[, score.focal.in.cmean  := data.table::fifelse(
          count.focal.in  > 0L, score.focal.in.sum  / count.focal.in,  NA_real_)]
        dir.wide[, score.focal.out.cmean := data.table::fifelse(
          count.focal.out > 0L, score.focal.out.sum / count.focal.out, NA_real_)]
        dir.wide[, score.cross.cmean     := data.table::fifelse(
          count.cross     > 0L, score.cross.sum     / count.cross,     NA_real_)]
      }
      if (!"sum" %in% methods.dir)
        dir.wide[, c("score.focal.in.sum",
                     "score.focal.out.sum",
                     "score.cross.sum") := NULL]

      result <- merge(total.agg, dir.wide, by = "cell", all = TRUE)
    } else {
      result <- total.agg
    }

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
      varlist = c(".nbr_centric_worker", "total.nbr.edges",
                  "methods.total", "methods.dir", "direction.todo"),
      envir   = environment()
    )
    results <- parallel::parLapply(
      cl, lrm.split,
      function(lrm.dt) .nbr_centric_worker(
        lrm.dt, total.nbr.edges, methods.total, methods.dir, direction.todo
      )
    )
  } else {
    results <- lapply(
      lrm.split,
      function(lrm.dt) .nbr_centric_worker(
        lrm.dt, total.nbr.edges, methods.total, methods.dir, direction.todo
      )
    )
  }

  ## Assemble final table ------------------------------------------------------
  # fill = TRUE handles LRMs where some direction columns were absent (all-zero
  # neighborhoods for that direction), which dcast may not produce.
  agg.dt <- rbindlist(results, fill = TRUE)
  setkeyv(agg.dt, c("cell", "LRM"))

  # Enforce column order: identity, total counts, total scores, dir counts,
  # dir scores grouped by direction.
  total.score.order <- c("score.total.sum",  "score.total.mean",  "score.total.cmean")
  dir.count.order   <- c("count.focal.in",   "count.focal.out",   "count.cross")
  dir.score.order   <- c(
    "score.focal.in.sum",   "score.focal.in.mean",   "score.focal.in.cmean",
    "score.focal.out.sum",  "score.focal.out.mean",  "score.focal.out.cmean",
    "score.cross.sum",      "score.cross.mean",       "score.cross.cmean"
  )
  col.order <- c(
    "cell", "LRM", "count.total",
    intersect(total.score.order, colnames(agg.dt)),
    intersect(dir.count.order,   colnames(agg.dt)),
    intersect(dir.score.order,   colnames(agg.dt))
  )
  setcolorder(agg.dt, col.order)

  ## Store result: full replace or incremental column join --------------------
  if (incremental && !is.null(existing.entry)) {
    # Add only columns not already present in the existing entry.
    # count.total is recomputed but not re-stored; existing value is
    # authoritative under overwrite = FALSE.
    new.cols <- setdiff(colnames(agg.dt), colnames(existing.entry))
    if (length(new.cols) > 0L) {
      join.dt        <- agg.dt[, c("cell", "LRM", new.cols), with = FALSE]
      existing.entry <- existing.entry[join.dt, on = c("cell", "LRM")]
    }
    niches.obj$aggregations$neighborhood.centric.edge.agg[[agg.key]] <-
      existing.entry
  } else {
    if (is.null(niches.obj$aggregations$neighborhood.centric.edge.agg))
      niches.obj$aggregations$neighborhood.centric.edge.agg <- list()
    niches.obj$aggregations$neighborhood.centric.edge.agg[[agg.key]] <- agg.dt
  }

  if (verbose)
    message(sprintf(
      "[add_neighborhood_centric_edge_agg] Done. %d rows | %d cells | %d LRM pairs | key: '%s'.",
      nrow(agg.dt),
      uniqueN(agg.dt$cell),
      uniqueN(agg.dt$LRM),
      agg.key
    ))

  niches.obj
}
