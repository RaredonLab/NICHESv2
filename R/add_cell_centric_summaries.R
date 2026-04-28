# add_cell_centric_summaries.R -----------------------------------------------
# Computes derived summary statistics from the cell-centric edge aggregation
# table ($aggregations$cell.centric.edge.agg) and appends them as new columns
# to the same table. No new data joins or parallel computation are needed --
# every statistic is a simple column-wise arithmetic transformation.
#
# Column naming convention: {base.metric}.{transformation}
#   e.g. score.sum.inout.ratio   = score.in.sum / score.out.sum
#        score.sum.inout.log2ratio = log2(score.in.sum / score.out.sum)
#        score.total.sum          = score.in.sum + score.out.sum  (cross-slot)
#
# Available summaries (11 total):
#   ratio.sum       -- score.sum.inout.ratio       score.in.sum / score.out.sum
#   ratio.mean      -- score.mean.inout.ratio       score.in.mean / score.out.mean
#   ratio.cmean     -- score.cmean.inout.ratio      score.in.cmean / score.out.cmean
#   log2ratio.sum   -- score.sum.inout.log2ratio    log2(score.in.sum / score.out.sum)
#   log2ratio.mean  -- score.mean.inout.log2ratio   log2(score.in.mean / score.out.mean)
#   log2ratio.cmean -- score.cmean.inout.log2ratio  log2(score.in.cmean / score.out.cmean)
#   total.sum       -- score.total.sum              score.in.sum + score.out.sum
#   total.mean      -- score.total.mean             score.in.mean + score.out.mean
#   total.cmean     -- score.total.cmean            score.in.cmean + score.out.cmean
#   asymmetry       -- asymmetry                    (count.in-count.out)/(count.in+count.out)
#   count.total     -- count.total                  count.in + count.out
#
# All undefined cases (division by zero, log of zero or negative) yield NA.
#
# Requires: $aggregations$cell.centric.edge.agg populated by
#           add_cell_centric_edge_agg(). Summaries involving a specific score
#           type additionally require that type to have been requested via
#           the methods argument of add_cell_centric_edge_agg().
#
# Source: $aggregations$cell.centric.edge.agg only -- no other slots read.
#
# DEPS: data.table


## CELL-CENTRIC SUMMARIES ----

#' Add Derived Summary Statistics to Cell-Centric Edge Aggregation
#'
#' Computes derived summary statistics from the cell-centric edge aggregation
#' table and appends them as new columns to
#' \code{$aggregations$cell.centric.edge.agg}. Three score types are covered
#' in parallel (\code{sum}, \code{mean}, \code{cmean}), plus count-based
#' statistics.
#'
#' \code{\link{add_cell_centric_edge_agg}} must be called first to populate
#' the source table. Summaries that require a specific score type will produce
#' a hard error if the corresponding base columns are absent; re-run
#' \code{add_cell_centric_edge_agg()} with the appropriate \code{methods} to
#' resolve this. In the standard workflow, \code{\link{aggregate_NICHESObject}}
#' calls \code{add_cell_centric_edge_agg} with all three methods before calling
#' this function, so the hard error is never triggered in practice.
#'
#' @section Column naming convention:
#' Summary columns follow the pattern \code{{base.metric}.{transformation}},
#' which groups related metrics alphabetically in any printout. For example,
#' all \code{score.sum.*} derived columns appear together. The exceptions are
#' \code{score.total.*} columns, which intentionally match the naming used in
#' \code{$aggregations$neighborhood.centric.edge.agg} for cross-slot
#' consistency.
#'
#' @section Summary definitions:
#' \strong{In/out ratio (receiver / sender):}
#' \describe{
#'   \item{\code{"ratio.sum"}, \code{"ratio.mean"}, \code{"ratio.cmean"}}{
#'     \code{score.in.{type} / score.out.{type}}. Values > 1 indicate the cell
#'     receives more of this signal than it sends; < 1 = sender-dominant.
#'     \code{NA} when denominator is zero or \code{NA}. Output columns:
#'     \code{score.sum.inout.ratio}, \code{score.mean.inout.ratio},
#'     \code{score.cmean.inout.ratio}.}
#' }
#' \strong{Log2 in/out ratio:}
#' \describe{
#'   \item{\code{"log2ratio.sum"}, \code{"log2ratio.mean"},
#'     \code{"log2ratio.cmean"}}{
#'     \code{log2(score.in.{type} / score.out.{type})}. Centred at 0
#'     (symmetric), positive = receiver-dominant, negative = sender-dominant.
#'     \code{NA} when either value is \eqn{\le 0} or \code{NA}. Output columns:
#'     \code{score.sum.inout.log2ratio}, \code{score.mean.inout.log2ratio},
#'     \code{score.cmean.inout.log2ratio}.}
#' }
#' \strong{Total signalling activity:}
#' \describe{
#'   \item{\code{"total.sum"}, \code{"total.mean"}, \code{"total.cmean"}}{
#'     \code{score.in.{type} + score.out.{type}}. Direction-agnostic total
#'     signalling involvement. For \code{total.cmean}: \code{NA} when either
#'     directional cmean is \code{NA} (i.e.\ that direction has no active
#'     edges for this LRM). Output columns: \code{score.total.sum},
#'     \code{score.total.mean}, \code{score.total.cmean}.}
#' }
#' \strong{Count-based summaries:}
#' \describe{
#'   \item{\code{"asymmetry"}}{
#'     \code{(count.in - count.out) / (count.in + count.out)}. Bounded in
#'     \eqn{[-1, 1]}: \code{+1} = pure receiver, \code{-1} = pure sender,
#'     \code{0} = symmetric. Based on edge counts (always available).
#'     \code{NA} when both counts are zero (cannot occur in a well-formed
#'     sparse table). Output column: \code{asymmetry}.}
#'   \item{\code{"count.total"}}{
#'     \code{count.in + count.out}. Total number of non-zero edges in either
#'     direction for this LRM. Always defined. Output column:
#'     \code{count.total}.}
#' }
#'
#' @section Incremental updates:
#' Re-calling with additional \code{summaries} and \code{overwrite = FALSE}
#' appends only the genuinely new columns; already-present columns are not
#' recomputed.
#'
#' @param niches.obj A \code{NICHESObject} with
#'   \code{$aggregations$cell.centric.edge.agg} populated by
#'   \code{\link{add_cell_centric_edge_agg}}.
#' @param summaries  Character vector. Any non-empty subset of
#'   \code{c("ratio.sum", "ratio.mean", "ratio.cmean",}
#'   \code{"log2ratio.sum", "log2ratio.mean", "log2ratio.cmean",}
#'   \code{"total.sum", "total.mean", "total.cmean",}
#'   \code{"asymmetry", "count.total")}. Default: all eleven.
#' @param overwrite  Logical. When \code{FALSE} (default), summary columns
#'   already present are skipped. When \code{TRUE}, all requested summaries
#'   are recomputed from scratch.
#' @param verbose    Logical. Emit progress messages. Default \code{TRUE}.
#'
#' @return The input \code{niches.obj} with new summary columns appended to
#'   \code{$aggregations$cell.centric.edge.agg}.
#'
#' @seealso \code{\link{add_cell_centric_edge_agg}},
#'   \code{\link{add_neighborhood_centric_summaries}},
#'   \code{\link{aggregate_NICHESObject}}
#'
#' @examples
#' \donttest{
#' data(NICHESv2_example)
#'
#' # Build aggregation table with all three methods first
#' obj <- add_cell_centric_edge_agg(
#'   NICHESv2_example,
#'   methods = c("sum", "mean", "cmean"),
#'   n.cores = 1L,
#'   verbose = FALSE
#' )
#'
#' # Add all eleven summary statistics
#' obj <- add_cell_centric_summaries(obj)
#' head(obj$aggregations$cell.centric.edge.agg)
#'
#' # Incremental: add only asymmetry and count.total
#' obj <- add_cell_centric_summaries(
#'   obj,
#'   summaries = c("asymmetry", "count.total"),
#'   overwrite = FALSE
#' )
#' }
#'
#' @export
add_cell_centric_summaries <- function(
    niches.obj,
    summaries = c(
      "ratio.sum",      "ratio.mean",      "ratio.cmean",
      "log2ratio.sum",  "log2ratio.mean",  "log2ratio.cmean",
      "total.sum",      "total.mean",      "total.cmean",
      "asymmetry",      "count.total"
    ),
    overwrite = FALSE,
    verbose   = TRUE
) {

  ## Input validation ----------------------------------------------------------
  if (!inherits(niches.obj, "NICHESObject"))
    stop("[add_cell_centric_summaries] 'niches.obj' must be a NICHESObject.",
         call. = FALSE)

  valid.summaries <- c(
    "ratio.sum",      "ratio.mean",      "ratio.cmean",
    "log2ratio.sum",  "log2ratio.mean",  "log2ratio.cmean",
    "total.sum",      "total.mean",      "total.cmean",
    "asymmetry",      "count.total"
  )
  summaries <- match.arg(summaries, valid.summaries, several.ok = TRUE)

  if (!is.logical(overwrite) || length(overwrite) != 1L)
    stop("[add_cell_centric_summaries] 'overwrite' must be a single logical.",
         call. = FALSE)

  ## Source slot guard ---------------------------------------------------------
  agg <- niches.obj$aggregations$cell.centric.edge.agg
  if (is.null(agg))
    stop(
      "[add_cell_centric_summaries] $aggregations$cell.centric.edge.agg is NULL. ",
      "Run add_cell_centric_edge_agg() first.",
      call. = FALSE
    )

  ## Base column requirement check --------------------------------------------
  # asymmetry and count.total only need count.in / count.out (always present).
  # All score-based summaries need the corresponding score type columns.
  agg.cols     <- colnames(agg)
  need.sum     <- c("ratio.sum",   "log2ratio.sum",   "total.sum")
  need.mean    <- c("ratio.mean",  "log2ratio.mean",  "total.mean")
  need.cmean   <- c("ratio.cmean", "log2ratio.cmean", "total.cmean")

  missing.methods <- character(0)
  if (any(summaries %in% need.sum) &&
      !all(c("score.in.sum",   "score.out.sum")   %in% agg.cols))
    missing.methods <- c(missing.methods, '"sum"')
  if (any(summaries %in% need.mean) &&
      !all(c("score.in.mean",  "score.out.mean")  %in% agg.cols))
    missing.methods <- c(missing.methods, '"mean"')
  if (any(summaries %in% need.cmean) &&
      !all(c("score.in.cmean", "score.out.cmean") %in% agg.cols))
    missing.methods <- c(missing.methods, '"cmean"')

  if (length(missing.methods) > 0L)
    stop(sprintf(
      "[add_cell_centric_summaries] Required base column(s) missing. ",
      "Re-run add_cell_centric_edge_agg() with methods including: %s.",
      paste(missing.methods, collapse = ", ")
    ), call. = FALSE)

  ## Map summary names to output column names ---------------------------------
  col.map <- c(
    ratio.sum        = "score.sum.inout.ratio",
    ratio.mean       = "score.mean.inout.ratio",
    ratio.cmean      = "score.cmean.inout.ratio",
    log2ratio.sum    = "score.sum.inout.log2ratio",
    log2ratio.mean   = "score.mean.inout.log2ratio",
    log2ratio.cmean  = "score.cmean.inout.log2ratio",
    total.sum        = "score.total.sum",
    total.mean       = "score.total.mean",
    total.cmean      = "score.total.cmean",
    asymmetry        = "asymmetry",
    count.total      = "count.total"
  )

  ## Resolve which summaries genuinely need computing -------------------------
  if (!overwrite) {
    summaries.todo <- summaries[!col.map[summaries] %in% agg.cols]
    if (length(summaries.todo) == 0L) {
      if (verbose)
        message("[add_cell_centric_summaries] All requested summaries already ",
                "present. Use overwrite = TRUE to recompute.")
      return(niches.obj)
    }
    if (verbose && length(summaries.todo) < length(summaries))
      message("[add_cell_centric_summaries] Skipping already-computed summaries; ",
              "adding: ", paste(summaries.todo, collapse = ", "), ".")
    summaries <- summaries.todo
  }

  if (verbose)
    message(sprintf(
      "[add_cell_centric_summaries] Computing: %s ...",
      paste(summaries, collapse = ", ")
    ))

  ## Compute summaries ---------------------------------------------------------

  # -- Ratio (receiver / sender) ----------------------------------------------
  # NA when denominator is zero (sum/mean can be zero-filled) or NA (cmean).
  # For cmean: values are either NA or strictly positive (sparse edge.data),
  # so NA propagates naturally without a separate fifelse guard.

  if ("ratio.sum" %in% summaries)
    agg[, score.sum.inout.ratio := fifelse(
      score.out.sum > 0, score.in.sum / score.out.sum, NA_real_
    )]

  if ("ratio.mean" %in% summaries)
    agg[, score.mean.inout.ratio := fifelse(
      score.out.mean > 0, score.in.mean / score.out.mean, NA_real_
    )]

  if ("ratio.cmean" %in% summaries)
    agg[, score.cmean.inout.ratio := score.in.cmean / score.out.cmean]

  # -- Log2 ratio (log2(receiver / sender)) -----------------------------------
  # NA when either value is <= 0 or NA.

  if ("log2ratio.sum" %in% summaries)
    agg[, score.sum.inout.log2ratio := fifelse(
      score.in.sum > 0 & score.out.sum > 0,
      log2(score.in.sum / score.out.sum),
      NA_real_
    )]

  if ("log2ratio.mean" %in% summaries)
    agg[, score.mean.inout.log2ratio := fifelse(
      score.in.mean > 0 & score.out.mean > 0,
      log2(score.in.mean / score.out.mean),
      NA_real_
    )]

  if ("log2ratio.cmean" %in% summaries)
    agg[, score.cmean.inout.log2ratio := fifelse(
      !is.na(score.in.cmean) & !is.na(score.out.cmean),
      log2(score.in.cmean / score.out.cmean),
      NA_real_
    )]

  # -- Total signalling activity (in + out) -----------------------------------
  # sum/mean are always defined when base columns exist.
  # total.cmean is NA when either directional cmean is NA (that direction
  # has no active edges for this LRM) -- NA propagates naturally.

  if ("total.sum" %in% summaries)
    agg[, score.total.sum   := score.in.sum   + score.out.sum]

  if ("total.mean" %in% summaries)
    agg[, score.total.mean  := score.in.mean  + score.out.mean]

  if ("total.cmean" %in% summaries)
    agg[, score.total.cmean := score.in.cmean + score.out.cmean]

  # -- Count-based summaries --------------------------------------------------

  # asymmetry: bounded [-1, 1]; NA when both counts are zero (cannot occur
  # in a well-formed sparse table).
  if ("asymmetry" %in% summaries)
    agg[, asymmetry := fifelse(
      (count.in + count.out) > 0L,
      (count.in - count.out) / (count.in + count.out),
      NA_real_
    )]

  # count.total: always defined.
  if ("count.total" %in% summaries)
    agg[, count.total := count.in + count.out]

  ## Write back to object ------------------------------------------------------
  niches.obj$aggregations$cell.centric.edge.agg <- agg

  if (verbose)
    message(sprintf(
      "[add_cell_centric_summaries] Done. Added %d summary column(s) to %d rows.",
      length(summaries),
      nrow(agg)
    ))

  niches.obj
}
