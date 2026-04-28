# add_neighborhood_centric_summaries.R ---------------------------------------
# Computes derived summary statistics from a neighborhood-centric edge
# aggregation entry ($aggregations$neighborhood.centric.edge.agg[[key]]) and
# appends them as new columns to that entry. No new data joins or parallel
# computation are needed -- every statistic is a simple column-wise arithmetic
# transformation.
#
# All summaries require by.direction = TRUE columns to have been computed by
# add_neighborhood_centric_edge_agg(). Score-based ratio summaries additionally
# require the corresponding score type columns.
#
# Column naming convention: {base.metric}.{transformation}
#   e.g. score.focal.sum.inout.ratio    = score.focal.in.sum / score.focal.out.sum
#        score.focal.sum.inout.log2ratio = log2(score.focal.in.sum / score.focal.out.sum)
#        count.focal.asymmetry           = (count.focal.in - count.focal.out) /
#                                          (count.focal.in + count.focal.out)
#        count.cross.fraction            = count.cross / count.total
#        count.focal.fraction            = (count.focal.in + count.focal.out) / count.total
#
# Available summaries (9 total):
#   ratio.focal.sum       -- score.focal.sum.inout.ratio
#   ratio.focal.mean      -- score.focal.mean.inout.ratio
#   ratio.focal.cmean     -- score.focal.cmean.inout.ratio
#   log2ratio.focal.sum   -- score.focal.sum.inout.log2ratio
#   log2ratio.focal.mean  -- score.focal.mean.inout.log2ratio
#   log2ratio.focal.cmean -- score.focal.cmean.inout.log2ratio
#   focal.asymmetry       -- count.focal.asymmetry
#   cross.fraction        -- count.cross.fraction
#   focal.fraction        -- count.focal.fraction
#
# All undefined cases (division by zero, log of zero or negative) yield NA.
#
# Source: $aggregations$neighborhood.centric.edge.agg[[key]] only.
#
# DEPS: data.table


## NEIGHBORHOOD-CENTRIC SUMMARIES ----

#' Add Derived Summary Statistics to Neighborhood-Centric Edge Aggregation
#'
#' Computes derived summary statistics from a neighborhood-centric edge
#' aggregation entry and appends them as new columns to that entry. All
#' summaries operate on the directional decomposition columns produced by
#' \code{\link{add_neighborhood_centric_edge_agg}} with
#' \code{by.direction = TRUE}.
#'
#' \code{\link{add_neighborhood_centric_edge_agg}} must be called first with
#' \code{by.direction = TRUE} to populate the directional source columns.
#' Score-based ratio summaries additionally require the corresponding score
#' type. In the standard workflow, \code{\link{aggregate_NICHESObject}} ensures
#' all prerequisites are met before calling this function.
#'
#' @section Key selection:
#' Results are stored in a named list keyed by neighborhood parameter
#' combination (see \code{\link{add_neighborhood_centric_edge_agg}} for the
#' key format). When \code{key = NULL} (default), the key is derived from
#' \code{attr(niches.obj, "edge.filter.mode")} and
#' \code{attr(niches.obj, "remove.autocrine")}. Supply \code{key} explicitly
#' to target a specific entry when multiple keys are present.
#'
#' @section Summary definitions:
#' \strong{Focal in/out ratio (receiver vs sender among focal edges):}
#' \describe{
#'   \item{\code{"ratio.focal.sum"}, \code{"ratio.focal.mean"},
#'     \code{"ratio.focal.cmean"}}{
#'     \code{score.focal.in.{type} / score.focal.out.{type}}. Focal-edge
#'     analogue of the cell-centric in/out ratio: values > 1 indicate the
#'     focal cell receives more signal via direct edges than it sends; < 1
#'     = sender-dominant. \code{NA} when denominator is zero or \code{NA}.
#'     Output columns: \code{score.focal.sum.inout.ratio},
#'     \code{score.focal.mean.inout.ratio},
#'     \code{score.focal.cmean.inout.ratio}.}
#' }
#' \strong{Log2 focal in/out ratio:}
#' \describe{
#'   \item{\code{"log2ratio.focal.sum"}, \code{"log2ratio.focal.mean"},
#'     \code{"log2ratio.focal.cmean"}}{
#'     \code{log2(score.focal.in.{type} / score.focal.out.{type})}. Log-scale
#'     version; centred at 0, positive = receiver-dominant. \code{NA} when
#'     either value is \eqn{\le 0} or \code{NA}. Output columns:
#'     \code{score.focal.sum.inout.log2ratio},
#'     \code{score.focal.mean.inout.log2ratio},
#'     \code{score.focal.cmean.inout.log2ratio}.}
#' }
#' \strong{Count-based neighborhood structure summaries:}
#' \describe{
#'   \item{\code{"focal.asymmetry"}}{
#'     \code{(count.focal.in - count.focal.out) / (count.focal.in + count.focal.out)}.
#'     Bounded in \eqn{[-1, 1]}: asymmetry among direct focal edges only.
#'     \code{NA} when the focal cell has no direct edges in its neighborhood
#'     for this LRM (i.e.\ all activity is cross-edge). Output column:
#'     \code{count.focal.asymmetry}.}
#'   \item{\code{"cross.fraction"}}{
#'     \code{count.cross / count.total}. Fraction of neighborhood edges that
#'     do not involve the focal cell — a measure of how much background
#'     signalling exists in the local environment. Bounded in \eqn{[0, 1]}.
#'     Output column: \code{count.cross.fraction}.}
#'   \item{\code{"focal.fraction"}}{
#'     \code{(count.focal.in + count.focal.out) / count.total}. Fraction of
#'     neighborhood edges that directly involve the focal cell. Bounded in
#'     \eqn{[0, 1]}. Note: \code{cross.fraction + focal.fraction = 1}.
#'     Output column: \code{count.focal.fraction}.}
#' }
#'
#' @section Incremental updates:
#' Re-calling with additional \code{summaries} and \code{overwrite = FALSE}
#' appends only the genuinely new columns; already-present columns are not
#' recomputed.
#'
#' @param niches.obj A \code{NICHESObject} with
#'   \code{$aggregations$neighborhood.centric.edge.agg} populated by
#'   \code{\link{add_neighborhood_centric_edge_agg}} with
#'   \code{by.direction = TRUE}.
#' @param summaries  Character vector. Any non-empty subset of
#'   \code{c("ratio.focal.sum", "ratio.focal.mean", "ratio.focal.cmean",}
#'   \code{"log2ratio.focal.sum", "log2ratio.focal.mean",}
#'   \code{"log2ratio.focal.cmean",}
#'   \code{"focal.asymmetry", "cross.fraction", "focal.fraction")}.
#'   Default: all nine.
#' @param key        Character or \code{NULL}. The named-list key identifying
#'   which entry of \code{$aggregations$neighborhood.centric.edge.agg} to
#'   operate on. When \code{NULL} (default), derived from the object's
#'   \code{"edge.filter.mode"} and \code{"remove.autocrine"} attributes.
#' @param overwrite  Logical. When \code{FALSE} (default), summary columns
#'   already present are skipped. When \code{TRUE}, all requested summaries
#'   are recomputed from scratch.
#' @param verbose    Logical. Emit progress messages. Default \code{TRUE}.
#'
#' @return The input \code{niches.obj} with new summary columns appended to
#'   the targeted entry of
#'   \code{$aggregations$neighborhood.centric.edge.agg}.
#'
#' @seealso \code{\link{add_neighborhood_centric_edge_agg}},
#'   \code{\link{add_cell_centric_summaries}},
#'   \code{\link{aggregate_NICHESObject}}
#'
#' @examples
#' \donttest{
#' data(NICHESv2_example)
#'
#' # Build aggregation table with all methods and directional decomposition
#' obj <- add_neighborhood_centric_edge_agg(
#'   NICHESv2_example,
#'   methods      = c("sum", "mean", "cmean"),
#'   by.direction = TRUE,
#'   n.cores      = 1L,
#'   verbose      = FALSE
#' )
#'
#' # Add all nine summary statistics for the default key
#' obj <- add_neighborhood_centric_summaries(obj)
#'
#' key <- paste0("both_in_neighborhood.ac=TRUE")
#' head(obj$aggregations$neighborhood.centric.edge.agg[[key]])
#'
#' # Target a specific key explicitly
#' obj2 <- add_neighborhoods(NICHESv2_example,
#'                           edge.filter.mode = "central_cell_only",
#'                           overwrite = TRUE)
#' obj2 <- add_neighborhood_centric_edge_agg(obj2,
#'                           methods = c("sum", "mean", "cmean"),
#'                           by.direction = TRUE)
#' obj2 <- add_neighborhood_centric_summaries(
#'   obj2,
#'   key = "central_cell_only.ac=TRUE"
#' )
#' }
#'
#' @export
add_neighborhood_centric_summaries <- function(
    niches.obj,
    summaries = c(
      "ratio.focal.sum",      "ratio.focal.mean",      "ratio.focal.cmean",
      "log2ratio.focal.sum",  "log2ratio.focal.mean",  "log2ratio.focal.cmean",
      "focal.asymmetry",      "cross.fraction",        "focal.fraction"
    ),
    key       = NULL,
    overwrite = FALSE,
    verbose   = TRUE
) {

  ## Input validation ----------------------------------------------------------
  if (!inherits(niches.obj, "NICHESObject"))
    stop("[add_neighborhood_centric_summaries] 'niches.obj' must be a NICHESObject.",
         call. = FALSE)

  valid.summaries <- c(
    "ratio.focal.sum",      "ratio.focal.mean",      "ratio.focal.cmean",
    "log2ratio.focal.sum",  "log2ratio.focal.mean",  "log2ratio.focal.cmean",
    "focal.asymmetry",      "cross.fraction",        "focal.fraction"
  )
  summaries <- match.arg(summaries, valid.summaries, several.ok = TRUE)

  if (!is.null(key) && (!is.character(key) || length(key) != 1L))
    stop("[add_neighborhood_centric_summaries] 'key' must be a single character string or NULL.",
         call. = FALSE)

  if (!is.logical(overwrite) || length(overwrite) != 1L)
    stop("[add_neighborhood_centric_summaries] 'overwrite' must be a single logical.",
         call. = FALSE)

  ## Resolve key ---------------------------------------------------------------
  if (is.null(key)) {
    efm      <- attr(niches.obj, "edge.filter.mode")
    rem.auto <- attr(niches.obj, "remove.autocrine")
    if (is.null(efm) || is.null(rem.auto))
      stop(
        "[add_neighborhood_centric_summaries] Object is missing 'edge.filter.mode' ",
        "or 'remove.autocrine' attributes. Re-run add_neighborhoods() to update, ",
        "or supply 'key' explicitly.",
        call. = FALSE
      )
    key <- paste0(efm, ".ac=", !rem.auto)
  }

  ## Source slot guard ---------------------------------------------------------
  agg.list <- niches.obj$aggregations$neighborhood.centric.edge.agg
  if (is.null(agg.list))
    stop(
      "[add_neighborhood_centric_summaries] $aggregations$neighborhood.centric.edge.agg ",
      "is NULL. Run add_neighborhood_centric_edge_agg() first.",
      call. = FALSE
    )
  agg <- agg.list[[key]]
  if (is.null(agg))
    stop(sprintf(
      "[add_neighborhood_centric_summaries] Key '%s' not found in ",
      "$aggregations$neighborhood.centric.edge.agg. Available keys: %s.",
      key,
      paste(names(agg.list), collapse = ", ")
    ), call. = FALSE)

  ## Base column requirement checks -------------------------------------------
  agg.cols <- colnames(agg)

  # All summaries in this function require the by.direction columns.
  # Presence of count.focal.in is the proxy for by.direction having been run.
  if (!"count.focal.in" %in% agg.cols)
    stop(
      "[add_neighborhood_centric_summaries] Directional columns absent from key '",
      key, "'. Re-run add_neighborhood_centric_edge_agg() with by.direction = TRUE.",
      call. = FALSE
    )

  # Score-type specific checks for ratio / log2ratio summaries.
  need.focal.sum   <- c("ratio.focal.sum",   "log2ratio.focal.sum")
  need.focal.mean  <- c("ratio.focal.mean",  "log2ratio.focal.mean")
  need.focal.cmean <- c("ratio.focal.cmean", "log2ratio.focal.cmean")

  missing.methods <- character(0)
  if (any(summaries %in% need.focal.sum) &&
      !all(c("score.focal.in.sum",   "score.focal.out.sum")   %in% agg.cols))
    missing.methods <- c(missing.methods, '"sum"')
  if (any(summaries %in% need.focal.mean) &&
      !all(c("score.focal.in.mean",  "score.focal.out.mean")  %in% agg.cols))
    missing.methods <- c(missing.methods, '"mean"')
  if (any(summaries %in% need.focal.cmean) &&
      !all(c("score.focal.in.cmean", "score.focal.out.cmean") %in% agg.cols))
    missing.methods <- c(missing.methods, '"cmean"')

  if (length(missing.methods) > 0L)
    stop(sprintf(
      "[add_neighborhood_centric_summaries] Required score columns missing for key '%s'. ",
      "Re-run add_neighborhood_centric_edge_agg() with by.direction = TRUE and ",
      "methods including: %s.",
      key,
      paste(missing.methods, collapse = ", ")
    ), call. = FALSE)

  ## Map summary names to output column names ---------------------------------
  col.map <- c(
    ratio.focal.sum       = "score.focal.sum.inout.ratio",
    ratio.focal.mean      = "score.focal.mean.inout.ratio",
    ratio.focal.cmean     = "score.focal.cmean.inout.ratio",
    log2ratio.focal.sum   = "score.focal.sum.inout.log2ratio",
    log2ratio.focal.mean  = "score.focal.mean.inout.log2ratio",
    log2ratio.focal.cmean = "score.focal.cmean.inout.log2ratio",
    focal.asymmetry       = "count.focal.asymmetry",
    cross.fraction        = "count.cross.fraction",
    focal.fraction        = "count.focal.fraction"
  )

  ## Resolve which summaries genuinely need computing -------------------------
  if (!overwrite) {
    summaries.todo <- summaries[!col.map[summaries] %in% agg.cols]
    if (length(summaries.todo) == 0L) {
      if (verbose)
        message("[add_neighborhood_centric_summaries] All requested summaries ",
                "already present for key '", key,
                "'. Use overwrite = TRUE to recompute.")
      return(niches.obj)
    }
    if (verbose && length(summaries.todo) < length(summaries))
      message("[add_neighborhood_centric_summaries] Skipping already-computed ",
              "summaries; adding: ",
              paste(summaries.todo, collapse = ", "), ".")
    summaries <- summaries.todo
  }

  if (verbose)
    message(sprintf(
      "[add_neighborhood_centric_summaries] Computing: %s | key: '%s' ...",
      paste(summaries, collapse = ", "), key
    ))

  ## Compute summaries ---------------------------------------------------------

  # -- Focal in/out ratio (receiver / sender among direct focal edges) --------
  # NA when denominator is zero or NA. cmean values are either NA or strictly
  # positive (sparse edge.data), so NA propagates naturally for cmean.

  if ("ratio.focal.sum" %in% summaries)
    agg[, score.focal.sum.inout.ratio := fifelse(
      score.focal.out.sum > 0, score.focal.in.sum / score.focal.out.sum, NA_real_
    )]

  if ("ratio.focal.mean" %in% summaries)
    agg[, score.focal.mean.inout.ratio := fifelse(
      score.focal.out.mean > 0, score.focal.in.mean / score.focal.out.mean, NA_real_
    )]

  if ("ratio.focal.cmean" %in% summaries)
    agg[, score.focal.cmean.inout.ratio := score.focal.in.cmean / score.focal.out.cmean]

  # -- Log2 focal in/out ratio ------------------------------------------------
  # NA when either value is <= 0 or NA.

  if ("log2ratio.focal.sum" %in% summaries)
    agg[, score.focal.sum.inout.log2ratio := fifelse(
      score.focal.in.sum > 0 & score.focal.out.sum > 0,
      log2(score.focal.in.sum / score.focal.out.sum),
      NA_real_
    )]

  if ("log2ratio.focal.mean" %in% summaries)
    agg[, score.focal.mean.inout.log2ratio := fifelse(
      score.focal.in.mean > 0 & score.focal.out.mean > 0,
      log2(score.focal.in.mean / score.focal.out.mean),
      NA_real_
    )]

  if ("log2ratio.focal.cmean" %in% summaries)
    agg[, score.focal.cmean.inout.log2ratio := fifelse(
      !is.na(score.focal.in.cmean) & !is.na(score.focal.out.cmean),
      log2(score.focal.in.cmean / score.focal.out.cmean),
      NA_real_
    )]

  # -- Count-based neighborhood structure summaries ---------------------------

  # focal.asymmetry: bounded [-1, 1] among direct focal edges.
  # NA when focal cell has no direct edges for this LRM (all activity is cross).
  if ("focal.asymmetry" %in% summaries)
    agg[, count.focal.asymmetry := fifelse(
      (count.focal.in + count.focal.out) > 0L,
      (count.focal.in - count.focal.out) / (count.focal.in + count.focal.out),
      NA_real_
    )]

  # cross.fraction: fraction of neighborhood edges not involving focal cell.
  # count.total is always > 0 in a well-formed sparse table, but guarded
  # with fifelse for robustness.
  if ("cross.fraction" %in% summaries)
    agg[, count.cross.fraction := fifelse(
      count.total > 0L,
      count.cross / count.total,
      NA_real_
    )]

  # focal.fraction: fraction of neighborhood edges directly involving focal cell.
  # Note: cross.fraction + focal.fraction = 1 by construction.
  if ("focal.fraction" %in% summaries)
    agg[, count.focal.fraction := fifelse(
      count.total > 0L,
      (count.focal.in + count.focal.out) / count.total,
      NA_real_
    )]

  ## Write back to object ------------------------------------------------------
  niches.obj$aggregations$neighborhood.centric.edge.agg[[key]] <- agg

  if (verbose)
    message(sprintf(
      "[add_neighborhood_centric_summaries] Done. Added %d summary column(s) | key: '%s'.",
      length(summaries), key
    ))

  niches.obj
}
