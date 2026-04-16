# CellToCell_InSitu_Single_LRM.R ------------------------------------------
# In situ visualisation of a single LRM pair's signalling across all
# spatial edges stored in a NICHESObject.
#
# FUNCTIONS:
#   CellToCell_InSitu_Single_LRM()  -- main plotting function
#
# DEPS: ggplot2, data.table, scales


# ==============================================================================
#' In Situ Visualisation of a Single Ligand-Receptor Pair
#'
#' @description
#' Draws a two-layer segment plot for one LRM pair across all spatial edges
#' in a \code{NICHESObject}:
#'
#' \enumerate{
#'   \item \strong{Background layer} -- all spatial edges, grey, thin, low
#'     alpha. Provides spatial context.
#'   \item \strong{Signal layer} -- edges with non-zero score for the
#'     requested LRM, coloured and alpha-scaled by score via
#'     \code{scale_color_viridis_c()}.
#' }
#'
#' Coordinates are read from \code{$cell.meta} (columns \code{x} and
#' \code{y}). Edge endpoints are resolved by joining \code{$edge.list}
#' against \code{$cell.meta}. Returns a \code{ggplot} object that can be
#' further customised with additional layers, themes, or labels before
#' being printed or saved.
#'
#' @section Direct vs batch use:
#' Call this function directly to visualise a single LRM pair and obtain
#' the live \code{ggplot} object -- useful when you want to add custom
#' layers (annotations, highlighting, faceting), modify the theme, or
#' embed the plot in a larger figure (e.g.\ via \pkg{patchwork} or
#' \pkg{cowplot}). Use \code{\link{CellToCell_InSitu}} when you need to
#' generate many PNGs in one batch -- it iterates over LRM pairs, writes
#' files to disk, and parallelises across cores.
#'
#' @param niches.obj     A \code{NICHESObject} with \code{x} and \code{y}
#'   columns in \code{$cell.meta}.
#' @param LRM            Character. LRM pair in \code{"ligand|receptor"}
#'   format, matching a value in \code{niches.obj$LRM.list}
#'   (e.g.\ \code{"MIF|CD74"}).
#' @param value.col      Character. Score column to map to the colour and
#'   alpha aesthetics. \code{"score"} (default) shows raw signal magnitude
#'   and is the right choice for absolute comparisons. \code{"score.norm"}
#'   shows per-edge proportional signal (sums to 1 across LRMs within an
#'   edge) and is the right choice for highlighting which LRM dominates
#'   the communication budget at each edge.
#' @param scale.upper    Numeric or \code{NULL}. Upper limit of the colour
#'   and alpha scale; values above are squished to the top colour.
#'   \code{NULL} (default) auto-computes the 90th percentile of non-zero
#'   values for this LRM, which produces a visually balanced individual
#'   plot but means colour values are \strong{not comparable} across
#'   plots. To compare across LRMs, set this manually to a fixed numeric
#'   value.
#' @param bg.color       Background edge colour. Default \code{"grey40"}.
#' @param bg.linewidth   Background edge line width. Default \code{0.05}.
#' @param bg.alpha       Background edge alpha. Default \code{0.1}.
#' @param sig.linewidth  Signal edge line width. Default \code{0.4}.
#' @param sig.alpha      Maximum signal edge alpha (applied to the
#'   highest-value edges). Lower-value edges scale down to near-zero
#'   alpha. Default \code{0.8}.
#' @param viridis.option Viridis palette letter (\code{"A"}-\code{"H"}).
#'   Default \code{"C"} (plasma).
#' @param title          Plot title. Default \code{NULL}, which produces
#'   \code{"ligand|receptor Connectivity"}.
#' @param precomputed.bg Optional pre-built background edge
#'   \code{data.table} (with columns \code{edge}, \code{SendingCell},
#'   \code{ReceivingCell}, \code{x.from}, \code{y.from}, \code{x.to},
#'   \code{y.to}). Used internally by \code{\link{CellToCell_InSitu}} to
#'   skip the per-call coordinate join during batch plotting. End users
#'   calling this function directly should leave this \code{NULL} (the
#'   default).
#'
#' @return A \code{ggplot} object. The plot is fully composed (data,
#'   scales, theme, fixed coordinate ratio) and is suitable for direct
#'   printing, saving via \code{ggplot2::ggsave()}, or further
#'   modification with the standard \code{+} layer-addition syntax.
#'
#' @seealso
#'   \code{\link{CellToCell_InSitu}} for batch plotting of many LRM pairs
#'     in one call;
#'   \code{\link{NICHESv2_example}} for an example object suitable for
#'     plotting;
#'   \code{\link{Aggregate_CellToCell_Edges}} for non-graphical per-cell
#'     summaries of the same edge data.
#'
#' @examples
#' \donttest{
#' data(NICHESv2_example)
#'
#' # Pick an LRM pair to visualise (the first one in $LRM.list)
#' lrm.to.plot <- NICHESv2_example$LRM.list[1]
#'
#' p <- CellToCell_InSitu_Single_LRM(
#'   niches.obj = NICHESv2_example,
#'   LRM        = lrm.to.plot,
#'   value.col  = "score"
#' )
#' print(p)
#'
#' # The returned ggplot can be customised further with the standard
#' # ggplot2 '+' syntax -- add layers, modify the title, change the theme:
#' p + ggplot2::labs(subtitle = "Customised after the call")
#' }
#'
#' @export
CellToCell_InSitu_Single_LRM <- function(niches.obj,
                                          LRM,
                                          value.col      = "score",
                                          scale.upper    = NULL,
                                          bg.color       = "grey40",
                                          bg.linewidth   = 0.05,
                                          bg.alpha       = 0.1,
                                          sig.linewidth  = 0.4,
                                          sig.alpha      = 0.8,
                                          viridis.option = "C",
                                          title          = NULL,
                                          precomputed.bg = NULL) {

  ## SECTION 1: Input validation ----

  if (!is_NICHESObject(niches.obj))
    stop("[CellToCell_InSitu_Single_LRM] 'niches.obj' must be a NICHESObject.")

  if (!all(c("x", "y") %in% colnames(niches.obj$cell.meta)))
    stop(
      "[CellToCell_InSitu_Single_LRM] '$cell.meta' must contain 'x' and 'y'.\n",
      "  Ensure 'meta.data' had 'x' and 'y' when create_NICHESObject() was called."
    )

  if (!is.character(LRM) || length(LRM) != 1L || !grepl("|", LRM, fixed = TRUE))
    stop("[CellToCell_InSitu_Single_LRM] 'LRM' must be a single string in ",
         "'ligand|receptor' format (e.g. 'MIF|CD74').")

  if (!LRM %in% niches.obj$LRM.list) {
    valid.lrms <- head(niches.obj$LRM.list, 5L)
    stop(sprintf(
      "[CellToCell_InSitu_Single_LRM] LRM '%s' not found in $LRM.list.\n  First 5 available: %s%s",
      LRM,
      paste(valid.lrms, collapse = ", "),
      if (length(niches.obj$LRM.list) > 5L)
        sprintf(" ... (%d total)", length(niches.obj$LRM.list)) else ""
    ))
  }

  if (!value.col %in% c("score", "score.norm"))
    stop("[CellToCell_InSitu_Single_LRM] 'value.col' must be \"score\" or \"score.norm\".")

  if (!is.null(scale.upper) &&
      (!is.numeric(scale.upper) || length(scale.upper) != 1L ||
       !is.finite(scale.upper) || scale.upper <= 0))
    stop("[CellToCell_InSitu_Single_LRM] 'scale.upper' must be a single positive finite number, or NULL.")

  ## SECTION 2: Build coordinate lookup table ----

  coords <- niches.obj$cell.meta[, .(cell, x, y)]

  ## SECTION 3: Build background edge table ----
  #
  # All unique inter-cell edges from $edge.list, with coordinates attached.
  # Autocrine edges (SendingCell == ReceivingCell) are excluded from plotting.
  #
  # When called from CellToCell_InSitu() in batch mode, precomputed.bg is
  # supplied and this section is skipped entirely -- the two coordinate joins
  # are O(n_edges) and identical across all LRMs, so pre-computing them once
  # in the caller and reusing here gives a large speedup for batch jobs.

  if (!is.null(precomputed.bg)) {
    all.edges <- precomputed.bg
  } else {
    all.edges <- niches.obj$edge.list[SendingCell != ReceivingCell]

    # Join sending coordinates
    all.edges <- coords[all.edges, on = .(cell = SendingCell)]
    setnames(all.edges, c("x", "y"), c("x.from", "y.from"))
    setnames(all.edges, "cell", "SendingCell")

    # Join receiving coordinates
    all.edges <- coords[all.edges, on = .(cell = ReceivingCell)]
    setnames(all.edges, c("x", "y"), c("x.to", "y.to"))
    setnames(all.edges, "cell", "ReceivingCell")
  }

  message(sprintf(
    "[CellToCell_InSitu_Single_LRM] %d total edges | LRM: %s",
    nrow(all.edges), LRM
  ))

  ## SECTION 4: Build signal edge table for this LRM ----
  #
  # Pull both score columns from $edge.data, left-join onto all.edges so every
  # edge appears, fill NA (absent for this LRM) with 0, then derive .plot.val
  # from the requested value.col.
  #
  # scale.upper is resolved here from non-zero .plot.val values if not supplied.

  lrm.val  <- LRM   # local capture for data.table scoping
  lr.edges <- niches.obj$edge.data[LRM == lrm.val, .(edge, score, score.norm)]

  sig.edges <- lr.edges[all.edges, on = "edge"]
  sig.edges[is.na(score),      score      := 0]
  sig.edges[is.na(score.norm), score.norm := 0]

  # Derive the plotting value column
  sig.edges[, .plot.val := get(value.col)]

  n.nonzero <- sum(sig.edges$.plot.val > 0)
  message(sprintf(
    "  Non-zero signal edges: %d / %d  (%.1f%%)",
    n.nonzero, nrow(sig.edges), 100 * n.nonzero / nrow(sig.edges)
  ))

  # Resolve upper scale limit from non-zero values only
  nonzero.vals <- sig.edges$.plot.val[sig.edges$.plot.val > 0]
  if (is.null(scale.upper)) {
    scale.upper <- if (length(nonzero.vals) > 0L) {
      as.numeric(quantile(nonzero.vals, 0.9))
    } else {
      1   # fallback: no signal at all for this LRM
    }
  }

  message(sprintf("  Scale upper limit: %.4g  (value.col = '%s')", scale.upper, value.col))

  ## SECTION 5: Build plot ----

  plot.title  <- if (!is.null(title)) title else sprintf("%s Connectivity", LRM)
  legend.name <- if (value.col == "score") "Score" else "Score\n(norm)"

  p <- ggplot() +

    # Layer 1: background -- all edges, fixed grey
    geom_segment(
      data      = all.edges,
      aes(x = x.from, y = y.from, xend = x.to, yend = y.to),
      color     = bg.color,
      linewidth = bg.linewidth,
      alpha     = bg.alpha,
      lineend   = "round"
    ) +

    # Layer 2: signal -- non-zero edges coloured AND alpha-scaled by .plot.val
    geom_segment(
      data      = sig.edges[.plot.val > 0],
      aes(x = x.from, y = y.from, xend = x.to, yend = y.to,
          color = .plot.val,
          alpha = .plot.val),
      linewidth = sig.linewidth,
      lineend   = "round"
    ) +

    scale_color_viridis_c(
      option = viridis.option,
      name   = legend.name,
      limits = c(0, scale.upper),
      oob    = scales::squish,
      guide  = guide_colorbar(
        barwidth       = 0.6,
        barheight      = 8,
        ticks          = FALSE,
        title.position = "top"
      )
    ) +

    # Alpha mirrors the colour scale; legend suppressed (redundant with colorbar)
    scale_alpha_continuous(
      range  = c(0.05, sig.alpha),
      limits = c(0, scale.upper),
      oob    = scales::squish,
      guide  = "none"
    ) +

    labs(title = plot.title, x = "x.Sending", y = "y.Sending") +

    theme_dark() +
    theme(
      plot.background   = element_rect(fill = "black", color = NA),
      panel.background  = element_rect(fill = "black", color = NA),
      panel.grid        = element_blank(),
      axis.line         = element_line(color = "white", linewidth = 0.3),
      axis.ticks        = element_line(color = "white", linewidth = 0.3),
      axis.text         = element_text(color = "white", size = 7),
      axis.title        = element_text(color = "white", size = 9),
      plot.title        = element_text(color = "white", size = 12,
                                       face = "bold", hjust = 0.5,
                                       margin = margin(b = 8)),
      legend.background = element_rect(fill = "black", color = NA),
      legend.title      = element_text(color = "white", size = 8, hjust = 0.5),
      legend.text       = element_text(color = "grey70", size = 7),
      legend.position   = "right",
      plot.margin       = margin(12, 12, 12, 12)
    ) +

    coord_equal()

  p
}
