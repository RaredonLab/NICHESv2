# ==============================================================================
# FILE:    CellToCell_InSitu_Single_LRM.R
# PURPOSE: In situ visualisation of a single LRM pair's signalling across all
#          spatial edges stored in a NICHESObject.
#
# FUNCTIONS:
#   CellToCell_InSitu_Single_LRM()  — main plotting function
#
# DEPS: ggplot2, data.table
# ==============================================================================


# ==============================================================================
#' In Situ Visualisation of a Single Ligand-Receptor Pair
#'
#' @description
#' Draws a two-layer segment plot for one LRM pair across all spatial edges
#' in a \code{NICHESObject}:
#'
#' \enumerate{
#'   \item \strong{Background layer} — all spatial edges, grey, thin, low alpha.
#'   \item \strong{Signal layer} — edges with non-zero score for the requested
#'     LRM, coloured by score via \code{scale_color_viridis_c()}.
#' }
#'
#' Coordinates are read from \code{$cell.meta} (columns \code{x} and \code{y}).
#' Edge endpoints are resolved by joining \code{$edge.list} against
#' \code{$cell.meta}.
#'
#' @param niches.obj     A \code{NICHESObject} with \code{x} and \code{y} in
#'   \code{$cell.meta}.
#' @param LRM            Character. LRM pair in \code{"ligand|receptor"} format,
#'   matching a value in \code{$LRM.list} (e.g. \code{"MIF|CD74"}).
#' @param bg.color       Background edge colour. Default \code{"grey40"}.
#' @param bg.linewidth   Background edge line width. Default \code{0.05}.
#' @param bg.alpha       Background edge alpha. Default \code{0.1}.
#' @param sig.linewidth  Signal edge line width. Default \code{0.4}.
#' @param sig.alpha      Signal edge alpha. Default \code{0.8}.
#' @param viridis.option Viridis palette letter (\code{"A"}–\code{"H"}).
#'   Default \code{"C"} (plasma).
#' @param title          Plot title. Defaults to \code{"ligand|receptor Connectivity"}.
#'
#' @return A \code{ggplot} object.
#'
#' @seealso \code{\link{CellToCell_InSitu}}
#'
#' @export
CellToCell_InSitu_Single_LRM <- function(niches.obj,
                                          LRM,
                                          bg.color       = "grey40",
                                          bg.linewidth   = 0.05,
                                          bg.alpha       = 0.1,
                                          sig.linewidth  = 0.4,
                                          sig.alpha      = 0.8,
                                          viridis.option = "C",
                                          title          = NULL) {

  # ════════════════════════════════════════════════════════════════════════════
  # SECTION 1 ── Input validation
  # ════════════════════════════════════════════════════════════════════════════

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

  # ════════════════════════════════════════════════════════════════════════════
  # SECTION 2 ── Build coordinate lookup table
  # ════════════════════════════════════════════════════════════════════════════

  coords <- niches.obj$cell.meta[, .(cell, x, y)]

  # ════════════════════════════════════════════════════════════════════════════
  # SECTION 3 ── Build background edge table
  #
  # All unique inter-cell edges from $edge.list, with coordinates attached.
  # Autocrine edges (SendingCell == ReceivingCell) are excluded from plotting.
  # ════════════════════════════════════════════════════════════════════════════

  all.edges <- niches.obj$edge.list[SendingCell != ReceivingCell]

  # Join sending coordinates
  all.edges <- coords[all.edges, on = .(cell = SendingCell)]
  setnames(all.edges, c("x", "y"), c("x.from", "y.from"))
  setnames(all.edges, "cell", "SendingCell")

  # Join receiving coordinates
  all.edges <- coords[all.edges, on = .(cell = ReceivingCell)]
  setnames(all.edges, c("x", "y"), c("x.to", "y.to"))
  setnames(all.edges, "cell", "ReceivingCell")

  message(sprintf(
    "[CellToCell_InSitu_Single_LRM] %d total edges | LRM: %s",
    nrow(all.edges), LRM
  ))

  # ════════════════════════════════════════════════════════════════════════════
  # SECTION 4 ── Build signal edge table for this LRM
  #
  # Filter $edge.data for the requested LRM, left-join onto all.edges so
  # every edge appears; edges absent for this LRM get score = 0.
  # ════════════════════════════════════════════════════════════════════════════

  lrm.val   <- LRM   # local capture for data.table scoping
  lr.edges  <- niches.obj$edge.data[LRM == lrm.val, .(edge, score)]

  sig.edges <- lr.edges[all.edges, on = "edge"]
  sig.edges[is.na(score), score := 0]

  n.nonzero <- sum(sig.edges$score > 0)
  message(sprintf(
    "  Non-zero signal edges: %d / %d  (%.1f%%)",
    n.nonzero, nrow(sig.edges), 100 * n.nonzero / nrow(sig.edges)
  ))

  # ════════════════════════════════════════════════════════════════════════════
  # SECTION 5 ── Build plot
  # ════════════════════════════════════════════════════════════════════════════

  plot.title <- if (!is.null(title)) title else
    sprintf("%s Connectivity", LRM)

  p <- ggplot() +

    # Layer 1: background — all edges, fixed grey
    geom_segment(
      data      = all.edges,
      aes(x = x.from, y = y.from, xend = x.to, yend = y.to),
      color     = bg.color,
      linewidth = bg.linewidth,
      alpha     = bg.alpha,
      lineend   = "round"
    ) +

    # Layer 2: signal — non-zero edges coloured by score
    geom_segment(
      data      = sig.edges[score > 0],
      aes(x = x.from, y = y.from, xend = x.to, yend = y.to,
          color = score),
      linewidth = sig.linewidth,
      alpha     = sig.alpha,
      lineend   = "round"
    ) +

    scale_color_viridis_c(
      option = viridis.option,
      name   = "Score",
      guide  = guide_colorbar(
        barwidth       = 0.6,
        barheight      = 8,
        ticks          = FALSE,
        title.position = "top"
      )
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
