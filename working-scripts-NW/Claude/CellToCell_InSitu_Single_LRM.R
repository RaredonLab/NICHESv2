# ==============================================================================
# FILE:    CellToCell_InSitu_Single_LRM.R
# PURPOSE: In situ visualization of a single ligand-receptor pair's signaling
#          across all spatial edges in an EdgeLRResult object.
#
# FUNCTIONS:
#   CellToCell_InSitu_Single_LRM()  — main plotting function
#
# DEPS: ggplot2, data.table
# ==============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(data.table)
})


# ==============================================================================
#' CellToCell_InSitu_Single_LRM
#'
#' @title In Situ Visualization of a Single Ligand-Receptor Pair
#'
#' @description
#' Draws a two-layer segment plot showing cell-cell signaling for one
#' ligand-receptor pair across all spatial edges:
#'
#' \enumerate{
#'   \item \strong{Background layer} — all spatial edges rendered as thin grey
#'     segments (low alpha) to show the full neighborhood connectivity.
#'   \item \strong{Signal layer} — the same edges colored by their
#'     ligand-receptor signal value using \code{scale_color_viridis_c()}.
#'     Edges absent from the signal table for this pair are treated as
#'     signal = 0 (darkest end of the viridis palette).
#' }
#'
#' Each segment runs from the sending cell (x, y) to the receiving cell
#' (x, y). Coordinates are read from \code{ctc$cell.meta}, which must
#' contain \code{x} and \code{y} columns before calling this function.
#'
#' @section Coordinates in cell_meta:
#' Spatial coordinates are read from \code{ctc$cell.meta$x} and
#' \code{ctc$cell.meta$y}. The easiest way to ensure they are present is
#' to pass your \code{coord} data.frame when calling
#' \code{computeCellToCellSpatial()}:
#' \preformatted{
#' ctc <- compute_CellToCell_Spatial(
#'   count.mtx = count.mtx,
#'   edge.list  = edge.list,
#'   LRM.list   = LRM.list
#' )
#' # Now ctc$cell.meta already has x and y — no Seurat object needed.
#' }
#' If you already have a \code{ctc} object without coordinates, add them with:
#' \preformatted{
#' library(data.table)
#' setDT(coord.df, keep.rownames = "cell")
#' ctc$cell.meta <- coord.df[ctc$cell.meta, on = "cell"]
#' }
#'
#' @param ctc       An \code{EdgeLRResult} object from
#'   \code{compute_CellToCell_Spatial()}. \code{$cell.meta} must contain
#'   \code{x} and \code{y} columns.
#' @param LR        Character. Ligand-receptor pair in \code{"Ligand.Receptor"}
#'   format, matching a value in \code{ctc$edges$LR} (e.g. \code{"MIF.CD74"}).
#' @param bg.color           Color for background edges. Default \code{"grey40"}.
#' @param bg.linewidth       Line width for background edges. Default \code{0.05}.
#' @param bg.alpha           Alpha for background edges. Default \code{0.1}.
#' @param sig.linewidth      Fixed line width for signal edges. Default \code{0.4}.
#' @param sig.alpha          Alpha for signal edges. Default \code{0.8}.
#' @param viridis.option     Viridis palette option. Default \code{"C"} (plasma).
#'   Any single letter A–H or named option accepted by
#'   \code{scale_color_viridis_c()} works.
#' @param title   Plot title. Defaults to \code{"Ligand.Receptor Connectivity"}.
#'
#' @return A \code{ggplot} object.
# ==============================================================================

CellToCell_InSitu_Single_LRM <- function(ctc,
                                          LR,
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

  if (!inherits(ctc, "EdgeLRResult"))
    stop("[CellToCell_InSitu_Single_LRM] 'ctc' must be an EdgeLRResult object.")

  if (!all(c("x", "y") %in% colnames(ctc$cell.meta))) {
    stop(
      "[CellToCell_InSitu_Single_LRM] 'ctc$cell.meta' must contain 'x' and 'y' columns.\n",
      "  Ensure create_Edgelist_Spatial() was used to build the edge list,\n",
      "  which embeds coordinates automatically into the ctc object."
    )
  }

  # Validate LR format: must be "Ligand.Receptor"
  if (!is.character(LR) || length(LR) != 1 || !grepl("\\.", LR))
    stop("[CellToCell_InSitu_Single_LRM] 'LR' must be a single string in 'Ligand.Receptor' format (e.g. 'MIF.CD74').")

  # Check the LR pair exists in the signals table
  if (!"LR" %in% colnames(ctc$edges))
    stop("[CellToCell_InSitu_Single_LRM] 'ctc$edges' does not have an 'LR' column.\n",
         "  Re-run compute_CellToCell_Spatial() to generate a compatible ctc object.")

  valid.lrs <- unique(ctc$edges$LR)
  if (!LR %in% valid.lrs)
    stop(sprintf(
      "[CellToCell_InSitu_Single_LRM] LR pair '%s' not found in ctc$edges$LR.\n  Available pairs: %s%s",
      LR,
      paste(head(valid.lrs, 5), collapse = ", "),
      if (length(valid.lrs) > 5) sprintf(" ... (%d total)", length(valid.lrs)) else ""
    ))

  # ════════════════════════════════════════════════════════════════════════════
  # SECTION 2 ── Build coordinate lookup table
  #
  # Use only the columns we need to keep joins fast.
  # ════════════════════════════════════════════════════════════════════════════

  coords <- ctc$cell.meta[, .(cell, x, y)]

  # ════════════════════════════════════════════════════════════════════════════
  # SECTION 3 ── Build background edge table
  #
  # All unique (from, to) pairs — the full spatial neighborhood graph.
  # Exclude autocrine edges (from == to) which would be zero-length segments.
  # ════════════════════════════════════════════════════════════════════════════

  all.edges <- unique(ctc$edges[from != to, .(from, to)])

  # Join sending cell coordinates
  all.edges <- coords[all.edges, on = .(cell = from)]
  setnames(all.edges, c("x", "y"), c("x.from", "y.from"))
  setnames(all.edges, "cell", "from")

  # Join receiving cell coordinates
  all.edges <- coords[all.edges, on = .(cell = to)]
  setnames(all.edges, c("x", "y"), c("x.to", "y.to"))
  setnames(all.edges, "cell", "to")

  message(sprintf(
    "[CellToCell_InSitu_Single_LRM] %d total edges | LR: %s",
    nrow(all.edges), LR
  ))

  # ════════════════════════════════════════════════════════════════════════════
  # SECTION 4 ── Build signal edge table (edge values for this LR)
  #
  # Filter edges for the requested LR pair, then left-join onto all.edges
  # so that every edge is represented. Edges absent for this pair get edge = 0.
  # ════════════════════════════════════════════════════════════════════════════

  # Capture LR in local variable to avoid data.table scoping issues with `==`
  lr.val <- LR

  lr.edges <- ctc$edges[LR == lr.val, .(from, to, edge)]

  # Left join: all.edges is the left table, so every edge is kept
  sig.edges <- lr.edges[all.edges, on = .(from, to)]
  sig.edges[is.na(edge), edge := 0]

  n.nonzero <- sum(sig.edges$edge > 0)
  message(sprintf(
    "  Non-zero signal edges: %d / %d  (%.1f%%)",
    n.nonzero, nrow(sig.edges), 100 * n.nonzero / nrow(sig.edges)
  ))

  # ════════════════════════════════════════════════════════════════════════════
  # SECTION 5 ── Build plot
  # ════════════════════════════════════════════════════════════════════════════

  # Title matches NICHES style: "Ligand.Receptor Connectivity"
  plot.title <- if (!is.null(title)) title else
    sprintf("%s Connectivity", LR)

  p <- ggplot() +

    # Layer 1: background — all edges, fixed grey, thin, low alpha
    geom_segment(
      data      = all.edges,
      aes(x = x.from, y = y.from, xend = x.to, yend = y.to),
      color     = bg.color,
      linewidth = bg.linewidth,
      alpha     = bg.alpha,
      lineend   = "round"
    ) +

    # Layer 2: edge signal — non-zero edges only, colored by edge value
    geom_segment(
      data      = sig.edges[edge > 0],
      aes(x = x.from, y = y.from, xend = x.to, yend = y.to,
          color = edge),
      linewidth = sig.linewidth,
      alpha     = sig.alpha,
      lineend   = "round"
    ) +

    # Color scale — single "Signal" legend
    scale_color_viridis_c(
      option = viridis.option,
      name   = "Edge",
      guide  = guide_colorbar(
        barwidth       = 0.6,
        barheight      = 8,
        ticks          = FALSE,
        title.position = "top"
      )
    ) +

    # Axis labels in sending-cell coordinate style
    labs(
      title = plot.title,
      x     = "x.Sending",
      y     = "y.Sending"
    ) +

    # Dark theme
    theme_dark() +
    theme(
      plot.background  = element_rect(fill = "black", color = NA),
      panel.background = element_rect(fill = "black", color = NA),
      panel.grid       = element_blank(),
      axis.line        = element_line(color = "white", linewidth = 0.3),
      axis.ticks       = element_line(color = "white", linewidth = 0.3),
      axis.text        = element_text(color = "white", size = 7),
      axis.title       = element_text(color = "white", size = 9),
      plot.title       = element_text(color = "white", size = 12,
                                      face = "bold", hjust = 0.5,
                                      margin = margin(b = 8)),
      legend.background = element_rect(fill = "black", color = NA),
      legend.title      = element_text(color = "white", size = 8, hjust = 0.5),
      legend.text       = element_text(color = "grey70", size = 7),
      legend.position   = "right",
      plot.margin       = margin(12, 12, 12, 12)
    ) +

    coord_equal()   # preserve spatial aspect ratio

  p
}
