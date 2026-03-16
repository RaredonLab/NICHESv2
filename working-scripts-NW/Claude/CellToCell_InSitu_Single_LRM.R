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
#' (x, y). Coordinates are read from \code{ctc$cell_meta}, which must
#' contain \code{x} and \code{y} columns before calling this function.
#'
#' @section Coordinates in cell_meta:
#' Spatial coordinates are read from \code{ctc$cell_meta$x} and
#' \code{ctc$cell_meta$y}. The easiest way to ensure they are present is
#' to pass your \code{coord} data.frame when calling
#' \code{computeCellToCellSpatial()}:
#' \preformatted{
#' ctc <- computeCellToCellSpatial(
#'   count.mtx = count_mtx,
#'   edge.list  = edge_list,
#'   LRM.list   = LRM.database,
#'   coord      = coord_df      # <-- same coord_df used in createEdgelistSpatial()
#' )
#' # Now ctc$cell_meta already has x and y — no Seurat object needed.
#' }
#' If you already have a \code{ctc} object without coordinates, add them with:
#' \preformatted{
#' library(data.table)
#' setDT(coord_df, keep.rownames = "cell")
#' ctc$cell_meta <- coord_df[ctc$cell_meta, on = "cell"]
#' }
#'
#' @param ctc       An \code{EdgeLRResult} object from
#'   \code{computeCellToCellSpatial()}. \code{$cell_meta} must contain
#'   \code{x} and \code{y} columns.
#' @param ligand    Character. Ligand gene symbol (must match a value in
#'   \code{ctc$signals$ligand}).
#' @param receptor  Character. Receptor gene symbol (must match a value in
#'   \code{ctc$signals$receptor}).
#' @param bg_color           Color for background edges. Default \code{"grey40"}.
#' @param bg_linewidth       Line width for background edges. Default \code{0.05}.
#' @param bg_alpha           Alpha for background edges. Default \code{0.1}.
#' @param sig_linewidth      Fixed line width for signal edges. Default \code{0.4}.
#' @param sig_alpha          Alpha for signal edges. Default \code{0.8}.
#' @param viridis_option     Viridis palette option. Default \code{"C"} (plasma).
#'   Any single letter A–H or named option accepted by
#'   \code{scale_color_viridis_c()} works.
#' @param title   Plot title. Defaults to \code{"Ligand.Receptor Connectivity"}.
#'
#' @return A \code{ggplot} object.
# ==============================================================================

CellToCell_InSitu_Single_LRM <- function(ctc,
                                          ligand,
                                          receptor,
                                          bg_color       = "grey40",
                                          bg_linewidth   = 0.05,
                                          bg_alpha       = 0.1,
                                          sig_linewidth  = 0.4,
                                          sig_alpha      = 0.8,
                                          viridis_option = "C",
                                          title          = NULL) {

  # ════════════════════════════════════════════════════════════════════════════
  # SECTION 1 ── Input validation
  # ════════════════════════════════════════════════════════════════════════════

  if (!inherits(ctc, "EdgeLRResult"))
    stop("[CellToCell_InSitu_Single_LRM] 'ctc' must be an EdgeLRResult object.")

  if (!all(c("x", "y") %in% colnames(ctc$cell_meta))) {
    stop(
      "[CellToCell_InSitu_Single_LRM] 'ctc$cell_meta' must contain 'x' and 'y' columns.\n",
      "  Pass 'coord' when computing ctc:\n",
      "    ctc <- computeCellToCellSpatial(..., coord = coord_df)\n",
      "  Or add coordinates to an existing ctc object:\n",
      "    setDT(coord_df, keep.rownames = 'cell')\n",
      "    ctc$cell_meta <- coord_df[ctc$cell_meta, on = 'cell']"
    )
  }

  if (!ligand %in% ctc$signals$ligand)
    stop(sprintf("[CellToCell_InSitu_Single_LRM] Ligand '%s' not found in ctc$signals.", ligand))

  if (!receptor %in% ctc$signals$receptor)
    stop(sprintf("[CellToCell_InSitu_Single_LRM] Receptor '%s' not found in ctc$signals.", receptor))

  # Check this specific pair exists
  pair_exists <- ctc$signals[ligand == get("ligand") & receptor == get("receptor"), .N] > 0
  if (!pair_exists)
    stop(sprintf(
      "[CellToCell_InSitu_Single_LRM] Pair '%s' → '%s' not found in ctc$signals.",
      ligand, receptor
    ))

  # ════════════════════════════════════════════════════════════════════════════
  # SECTION 2 ── Build coordinate lookup table
  #
  # Use only the columns we need to keep joins fast.
  # ════════════════════════════════════════════════════════════════════════════

  coords <- ctc$cell_meta[, .(cell, x, y)]

  # ════════════════════════════════════════════════════════════════════════════
  # SECTION 3 ── Build background edge table
  #
  # All unique (from, to) pairs — the full spatial neighborhood graph.
  # Exclude autocrine edges (from == to) which would be zero-length segments.
  # ════════════════════════════════════════════════════════════════════════════

  all_edges <- unique(ctc$signals[from != to, .(from, to)])

  # Join sending cell coordinates
  all_edges <- coords[all_edges, on = .(cell = from)]
  setnames(all_edges, c("x", "y"), c("x_from", "y_from"))
  setnames(all_edges, "cell", "from")

  # Join receiving cell coordinates
  all_edges <- coords[all_edges, on = .(cell = to)]
  setnames(all_edges, c("x", "y"), c("x_to", "y_to"))
  setnames(all_edges, "cell", "to")

  message(sprintf(
    "[CellToCell_InSitu_Single_LRM] %d total edges | Ligand: %s | Receptor: %s",
    nrow(all_edges), ligand, receptor
  ))

  # ════════════════════════════════════════════════════════════════════════════
  # SECTION 4 ── Build signal edge table
  #
  # Filter signals for the requested LR pair, then left-join onto all_edges
  # so that every edge is represented. Edges absent from the signal table
  # for this pair get signal = 0.
  # ════════════════════════════════════════════════════════════════════════════

  # Capture in local variables to avoid data.table scoping issues with `==`
  lig <- ligand
  rec <- receptor

  lr_signals <- ctc$signals[ligand == lig & receptor == rec, .(from, to, signal)]

  # Left join: all_edges is the left table, so every edge is kept
  sig_edges <- lr_signals[all_edges, on = .(from, to)]
  sig_edges[is.na(signal), signal := 0]

  n_nonzero <- sum(sig_edges$signal > 0)
  message(sprintf(
    "  Non-zero signal edges: %d / %d  (%.1f%%)",
    n_nonzero, nrow(sig_edges), 100 * n_nonzero / nrow(sig_edges)
  ))

  # ════════════════════════════════════════════════════════════════════════════
  # SECTION 5 ── Build plot
  # ════════════════════════════════════════════════════════════════════════════

  # Title matches NICHES style: "Ligand.Receptor Connectivity"
  plot_title <- if (!is.null(title)) title else
    sprintf("%s.%s Connectivity", ligand, receptor)

  p <- ggplot() +

    # Layer 1: background — all edges, fixed grey, thin, low alpha
    geom_segment(
      data      = all_edges,
      aes(x = x_from, y = y_from, xend = x_to, yend = y_to),
      color     = bg_color,
      linewidth = bg_linewidth,
      alpha     = bg_alpha,
      lineend   = "round"
    ) +

    # Layer 2: signal — non-zero edges only, colored by LR signal value
    geom_segment(
      data      = sig_edges[signal > 0],
      aes(x = x_from, y = y_from, xend = x_to, yend = y_to,
          color = signal),
      linewidth = sig_linewidth,
      alpha     = sig_alpha,
      lineend   = "round"
    ) +

    # Color scale — single "Signal" legend
    scale_color_viridis_c(
      option = viridis_option,
      name   = "Signal",
      guide  = guide_colorbar(
        barwidth       = 0.6,
        barheight      = 8,
        ticks          = FALSE,
        title.position = "top"
      )
    ) +

    # Axis labels in sending-cell coordinate style
    labs(
      title = plot_title,
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
