# ==============================================================================
# FILE:    CellToCell_InSitu.R
# PURPOSE: Batch in-situ plotting of ligand-receptor signals across all or
#          a selected set of LR pairs in an EdgeLRResult object.
#          Calls CellToCell_InSitu_Single_LRM() for each pair.
#
# FUNCTIONS:
#   CellToCell_InSitu()  — main batch plotting function
#
# DEPS: ggplot2, data.table, parallel
# ==============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(data.table)
  library(parallel)
})

# Source CellToCell_InSitu_Single_LRM if not already loaded
if (!exists("CellToCell_InSitu_Single_LRM", mode = "function")) {
  source("~/Documents/GitHub/NICHESv2-NW/working-scripts-NW/Claude/CellToCell_InSitu_Single_LRM.R")
}


# ==============================================================================
#' CellToCell_InSitu
#'
#' @title Batch In-Situ Visualization of Ligand-Receptor Signals
#'
#' @description
#' Iterates over LR pairs and calls \code{CellToCell_InSitu_Single_LRM()} for
#' each, saving one PNG per pair. Two modes:
#'
#' \itemize{
#'   \item \strong{LRM.list = NULL} — plots all LR pairs present in
#'     \code{ctc$signals}. Files are saved to \code{folder} (or the working
#'     directory if \code{folder = NULL}).
#'   \item \strong{LRM.list = data.frame} — plots only the supplied pairs.
#'     Files are saved to a subfolder named \code{"selected_LRMs"} inside
#'     \code{folder} (or inside the working directory if \code{folder = NULL}).
#' }
#'
#' @section LRM.list format:
#' A \code{data.frame} with two columns: \code{ligand} and \code{receptor}.
#' One row per pair. Convert from a character vector with em-dash separators:
#' \preformatted{
#' vec <- c("MIF—CD74", "VEGFA—KDR")
#' lrm_df <- data.frame(
#'   ligand   = sub("\u2014.*", "", vec),
#'   receptor = sub(".*\u2014", "", vec)
#' )
#' }
#'
#' @param ctc          An \code{EdgeLRResult} object from
#'   \code{compute_CellToCell_Spatial()}.
#' @param LRM.list     \code{NULL} (plot all) or a \code{character} vector of
#'   LR pair strings in \code{"Ligand.Receptor"} format (e.g.
#'   \code{c("MIF.CD74", "VEGFA.FLT1")}). When supplied, only these pairs
#'   are plotted and saved to a \code{"selected_LRMs"} subfolder.
#' @param folder       Output directory name. Created if it does not exist.
#'   Set to \code{NULL} to save directly in the working directory.
#' @param bg.color     Background edge color. Default \code{"grey50"}.
#' @param bg.linewidth Background edge line width. Default \code{0.1}.
#' @param bg.alpha     Background edge alpha. Default \code{0.2}.
#' @param sig.linewidth Signal edge line width. Default \code{0.4}.
#' @param sig.alpha    Signal edge alpha. Default \code{0.8}.
#' @param viridis.option Viridis palette option. Default \code{"C"}.
#' @param width        Plot width in inches. Default \code{7}.
#' @param height       Plot height in inches. Default \code{6}.
#' @param dpi          Plot resolution. Default \code{150}.
#' @param n.cores      Number of cores for parallel saving. Default \code{1}.
#'   Fork-based on Unix/macOS; falls back to serial on Windows.
#'
#' @return Invisibly returns a data.frame summarising each pair:
#'   \code{LR}, \code{status} ("saved" or "skipped"), \code{path}.
# ==============================================================================

CellToCell_InSitu <- function(ctc,
                               LRM.list       = NULL,
                               folder         = NULL,
                               bg.color       = "grey50",
                               bg.linewidth   = 0.1,
                               bg.alpha       = 0.2,
                               sig.linewidth  = 0.4,
                               sig.alpha      = 0.8,
                               viridis.option = "C",
                               width          = 7,
                               height         = 6,
                               dpi            = 150,
                               n.cores        = 1L) {

  if (!inherits(ctc, "EdgeLRResult"))
    stop("[CellToCell_InSitu] 'ctc' must be an EdgeLRResult object.")

  n.cores <- as.integer(n.cores)

  # ════════════════════════════════════════════════════════════════════════════
  # SECTION 1 ── Resolve output directory and LR pair list
  # ════════════════════════════════════════════════════════════════════════════

  selected.mode <- !is.null(LRM.list)

  # Determine output directory
  if (selected.mode) {
    # Selected LRMs → always go into "selected_LRMs" subfolder
    out.dir <- if (!is.null(folder)) file.path(folder, "selected_LRMs")
               else "selected_LRMs"
  } else {
    # All LRMs → use folder as-is, or working dir
    out.dir <- if (!is.null(folder)) folder else "."
  }

  if (out.dir != ".") {
    dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)
    message(sprintf("[CellToCell_InSitu] Output directory: %s", out.dir))
  } else {
    message("[CellToCell_InSitu] Saving to working directory.")
  }

  # Resolve which LR pairs to plot
  # Resolve which LR pairs to plot
  if (selected.mode) {

    if (!is.character(LRM.list) || length(LRM.list) == 0)
      stop("[CellToCell_InSitu] 'LRM.list' must be a non-empty character vector of ",
           "'Ligand.Receptor' strings (e.g. c(\"MIF.CD74\", \"VEGFA.FLT1\")).")

    pairs <- data.table(LR = LRM.list)

  } else {

    # All unique LR pairs present in the edge table
    pairs <- unique(ctc$edges[, .(LR)])
  }

  n.pairs <- nrow(pairs)
  message(sprintf("[CellToCell_InSitu] %d LR pair(s) to plot.", n.pairs))

  # ════════════════════════════════════════════════════════════════════════════
  # SECTION 2 ── Define per-pair plotting job
  # ════════════════════════════════════════════════════════════════════════════

  .plot_one <- function(i) {

    lr    <- pairs$LR[i]
    fname <- sprintf("%s.png", lr)
    fpath <- file.path(out.dir, fname)

    result <- tryCatch({

      p <- CellToCell_InSitu_Single_LRM(
        ctc            = ctc,
        LR             = lr,
        bg.color       = bg.color,
        bg.linewidth   = bg.linewidth,
        bg.alpha       = bg.alpha,
        sig.linewidth  = sig.linewidth,
        sig.alpha      = sig.alpha,
        viridis.option = viridis.option
      )

      ggsave(
        filename = fpath,
        plot     = p,
        width    = width,
        height   = height,
        dpi      = dpi,
        bg       = "black"
      )

      list(LR = lr, status = "saved", path = fpath)

    }, error = function(e) {
      message(sprintf("  [SKIP] %s : %s", lr, conditionMessage(e)))
      list(LR = lr, status = "skipped", path = NA_character_)
    })

    result
  }

  # ════════════════════════════════════════════════════════════════════════════
  # SECTION 3 ── Execute (parallel or serial)
  # ════════════════════════════════════════════════════════════════════════════

  message(sprintf("[CellToCell_InSitu] Plotting with %d core(s) ...", n.cores))

  # macOS + ggplot2 crash: fork() is unsafe after Obj-C frameworks initialize,
  # which ggplot2 triggers for font rendering. Sys.setenv() cannot fix this
  # at runtime — the env var must be set before R starts. Serial execution
  # is the reliable solution on macOS; disk I/O dominates anyway for ggsave.
  is.mac   <- Sys.info()[["sysname"]] == "Darwin"
  use.fork <- n.cores > 1L && .Platform$OS.type == "unix" && !is.mac

  if (n.cores > 1L && is.mac) {
    message("  macOS detected: falling back to serial execution for ggsave safety.")
    message("  (fork() + ggplot2 font rendering crashes on macOS at runtime.)")
  }

  if (use.fork) {
    job.results <- parallel::mclapply(
      seq_len(n.pairs), .plot_one,
      mc.cores = n.cores
    )
  } else {
    if (n.cores > 1L && !is.mac && .Platform$OS.type != "unix")
      message("  Fork parallelism unavailable on Windows — running serially.")
    job.results <- lapply(seq_len(n.pairs), .plot_one)
  }

  # ════════════════════════════════════════════════════════════════════════════
  # SECTION 4 ── Summarise and report
  # ════════════════════════════════════════════════════════════════════════════

  summary.df <- rbindlist(lapply(job.results, as.data.table))

  n.saved   <- sum(summary.df$status == "saved")
  n.skipped <- sum(summary.df$status == "skipped")

  message(sprintf(
    "[CellToCell_InSitu] Done. Saved: %d | Skipped: %d | Directory: %s",
    n.saved, n.skipped, normalizePath(out.dir, mustWork = FALSE)
  ))

  invisible(as.data.frame(summary.df))
}
