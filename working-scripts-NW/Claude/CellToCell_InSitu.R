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
#'   \code{computeCellToCellSpatial()}.
#' @param LRM.list     \code{NULL} (plot all) or a \code{data.frame} with
#'   columns \code{ligand} and \code{receptor} (plot selected pairs only).
#' @param folder       Output directory name. Created if it does not exist.
#'   Set to \code{NULL} to save directly in the working directory.
#' @param bg_color     Background edge color. Default \code{"grey50"}.
#' @param bg_linewidth Background edge line width. Default \code{0.1}.
#' @param bg_alpha     Background edge alpha. Default \code{0.2}.
#' @param sig_linewidth Signal edge line width. Default \code{0.4}.
#' @param sig_alpha    Signal edge alpha. Default \code{0.8}.
#' @param viridis_option Viridis palette option. Default \code{"C"}.
#' @param width        Plot width in inches. Default \code{7}.
#' @param height       Plot height in inches. Default \code{6}.
#' @param dpi          Plot resolution. Default \code{150}.
#' @param n_cores      Number of cores for parallel saving. Default \code{1}.
#'   Fork-based on Unix/macOS; falls back to serial on Windows.
#'
#' @return Invisibly returns a data.frame summarising each pair:
#'   \code{ligand}, \code{receptor}, \code{status} ("saved" or "skipped"),
#'   \code{path}.
# ==============================================================================

CellToCell_InSitu <- function(ctc,
                               LRM.list      = NULL,
                               folder        = NULL,
                               bg_color      = "grey50",
                               bg_linewidth  = 0.1,
                               bg_alpha      = 0.2,
                               sig_linewidth = 0.4,
                               sig_alpha     = 0.8,
                               viridis_option = "C",
                               width         = 7,
                               height        = 6,
                               dpi           = 150,
                               n_cores       = 1L) {

  if (!inherits(ctc, "EdgeLRResult"))
    stop("[CellToCell_InSitu] 'ctc' must be an EdgeLRResult object.")

  n_cores <- as.integer(n_cores)

  # ════════════════════════════════════════════════════════════════════════════
  # SECTION 1 ── Resolve output directory and LR pair list
  # ════════════════════════════════════════════════════════════════════════════

  selected_mode <- !is.null(LRM.list)

  # Determine output directory
  if (selected_mode) {
    # Selected LRMs → always go into "selected_LRMs" subfolder
    out_dir <- if (!is.null(folder)) file.path(folder, "selected_LRMs")
               else "selected_LRMs"
  } else {
    # All LRMs → use folder as-is, or working dir
    out_dir <- if (!is.null(folder)) folder else "."
  }

  if (out_dir != ".") {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    message(sprintf("[CellToCell_InSitu] Output directory: %s", out_dir))
  } else {
    message("[CellToCell_InSitu] Saving to working directory.")
  }

  # Resolve which LR pairs to plot
  if (selected_mode) {

    if (!is.data.frame(LRM.list) ||
        !all(c("ligand", "receptor") %in% colnames(LRM.list)))
      stop("[CellToCell_InSitu] 'LRM.list' must be a data.frame with columns ",
           "'ligand' and 'receptor'.")

    pairs <- as.data.table(LRM.list)[, .(ligand, receptor)]

  } else {

    # All unique LR pairs that have at least one non-zero signal
    pairs <- unique(ctc$signals[, .(ligand, receptor)])
  }

  n_pairs <- nrow(pairs)
  message(sprintf("[CellToCell_InSitu] %d LR pair(s) to plot.", n_pairs))

  # ════════════════════════════════════════════════════════════════════════════
  # SECTION 2 ── Define per-pair plotting job
  # ════════════════════════════════════════════════════════════════════════════

  .plot_one <- function(i) {

    lig <- pairs$ligand[i]
    rec <- pairs$receptor[i]
    fname <- sprintf("%s_%s.png", lig, rec)
    fpath <- file.path(out_dir, fname)

    result <- tryCatch({

      p <- CellToCell_InSitu_Single_LRM(
        ctc            = ctc,
        ligand         = lig,
        receptor       = rec,
        bg_color       = bg_color,
        bg_linewidth   = bg_linewidth,
        bg_alpha       = bg_alpha,
        sig_linewidth  = sig_linewidth,
        sig_alpha      = sig_alpha,
        viridis_option = viridis_option
      )

      ggsave(
        filename = fpath,
        plot     = p,
        width    = width,
        height   = height,
        dpi      = dpi,
        bg       = "black"
      )

      list(ligand = lig, receptor = rec, status = "saved", path = fpath)

    }, error = function(e) {
      message(sprintf("  [SKIP] %s_%s : %s", lig, rec, conditionMessage(e)))
      list(ligand = lig, receptor = rec, status = "skipped", path = NA_character_)
    })

    result
  }

  # ════════════════════════════════════════════════════════════════════════════
  # SECTION 3 ── Execute (parallel or serial)
  # ════════════════════════════════════════════════════════════════════════════

  message(sprintf("[CellToCell_InSitu] Plotting with %d core(s) ...", n_cores))

  if (n_cores > 1L && .Platform$OS.type == "unix") {
    # macOS: disable Objective-C fork safety check that crashes forked workers
    # when GUI frameworks (e.g. ggplot2 font system) are loaded in the parent.
    old_env <- Sys.getenv("OBJC_DISABLE_INITIALIZE_FORK_SAFETY")
    Sys.setenv(OBJC_DISABLE_INITIALIZE_FORK_SAFETY = "YES")
    on.exit(Sys.setenv(OBJC_DISABLE_INITIALIZE_FORK_SAFETY = old_env), add = TRUE)

    job_results <- parallel::mclapply(
      seq_len(n_pairs), .plot_one,
      mc.cores = n_cores
    )
  } else {
    if (n_cores > 1L && .Platform$OS.type != "unix")
      message("  Fork parallelism unavailable on Windows — running serially.")
    job_results <- lapply(seq_len(n_pairs), .plot_one)
  }

  # ════════════════════════════════════════════════════════════════════════════
  # SECTION 4 ── Summarise and report
  # ════════════════════════════════════════════════════════════════════════════

  summary_df <- rbindlist(lapply(job_results, as.data.table))

  n_saved   <- sum(summary_df$status == "saved")
  n_skipped <- sum(summary_df$status == "skipped")

  message(sprintf(
    "[CellToCell_InSitu] Done. Saved: %d | Skipped: %d | Directory: %s",
    n_saved, n_skipped, normalizePath(out_dir, mustWork = FALSE)
  ))

  invisible(as.data.frame(summary_df))
}
