# ==============================================================================
# FILE:    CellToCell_InSitu.R
# PURPOSE: Batch in situ plotting of LRM signals across all or a selected
#          subset of pairs in a NICHESObject. Calls
#          CellToCell_InSitu_Single_LRM() for each pair and saves PNGs.
#
# FUNCTIONS:
#   CellToCell_InSitu()  — main batch plotting function
#
# DEPS: ggplot2, data.table, parallel
# ==============================================================================


# ==============================================================================
#' Batch In Situ Visualisation of Ligand-Receptor Signals
#'
#' @description
#' Iterates over LRM pairs in a \code{NICHESObject} and calls
#' \code{CellToCell_InSitu_Single_LRM()} for each, saving one PNG per pair.
#'
#' Two output modes:
#' \itemize{
#'   \item \strong{\code{LRM.list = NULL}} — plots all pairs in
#'     \code{$LRM.list}; saved to \code{folder}.
#'   \item \strong{\code{LRM.list = character vector}} — plots only the named
#'     pairs; saved to \code{folder/selected_LRMs/}.
#' }
#'
#' @param niches.obj     A \code{NICHESObject}.
#' @param LRM.list       \code{NULL} (all pairs) or a character vector of
#'   \code{"ligand|receptor"} strings to plot selectively.
#' @param folder         Output directory. Created if absent. Default
#'   \code{"niches_plots"}.
#' @param bg.color       Background edge colour. Default \code{"grey50"}.
#' @param bg.linewidth   Background edge line width. Default \code{0.1}.
#' @param bg.alpha       Background edge alpha. Default \code{0.2}.
#' @param sig.linewidth  Signal edge line width. Default \code{0.4}.
#' @param sig.alpha      Signal edge alpha. Default \code{0.8}.
#' @param viridis.option Viridis palette letter. Default \code{"C"}.
#' @param width          Plot width in inches. Default \code{7}.
#' @param height         Plot height in inches. Default \code{6}.
#' @param dpi            Plot resolution. Default \code{150}.
#' @param n.cores        Parallel cores. Default \code{1}. macOS forces serial
#'   (fork + ggplot2 font-rendering incompatibility).
#'
#' @return Invisibly, a \code{data.frame} with columns \code{LRM},
#'   \code{status} (\code{"saved"} or \code{"skipped"}), and \code{path}.
#'
#' @seealso \code{\link{CellToCell_InSitu_Single_LRM}}
#'
#' @export
CellToCell_InSitu <- function(niches.obj,
                               LRM.list       = NULL,
                               folder         = "niches_plots",
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

  if (!is_NICHESObject(niches.obj))
    stop("[CellToCell_InSitu] 'niches.obj' must be a NICHESObject.")

  n.cores <- as.integer(n.cores)

  # ════════════════════════════════════════════════════════════════════════════
  # SECTION 1 ── Resolve output directory and LRM pair list
  # ════════════════════════════════════════════════════════════════════════════

  selected.mode <- !is.null(LRM.list)

  out.dir <- if (selected.mode) {
    file.path(folder, "selected_LRMs")
  } else {
    folder
  }

  dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)
  message(sprintf("[CellToCell_InSitu] Output directory: %s", normalizePath(out.dir, mustWork = FALSE)))

  if (selected.mode) {
    if (!is.character(LRM.list) || length(LRM.list) == 0L)
      stop("[CellToCell_InSitu] 'LRM.list' must be a non-empty character vector ",
           "of 'ligand|receptor' strings.")

    # Warn about any pairs not found in the object
    missing.lrms <- setdiff(LRM.list, niches.obj$LRM.list)
    if (length(missing.lrms) > 0L)
      warning(sprintf(
        "[CellToCell_InSitu] %d LRM(s) not found in $LRM.list and will be skipped: %s",
        length(missing.lrms), paste(head(missing.lrms, 3L), collapse = ", ")
      ))

    pairs <- data.table(LRM = LRM.list)
  } else {
    pairs <- data.table(LRM = niches.obj$LRM.list)
  }

  n.pairs <- nrow(pairs)
  message(sprintf("[CellToCell_InSitu] %d LRM pair(s) to plot.", n.pairs))

  # ════════════════════════════════════════════════════════════════════════════
  # SECTION 2 ── Per-pair plotting job
  # ════════════════════════════════════════════════════════════════════════════

  .plot_one <- function(i) {
    lrm   <- pairs$LRM[i]
    fpath <- file.path(out.dir, sprintf("%s.png", gsub("|", "_", lrm, fixed = TRUE)))

    tryCatch({
      p <- CellToCell_InSitu_Single_LRM(
        niches.obj     = niches.obj,
        LRM            = lrm,
        bg.color       = bg.color,
        bg.linewidth   = bg.linewidth,
        bg.alpha       = bg.alpha,
        sig.linewidth  = sig.linewidth,
        sig.alpha      = sig.alpha,
        viridis.option = viridis.option
      )
      ggsave(filename = fpath, plot = p,
             width = width, height = height, dpi = dpi, bg = "black")
      list(LRM = lrm, status = "saved", path = fpath)
    }, error = function(e) {
      message(sprintf("  [SKIP] %s : %s", lrm, conditionMessage(e)))
      list(LRM = lrm, status = "skipped", path = NA_character_)
    })
  }

  # ════════════════════════════════════════════════════════════════════════════
  # SECTION 3 ── Execute (parallel on Linux, serial on macOS/Windows)
  #
  # macOS: fork() + ggplot2 font rendering crashes at runtime.
  # Windows: no fork-based parallelism.
  # Both fall back to serial automatically.
  # ════════════════════════════════════════════════════════════════════════════

  is.mac   <- Sys.info()[["sysname"]] == "Darwin"
  use.fork <- n.cores > 1L && .Platform$OS.type == "unix" && !is.mac

  if (n.cores > 1L && is.mac) {
    message("  macOS: falling back to serial execution (fork + ggplot2 incompatibility).")
  }

  message(sprintf("[CellToCell_InSitu] Plotting with %d core(s) ...", if (use.fork) n.cores else 1L))

  job.results <- if (use.fork) {
    parallel::mclapply(seq_len(n.pairs), .plot_one, mc.cores = n.cores)
  } else {
    lapply(seq_len(n.pairs), .plot_one)
  }

  # ════════════════════════════════════════════════════════════════════════════
  # SECTION 4 ── Summarise and report
  # ════════════════════════════════════════════════════════════════════════════

  summary.df <- rbindlist(lapply(job.results, as.data.table))
  n.saved    <- sum(summary.df$status == "saved")
  n.skipped  <- sum(summary.df$status == "skipped")

  message(sprintf(
    "[CellToCell_InSitu] Done. Saved: %d | Skipped: %d | Directory: %s",
    n.saved, n.skipped, normalizePath(out.dir, mustWork = FALSE)
  ))

  invisible(as.data.frame(summary.df))
}
