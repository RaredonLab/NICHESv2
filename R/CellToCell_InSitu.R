# ==============================================================================
# FILE:    CellToCell_InSitu.R
# PURPOSE: Batch in situ plotting of LRM signals across all or a selected
#          subset of pairs in a NICHESObject. Calls
#          CellToCell_InSitu_Single_LRM() for each pair and saves PNGs.
#
# FUNCTIONS:
#   CellToCell_InSitu()  — main batch plotting function
#
# DEPS: ggplot2, data.table, scales, parallel
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
#' @param value.col      Character. Score column to map to colour and alpha.
#'   \code{"score"} (raw, default) or \code{"score.norm"}. Passed to
#'   \code{CellToCell_InSitu_Single_LRM()}.
#' @param scale.upper    Numeric or \code{NULL}. Upper colour-scale limit.
#'   \code{NULL} (default) auto-computes the 90th percentile of non-zero values
#'   per pair. Passed to \code{CellToCell_InSitu_Single_LRM()}.
#' @param bg.color       Background edge colour. Default \code{"grey50"}.
#' @param bg.linewidth   Background edge line width. Default \code{0.1}.
#' @param bg.alpha       Background edge alpha. Default \code{0.2}.
#' @param sig.linewidth  Signal edge line width. Default \code{0.4}.
#' @param sig.alpha      Maximum signal edge alpha. Default \code{0.8}.
#' @param viridis.option Viridis palette letter. Default \code{"C"}.
#' @param width          Plot width in inches. Default \code{7}.
#' @param height         Plot height in inches. Default \code{6}.
#' @param dpi            Plot resolution. Default \code{150}.
#' @param n.cores        Parallel cores. Default \code{1}. Uses socket-based
#'   clusters (\code{parallel::makeCluster}) which work on all platforms
#'   (Windows, macOS, Linux).
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
                               value.col      = "score",
                               scale.upper    = NULL,
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
        value.col      = value.col,
        scale.upper    = scale.upper,
        bg.color       = bg.color,
        bg.linewidth   = bg.linewidth,
        bg.alpha       = bg.alpha,
        sig.linewidth  = sig.linewidth,
        sig.alpha      = sig.alpha,
        viridis.option = viridis.option,
        precomputed.bg = all.edges
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
  # SECTION 3 ── Slim the NICHESObject before export
  #
  # Socket clusters serialise every exported object and send it over a local
  # socket to each worker. The full niches.obj is expensive to ship because:
  #   $cell.data  — large sparse table, completely unused by plotting
  #   $edge.data  — contains ALL LRMs; only the pairs being plotted are needed
  #
  # We overwrite niches.obj in this local environment with a stripped version:
  #   $cell.data  → NULL
  #   $edge.data  → filtered to rows whose LRM is in the current pair list
  #
  # .plot_one() looks up niches.obj at call time (R lexical scoping), so it
  # automatically uses the slim version with no other code changes required.
  # ════════════════════════════════════════════════════════════════════════════

  niches.obj$cell.data <- NULL
  niches.obj$edge.data <- niches.obj$edge.data[LRM %in% pairs$LRM]

  # ════════════════════════════════════════════════════════════════════════════
  # SECTION 4 ── Pre-compute background edge table once
  #
  # The background (all inter-cell edges with coordinates attached) is identical
  # for every LRM — building it inside Single_LRM would repeat two O(n_edges)
  # joins for every plot. We build it once here and pass it to every worker
  # via precomputed.bg, so Single_LRM skips Section 3 entirely in batch mode.
  # ════════════════════════════════════════════════════════════════════════════

  coords    <- niches.obj$cell.meta[, .(cell, x, y)]
  all.edges <- niches.obj$edge.list[SendingCell != ReceivingCell]

  all.edges <- coords[all.edges, on = .(cell = SendingCell)]
  setnames(all.edges, c("x", "y"), c("x.from", "y.from"))
  setnames(all.edges, "cell", "SendingCell")

  all.edges <- coords[all.edges, on = .(cell = ReceivingCell)]
  setnames(all.edges, c("x", "y"), c("x.to", "y.to"))
  setnames(all.edges, "cell", "ReceivingCell")

  # ════════════════════════════════════════════════════════════════════════════
  # SECTION 5 ── Execute (socket cluster on all platforms)
  #
  # Socket clusters (parallel::makeCluster + parLapply) work on Windows,
  # macOS, and Linux. Unlike fork-based mclapply, they avoid the CoreText
  # crash on macOS that occurs when ggplot2 font rendering state is inherited
  # by forked child processes.
  #
  # Workers receive:
  #   - CellToCell_InSitu_Single_LRM from the NICHESv2 namespace
  #   - slim niches.obj, pre-computed all.edges, and all plot parameters
  #   - ggplot2, data.table, scales loaded; data.table threads set to 1
  #     to prevent OpenMP contention across socket-cluster workers
  #
  # on.exit() guarantees stopCluster() even if an error occurs mid-batch.
  # ════════════════════════════════════════════════════════════════════════════

  message(sprintf("[CellToCell_InSitu] Plotting with %d core(s) ...", n.cores))

  if (n.cores > 1L) {
    cl <- parallel::makeCluster(n.cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterEvalQ(cl, {
      library(ggplot2)
      library(data.table)
      library(scales)
      data.table::setDTthreads(1)
    })
    # Ship the plotting function from the package namespace.
    # asNamespace() works with both devtools::load_all() and installed builds.
    parallel::clusterExport(cl,
      varlist = c("CellToCell_InSitu_Single_LRM"),
      envir   = asNamespace("NICHESv2"))
    # Ship the slim niches.obj, the pre-computed background, and all plot
    # parameters. $cell.data is NULL and $edge.data is pre-filtered, so
    # serialisation cost is minimal. all.edges is shipped once and reused
    # across every plot job, eliminating repeated coordinate joins in workers.
    parallel::clusterExport(cl,
      varlist = c("pairs", "out.dir", "niches.obj", "all.edges",
                  "value.col", "scale.upper",
                  "bg.color", "bg.linewidth", "bg.alpha", "sig.linewidth",
                  "sig.alpha", "viridis.option", "width", "height", "dpi"),
      envir = environment())
    job.results <- parallel::parLapply(cl, seq_len(n.pairs), .plot_one)
  } else {
    job.results <- lapply(seq_len(n.pairs), .plot_one)
  }

  # ════════════════════════════════════════════════════════════════════════════
  # SECTION 5 ── Summarise and report
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
