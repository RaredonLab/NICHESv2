# CellToCell_InSitu.R -----------------------------------------------------
# Batch in situ plotting of LRM signals across all or a selected
# subset of pairs in a NICHESObject. Calls
# CellToCell_InSitu_Single_LRM() for each pair and saves PNGs.
#
# FUNCTIONS:
#   CellToCell_InSitu()  -- main batch plotting function
#
# DEPS: ggplot2, data.table, scales, parallel


# ==============================================================================
#' Batch In Situ Visualisation of Ligand-Receptor Signals
#'
#' @description
#' Iterates over LRM pairs in a \code{NICHESObject} and calls
#' \code{\link{CellToCell_InSitu_Single_LRM}} for each, saving one PNG per
#' pair to disk. Designed for batch reporting -- generating dozens or
#' hundreds of in situ visualisations in one call without holding all
#' figure objects in memory.
#'
#' Two output modes:
#' \itemize{
#'   \item \strong{\code{LRM.list = NULL}} -- plots all pairs in
#'     \code{$LRM.list}; PNGs saved directly into \code{folder}.
#'   \item \strong{\code{LRM.list = character vector}} -- plots only the
#'     named pairs; PNGs saved into \code{folder/selected_LRMs/}.
#' }
#'
#' Each plot is rendered against a shared background of all inter-cell
#' edges in the object's \code{$edge.list}, with the LRM-specific signal
#' overlaid in viridis colours.
#'
#' @section Files written:
#' For each LRM pair, one PNG file is created. The filename is the LRM
#' string with the \code{"|"} separator replaced by \code{"_"} -- so
#' LRM \code{"VEGFA|KDR"} produces file \code{VEGFA_KDR.png}. Pairs that
#' fail to plot (e.g.\ all-zero scores after filtering) are skipped with
#' a message; their absence from the output directory is intentional.
#' The \code{folder} argument is created (recursively) if it does not
#' already exist.
#'
#' @section Performance:
#' For large \code{NICHESObject}s (tens of thousands of cells, hundreds of
#' LRM pairs) two optimisations make this function tractable:
#' \enumerate{
#'   \item \strong{Slim object before serialisation.} \code{$cell.data} is
#'         set to \code{NULL} and \code{$edge.data} is pre-filtered to the
#'         pairs being plotted before being shipped to workers, eliminating
#'         the bulk of the per-worker payload.
#'   \item \strong{Pre-computed background edge table.} The inter-cell edge
#'         coordinates -- identical for every plot -- are joined once in
#'         the parent process and shipped to every worker via
#'         \code{precomputed.bg}, avoiding two large coordinate joins per
#'         plot.
#' }
#' Both optimisations are transparent to the caller; they are documented
#' here only so power users understand the intended scaling behaviour.
#'
#' @section Parallelization:
#' \code{n.cores > 1} uses a socket cluster
#' (\code{parallel::makeCluster}). Socket clusters are required here, not
#' merely preferred: fork-based parallelism (\code{parallel::mclapply})
#' inherits the parent's CoreText state on macOS, which crashes when
#' ggplot2 renders fonts inside a forked child. Socket clusters work
#' identically on Windows, macOS, and Linux. Each worker has its
#' \code{data.table} thread count throttled to 1 to prevent OpenMP
#' contention with the outer parallelism. The cluster is torn down via
#' \code{on.exit(stopCluster(cl))} so it is released even if a render
#' errors midway.
#'
#' @param niches.obj     A \code{NICHESObject} produced by
#'   \code{\link{create_NICHESObject}}. Must have \code{$cell.meta} with
#'   \code{x}/\code{y} coordinates and a populated \code{$edge.data}.
#' @param LRM.list       \code{NULL} (plot all pairs in
#'   \code{niches.obj$LRM.list}, default) or a character vector of
#'   \code{"ligand|receptor"} strings to plot selectively. Names not
#'   present in \code{$LRM.list} produce a warning and are skipped.
#' @param folder         Output directory. Created (recursively) if it
#'   does not exist. Default \code{"niches_plots"} -- in non-interactive
#'   contexts, supply an absolute path or one rooted under
#'   \code{tempdir()}.
#' @param value.col      Character. Score column to map to the colour and
#'   alpha aesthetics. \code{"score"} (default) shows raw signal magnitude
#'   and is the right choice for absolute comparisons across LRMs.
#'   \code{"score.norm"} shows per-edge proportional signal (sums to 1
#'   across LRMs within an edge) and is the right choice for highlighting
#'   which LRM dominates the communication budget at each edge.
#' @param scale.upper    Numeric or \code{NULL}. Upper limit of the colour
#'   scale. \code{NULL} (default) auto-computes the 90th percentile of
#'   non-zero values \emph{per pair}, which produces visually balanced
#'   individual plots but means colour values are \strong{not comparable}
#'   across PNGs. To compare across LRMs, set this manually to a fixed
#'   numeric value.
#' @param bg.color       Background edge colour. Default \code{"grey50"}.
#' @param bg.linewidth   Background edge line width. Default \code{0.1}.
#' @param bg.alpha       Background edge alpha. Default \code{0.2}.
#' @param sig.linewidth  Signal edge line width. Default \code{0.4}.
#' @param sig.alpha      Maximum signal edge alpha. Default \code{0.8}.
#' @param viridis.option Viridis palette letter
#'   (\code{"A"}--\code{"H"}). Default \code{"C"} (plasma).
#' @param width          Plot width in inches. Default \code{7}.
#' @param height         Plot height in inches. Default \code{6}.
#' @param dpi            Plot resolution. Default \code{150}.
#' @param n.cores        Integer. Parallel cores for batch rendering.
#'   Default \code{1L}. See the Parallelization section.
#'
#' @return Invisibly, a \code{data.frame} with one row per LRM pair and
#'   columns:
#' \describe{
#'   \item{\code{LRM}}{The \code{"ligand|receptor"} string.}
#'   \item{\code{status}}{\code{"saved"} when the PNG was written
#'     successfully, \code{"skipped"} when the per-pair plotter errored
#'     (the underlying error message is sent to the console as it occurs).}
#'   \item{\code{path}}{Absolute path of the PNG, or \code{NA} when
#'     skipped.}
#' }
#' Useful for programmatic auditing of large batch runs.
#'
#' @seealso
#'   \code{\link{CellToCell_InSitu_Single_LRM}} for the per-pair plotter
#'     called internally;
#'   \code{\link{NICHESv2_example}} for an example object suitable for
#'     batch plotting;
#'   \code{\link{aggregate_NICHESObject}} for non-graphical per-cell
#'     summaries of the same edge data.
#'
#' @examples
#' \donttest{
#' data(NICHESv2_example)
#'
#' # Use tempdir() for the output -- required for CRAN-safe examples.
#' tmp.dir <- file.path(tempdir(), "niches_plots_demo")
#'
#' # Plot just 2 of the 15 LRM pairs to keep the example fast.
#' # Selective mode writes into folder/selected_LRMs/.
#' pairs.to.plot <- head(NICHESv2_example$LRM.list, 2L)
#'
#' summary.df <- CellToCell_InSitu(
#'   niches.obj = NICHESv2_example,
#'   LRM.list   = pairs.to.plot,
#'   folder     = tmp.dir,
#'   n.cores    = 1L
#' )
#'
#' summary.df              # one row per pair: LRM, status, path
#' list.files(file.path(tmp.dir, "selected_LRMs"))   # generated PNGs
#'
#' # Clean up
#' unlink(tmp.dir, recursive = TRUE)
#' }
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

  ## SECTION 1: Resolve output directory and LRM pair list ----

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

  ## SECTION 2: Per-pair plotting job ----

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

  ## SECTION 3: Slim the NICHESObject before export ----
  #
  # Socket clusters serialise every exported object and send it over a local
  # socket to each worker. The full niches.obj is expensive to ship because:
  #   $cell.data  -- large sparse table, completely unused by plotting
  #   $edge.data  -- contains ALL LRMs; only the pairs being plotted are needed
  #
  # We overwrite niches.obj in this local environment with a stripped version:
  #   $cell.data  -> NULL
  #   $edge.data  -> filtered to rows whose LRM is in the current pair list
  #
  # .plot_one() looks up niches.obj at call time (R lexical scoping), so it
  # automatically uses the slim version with no other code changes required.

  niches.obj$cell.data <- NULL
  niches.obj$edge.data <- niches.obj$edge.data[LRM %in% pairs$LRM]

  ## SECTION 4: Pre-compute background edge table once ----
  #
  # The background (all inter-cell edges with coordinates attached) is identical
  # for every LRM -- building it inside Single_LRM would repeat two O(n_edges)
  # joins for every plot. We build it once here and pass it to every worker
  # via precomputed.bg, so Single_LRM skips Section 3 entirely in batch mode.

  coords    <- niches.obj$cell.meta[, .(cell, x, y)]
  all.edges <- niches.obj$edge.list[SendingCell != ReceivingCell]

  all.edges <- coords[all.edges, on = .(cell = SendingCell)]
  setnames(all.edges, c("x", "y"), c("x.from", "y.from"))
  setnames(all.edges, "cell", "SendingCell")

  all.edges <- coords[all.edges, on = .(cell = ReceivingCell)]
  setnames(all.edges, c("x", "y"), c("x.to", "y.to"))
  setnames(all.edges, "cell", "ReceivingCell")

  ## SECTION 5: Execute (socket cluster on all platforms) ----
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

  ## SECTION 6: Summarise and report ----

  summary.df <- rbindlist(lapply(job.results, as.data.table))
  n.saved    <- sum(summary.df$status == "saved")
  n.skipped  <- sum(summary.df$status == "skipped")

  message(sprintf(
    "[CellToCell_InSitu] Done. Saved: %d | Skipped: %d | Directory: %s",
    n.saved, n.skipped, normalizePath(out.dir, mustWork = FALSE)
  ))

  invisible(as.data.frame(summary.df))
}
