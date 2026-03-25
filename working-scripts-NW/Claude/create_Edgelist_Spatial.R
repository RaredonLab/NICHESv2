# ==============================================================================
# FILE:    create_Edgelist_Spatial.R
# PURPOSE: Build a spatial neighborhood edge list from 2D cell coordinates.
#          Supports k-nearest neighbor (mutual kNN) and fixed-radius modes.
#          Uses kd-tree indexing to avoid O(n²) full distance matrix computation.
#
# DEPS:    RANN (>= 0.9), dbscan (>= 1.1), data.table (>= 1.14), parallel
# AUTHOR:  [Your Name]
# DATE:    2026-03-11
# ==============================================================================


# ── Package loading ─────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  for (pkg in c("RANN", "dbscan", "data.table")) {
    if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
  }
  library(RANN)        # Fast kd-tree kNN:         nn2()
  library(dbscan)      # Fast kd-tree radius NN:   frNN()
  library(data.table)  # Efficient tabular storage: data.table, rbindlist, fintersect
  library(parallel)    # Parallelization:           mclapply / parLapply
})


# ==============================================================================
#' create_Edgelist_Spatial
#'
#' @title Build a Spatial Neighborhood Edge List from 2D Cell Coordinates
#'
#' @description
#' Constructs a directed edge list representing spatial cell neighborhoods.
#' \strong{Autocrine (self-loop) edges A→A with \code{weight = 0} are always
#' appended} for every cell, in both modes.
#' Two computation modes are available for inter-cell edges:
#'
#' \itemize{
#'   \item \strong{kNN mode}    — Finds k nearest neighbors per cell via a kd-tree
#'         (RANN::nn2). Only \emph{mutual} edges are retained: edge A→B is kept
#'         only if B→A also exists. Non-mutual edges are silently discarded.
#'   \item \strong{Radius mode} — Finds all cells within a fixed Euclidean radius
#'         via a kd-tree (dbscan::frNN). Edges are inherently symmetric.
#' }
#'
#' The \code{weight} column (Euclidean distance) is always included.
#' Both modes avoid the O(n²) full distance matrix via kd-tree spatial indexing.
#'
#' @param coord   \code{data.frame}. Two-column spatial table. The first column
#'   is interpreted as x and the second as y. \strong{Row names must be unique
#'   cell identifiers.}
#' @param k       \code{integer} or \code{NULL}. Number of nearest neighbors in
#'   kNN mode. Must be a positive integer ≥ 1.
#' @param rad     \code{numeric} or \code{NULL}. Neighborhood radius in radius
#'   mode. Must be a positive finite number. Units match those of \code{coord}.
#' @param n.cores \code{integer}. Number of CPU cores for parallelization.
#'   Default = 1 (serial). Fork-based on Unix/macOS (\code{mclapply}); falls back
#'   to serial on Windows with a message (see Details).
#'
#' @return A \code{data.table} with columns:
#' \describe{
#'   \item{\code{from}}{Source cell name (character).}
#'   \item{\code{to}}{Target cell name (character).}
#'   \item{\code{weight}}{Euclidean distance (numeric). Self-loops (autocrine
#'     edges) carry \code{weight = 0}.}
#'   \item{\code{x.from}, \code{y.from}}{Spatial coordinates of the sending cell.}
#'   \item{\code{x.to},   \code{y.to}  }{Spatial coordinates of the receiving cell.}
#' }
#' Keyed on \code{(from, to)} for O(log n) downstream lookups. Carries
#' metadata attributes: \code{mode}, \code{param}, \code{n_cells},
#' \code{cell_names}.
#'
#' @details
#' \strong{Parameter resolution (priority order):}
#' \enumerate{
#'   \item Both \code{k} and \code{rad} NULL / invalid →
#'         default \code{k = 4}, user notified via \code{message()}.
#'   \item Both \code{k} and \code{rad} valid →
#'         radius mode is used, \code{k} overwritten, user notified.
#'   \item Only \code{rad} valid → radius mode.
#'   \item Only \code{k} valid  → kNN mode (mutual edges only).
#' }
#'
#' \strong{Mutual kNN filtering:}
#' After building the full directed kNN graph, a set intersection via
#' \code{data.table::fintersect()} retains only edge pairs where both
#' directions exist. Cost is O(E log E) — much cheaper than O(n²).
#'
#' \strong{Parallelization note (Windows):}
#' Windows does not support fork-based parallelism. On Windows, \code{n_cores > 1}
#' emits a message and falls back to serial execution. To use parallelism on
#' Windows, consider wrapping this function with a PSOCK cluster via
#' \code{parallel::makeCluster()} or using the \code{future.apply} package.
#'
#' @seealso \code{\link[RANN]{nn2}}, \code{\link[dbscan]{frNN}},
#'          \code{\link[data.table]{fintersect}}
#'
#' @examples
#' set.seed(42)
#' n  <- 200
#' co <- data.frame(x = runif(n), y = runif(n), row.names = paste0("cell_", seq_len(n)))
#' el.knn <- create_Edgelist_Spatial(co, k = 5)
#' el.rad <- create_Edgelist_Spatial(co, rad = 0.1)
# ==============================================================================

create_Edgelist_Spatial <- function(coord,
                                    k       = NULL,
                                    rad     = NULL,
                                    n.cores = 1L) {

  # ════════════════════════════════════════════════════════════════════════════
  # SECTION 1 ── Validate core inputs
  # ════════════════════════════════════════════════════════════════════════════

  if (!is.data.frame(coord))
    stop("[create_Edgelist_Spatial] 'coord' must be a data.frame.")
  if (ncol(coord) < 2)
    stop("[create_Edgelist_Spatial] 'coord' must have at least 2 columns (x, y).")
  if (is.null(rownames(coord)) || any(rownames(coord) == ""))
    stop("[create_Edgelist_Spatial] 'coord' must have unique, non-empty row names.")
  if (anyDuplicated(rownames(coord)))
    stop("[create_Edgelist_Spatial] Row names of 'coord' must be unique.")
  if (!is.numeric(n.cores) || length(n.cores) != 1 || n.cores < 1)
    stop("[create_Edgelist_Spatial] 'n.cores' must be a single positive integer.")

  n.cores    <- as.integer(n.cores)
  cell.names <- rownames(coord)
  n.cells    <- nrow(coord)
  xy         <- as.matrix(coord[, 1:2, drop = FALSE])  # n x 2 numeric matrix
  colnames(xy) <- NULL                                  # clean for kd-tree calls

  # Coordinate lookup vectors (named by cell) — used to attach x/y to edges
  x.coords <- setNames(xy[, 1], cell.names)
  y.coords <- setNames(xy[, 2], cell.names)

  if (n.cells < 2)
    stop("[create_Edgelist_Spatial] 'coord' must contain at least 2 cells.")

  # ════════════════════════════════════════════════════════════════════════════
  # SECTION 2 ── Parameter helpers & mode resolution
  # ════════════════════════════════════════════════════════════════════════════

  # Returns TRUE only for a single, finite, positive integer value
  .is_valid_k <- function(x)
    !is.null(x) && is.numeric(x) && length(x) == 1L &&
    !is.na(x) && is.finite(x) && x >= 1 && x == floor(x)

  # Returns TRUE only for a single, finite, strictly-positive numeric value
  .is_valid_rad <- function(x)
    !is.null(x) && is.numeric(x) && length(x) == 1L &&
    !is.na(x) && is.finite(x) && x > 0

  k.valid   <- .is_valid_k(k)
  rad.valid <- .is_valid_rad(rad)

  # Resolve mode according to the four-case priority table (see @details)
  if (!k.valid && !rad.valid) {
    message("[create_Edgelist_Spatial] Neither 'k' nor 'rad' is valid or provided. ",
            "Defaulting to k = 4.")
    k          <- 4L
    use.radius <- FALSE

  } else if (k.valid && rad.valid) {
    message("[create_Edgelist_Spatial] Both 'k' and 'rad' are valid. ",
            "'k' is overwritten — switching to radius mode.")
    use.radius <- TRUE

  } else if (rad.valid) {
    use.radius <- TRUE

  } else {
    k          <- as.integer(k)
    use.radius <- FALSE
  }

  # ════════════════════════════════════════════════════════════════════════════
  # SECTION 3 ── Cross-platform parallel helper
  #
  # Fork-based mclapply (Unix/macOS) transparently shares the parent's memory,
  # so closures referencing local variables work without explicit export.
  # Windows does not support fork; we fall back to serial with a message.
  # ════════════════════════════════════════════════════════════════════════════

  .par_lapply <- function(X, FUN) {
    if (n.cores > 1L) {
      if (.Platform$OS.type == "unix") {
        parallel::mclapply(X, FUN, mc.cores = n.cores)
      } else {
        message("[create_Edgelist_Spatial] Fork-based parallelism is unavailable ",
                "on Windows. Falling back to serial execution.\n",
                "  Tip: Use future.apply::future_lapply() with a PSOCK plan ",
                "for cross-platform parallelism.")
        lapply(X, FUN)
      }
    } else {
      lapply(X, FUN)
    }
  }

  # ════════════════════════════════════════════════════════════════════════════
  # SECTION 4A ── Radius mode  (dbscan::frNN, kd-tree)
  #
  # frNN builds a kd-tree once over all n points (O(n log n)) and then
  # queries each point's ε-neighborhood without computing the n×n matrix.
  # Returned distances are Euclidean and already computed — we reuse them
  # instead of recalculating.
  # ════════════════════════════════════════════════════════════════════════════

  if (use.radius) {

    fnn       <- dbscan::frNN(xy, eps = rad, sort = FALSE)
    nb.list   <- fnn$id    # list[[i]]: integer indices of cell i's neighbors
    dist.list <- fnn$dist  # list[[i]]: matching Euclidean distances

    # Build one data.table per cell, then rbindlist (faster than growing a df)
    .build_radius_edges <- function(i) {
      nb    <- nb.list[[i]]
      dists <- dist.list[[i]]

      # frNN may include self (dist = 0 < eps), remove defensively
      keep  <- nb != i
      nb    <- nb[keep]
      dists <- dists[keep]

      if (length(nb) == 0L) return(NULL)

      data.table(from = cell.names[i], to = cell.names[nb], weight = dists)
    }

    edge.lst <- .par_lapply(seq_len(n.cells), .build_radius_edges)
    edge.dt   <- rbindlist(Filter(Negate(is.null), edge.lst))

  # ════════════════════════════════════════════════════════════════════════════
  # SECTION 4B ── kNN mode  (RANN::nn2, kd-tree)
  #
  # nn2 queries all n points at once for their k+1 nearest neighbors
  # (including self as the closest) in O(n log n).
  # We then work entirely with integer indices for speed, translating to
  # cell names only after the mutual-edge filter is applied.
  # ════════════════════════════════════════════════════════════════════════════

  } else {

    # k+1 because nn2 always returns the query point itself as neighbor 1
    k.query  <- min(k + 1L, n.cells)
    nn       <- RANN::nn2(data = xy, k = k.query)

    # Drop self (column 1 is always the query point at distance 0)
    nn.idx   <- nn$nn.idx[,   -1L, drop = FALSE]   # n.cells × k  integer matrix
    nn.dists <- nn$nn.dists[, -1L, drop = FALSE]   # n.cells × k  numeric matrix

    # Build directed edge table using integer indices (cheaper for set ops)
    .build_knn_edges <- function(i) {
      nb    <- nn.idx[i, ]
      valid <- nb > 0L          # nn2 pads with 0 if k.query > n.cells
      nb    <- nb[valid]
      if (length(nb) == 0L) return(NULL)

      data.table(from = i, to = nb, weight = nn.dists[i, valid])
    }

    edge.lst    <- .par_lapply(seq_len(n.cells), .build_knn_edges)
    edge.dt.idx <- rbindlist(Filter(Negate(is.null), edge.lst))

    # ── Mutual edge filter ──────────────────────────────────────────────────
    # Strategy: build the reverse edge table; keep only pairs present in both.
    #
    #   fwd.pairs: all directed edges as (from, to) integer pairs
    #   rev.pairs: same edges with from/to swapped
    #   fintersect: set intersection — O(E log E), no O(n²) comparison
    #
    # Result: only (i→j) pairs where (j→i) also exists are retained.
    # Both directions are preserved (A→B AND B→A), matching undirected
    # graph convention while using a directed data structure.

    fwd.pairs <- edge.dt.idx[, .(from, to)]
    rev.pairs <- edge.dt.idx[, .(from = to, to = from)]  # flip all edges
    mutual.pairs <- fintersect(fwd.pairs, rev.pairs)      # intersection

    # Re-join to recover edge weights from the original directed table
    setkey(edge.dt.idx,  from, to)
    setkey(mutual.pairs, from, to)
    edge.dt.idx <- edge.dt.idx[mutual.pairs, nomatch = 0L]

    # Translate integer indices → cell name strings
    edge.dt <- copy(edge.dt.idx)
    edge.dt[, from := cell.names[from]]
    edge.dt[, to   := cell.names[to]]
  }

  # ════════════════════════════════════════════════════════════════════════════
  # SECTION 5 ── Append autocrine (self-loop) edges, then finalize
  #
  # Every cell gets an A→A edge with weight = 0, representing autocrine
  # signaling / self-influence. These are appended after inter-cell edges
  # so the inter-cell logic above stays unaffected.
  # ════════════════════════════════════════════════════════════════════════════

  autocrine.dt <- data.table(
    from   = cell.names,
    to     = cell.names,
    weight = 0
  )

  edge.dt <- rbindlist(list(edge.dt, autocrine.dt), use.names = TRUE)

  # Warn if no inter-cell edges were produced (helps diagnose tight parameters)
  n.intercell <- nrow(edge.dt) - n.cells
  if (n.intercell == 0L)
    warning("[create_Edgelist_Spatial] No inter-cell edges found. ",
            "Consider increasing 'k' or 'rad'. Only autocrine edges are present.")

  # ── Attach sending and receiving cell coordinates ───────────────────────────
  # Looked up from the named vectors built at the top of the function.
  # These travel with the edge list so downstream functions (computeCellToCellSpatial,
  # plotting) never need to re-join the original coord data.frame.
  edge.dt[, x.from := x.coords[from]]
  edge.dt[, y.from := y.coords[from]]
  edge.dt[, x.to   := x.coords[to]]
  edge.dt[, y.to   := y.coords[to]]

  # Canonical column ordering: from | to | weight | x.from | y.from | x.to | y.to
  setcolorder(edge.dt, c("from", "to", "weight", "x.from", "y.from", "x.to", "y.to"))

  # Key on (from, to) → O(log n) lookup of all edges from a given cell
  setkeyv(edge.dt, c("from", "to"))

  # Attach run metadata as attributes for downstream functions
  attr(edge.dt, "mode")       <- if (use.radius) "radius" else "kNN"
  attr(edge.dt, "param")      <- if (use.radius) rad else k
  attr(edge.dt, "n.cells")    <- n.cells
  attr(edge.dt, "cell.names") <- cell.names

  message(sprintf(
    "[create_Edgelist_Spatial] Done. Mode: %-6s | Param: %s | Cells: %d | Inter-cell edges: %d | Autocrine edges: %d",
    attr(edge.dt, "mode"),
    attr(edge.dt, "param"),
    n.cells,
    n.intercell,
    n.cells
  ))

  return(invisible(edge.dt))
}
