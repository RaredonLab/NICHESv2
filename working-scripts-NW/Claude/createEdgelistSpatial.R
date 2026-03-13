# ==============================================================================
# FILE:    createEdgelistSpatial.R
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
#' createEdgelistSpatial
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
#' @param n_cores \code{integer}. Number of CPU cores for parallelization.
#'   Default = 1 (serial). Fork-based on Unix/macOS (\code{mclapply}); falls back
#'   to serial on Windows with a message (see Details).
#'
#' @return A \code{data.table} with columns:
#' \describe{
#'   \item{\code{from}}{Source cell name (character).}
#'   \item{\code{to}}{Target cell name (character).}
#'   \item{\code{weight}}{Euclidean distance (numeric). Self-loops (autocrine
#'     edges) carry \code{weight = 0}.}
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
#' el_knn <- createEdgelistSpatial(co, k = 5)
#' el_rad <- createEdgelistSpatial(co, rad = 0.1)
# ==============================================================================

createEdgelistSpatial <- function(coord,
                                   k       = NULL,
                                   rad     = NULL,
                                   n_cores = 1L) {

  # ════════════════════════════════════════════════════════════════════════════
  # SECTION 1 ── Validate core inputs
  # ════════════════════════════════════════════════════════════════════════════

  if (!is.data.frame(coord))
    stop("[createEdgelistSpatial] 'coord' must be a data.frame.")
  if (ncol(coord) < 2)
    stop("[createEdgelistSpatial] 'coord' must have at least 2 columns (x, y).")
  if (is.null(rownames(coord)) || any(rownames(coord) == ""))
    stop("[createEdgelistSpatial] 'coord' must have unique, non-empty row names.")
  if (anyDuplicated(rownames(coord)))
    stop("[createEdgelistSpatial] Row names of 'coord' must be unique.")
  if (!is.numeric(n_cores) || length(n_cores) != 1 || n_cores < 1)
    stop("[createEdgelistSpatial] 'n_cores' must be a single positive integer.")

  n_cores    <- as.integer(n_cores)
  cell_names <- rownames(coord)
  n_cells    <- nrow(coord)
  xy         <- as.matrix(coord[, 1:2, drop = FALSE])  # n x 2 numeric matrix
  colnames(xy) <- NULL                                  # clean for kd-tree calls

  if (n_cells < 2)
    stop("[createEdgelistSpatial] 'coord' must contain at least 2 cells.")

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

  k_valid   <- .is_valid_k(k)
  rad_valid <- .is_valid_rad(rad)

  # Resolve mode according to the four-case priority table (see @details)
  if (!k_valid && !rad_valid) {
    message("[createEdgelistSpatial] Neither 'k' nor 'rad' is valid or provided. ",
            "Defaulting to k = 4.")
    k          <- 4L
    use_radius <- FALSE

  } else if (k_valid && rad_valid) {
    message("[createEdgelistSpatial] Both 'k' and 'rad' are valid. ",
            "'k' is overwritten — switching to radius mode.")
    use_radius <- TRUE

  } else if (rad_valid) {
    use_radius <- TRUE

  } else {
    k          <- as.integer(k)
    use_radius <- FALSE
  }

  # ════════════════════════════════════════════════════════════════════════════
  # SECTION 3 ── Cross-platform parallel helper
  #
  # Fork-based mclapply (Unix/macOS) transparently shares the parent's memory,
  # so closures referencing local variables work without explicit export.
  # Windows does not support fork; we fall back to serial with a message.
  # ════════════════════════════════════════════════════════════════════════════

  .par_lapply <- function(X, FUN) {
    if (n_cores > 1L) {
      if (.Platform$OS.type == "unix") {
        parallel::mclapply(X, FUN, mc.cores = n_cores)
      } else {
        message("[createEdgelistSpatial] Fork-based parallelism is unavailable ",
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

  if (use_radius) {

    fnn       <- dbscan::frNN(xy, eps = rad, sort = FALSE)
    nb_list   <- fnn$id    # list[[i]]: integer indices of cell i's neighbors
    dist_list <- fnn$dist  # list[[i]]: matching Euclidean distances

    # Build one data.table per cell, then rbindlist (faster than growing a df)
    .build_radius_edges <- function(i) {
      nb    <- nb_list[[i]]
      dists <- dist_list[[i]]

      # frNN may include self (dist = 0 < eps), remove defensively
      keep  <- nb != i
      nb    <- nb[keep]
      dists <- dists[keep]

      if (length(nb) == 0L) return(NULL)

      data.table(from = cell_names[i], to = cell_names[nb], weight = dists)
    }

    edge_list <- .par_lapply(seq_len(n_cells), .build_radius_edges)
    edge_dt   <- rbindlist(Filter(Negate(is.null), edge_list))

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
    k_query  <- min(k + 1L, n_cells)
    nn       <- RANN::nn2(data = xy, k = k_query)

    # Drop self (column 1 is always the query point at distance 0)
    nn_idx   <- nn$nn.idx[,   -1L, drop = FALSE]   # n_cells × k  integer matrix
    nn_dists <- nn$nn.dists[, -1L, drop = FALSE]   # n_cells × k  numeric matrix

    # Build directed edge table using integer indices (cheaper for set ops)
    .build_knn_edges <- function(i) {
      nb    <- nn_idx[i, ]
      valid <- nb > 0L          # nn2 pads with 0 if k_query > n_cells
      nb    <- nb[valid]
      if (length(nb) == 0L) return(NULL)

      data.table(from = i, to = nb, weight = nn_dists[i, valid])
    }

    edge_list   <- .par_lapply(seq_len(n_cells), .build_knn_edges)
    edge_dt_idx <- rbindlist(Filter(Negate(is.null), edge_list))

    # ── Mutual edge filter ──────────────────────────────────────────────────
    # Strategy: build the reverse edge table; keep only pairs present in both.
    #
    #   fwd_pairs: all directed edges as (from, to) integer pairs
    #   rev_pairs: same edges with from/to swapped
    #   fintersect: set intersection — O(E log E), no O(n²) comparison
    #
    # Result: only (i→j) pairs where (j→i) also exists are retained.
    # Both directions are preserved (A→B AND B→A), matching undirected
    # graph convention while using a directed data structure.

    fwd_pairs <- edge_dt_idx[, .(from, to)]
    rev_pairs <- edge_dt_idx[, .(from = to, to = from)]  # flip all edges
    mutual_pairs <- fintersect(fwd_pairs, rev_pairs)      # intersection

    # Re-join to recover edge weights from the original directed table
    setkey(edge_dt_idx,  from, to)
    setkey(mutual_pairs, from, to)
    edge_dt_idx <- edge_dt_idx[mutual_pairs, nomatch = 0L]

    # Translate integer indices → cell name strings
    edge_dt <- copy(edge_dt_idx)
    edge_dt[, from := cell_names[from]]
    edge_dt[, to   := cell_names[to]]
  }

  # ════════════════════════════════════════════════════════════════════════════
  # SECTION 5 ── Append autocrine (self-loop) edges, then finalize
  #
  # Every cell gets an A→A edge with weight = 0, representing autocrine
  # signaling / self-influence. These are appended after inter-cell edges
  # so the inter-cell logic above stays unaffected.
  # ════════════════════════════════════════════════════════════════════════════

  autocrine_dt <- data.table(
    from   = cell_names,
    to     = cell_names,
    weight = 0
  )

  edge_dt <- rbindlist(list(edge_dt, autocrine_dt), use.names = TRUE)

  # Warn if no inter-cell edges were produced (helps diagnose tight parameters)
  n_intercell <- nrow(edge_dt) - n_cells
  if (n_intercell == 0L)
    warning("[createEdgelistSpatial] No inter-cell edges found. ",
            "Consider increasing 'k' or 'rad'. Only autocrine edges are present.")

  # Canonical column ordering: from | to | weight
  setcolorder(edge_dt, c("from", "to", "weight"))

  # Key on (from, to) → O(log n) lookup of all edges from a given cell
  setkeyv(edge_dt, c("from", "to"))

  # Attach run metadata as attributes for downstream functions
  attr(edge_dt, "mode")       <- if (use_radius) "radius" else "kNN"
  attr(edge_dt, "param")      <- if (use_radius) rad else k
  attr(edge_dt, "n_cells")    <- n_cells
  attr(edge_dt, "cell_names") <- cell_names

  message(sprintf(
    "[createEdgelistSpatial] Done. Mode: %-6s | Param: %s | Cells: %d | Inter-cell edges: %d | Autocrine edges: %d",
    attr(edge_dt, "mode"),
    attr(edge_dt, "param"),
    n_cells,
    n_intercell,
    n_cells
  ))

  return(invisible(edge_dt))
}
