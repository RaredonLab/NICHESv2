# create_Edgelist_Spatial.R -----------------------------------------------
# Build a spatial neighbourhood edge list from 2-D cell coordinates.
# Supports mutual k-nearest-neighbour (kNN) and fixed-radius modes.
# Uses kd-tree indexing to avoid O(n^2) full distance matrix computation.
# Autocrine (self-loop) edges are always appended for every cell.
#
# FUNCTIONS:
#   create_Edgelist_Spatial()  -- main function
#
# DEPS: RANN (>= 2.6), dbscan (>= 1.1), data.table (>= 1.14), parallel


# ==============================================================================
#' Build a Spatial Neighborhood Edge List from 2-D Cell Coordinates
#'
#' @description
#' Constructs a directed edge list representing spatial cell neighborhoods.
#' Two computation modes are available for inter-cell edges:
#'
#' \itemize{
#'   \item \strong{kNN mode} -- finds \code{k} nearest neighbors per cell via a
#'     kd-tree (\code{RANN::nn2}). Only \emph{mutual} edges are retained (A->B
#'     kept only if B->A also exists). Non-mutual edges are discarded.
#'   \item \strong{Radius mode} -- finds all cells within Euclidean distance
#'     \code{rad} via a kd-tree (\code{dbscan::frNN}). Edges are inherently
#'     symmetric.
#' }
#'
#' Autocrine self-loop edges (A->A, \code{distance = 0}) are always appended for
#' every cell in both modes.
#'
#' @param coord   \code{data.frame}. First column = x, second column = y.
#'   \strong{Row names must be unique cell identifiers.}
#' @param k       \code{integer} or \code{NULL}. Number of nearest neighbors
#'   (kNN mode). Positive integer.
#' @param rad     \code{numeric} or \code{NULL}. Neighborhood radius (radius
#'   mode). Positive finite number in the same units as \code{coord}.
#' @param n.cores \code{integer}. CPU cores. Default 1. Uses socket-based
#'   clusters (\code{parallel::makeCluster}) which work on all platforms
#'   (Windows, macOS, Linux).
#'
#' @return A \code{data.table} keyed on \code{(from, to)} with columns:
#' \describe{
#'   \item{\code{from}}{Sending cell barcode (character).}
#'   \item{\code{to}}{Receiving cell barcode (character).}
#'   \item{\code{distance}}{Euclidean distance; 0 for self-loops (numeric).}
#'   \item{\code{x.from}, \code{y.from}}{Sending cell coordinates (numeric).}
#'   \item{\code{x.to},   \code{y.to}  }{Receiving cell coordinates (numeric).}
#' }
#' Carries metadata attributes: \code{mode}, \code{param}, \code{n.cells},
#' \code{cell.names}.
#'
#' @details
#' This function is a low-level primitive: it builds the spatial edge
#' list that \code{\link{create_NICHESObject}} consumes internally. Most
#' users should call \code{create_NICHESObject(mode = "spatial", ...)},
#' which invokes this function and assembles the result into a full
#' \code{NICHESObject} in one step. \code{create_Edgelist_Spatial()} is
#' exposed for power users who want to inspect or post-process the edge
#' list, benchmark alternative \code{k} or \code{rad} values without
#' re-running the scoring step, or feed the edges into a non-NICHESv2
#' scoring pipeline.
#'
#' \strong{Parameter resolution (priority):}
#' \enumerate{
#'   \item Both NULL / invalid -> default \code{k = 4}.
#'   \item Both valid -> radius mode; user notified.
#'   \item Only \code{rad} valid -> radius mode.
#'   \item Only \code{k} valid -> kNN mode (mutual edges only).
#' }
#'
#' @examples
#' set.seed(1)
#' co <- data.frame(x = runif(50), y = runif(50),
#'                  row.names = paste0("cell_", 1:50))
#' el <- create_Edgelist_Spatial(co, k = 5)
#' el.rad <- create_Edgelist_Spatial(co, rad = 0.2)
#'
#' @keywords internal
#' @export
create_Edgelist_Spatial <- function(coord,
                                     k       = NULL,
                                     rad     = NULL,
                                     n.cores = 1L) {

  ## SECTION 1: Input validation ----

  if (!is.data.frame(coord))
    stop("[create_Edgelist_Spatial] 'coord' must be a data.frame.")
  if (ncol(coord) < 2L)
    stop("[create_Edgelist_Spatial] 'coord' must have at least 2 columns (x, y).")
  if (is.null(rownames(coord)) || any(rownames(coord) == ""))
    stop("[create_Edgelist_Spatial] 'coord' must have unique, non-empty row names.")
  if (anyDuplicated(rownames(coord))) {
    dupes <- unique(rownames(coord)[duplicated(rownames(coord))])
    stop(paste0(
      "[create_Edgelist_Spatial] Row names of 'coord' must be unique cell barcodes.\n",
      "  Duplicate(s) found: ", paste(head(dupes, 5), collapse = ", "), "\n",
      "  Each row must correspond to exactly one cell."
    ))
  }
  if (!is.numeric(n.cores) || length(n.cores) != 1L || n.cores < 1L)
    stop("[create_Edgelist_Spatial] 'n.cores' must be a single positive integer.")

  n.cores    <- as.integer(n.cores)
  cell.names <- rownames(coord)
  n.cells    <- nrow(coord)
  xy         <- as.matrix(coord[, 1:2, drop = FALSE])
  colnames(xy) <- NULL

  x.coords <- setNames(xy[, 1L], cell.names)
  y.coords <- setNames(xy[, 2L], cell.names)

  if (n.cells < 2L)
    stop("[create_Edgelist_Spatial] 'coord' must contain at least 2 cells.")

  ## SECTION 2: Mode resolution ----

  .is_valid_k <- function(x)
    !is.null(x) && is.numeric(x) && length(x) == 1L &&
    !is.na(x) && is.finite(x) && x >= 1L && x == floor(x)

  .is_valid_rad <- function(x)
    !is.null(x) && is.numeric(x) && length(x) == 1L &&
    !is.na(x) && is.finite(x) && x > 0

  k.valid   <- .is_valid_k(k)
  rad.valid <- .is_valid_rad(rad)

  if (!k.valid && !rad.valid) {
    message("[create_Edgelist_Spatial] Neither 'k' nor 'rad' provided. Defaulting to k = 4.")
    k          <- 4L
    use.radius <- FALSE
  } else if (k.valid && rad.valid) {
    message("[create_Edgelist_Spatial] Both 'k' and 'rad' valid -- switching to radius mode.")
    use.radius <- TRUE
  } else if (rad.valid) {
    use.radius <- TRUE
  } else {
    k          <- as.integer(k)
    use.radius <- FALSE
  }

  ## SECTION 3: Socket cluster setup (all platforms) ----
  #
  # parallel::makeCluster spawns independent worker processes via sockets.
  # Unlike mclapply (fork-based, Unix only), socket clusters work on Windows,
  # macOS, and Linux without OS-specific guards or font-rendering instability.
  # on.exit() guarantees stopCluster() is called even if the function errors.

  if (n.cores > 1L) {
    cl <- parallel::makeCluster(n.cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterEvalQ(cl, library(data.table))
  }

  ## SECTION 4A -- Radius mode (dbscan::frNN) ----

  if (use.radius) {

    fnn       <- dbscan::frNN(xy, eps = rad, sort = FALSE)
    nb.list   <- fnn$id
    dist.list <- fnn$dist

    .build_radius_edges <- function(i) {
      nb    <- nb.list[[i]]
      dists <- dist.list[[i]]
      keep  <- nb != i
      nb    <- nb[keep]; dists <- dists[keep]
      if (length(nb) == 0L) return(NULL)
      data.table(from = cell.names[i], to = cell.names[nb], distance = dists)
    }

    if (n.cores > 1L) {
      parallel::clusterExport(cl, c("nb.list", "dist.list", "cell.names"),
                              envir = environment())
      edge.lst <- parallel::parLapply(cl, seq_len(n.cells), .build_radius_edges)
    } else {
      edge.lst <- lapply(seq_len(n.cells), .build_radius_edges)
    }
    edge.dt <- rbindlist(Filter(Negate(is.null), edge.lst))

  ## SECTION 4B -- kNN mode (RANN::nn2, mutual edges only) ----

  } else {

    k.query  <- min(k + 1L, n.cells)
    nn       <- RANN::nn2(data = xy, k = k.query)
    nn.idx   <- nn$nn.idx[,   -1L, drop = FALSE]
    nn.dists <- nn$nn.dists[, -1L, drop = FALSE]

    .build_knn_edges <- function(i) {
      nb    <- nn.idx[i, ]
      valid <- nb > 0L
      nb    <- nb[valid]
      if (length(nb) == 0L) return(NULL)
      data.table(from = i, to = nb, distance = nn.dists[i, valid])
    }

    if (n.cores > 1L) {
      parallel::clusterExport(cl, c("nn.idx", "nn.dists", "cell.names"),
                              envir = environment())
      edge.lst <- parallel::parLapply(cl, seq_len(n.cells), .build_knn_edges)
    } else {
      edge.lst <- lapply(seq_len(n.cells), .build_knn_edges)
    }
    edge.dt.idx <- rbindlist(Filter(Negate(is.null), edge.lst))

    # Mutual edge filter via set intersection (O(E log E))
    fwd.pairs    <- edge.dt.idx[, .(from, to)]
    rev.pairs    <- edge.dt.idx[, .(from = to, to = from)]
    mutual.pairs <- fintersect(fwd.pairs, rev.pairs)
    setkeyv(edge.dt.idx, c("from", "to"))
    setkeyv(mutual.pairs, c("from", "to"))
    edge.dt.idx <- edge.dt.idx[mutual.pairs, nomatch = 0L]

    edge.dt <- copy(edge.dt.idx)
    edge.dt[, from := cell.names[from]]
    edge.dt[, to   := cell.names[to]]
  }

  ## SECTION 5: Append autocrine edges & attach coordinates ----

  autocrine.dt <- data.table(from = cell.names, to = cell.names, distance = 0)
  edge.dt      <- rbindlist(list(edge.dt, autocrine.dt), use.names = TRUE)

  n.intercell <- nrow(edge.dt) - n.cells
  if (n.intercell == 0L)
    warning("[create_Edgelist_Spatial] No inter-cell edges found. ",
            "Consider increasing 'k' or 'rad'.")

  edge.dt[, x.from := x.coords[from]]
  edge.dt[, y.from := y.coords[from]]
  edge.dt[, x.to   := x.coords[to]]
  edge.dt[, y.to   := y.coords[to]]

  setcolorder(edge.dt, c("from", "to", "distance", "x.from", "y.from", "x.to", "y.to"))
  setkeyv(edge.dt, c("from", "to"))

  attr(edge.dt, "mode")       <- if (use.radius) "radius" else "kNN"
  attr(edge.dt, "param")      <- if (use.radius) rad else k
  attr(edge.dt, "n.cells")    <- n.cells
  attr(edge.dt, "cell.names") <- cell.names

  message(sprintf(
    "[create_Edgelist_Spatial] Mode: %-6s | Param: %s | Cells: %d | Inter-cell: %d | Autocrine: %d",
    attr(edge.dt, "mode"), attr(edge.dt, "param"),
    n.cells, n.intercell, n.cells
  ))

  invisible(edge.dt)
}
