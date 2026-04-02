#' Explore the relationship between search radius and neighbourhood size
#'
#' Estimates the mean number of neighbouring cells as a function of search
#' radius, using a randomly sampled subset of cells for efficiency. Designed
#' to help users choose a radius parameter before calling
#' \code{create_NICHESObject()}.
#'
#' A single kd-tree query via \code{dbscan::frNN()} is run at \code{r.max};
#' all smaller radii are evaluated by thresholding the already-retrieved
#' distance lists, so the kd-tree is built only once regardless of how many
#' radius steps are requested.
#'
#' @param coord data.frame. First two columns must be numeric x and y
#'   coordinates. Row names must be unique cell barcodes.
#' @param r.min numeric(1). Minimum radius to evaluate. Default: 1\% of the
#'   bounding-box diagonal.
#' @param r.max numeric(1). Maximum radius to evaluate. Default: 20\% of the
#'   bounding-box diagonal.
#' @param r.step numeric(1). Step size between evaluated radii. Default:
#'   \code{(r.max - r.min) / 50} (50 evenly spaced evaluation points).
#' @param n.sample integer(1) or NULL. Number of cells to randomly sample for
#'   the neighbourhood estimate. Default 2000. Set to \code{NULL} to use all
#'   cells — not recommended for large datasets.
#' @param seed integer(1) or NULL. Random seed for reproducible sampling.
#'   Default 42.
#' @param n.cores integer(1). Number of parallel workers (socket cluster) for
#'   distributing the radius evaluation. Default 1.
#' @param show.plot logical. Whether to print the summary ggplot before
#'   returning. Default \code{TRUE}.
#'
#' @return A \code{data.table} with columns:
#'   \describe{
#'     \item{radius}{evaluated radius value}
#'     \item{mean.neighbors}{mean neighbour count across sampled cells}
#'     \item{sd.neighbors}{standard deviation of neighbour count}
#'     \item{median.neighbors}{median neighbour count (robust to density outliers)}
#'   }
#'   Returned invisibly when \code{show.plot = TRUE}.
#'
#' @export
Explore_Radius_Neighborhood <- function(
  coord,
  r.min     = NULL,
  r.max     = NULL,
  r.step    = NULL,
  n.sample  = 2000L,
  seed      = 42L,
  n.cores   = 1L,
  show.plot = TRUE
) {

  # ── 1. Validate coord ───────────────────────────────────────────────────────
  if (!is.data.frame(coord)) {
    stop(
      "coord must be a data.frame with row names as cell barcodes and ",
      "the first two columns as numeric x and y coordinates."
    )
  }
  if (ncol(coord) < 2L) {
    stop("coord must have at least two columns (x and y).")
  }
  if (!is.numeric(coord[[1L]]) || !is.numeric(coord[[2L]])) {
    stop("The first two columns of coord must be numeric (x and y coordinates).")
  }
  n.cells <- nrow(coord)
  if (n.cells < 2L) {
    stop("coord must contain at least 2 cells.")
  }

  # ── 2. Sampling ─────────────────────────────────────────────────────────────
  if (is.null(n.sample)) {
    warning(
      "n.sample = NULL: neighbourhood sizes will be estimated from ALL ",
      n.cells, " cells. This is not recommended for large datasets ",
      "and may be very slow."
    )
    sample.idx <- seq_len(n.cells)

  } else {
    n.sample <- suppressWarnings(as.integer(n.sample))
    if (is.na(n.sample) || n.sample <= 0L) {
      stop("n.sample must be a positive integer.")
    }
    if (n.sample >= n.cells) {
      # Requested more than available — silently use all cells.
      sample.idx <- seq_len(n.cells)
    } else {
      if (!is.null(seed)) set.seed(as.integer(seed))
      sample.idx <- sort(sample.int(n.cells, n.sample))
    }
  }
  n.sample.actual <- length(sample.idx)

  # ── 3. Coordinate matrices ──────────────────────────────────────────────────
  coords.mat    <- as.matrix(coord[, 1:2, drop = FALSE])
  sample.coords <- coords.mat[sample.idx, , drop = FALSE]

  # ── 4. Auto-compute radius defaults from bounding box ──────────────────────
  x.range   <- range(coords.mat[, 1L])
  y.range   <- range(coords.mat[, 2L])
  bbox.diag <- sqrt(diff(x.range)^2 + diff(y.range)^2)

  if (is.null(r.min)) r.min <- bbox.diag * 0.01
  if (is.null(r.max)) r.max <- bbox.diag * 0.20

  # ── 5. Validate radius parameters ──────────────────────────────────────────
  if (!is.numeric(r.min) || length(r.min) != 1L || is.na(r.min) || r.min < 0) {
    stop("r.min must be a single non-negative numeric value.")
  }
  if (!is.numeric(r.max) || length(r.max) != 1L || is.na(r.max) || r.max <= 0) {
    stop("r.max must be a single positive numeric value.")
  }
  if (r.min >= r.max) {
    stop(
      "r.min (", round(r.min, 4L), ") must be strictly less than ",
      "r.max (", round(r.max, 4L), ")."
    )
  }

  if (is.null(r.step)) r.step <- (r.max - r.min) / 49

  if (!is.numeric(r.step) || length(r.step) != 1L || is.na(r.step) || r.step <= 0) {
    stop("r.step must be a single positive numeric value.")
  }
  if (r.step >= (r.max - r.min)) {
    stop(
      "r.step (", round(r.step, 4L), ") must be smaller than ",
      "r.max - r.min (", round(r.max - r.min, 4L), ")."
    )
  }

  # ── 6. Build radius sequence ────────────────────────────────────────────────
  r.seq   <- seq(r.min, r.max, by = r.step)
  n.radii <- length(r.seq)

  if (n.radii > 10000L) {
    warning(
      "The radius sequence contains ", n.radii, " evaluation points. ",
      "Consider increasing r.step to reduce computation time."
    )
  }

  # ── 7. Single frNN query at r.max ───────────────────────────────────────────
  # Query sampled cells against ALL cells within r.max.
  # frNN returns distances sorted ascending (sort = TRUE), enabling
  # fast O(log k) threshold counting via findInterval() below.
  frnn.result <- dbscan::frNN(
    x     = coords.mat,
    eps   = r.max,
    query = sample.coords,
    sort  = TRUE
  )

  # Exclude self-neighbor: sample cell i is row sample.idx[i] in coords.mat.
  dists.clean <- vector("list", n.sample.actual)
  for (i in seq_len(n.sample.actual)) {
    ids              <- frnn.result$id[[i]]
    ds               <- frnn.result$dist[[i]]
    dists.clean[[i]] <- ds[ids != sample.idx[i]]
  }

  # ── 8. Count neighbors at each radius ──────────────────────────────────────
  # For each sampled cell and radius r: count of distances <= r.
  # findInterval(r, d) on a sorted d gives exactly this in O(log k).
  # Result: n.sample × n.radii integer matrix.

  n.workers <- min(as.integer(n.cores), n.radii)

  if (n.workers > 1L) {
    cl <- parallel::makeCluster(n.workers)
    on.exit(parallel::stopCluster(cl), add = TRUE)

    parallel::clusterEvalQ(cl, {
      data.table::setDTthreads(1L)
      NULL
    })
    parallel::clusterExport(cl, varlist = "dists.clean", envir = environment())

    r.chunks <- split(
      r.seq,
      cut(seq_along(r.seq), breaks = n.workers, labels = FALSE)
    )

    chunk.mats <- parallel::parLapply(cl, r.chunks, function(r.chunk) {
      # Returns n.sample × length(r.chunk) integer matrix
      t(vapply(
        dists.clean,
        function(d) vapply(r.chunk, function(r) findInterval(r, d), integer(1L)),
        integer(length(r.chunk))
      ))
    })

    count.mat <- do.call(cbind, chunk.mats)

  } else {
    count.mat <- t(vapply(
      dists.clean,
      function(d) vapply(r.seq, function(r) findInterval(r, d), integer(1L)),
      integer(n.radii)
    ))
  }

  # ── 9. Summarise across sampled cells ──────────────────────────────────────
  result.dt <- data.table::data.table(
    radius           = r.seq,
    mean.neighbors   = colMeans(count.mat),
    sd.neighbors     = apply(count.mat, 2L, sd),
    median.neighbors = apply(count.mat, 2L, median)
  )

  # ── 10. Plot ────────────────────────────────────────────────────────────────
  if (show.plot) {
    p <- ggplot2::ggplot(
      result.dt,
      ggplot2::aes(x = radius, y = mean.neighbors)
    ) +
      ggplot2::geom_ribbon(
        ggplot2::aes(
          ymin = pmax(0, mean.neighbors - sd.neighbors),
          ymax = mean.neighbors + sd.neighbors
        ),
        fill  = "#4C72B0",
        alpha = 0.15
      ) +
      ggplot2::geom_line(color = "#4C72B0", linewidth = 0.9) +
      ggplot2::labs(
        x        = "Radius",
        y        = "Mean Neighbour Count",
        title    = "Radius vs. Neighbourhood Size",
        subtitle = paste0(
          "Estimated from ", n.sample.actual, " sampled cells  \u00b7  ",
          n.radii, " radius steps"
        )
      ) +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::theme(
        panel.grid.minor = ggplot2::element_blank(),
        plot.subtitle    = ggplot2::element_text(color = "grey50", size = 9)
      )
    print(p)
    return(invisible(result.dt))
  }

  result.dt
}
