# ==============================================================================
# FILE:    tests/testthat/test-Explore_Radius_Neighborhood.R
# PURPOSE: Unit tests for Explore_Radius_Neighborhood().
#          All inputs synthesised via make_grid_coord() from helper-testutils.R
#          — no shared fixtures needed.
# ==============================================================================

# ── Output structure ───────────────────────────────────────────────────────────

test_that("Explore_Radius_Neighborhood: returns data.table with four required columns", {
  coord  <- make_grid_coord(5L)
  result <- Explore_Radius_Neighborhood(
    coord, r.min = 0.5, r.max = 3.0, r.step = 0.5,
    n.sample = 10L, seed = 1L, show.plot = FALSE
  )
  expect_s3_class(result, "data.table")
  expect_named(result, c("radius", "mean.neighbors",
                         "sd.neighbors", "median.neighbors"))
})

test_that("Explore_Radius_Neighborhood: number of rows matches radius sequence length", {
  coord  <- make_grid_coord(5L)
  r.seq  <- seq(0.5, 3.0, by = 0.5)
  result <- Explore_Radius_Neighborhood(
    coord, r.min = 0.5, r.max = 3.0, r.step = 0.5,
    n.sample = 10L, seed = 1L, show.plot = FALSE
  )
  expect_equal(nrow(result), length(r.seq))
  expect_equal(result$radius, r.seq)
})

test_that("Explore_Radius_Neighborhood: mean.neighbors is non-decreasing as radius grows", {
  coord  <- make_grid_coord(8L)
  result <- Explore_Radius_Neighborhood(
    coord, r.min = 0.5, r.max = 4.0, r.step = 0.25,
    n.sample = 20L, seed = 1L, show.plot = FALSE
  )
  expect_true(all(diff(result$mean.neighbors) >= -1e-10))
})

test_that("Explore_Radius_Neighborhood: mean.neighbors at r.max exceeds r.min on non-trivial grid", {
  coord  <- make_grid_coord(6L)
  result <- suppressWarnings(Explore_Radius_Neighborhood(
    coord, r.min = 0.5, r.max = 4.0, r.step = 0.5,
    n.sample = NULL, seed = 1L, show.plot = FALSE
  ))
  expect_gt(result$mean.neighbors[nrow(result)], result$mean.neighbors[1L])
})

test_that("Explore_Radius_Neighborhood: all neighbor statistics are non-negative", {
  coord  <- make_grid_coord(5L)
  result <- Explore_Radius_Neighborhood(
    coord, r.min = 0.5, r.max = 2.0, r.step = 0.5,
    n.sample = 15L, seed = 1L, show.plot = FALSE
  )
  expect_true(all(result$mean.neighbors >= 0))
  expect_true(all(result$sd.neighbors   >= 0 | is.na(result$sd.neighbors)))
})

# ── Input validation ───────────────────────────────────────────────────────────

test_that("Explore_Radius_Neighborhood: non-data.frame coord gives informative error", {
  expect_error(
    Explore_Radius_Neighborhood(matrix(1:10, nrow = 5L), show.plot = FALSE),
    regexp = "data.frame"
  )
})

test_that("Explore_Radius_Neighborhood: coord with fewer than 2 columns gives informative error", {
  coord <- data.frame(x = 1:5)
  rownames(coord) <- paste0("cell_", 1:5)
  expect_error(
    Explore_Radius_Neighborhood(coord, show.plot = FALSE),
    regexp = "two columns"
  )
})

test_that("Explore_Radius_Neighborhood: non-numeric coordinate columns give informative error", {
  coord <- data.frame(x = letters[1:5], y = 1:5, stringsAsFactors = FALSE)
  rownames(coord) <- paste0("cell_", 1:5)
  expect_error(
    Explore_Radius_Neighborhood(coord, show.plot = FALSE),
    regexp = "numeric"
  )
})

test_that("Explore_Radius_Neighborhood: r.min >= r.max gives informative error", {
  coord <- make_grid_coord(5L)
  expect_error(
    Explore_Radius_Neighborhood(coord, r.min = 5, r.max = 1, show.plot = FALSE),
    regexp = "r.min"
  )
  expect_error(
    Explore_Radius_Neighborhood(coord, r.min = 3, r.max = 3, show.plot = FALSE),
    regexp = "r.min"
  )
})

test_that("Explore_Radius_Neighborhood: r.step >= (r.max - r.min) gives informative error", {
  coord <- make_grid_coord(5L)
  expect_error(
    Explore_Radius_Neighborhood(coord, r.min = 1, r.max = 3, r.step = 10,
                                show.plot = FALSE),
    regexp = "r.step"
  )
  expect_error(
    Explore_Radius_Neighborhood(coord, r.min = 1, r.max = 3, r.step = 2,
                                show.plot = FALSE),
    regexp = "r.step"
  )
})

test_that("Explore_Radius_Neighborhood: invalid n.sample values give informative error", {
  coord <- make_grid_coord(5L)
  expect_error(
    Explore_Radius_Neighborhood(coord, n.sample = -5L, show.plot = FALSE),
    regexp = "positive integer"
  )
  expect_error(
    Explore_Radius_Neighborhood(coord, n.sample = 0L, show.plot = FALSE),
    regexp = "positive integer"
  )
})

test_that("Explore_Radius_Neighborhood: n.sample = NULL triggers a warning mentioning NULL", {
  coord <- make_grid_coord(5L)
  expect_warning(
    Explore_Radius_Neighborhood(
      coord, r.min = 0.5, r.max = 2.0, r.step = 0.5,
      n.sample = NULL, seed = 1L, show.plot = FALSE
    ),
    regexp = "NULL"
  )
})

test_that("Explore_Radius_Neighborhood: n.sample > n.cells uses all cells without warning", {
  coord <- make_grid_coord(4L)   # 16 cells
  expect_no_warning(
    Explore_Radius_Neighborhood(
      coord, r.min = 0.5, r.max = 2.0, r.step = 0.5,
      n.sample = 10000L, seed = 1L, show.plot = FALSE
    )
  )
})

# ── Reproducibility ────────────────────────────────────────────────────────────

test_that("Explore_Radius_Neighborhood: same seed gives identical results", {
  coord <- make_grid_coord(10L)
  run.1 <- Explore_Radius_Neighborhood(
    coord, r.min = 0.5, r.max = 3.0, r.step = 0.5,
    n.sample = 30L, seed = 7L, show.plot = FALSE
  )
  run.2 <- Explore_Radius_Neighborhood(
    coord, r.min = 0.5, r.max = 3.0, r.step = 0.5,
    n.sample = 30L, seed = 7L, show.plot = FALSE
  )
  expect_identical(run.1, run.2)
})

# ── show.plot visibility contract ─────────────────────────────────────────────

test_that("Explore_Radius_Neighborhood: show.plot=TRUE returns data.table invisibly", {
  coord <- make_grid_coord(4L)
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  out <- withVisible(
    Explore_Radius_Neighborhood(
      coord, r.min = 0.5, r.max = 2.0, r.step = 0.5,
      n.sample = 8L, seed = 1L, show.plot = TRUE
    )
  )
  expect_false(out$visible)
  expect_s3_class(out$value, "data.table")
})

test_that("Explore_Radius_Neighborhood: show.plot=FALSE returns data.table visibly", {
  coord <- make_grid_coord(4L)
  out <- withVisible(
    Explore_Radius_Neighborhood(
      coord, r.min = 0.5, r.max = 2.0, r.step = 0.5,
      n.sample = 8L, seed = 1L, show.plot = FALSE
    )
  )
  expect_true(out$visible)
  expect_s3_class(out$value, "data.table")
})

# ── Auto-default radius range ──────────────────────────────────────────────────

test_that("Explore_Radius_Neighborhood: auto-defaults produce 50-row result with positive neighbors", {
  coord  <- make_grid_coord(10L)  # 10x10 grid, spacing 1 → bbox diagonal ≈ 12.73
  result <- Explore_Radius_Neighborhood(
    coord, n.sample = 20L, seed = 1L, show.plot = FALSE
  )
  expect_equal(nrow(result), 50L)
  expect_gt(max(result$mean.neighbors), 0)
})

# ── Parallel path (skip on CRAN) ───────────────────────────────────────────────

test_that("Explore_Radius_Neighborhood: n.cores=2 gives same result as n.cores=1", {
  skip_on_cran()
  if (parallel::detectCores() < 2L) skip("only one core available")
  
  coord    <- make_grid_coord(7L)
  result.s <- Explore_Radius_Neighborhood(
    coord, r.min = 0.5, r.max = 3.0, r.step = 0.5,
    n.sample = 20L, seed = 42L, n.cores = 1L, show.plot = FALSE
  )
  result.p <- Explore_Radius_Neighborhood(
    coord, r.min = 0.5, r.max = 3.0, r.step = 0.5,
    n.sample = 20L, seed = 42L, n.cores = 2L, show.plot = FALSE
  )
  expect_equal(result.s$mean.neighbors,   result.p$mean.neighbors,   tolerance = 1e-10)
  expect_equal(result.s$sd.neighbors,     result.p$sd.neighbors,     tolerance = 1e-10)
  expect_equal(result.s$median.neighbors, result.p$median.neighbors, tolerance = 1e-10)
})