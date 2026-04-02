test_that("returns a data.table with the correct four columns", {
  coord <- make_grid_coord(5)
  result <- Explore_Radius_Neighborhood(
    coord, r.min = 0.5, r.max = 3.0, r.step = 0.5,
    n.sample = 10L, seed = 1L, show.plot = FALSE
  )
  expect_s3_class(result, "data.table")
  expect_named(result, c("radius", "mean.neighbors", "sd.neighbors", "median.neighbors"))
})

test_that("number of rows matches the length of the radius sequence", {
  coord  <- make_grid_coord(5)
  r.seq  <- seq(0.5, 3.0, by = 0.5)
  result <- Explore_Radius_Neighborhood(
    coord, r.min = 0.5, r.max = 3.0, r.step = 0.5,
    n.sample = 10L, seed = 1L, show.plot = FALSE
  )
  expect_equal(nrow(result), length(r.seq))
  expect_equal(result$radius, r.seq)
})

test_that("mean.neighbors is non-decreasing as radius increases", {
  coord  <- make_grid_coord(8)
  result <- Explore_Radius_Neighborhood(
    coord, r.min = 0.5, r.max = 4.0, r.step = 0.25,
    n.sample = 20L, seed = 1L, show.plot = FALSE
  )
  diffs <- diff(result$mean.neighbors)
  expect_true(all(diffs >= -1e-10))          # allow floating-point noise
})

test_that("mean.neighbors at r.max is greater than at r.min for a non-trivial grid", {
  coord  <- make_grid_coord(6)
  result <- suppressWarnings(Explore_Radius_Neighborhood(
    coord, r.min = 0.5, r.max = 4.0, r.step = 0.5,
    n.sample = NULL, seed = 1L, show.plot = FALSE
  ))
  expect_gt(result$mean.neighbors[nrow(result)], result$mean.neighbors[1L])
})

test_that("mean.neighbors values are non-negative", {
  coord  <- make_grid_coord(5)
  result <- Explore_Radius_Neighborhood(
    coord, r.min = 0.5, r.max = 2.0, r.step = 0.5,
    n.sample = 15L, seed = 1L, show.plot = FALSE
  )
  expect_true(all(result$mean.neighbors >= 0))
  expect_true(all(result$sd.neighbors   >= 0 | is.na(result$sd.neighbors)))
})

test_that("n.sample = NULL triggers a warning mentioning NULL", {
  coord <- make_grid_coord(5)
  expect_warning(
    Explore_Radius_Neighborhood(
      coord, r.min = 0.5, r.max = 2.0, r.step = 0.5,
      n.sample = NULL, seed = 1L, show.plot = FALSE
    ),
    regexp = "NULL"
  )
})

test_that("n.sample > n.cells silently uses all cells without error", {
  coord <- make_grid_coord(4)   # 16 cells
  expect_no_warning(
    Explore_Radius_Neighborhood(
      coord, r.min = 0.5, r.max = 2.0, r.step = 0.5,
      n.sample = 10000L, seed = 1L, show.plot = FALSE
    )
  )
})

test_that("r.min >= r.max gives an informative error", {
  coord <- make_grid_coord(5)
  expect_error(
    Explore_Radius_Neighborhood(coord, r.min = 5, r.max = 1, show.plot = FALSE),
    regexp = "r.min"
  )
  expect_error(
    Explore_Radius_Neighborhood(coord, r.min = 3, r.max = 3, show.plot = FALSE),
    regexp = "r.min"
  )
})

test_that("r.step >= (r.max - r.min) gives an informative error", {
  coord <- make_grid_coord(5)
  expect_error(
    Explore_Radius_Neighborhood(
      coord, r.min = 1, r.max = 3, r.step = 10, show.plot = FALSE
    ),
    regexp = "r.step"
  )
  expect_error(
    Explore_Radius_Neighborhood(
      coord, r.min = 1, r.max = 3, r.step = 2, show.plot = FALSE
    ),
    regexp = "r.step"
  )
})

test_that("non-numeric coordinate columns give an informative error", {
  coord <- data.frame(
    x = letters[1:5], y = 1:5, stringsAsFactors = FALSE
  )
  rownames(coord) <- paste0("cell_", 1:5)
  expect_error(
    Explore_Radius_Neighborhood(coord, show.plot = FALSE),
    regexp = "numeric"
  )
})

test_that("coord with fewer than 2 columns gives an informative error", {
  coord <- data.frame(x = 1:5)
  rownames(coord) <- paste0("cell_", 1:5)
  expect_error(
    Explore_Radius_Neighborhood(coord, show.plot = FALSE),
    regexp = "two columns"
  )
})

test_that("non-data.frame coord gives an informative error", {
  expect_error(
    Explore_Radius_Neighborhood(matrix(1:10, nrow = 5), show.plot = FALSE),
    regexp = "data.frame"
  )
})

test_that("invalid n.sample values give an error", {
  coord <- make_grid_coord(5)
  expect_error(
    Explore_Radius_Neighborhood(coord, n.sample = -5L, show.plot = FALSE),
    regexp = "positive integer"
  )
  expect_error(
    Explore_Radius_Neighborhood(coord, n.sample = 0L, show.plot = FALSE),
    regexp = "positive integer"
  )
})

test_that("n.cores > 1 gives the same result as serial (n.cores = 1)", {
  coord <- make_grid_coord(7)
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

test_that("show.plot = TRUE returns the data.table invisibly", {
  coord <- make_grid_coord(4)
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

test_that("show.plot = FALSE returns the data.table visibly", {
  coord <- make_grid_coord(4)
  out <- withVisible(
    Explore_Radius_Neighborhood(
      coord, r.min = 0.5, r.max = 2.0, r.step = 0.5,
      n.sample = 8L, seed = 1L, show.plot = FALSE
    )
  )
  expect_true(out$visible)
  expect_s3_class(out$value, "data.table")
})

test_that("results are reproducible with the same seed", {
  coord  <- make_grid_coord(10)
  run.1  <- Explore_Radius_Neighborhood(
    coord, r.min = 0.5, r.max = 3.0, r.step = 0.5,
    n.sample = 30L, seed = 7L, show.plot = FALSE
  )
  run.2  <- Explore_Radius_Neighborhood(
    coord, r.min = 0.5, r.max = 3.0, r.step = 0.5,
    n.sample = 30L, seed = 7L, show.plot = FALSE
  )
  expect_identical(run.1, run.2)
})

test_that("auto-defaults produce a sensible radius range from bounding box", {
  coord  <- make_grid_coord(10)   # 10×10 grid, spacing 1 → bbox diagonal ≈ 12.73
  result <- Explore_Radius_Neighborhood(
    coord, n.sample = 20L, seed = 1L, show.plot = FALSE
  )
  # r.min default = 1% of diag ≈ 0.127; r.max = 20% ≈ 2.55
  # At r.max, interior cells definitely have neighbours, so max mean > 0
  expect_gt(max(result$mean.neighbors), 0)
  # Radius sequence should contain exactly 50 steps (default)
  expect_equal(nrow(result), 50L)
})


