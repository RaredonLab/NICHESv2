test_that("smooth_expression: returns matrix of same dimensions as input", {
  set.seed(20L)
  count.mtx <- Matrix::rsparsematrix(50L, 30L, density = 0.4)
  count.mtx <- abs(count.mtx) + 1
  rownames(count.mtx) <- paste0("g", 1:50)
  colnames(count.mtx) <- paste0("c", 1:30)

  out <- smooth_expression(count.mtx,
                           barcodes = c("c1","c2","c3"),
                           k = 5L, n.pcs = 10L)
  expect_equal(dim(out), dim(count.mtx))
  expect_equal(rownames(out), rownames(count.mtx))
  expect_equal(colnames(out), colnames(count.mtx))
})

test_that("smooth_expression: non-target columns unchanged", {
  set.seed(21L)
  count.mtx <- Matrix::rsparsematrix(40L, 25L, density = 0.5)
  count.mtx <- abs(count.mtx) + 1
  rownames(count.mtx) <- paste0("g", 1:40)
  colnames(count.mtx) <- paste0("c", 1:25)

  target.barcodes <- c("c1","c2")
  out <- smooth_expression(count.mtx, barcodes = target.barcodes,
                           k = 3L, n.pcs = 8L)

  untouched <- setdiff(colnames(count.mtx), target.barcodes)
  expect_equal(
    as.matrix(out[, untouched]),
    as.matrix(count.mtx[, untouched])
  )
})

test_that("smooth_expression: errors when matrix too small for any PCA component", {
  count.mtx <- Matrix::Matrix(matrix(1:4, nrow = 2L), sparse = TRUE)
  rownames(count.mtx) <- c("g1","g2")
  colnames(count.mtx) <- c("c1","c2")
  # n.pcs cap is min(n.pcs, n.cells - 1, n.genes - 1) = min(_, 1, 1) = 1
  # So PCA still runs with 1 PC. To force the error we need a 1x1 matrix.
  tiny <- Matrix::Matrix(matrix(5, nrow = 1L), sparse = TRUE)
  rownames(tiny) <- "g1"; colnames(tiny) <- "c1"
  expect_error(
    smooth_expression(tiny, barcodes = "c1", k = 1L, n.pcs = 5L),
    "Not enough cells or genes"
  )
})

test_that("smooth_expression: handles cells with zero library size without crashing", {
  set.seed(22L)
  count.mtx <- Matrix::rsparsematrix(20L, 15L, density = 0.4)
  count.mtx <- abs(count.mtx)
  rownames(count.mtx) <- paste0("g", 1:20)
  colnames(count.mtx) <- paste0("c", 1:15)
  # zero out one cell entirely
  count.mtx[, "c5"] <- 0

  expect_no_error({
    out <- smooth_expression(count.mtx, barcodes = c("c1","c2"),
                             k = 3L, n.pcs = 5L)
  })
  expect_equal(dim(out), dim(count.mtx))
})

test_that("smooth_expression: smoothed values are bounded by row min/max of input", {
  set.seed(23L)
  count.mtx <- Matrix::rsparsematrix(15L, 20L, density = 0.6)
  count.mtx <- abs(count.mtx) + 1
  rownames(count.mtx) <- paste0("g", 1:15)
  colnames(count.mtx) <- paste0("c", 1:20)

  out <- smooth_expression(count.mtx, barcodes = "c1",
                           k = 4L, n.pcs = 5L)

  # The smoothed column for c1 should equal the rowMeans of some 5 columns
  # of the original matrix (c1 itself + 4 neighbors). We don't know which 5
  # without recomputing PCA — but we can assert it equals SOME 5-column mean.
  # Weaker but still meaningful check: smoothed col is in the convex hull
  # of original column ranges row-wise.
  smoothed.col <- as.numeric(out[, "c1"])
  row.mins <- apply(as.matrix(count.mtx), 1L, min)
  row.maxs <- apply(as.matrix(count.mtx), 1L, max)
  expect_true(all(smoothed.col >= row.mins - 1e-9))
  expect_true(all(smoothed.col <= row.maxs + 1e-9))
})

test_that("smooth_expression: defaults k=9, n.pcs=30 are accepted", {
  set.seed(24L)
  count.mtx <- Matrix::rsparsematrix(50L, 40L, density = 0.4)
  count.mtx <- abs(count.mtx) + 1
  rownames(count.mtx) <- paste0("g", 1:50)
  colnames(count.mtx) <- paste0("c", 1:40)
  # irlba warns on small matrices when n.pcs is a large fraction of min(dim);
  # harmless here — this test only asserts that the defaults are accepted.
  expect_no_error(
    suppressWarnings(smooth_expression(count.mtx, barcodes = c("c1","c2")))
  )
})
