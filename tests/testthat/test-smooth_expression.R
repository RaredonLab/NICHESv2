# ==============================================================================
# FILE:    tests/testthat/test-smooth_expression.R
# PURPOSE: Unit tests for the internal helper .smooth_expression().
#          All inputs are synthesised inline — no shared fixtures needed.
#          Access via NICHESv2::: because the function is unexported.
# ==============================================================================

test_that("smooth_expression: returns matrix of same dimensions as input", {
  set.seed(20L)
  count.mtx <- Matrix::rsparsematrix(50L, 30L, density = 0.4)
  count.mtx <- abs(count.mtx) + 1
  rownames(count.mtx) <- paste0("g", 1:50)
  colnames(count.mtx) <- paste0("c", 1:30)
  
  out <- NICHESv2:::.smooth_expression(count.mtx,
                                       barcodes = c("c1", "c2", "c3"),
                                       k = 5L, n.pcs = 10L)
  expect_equal(dim(out), dim(count.mtx))
  expect_equal(rownames(out), rownames(count.mtx))
  expect_equal(colnames(out), colnames(count.mtx))
})

test_that("smooth_expression: non-target columns are unchanged", {
  set.seed(21L)
  count.mtx <- Matrix::rsparsematrix(40L, 25L, density = 0.5)
  count.mtx <- abs(count.mtx) + 1
  rownames(count.mtx) <- paste0("g", 1:40)
  colnames(count.mtx) <- paste0("c", 1:25)
  
  target.barcodes <- c("c1", "c2")
  out <- NICHESv2:::.smooth_expression(count.mtx, barcodes = target.barcodes,
                                       k = 3L, n.pcs = 8L)
  
  untouched <- setdiff(colnames(count.mtx), target.barcodes)
  expect_equal(
    as.matrix(out[, untouched]),
    as.matrix(count.mtx[, untouched])
  )
})

test_that("smooth_expression: errors when matrix is too small for any PCA component", {
  # n.pcs is capped at min(n.pcs, n.cells - 1, n.genes - 1).
  # A 1x1 matrix gives cap = min(_, 0, 0) = 0, which must error.
  tiny <- Matrix::Matrix(matrix(5, nrow = 1L), sparse = TRUE)
  rownames(tiny) <- "g1"
  colnames(tiny) <- "c1"
  expect_error(
    NICHESv2:::.smooth_expression(tiny, barcodes = "c1", k = 1L, n.pcs = 5L),
    "Not enough cells or genes"
  )
})

test_that("smooth_expression: handles cells with zero library size without crashing", {
  set.seed(22L)
  count.mtx <- Matrix::rsparsematrix(20L, 15L, density = 0.4)
  count.mtx <- abs(count.mtx)
  rownames(count.mtx) <- paste0("g", 1:20)
  colnames(count.mtx) <- paste0("c", 1:15)
  count.mtx[, "c5"] <- 0  # zero out one cell entirely
  
  expect_no_error({
    out <- NICHESv2:::.smooth_expression(count.mtx, barcodes = c("c1", "c2"),
                                         k = 3L, n.pcs = 5L)
  })
  expect_equal(dim(out), dim(count.mtx))
})

test_that("smooth_expression: smoothed values are within the row-wise range of input", {
  set.seed(23L)
  count.mtx <- Matrix::rsparsematrix(15L, 20L, density = 0.6)
  count.mtx <- abs(count.mtx) + 1
  rownames(count.mtx) <- paste0("g", 1:15)
  colnames(count.mtx) <- paste0("c", 1:20)
  
  out <- NICHESv2:::.smooth_expression(count.mtx, barcodes = "c1",
                                       k = 4L, n.pcs = 5L)
  
  # Smoothing averages a cell's neighborhood in PCA space, so the result must
  # lie within the convex hull of the original values row-wise.
  smoothed.col <- as.numeric(out[, "c1"])
  row.mins <- apply(as.matrix(count.mtx), 1L, min)
  row.maxs <- apply(as.matrix(count.mtx), 1L, max)
  expect_true(all(smoothed.col >= row.mins - 1e-9))
  expect_true(all(smoothed.col <= row.maxs + 1e-9))
})

test_that("smooth_expression: default arguments k=9 and n.pcs=30 are accepted", {
  set.seed(24L)
  count.mtx <- Matrix::rsparsematrix(50L, 40L, density = 0.4)
  count.mtx <- abs(count.mtx) + 1
  rownames(count.mtx) <- paste0("g", 1:50)
  colnames(count.mtx) <- paste0("c", 1:40)
  
  # irlba warns when n.pcs is a large fraction of min(nrow, ncol) on a small
  # matrix — this is a known benign side-effect, not the behaviour under test.
  expect_no_error(
    suppressWarnings(
      NICHESv2:::.smooth_expression(count.mtx, barcodes = c("c1", "c2"))
    )
  )
})