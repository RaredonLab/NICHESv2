# ==============================================================================
# FILE:    tests/testthat/test-compute_LRM_Score.R
# PURPOSE: Unit tests for the internal helper .compute_LRM_Score().
#          All inputs are synthesised inline — no shared fixtures needed.
#          Access via NICHESv2::: because the function is unexported.
# ==============================================================================

test_that("compute_LRM_Score: simple single-gene ligand and receptor (product)", {
  set.seed(1L)
  count.mtx <- Matrix::rsparsematrix(20L, 10L, density = 0.5)
  count.mtx <- abs(count.mtx) + 1
  rownames(count.mtx) <- paste0("g", 1:20)
  colnames(count.mtx) <- paste0("c", 1:10)
  
  scores <- NICHESv2:::.compute_LRM_Score(
    lig.str       = "g1",
    rec.str       = "g2",
    count.mtx     = count.mtx,
    sending.cells = c("c1", "c2", "c3"),
    receiv.cells  = c("c4", "c5", "c6"),
    method        = "product"
  )
  
  expected <- as.numeric(count.mtx["g1", c("c1","c2","c3")]) *
    as.numeric(count.mtx["g2", c("c4","c5","c6")])
  
  expect_length(scores, 3L)
  expect_equal(scores, expected)
})

test_that("compute_LRM_Score: minimum method takes element-wise pmin", {
  set.seed(2L)
  count.mtx <- Matrix::rsparsematrix(10L, 6L, density = 0.7)
  count.mtx <- abs(count.mtx) + 0.5
  rownames(count.mtx) <- paste0("g", 1:10)
  colnames(count.mtx) <- paste0("c", 1:6)
  
  scores <- NICHESv2:::.compute_LRM_Score(
    lig.str       = "g1",
    rec.str       = "g2",
    count.mtx     = count.mtx,
    sending.cells = c("c1", "c2"),
    receiv.cells  = c("c3", "c4"),
    method        = "minimum"
  )
  
  expected <- pmin(
    as.numeric(count.mtx["g1", c("c1","c2")]),
    as.numeric(count.mtx["g2", c("c3","c4")])
  )
  expect_equal(scores, expected)
})

test_that("compute_LRM_Score: complex ligand (underscore-joined) combines components via product", {
  set.seed(3L)
  count.mtx <- matrix(1:30, nrow = 6L, ncol = 5L)
  rownames(count.mtx) <- paste0("g", 1:6)
  colnames(count.mtx) <- paste0("c", 1:5)
  
  scores <- NICHESv2:::.compute_LRM_Score(
    lig.str       = "g1_g2",
    rec.str       = "g3",
    count.mtx     = count.mtx,
    sending.cells = c("c1", "c2"),
    receiv.cells  = c("c4", "c5"),
    method        = "product"
  )
  
  expected <- count.mtx["g1", c("c1","c2")] *
    count.mtx["g2", c("c1","c2")] *
    count.mtx["g3", c("c4","c5")]
  expect_equal(as.numeric(scores), as.numeric(expected))
})

test_that("compute_LRM_Score: complex receptor (underscore-joined) combines components via product", {
  set.seed(4L)
  count.mtx <- matrix(seq(2, 24, by = 2), nrow = 6L, ncol = 2L)
  rownames(count.mtx) <- paste0("g", 1:6)
  colnames(count.mtx) <- c("c1", "c2")
  
  scores <- NICHESv2:::.compute_LRM_Score(
    lig.str       = "g1",
    rec.str       = "g3_g4",
    count.mtx     = count.mtx,
    sending.cells = "c1",
    receiv.cells  = "c2",
    method        = "product"
  )
  
  # product = g1[c1] * g3[c2] * g4[c2]
  expected <- count.mtx["g1", "c1"] *
    count.mtx["g3", "c2"] *
    count.mtx["g4", "c2"]
  expect_equal(as.numeric(scores), as.numeric(expected))
})

test_that("compute_LRM_Score: complex receptor (underscore-joined) combines components via minimum", {
  set.seed(5L)
  count.mtx <- matrix(c(10, 3, 7, 2, 5, 9), nrow = 3L, ncol = 2L)
  rownames(count.mtx) <- c("g1", "g2", "g3")
  colnames(count.mtx) <- c("c1", "c2")
  
  scores <- NICHESv2:::.compute_LRM_Score(
    lig.str       = "g1",
    rec.str       = "g2_g3",
    count.mtx     = count.mtx,
    sending.cells = "c1",
    receiv.cells  = "c2",
    method        = "minimum"
  )
  
  # minimum of: g1[c1]=10, g2[c2]=5, g3[c2]=9  → 5
  expect_equal(as.numeric(scores), 5)
})

test_that("compute_LRM_Score: zeros in input produce zero scores under product", {
  count.mtx <- matrix(c(0, 1, 1, 1), nrow = 2L)
  rownames(count.mtx) <- c("g1", "g2")
  colnames(count.mtx) <- c("c1", "c2")
  
  scores <- NICHESv2:::.compute_LRM_Score(
    lig.str       = "g1",
    rec.str       = "g2",
    count.mtx     = count.mtx,
    sending.cells = c("c1", "c2"),
    receiv.cells  = c("c1", "c2"),
    method        = "product"
  )
  expect_equal(scores, c(0, 1))
})

test_that("compute_LRM_Score: default method is 'product'", {
  count.mtx <- matrix(c(2, 3, 4, 5), nrow = 2L)
  rownames(count.mtx) <- c("g1", "g2")
  colnames(count.mtx) <- c("c1", "c2")
  
  default <- NICHESv2:::.compute_LRM_Score("g1", "g2", count.mtx, "c1", "c2")
  product <- NICHESv2:::.compute_LRM_Score("g1", "g2", count.mtx, "c1", "c2",
                                           method = "product")
  expect_equal(default, product)
})

test_that("compute_LRM_Score: invalid method errors via match.arg", {
  count.mtx <- matrix(1:4, nrow = 2L,
                      dimnames = list(c("g1", "g2"), c("c1", "c2")))
  expect_error(
    NICHESv2:::.compute_LRM_Score("g1", "g2", count.mtx, "c1", "c2",
                                  method = "garbage")
  )
})