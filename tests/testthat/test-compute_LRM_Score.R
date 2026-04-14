test_that("compute_LRM_Score: simple single-gene ligand and receptor (product)", {
  set.seed(1L)
  count.mtx <- Matrix::rsparsematrix(20L, 10L, density = 0.5)
  count.mtx <- abs(count.mtx) + 1  # ensure no zeros so product is non-trivial
  rownames(count.mtx) <- paste0("g", 1:20)
  colnames(count.mtx) <- paste0("c", 1:10)

  scores <- compute_LRM_Score(
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

  scores <- compute_LRM_Score(
    lig.str = "g1", rec.str = "g2",
    count.mtx = count.mtx,
    sending.cells = c("c1","c2"),
    receiv.cells  = c("c3","c4"),
    method = "minimum"
  )

  expected <- pmin(
    as.numeric(count.mtx["g1", c("c1","c2")]),
    as.numeric(count.mtx["g2", c("c3","c4")])
  )
  expect_equal(scores, expected)
})

test_that("compute_LRM_Score: complex ligand (underscore-joined) combines components", {
  set.seed(3L)
  count.mtx <- matrix(1:30, nrow = 6L, ncol = 5L)
  rownames(count.mtx) <- paste0("g", 1:6)
  colnames(count.mtx) <- paste0("c", 1:5)

  scores <- compute_LRM_Score(
    lig.str = "g1_g2",   # complex ligand
    rec.str = "g3",
    count.mtx = count.mtx,
    sending.cells = c("c1","c2"),
    receiv.cells  = c("c4","c5"),
    method = "product"
  )

  expected <- count.mtx["g1", c("c1","c2")] *
              count.mtx["g2", c("c1","c2")] *
              count.mtx["g3", c("c4","c5")]
  expect_equal(as.numeric(scores), as.numeric(expected))
})

test_that("compute_LRM_Score: zeros in input produce zero scores under product", {
  count.mtx <- matrix(c(0, 1, 1, 1), nrow = 2L)
  rownames(count.mtx) <- c("g1","g2")
  colnames(count.mtx) <- c("c1","c2")

  scores <- compute_LRM_Score(
    lig.str = "g1", rec.str = "g2",
    count.mtx = count.mtx,
    sending.cells = c("c1","c2"),
    receiv.cells  = c("c1","c2"),
    method = "product"
  )
  expect_equal(scores, c(0, 1))
})

test_that("compute_LRM_Score: default method is 'product'", {
  count.mtx <- matrix(c(2, 3, 4, 5), nrow = 2L)
  rownames(count.mtx) <- c("g1","g2")
  colnames(count.mtx) <- c("c1","c2")

  default <- compute_LRM_Score("g1","g2", count.mtx, "c1","c2")
  product <- compute_LRM_Score("g1","g2", count.mtx, "c1","c2", method = "product")
  expect_equal(default, product)
})

test_that("compute_LRM_Score: invalid method errors via match.arg", {
  count.mtx <- matrix(1:4, nrow = 2L,
                      dimnames = list(c("g1","g2"), c("c1","c2")))
  expect_error(
    compute_LRM_Score("g1","g2", count.mtx, "c1","c2", method = "garbage")
  )
})
