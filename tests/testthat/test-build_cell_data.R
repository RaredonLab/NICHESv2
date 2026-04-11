test_that("build_cell_data: basic structure (cell, gene, score, score.norm columns)", {
  set.seed(10L)
  count.mtx <- Matrix::rsparsematrix(8L, 5L, density = 0.5)
  count.mtx <- abs(count.mtx) + 1
  rownames(count.mtx) <- paste0("g", 1:8)
  colnames(count.mtx) <- paste0("c", 1:5)

  cd <- build_cell_data(count.mtx, normalize.method = "prop")

  expect_s3_class(cd, "data.table")
  expect_equal(sort(names(cd)), sort(c("cell","gene","score","score.norm")))
  expect_equal(data.table::key(cd), c("cell","gene"))
  expect_true(all(cd$score > 0))  # zero rows must be dropped
})

test_that("build_cell_data: prop normalization sums to 1 within each cell", {
  count.mtx <- Matrix::Matrix(
    c(1, 2, 3,
      4, 5, 6,
      7, 8, 9),
    nrow = 3L, ncol = 3L, sparse = TRUE
  )
  rownames(count.mtx) <- c("g1","g2","g3")
  colnames(count.mtx) <- c("c1","c2","c3")

  cd <- build_cell_data(count.mtx, normalize.method = "prop")
  per.cell.sums <- cd[, sum(score.norm), by = cell]
  expect_equal(per.cell.sums$V1, rep(1, 3L))
})

test_that("build_cell_data: non-prop normalization leaves score.norm == score", {
  count.mtx <- Matrix::Matrix(c(1,2,3,4), nrow = 2L,
                              dimnames = list(c("g1","g2"), c("c1","c2")),
                              sparse = TRUE)
  cd <- build_cell_data(count.mtx, normalize.method = "none")
  expect_equal(cd$score, cd$score.norm)
})

test_that("build_cell_data: zero-only rows dropped, sparse-zero entries absent", {
  count.mtx <- Matrix::Matrix(
    matrix(c(0, 0, 0,
             1, 0, 2),
           nrow = 2L, ncol = 3L, byrow = TRUE),
    sparse = TRUE
  )
  rownames(count.mtx) <- c("g1","g2")
  colnames(count.mtx) <- c("c1","c2","c3")
  
  cd <- build_cell_data(count.mtx, normalize.method = "prop")
  # g1 should not appear at all (all zeros); only nonzero g2 entries present
  expect_false("g1" %in% cd$gene)
  expect_equal(nrow(cd), 2L)  # g2 in c1 and c3 only
  expect_setequal(cd$cell, c("c1","c3"))
})

test_that("build_cell_data: default normalize.method is 'prop'", {
  count.mtx <- Matrix::Matrix(c(1,2,3,4), nrow = 2L,
                              dimnames = list(c("g1","g2"), c("c1","c2")),
                              sparse = TRUE)
  cd.default <- build_cell_data(count.mtx)
  cd.prop    <- build_cell_data(count.mtx, normalize.method = "prop")
  expect_equal(cd.default, cd.prop)
})
