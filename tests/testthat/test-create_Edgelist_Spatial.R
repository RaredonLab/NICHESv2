# ==============================================================================
# FILE:    tests/testthat/test-create_Edgelist_Spatial.R
# ==============================================================================

test_that("create_Edgelist_Spatial: kNN mode returns correct structure", {
  .setup_fixtures()
  el <- create_Edgelist_Spatial(
    coord = data.frame(x = meta.data$x, y = meta.data$y,
                       row.names = rownames(meta.data)),
    k = 4L
  )
  expect_s3_class(el, "data.table")
  expect_true(all(c("from", "to", "distance", "x.from", "y.from", "x.to", "y.to") %in% colnames(el)))
  expect_equal(attr(el, "mode"), "kNN")
  expect_equal(attr(el, "param"), 4L)
})

test_that("create_Edgelist_Spatial: radius mode returns correct structure", {
  .setup_fixtures()
  el <- create_Edgelist_Spatial(
    coord = data.frame(x = meta.data$x, y = meta.data$y,
                       row.names = rownames(meta.data)),
    rad = 0.3
  )
  expect_equal(attr(el, "mode"), "radius")
})

test_that("create_Edgelist_Spatial: autocrine self-loops always present", {
  .setup_fixtures()
  el <- create_Edgelist_Spatial(
    coord = data.frame(x = meta.data$x, y = meta.data$y,
                       row.names = rownames(meta.data)),
    k = 4L
  )
  n.cells <- nrow(meta.data)
  autocrine <- el[el$from == el$to, ]
  expect_equal(nrow(autocrine), n.cells)
  expect_true(all(autocrine$distance == 0))
})

test_that("create_Edgelist_Spatial: kNN edges are mutual", {
  .setup_fixtures()
  el <- create_Edgelist_Spatial(
    coord = data.frame(x = meta.data$x, y = meta.data$y,
                       row.names = rownames(meta.data)),
    k = 4L
  )
  inter <- el[el$from != el$to, ]
  # Every edge A→B should have a corresponding B→A
  fwd <- paste(inter$from, inter$to)
  rev <- paste(inter$to,   inter$from)
  expect_true(all(fwd %in% rev))
})

test_that("create_Edgelist_Spatial: defaults to k=4 when no params given", {
  .setup_fixtures()
  coord.df <- data.frame(x = meta.data$x[1:10], y = meta.data$y[1:10],
                          row.names = rownames(meta.data)[1:10])
  expect_message(
    el <- create_Edgelist_Spatial(coord.df),
    regexp = "Defaulting to k = 4"
  )
  expect_equal(attr(el, "param"), 4L)
})

test_that("create_Edgelist_Spatial: errors on non-data.frame input", {
  expect_error(create_Edgelist_Spatial(matrix(1:10, 5)), regexp = "data.frame")
})

test_that("create_Edgelist_Spatial: errors on duplicate row names", {
  # data.frame() itself rejects duplicate row.names, so we use structure()
  # to attach them after construction — this delivers the bad input directly
  # to our function, which then catches it and emits a clear barcode error.
  bad.coord <- structure(data.frame(x = 1:3, y = 1:3), row.names = c("A", "A", "B"))
  expect_error(create_Edgelist_Spatial(bad.coord), regexp = "unique cell barcodes")
})

test_that("create_Edgelist_Spatial: coordinates attached correctly", {
  .setup_fixtures()
  coord.df <- data.frame(x = meta.data$x, y = meta.data$y,
                          row.names = rownames(meta.data))
  el <- create_Edgelist_Spatial(coord.df, k = 4L)
  # x.from should match coord x for the 'from' cell
  test.row <- el[el$from != el$to, ][1, ]
  expect_equal(test.row$x.from, coord.df[test.row$from, "x"])
  expect_equal(test.row$y.to,   coord.df[test.row$to,   "y"])
})
