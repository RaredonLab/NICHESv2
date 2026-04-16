# ==============================================================================
# FILE:    tests/testthat/test-create_Edgelist_Spatial.R
# PURPOSE: Unit tests for create_Edgelist_Spatial().
#          Coordinates synthesised via make_grid_coord() from helper-testutils.R
#          — no shared fixtures needed.
# ==============================================================================

# ── Structure and attributes ───────────────────────────────────────────────────

test_that("create_Edgelist_Spatial: kNN mode returns correct structure", {
  coord <- make_grid_coord(5L)   # 25-cell regular grid
  el    <- create_Edgelist_Spatial(coord, k = 4L)
  
  expect_s3_class(el, "data.table")
  expect_true(all(c("from", "to", "distance",
                    "x.from", "y.from", "x.to", "y.to") %in% colnames(el)))
  expect_equal(attr(el, "mode"),  "kNN")
  expect_equal(attr(el, "param"), 4L)
})

test_that("create_Edgelist_Spatial: radius mode returns correct structure", {
  coord <- make_grid_coord(5L)
  el    <- create_Edgelist_Spatial(coord, rad = 1.5)
  
  expect_s3_class(el, "data.table")
  expect_equal(attr(el, "mode"), "radius")
})

test_that("create_Edgelist_Spatial: autocrine self-loops present for every cell", {
  coord   <- make_grid_coord(5L)
  el      <- create_Edgelist_Spatial(coord, k = 4L)
  n.cells <- nrow(coord)
  
  autocrine <- el[el$from == el$to, ]
  expect_equal(nrow(autocrine), n.cells)
  expect_true(all(autocrine$distance == 0))
})

test_that("create_Edgelist_Spatial: kNN edges are mutual (A->B implies B->A)", {
  coord <- make_grid_coord(5L)
  el    <- create_Edgelist_Spatial(coord, k = 4L)
  inter <- el[el$from != el$to, ]
  
  fwd <- paste(inter$from, inter$to)
  rev <- paste(inter$to,   inter$from)
  expect_true(all(fwd %in% rev))
})

test_that("create_Edgelist_Spatial: defaults to k=4 when no params given", {
  coord <- make_grid_coord(3L)   # 9-cell grid — small enough to be fast
  expect_message(
    el <- create_Edgelist_Spatial(coord),
    regexp = "Defaulting to k = 4"
  )
  expect_equal(attr(el, "param"), 4L)
})

test_that("create_Edgelist_Spatial: coordinates attached correctly to edges", {
  coord <- make_grid_coord(5L)
  el    <- create_Edgelist_Spatial(coord, k = 4L)
  
  # Pick any inter-cell edge and verify x.from / y.to match the coord lookup
  test.row <- el[el$from != el$to, ][1L, ]
  expect_equal(test.row$x.from, coord[test.row$from, "x"])
  expect_equal(test.row$y.to,   coord[test.row$to,   "y"])
})

# ── Input validation ───────────────────────────────────────────────────────────

test_that("create_Edgelist_Spatial: errors on non-data.frame input", {
  expect_error(
    create_Edgelist_Spatial(matrix(1:10, 5L)),
    regexp = "data.frame"
  )
})

test_that("create_Edgelist_Spatial: errors on duplicate row names", {
  # data.frame() rejects duplicate rownames itself, so we attach them via
  # structure() to deliver the bad input directly to the function under test.
  bad.coord <- structure(
    data.frame(x = 1:3, y = 1:3),
    row.names = c("A", "A", "B")
  )
  expect_error(
    create_Edgelist_Spatial(bad.coord),
    regexp = "unique cell barcodes"
  )
})