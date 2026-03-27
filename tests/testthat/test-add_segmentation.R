# ==============================================================================
# FILE:    tests/testthat/test-add_segmentation.R
#
# All tests skip automatically if 'sf' is not installed, since sf is in
# Suggests (optional). Tests build minimal sf objects inline — no changes to
# generate_testdata.R or re-running it are needed.
# ==============================================================================


# ------------------------------------------------------------------------------
# Helper: build a minimal sf object for a given character vector of cell names.
# Each cell gets a tiny unit square polygon — the geometry is irrelevant for
# these tests; we only care that the slot is correctly stored and filtered.
# ------------------------------------------------------------------------------

.make_seg_sf <- function(cells) {
  polys <- lapply(seq_along(cells), function(i) {
    sf::st_polygon(list(matrix(
      c(i,   i,
        i+1, i,
        i+1, i+1,
        i,   i+1,
        i,   i),
      ncol = 2, byrow = TRUE
    )))
  })
  sf::st_sf(
    cell     = cells,
    geometry = sf::st_sfc(polys)
  )
}


# ------------------------------------------------------------------------------
# add_segmentation: basic attachment
# ------------------------------------------------------------------------------

test_that("add_segmentation: returns a NICHESObject", {
  skip_if_not_installed("sf")
  .setup_fixtures()
  seg.sf  <- .make_seg_sf(niches.obj$cell.list)
  result  <- add_segmentation(niches.obj, seg.sf)
  expect_true(is_NICHESObject(result))
})

test_that("add_segmentation: $segmentation is an sf object after attachment", {
  skip_if_not_installed("sf")
  .setup_fixtures()
  seg.sf <- .make_seg_sf(niches.obj$cell.list)
  result <- add_segmentation(niches.obj, seg.sf)
  expect_true(inherits(result$segmentation, "sf"))
})

test_that("add_segmentation: $segmentation has correct number of rows", {
  skip_if_not_installed("sf")
  .setup_fixtures()
  seg.sf <- .make_seg_sf(niches.obj$cell.list)
  result <- add_segmentation(niches.obj, seg.sf)
  expect_equal(nrow(result$segmentation), length(niches.obj$cell.list))
})

test_that("add_segmentation: $segmentation is NULL before attachment", {
  .setup_fixtures()
  expect_null(niches.obj$segmentation)
})


# ------------------------------------------------------------------------------
# add_segmentation: input validation errors
# ------------------------------------------------------------------------------

test_that("add_segmentation: errors on non-NICHESObject input", {
  skip_if_not_installed("sf")
  seg.sf <- .make_seg_sf(c("CellA", "CellB"))
  expect_error(
    add_segmentation(list(), seg.sf),
    regexp = "NICHESObject"
  )
})

test_that("add_segmentation: errors if seg.sf is not an sf object", {
  skip_if_not_installed("sf")
  .setup_fixtures()
  expect_error(
    add_segmentation(niches.obj, data.frame(cell = niches.obj$cell.list)),
    regexp = "sf object"
  )
})

test_that("add_segmentation: errors if seg.sf has no 'cell' column", {
  skip_if_not_installed("sf")
  .setup_fixtures()
  # Build sf with wrong column name
  poly  <- sf::st_polygon(list(matrix(c(0,0,1,0,1,1,0,1,0,0), ncol=2, byrow=TRUE)))
  bad.sf <- sf::st_sf(barcode = "x", geometry = sf::st_sfc(poly))
  expect_error(
    add_segmentation(niches.obj, bad.sf),
    regexp = "'cell' column"
  )
})


# ------------------------------------------------------------------------------
# add_segmentation: partial coverage warning
# ------------------------------------------------------------------------------

test_that("add_segmentation: errors when some cells have no polygon", {
  skip_if_not_installed("sf")
  .setup_fixtures()
  # Only provide polygons for the first half of cells — this is a hard error
  half.cells <- niches.obj$cell.list[seq_len(floor(length(niches.obj$cell.list) / 2))]
  seg.sf     <- .make_seg_sf(half.cells)
  expect_error(
    add_segmentation(niches.obj, seg.sf),
    regexp = "no polygon"
  )
})

test_that("add_segmentation: warns when seg.sf contains extra cells not in object", {
  skip_if_not_installed("sf")
  .setup_fixtures()
  extra.cells <- c(niches.obj$cell.list, "PHANTOM_CELL_1", "PHANTOM_CELL_2")
  seg.sf      <- .make_seg_sf(extra.cells)
  expect_warning(
    add_segmentation(niches.obj, seg.sf),
    regexp = "not in the NICHESObject"
  )
})


# ------------------------------------------------------------------------------
# Subsetting: $segmentation is filtered correctly by [.NICHESObject
# ------------------------------------------------------------------------------

test_that("subset [: $segmentation is filtered to retained cells only", {
  skip_if_not_installed("sf")
  .setup_fixtures()
  seg.sf  <- .make_seg_sf(niches.obj$cell.list)
  obj     <- add_segmentation(niches.obj, seg.sf)

  keep    <- niches.obj$cell.list[1:10]
  sub.obj <- obj[keep]

  expect_equal(nrow(sub.obj$segmentation), 10L)
  expect_true(all(sub.obj$segmentation$cell %in% keep))
})

test_that("subset [: $segmentation remains NULL when no segmentation attached", {
  .setup_fixtures()
  keep    <- niches.obj$cell.list[1:10]
  sub.obj <- niches.obj[keep]
  expect_null(sub.obj$segmentation)
})
