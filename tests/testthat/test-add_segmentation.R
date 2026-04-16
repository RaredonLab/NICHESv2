# ==============================================================================
# FILE:    tests/testthat/test-add_segmentation.R
# PURPOSE: Unit tests for add_segmentation().
#          Uses data(NICHESv2_example) for the full NICHESObject — loaded
#          inline inside each block to preserve block independence.
#          All blocks that touch sf objects skip when sf is not installed.
# ==============================================================================

# File-level helper: builds a minimal sf object (unit-square polygons) for a
# given character vector of cell names. Geometry values are irrelevant for
# these tests — only slot population and filtering logic is being tested.
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

# ── Basic attachment ───────────────────────────────────────────────────────────

test_that("add_segmentation: $segmentation is NULL before attachment", {
  data(NICHESv2_example)
  expect_null(NICHESv2_example$segmentation)
})

test_that("add_segmentation: returns a NICHESObject", {
  skip_if_not_installed("sf")
  data(NICHESv2_example)
  seg.sf <- .make_seg_sf(NICHESv2_example$cell.list)
  result <- add_segmentation(NICHESv2_example, seg.sf)
  expect_true(is_NICHESObject(result))
})

test_that("add_segmentation: $segmentation is an sf object after attachment", {
  skip_if_not_installed("sf")
  data(NICHESv2_example)
  seg.sf <- .make_seg_sf(NICHESv2_example$cell.list)
  result <- add_segmentation(NICHESv2_example, seg.sf)
  expect_true(inherits(result$segmentation, "sf"))
})

test_that("add_segmentation: $segmentation has one row per cell", {
  skip_if_not_installed("sf")
  data(NICHESv2_example)
  seg.sf <- .make_seg_sf(NICHESv2_example$cell.list)
  result <- add_segmentation(NICHESv2_example, seg.sf)
  expect_equal(nrow(result$segmentation), length(NICHESv2_example$cell.list))
})

# ── Input validation ───────────────────────────────────────────────────────────

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
  data(NICHESv2_example)
  expect_error(
    add_segmentation(NICHESv2_example,
                     data.frame(cell = NICHESv2_example$cell.list)),
    regexp = "sf object"
  )
})

test_that("add_segmentation: errors if seg.sf has no 'cell' column", {
  skip_if_not_installed("sf")
  data(NICHESv2_example)
  poly   <- sf::st_polygon(list(matrix(c(0,0,1,0,1,1,0,1,0,0),
                                       ncol = 2, byrow = TRUE)))
  bad.sf <- sf::st_sf(barcode = "x", geometry = sf::st_sfc(poly))
  expect_error(
    add_segmentation(NICHESv2_example, bad.sf),
    regexp = "'cell' column"
  )
})

# ── Coverage validation ────────────────────────────────────────────────────────

test_that("add_segmentation: errors when some cells have no polygon", {
  skip_if_not_installed("sf")
  data(NICHESv2_example)
  half.cells <- NICHESv2_example$cell.list[
    seq_len(floor(length(NICHESv2_example$cell.list) / 2))
  ]
  seg.sf <- .make_seg_sf(half.cells)
  expect_error(
    add_segmentation(NICHESv2_example, seg.sf),
    regexp = "no polygon"
  )
})

test_that("add_segmentation: warns when seg.sf contains extra cells not in object", {
  skip_if_not_installed("sf")
  data(NICHESv2_example)
  extra.cells <- c(NICHESv2_example$cell.list, "PHANTOM_CELL_1", "PHANTOM_CELL_2")
  seg.sf      <- .make_seg_sf(extra.cells)
  expect_warning(
    add_segmentation(NICHESv2_example, seg.sf),
    regexp = "not in the NICHESObject"
  )
})

# ── Subsetting with segmentation ───────────────────────────────────────────────

test_that("add_segmentation: subset [ filters $segmentation to retained cells", {
  skip_if_not_installed("sf")
  data(NICHESv2_example)
  seg.sf  <- .make_seg_sf(NICHESv2_example$cell.list)
  obj     <- add_segmentation(NICHESv2_example, seg.sf)
  keep    <- NICHESv2_example$cell.list[1:10]
  sub.obj <- obj[keep]
  
  expect_equal(nrow(sub.obj$segmentation), 10L)
  expect_true(all(sub.obj$segmentation$cell %in% keep))
})

test_that("add_segmentation: subset [ leaves $segmentation NULL when none attached", {
  data(NICHESv2_example)
  keep    <- NICHESv2_example$cell.list[1:10]
  sub.obj <- NICHESv2_example[keep]
  expect_null(sub.obj$segmentation)
})