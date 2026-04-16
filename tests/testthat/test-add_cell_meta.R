# ==============================================================================
# FILE:    tests/testthat/test-add_cell_meta.R
# PURPOSE: Unit tests for add_cell_meta().
#          Uses data(NICHESv2_example) for the full NICHESObject — loaded
#          inline inside each block to preserve block independence.
# ==============================================================================

test_that("add_cell_meta: adds a new column to $cell.meta", {
  data(NICHESv2_example)
  set.seed(1L)
  new.col <- data.frame(
    cell    = NICHESv2_example$cell.list,
    cluster = sample(1:5, length(NICHESv2_example$cell.list), replace = TRUE),
    stringsAsFactors = FALSE
  )
  updated <- add_cell_meta(NICHESv2_example, new.col)
  expect_true("cluster" %in% colnames(updated$cell.meta))
})

test_that("add_cell_meta: overwrites an existing column", {
  data(NICHESv2_example)
  override <- data.frame(
    cell     = NICHESv2_example$cell.list,
    celltype = "overridden",
    stringsAsFactors = FALSE
  )
  updated <- add_cell_meta(NICHESv2_example, override)
  expect_true(all(updated$cell.meta$celltype == "overridden"))
})

test_that("add_cell_meta: rebuilds $edge.meta with new Joint column", {
  data(NICHESv2_example)
  set.seed(2L)
  new.col <- data.frame(
    cell    = NICHESv2_example$cell.list,
    cluster = paste0("C", sample(1:3, length(NICHESv2_example$cell.list),
                                 replace = TRUE)),
    stringsAsFactors = FALSE
  )
  updated <- add_cell_meta(NICHESv2_example, new.col)
  expect_true("Joint.cluster" %in% colnames(updated$edge.meta))
})

test_that("add_cell_meta: returns a NICHESObject", {
  data(NICHESv2_example)
  dummy   <- data.frame(cell = NICHESv2_example$cell.list[1L], extra = "x")
  updated <- add_cell_meta(NICHESv2_example, dummy)
  expect_s3_class(updated, "NICHESObject")
})

test_that("add_cell_meta: errors when 'cell' column missing from new.meta", {
  data(NICHESv2_example)
  bad.meta <- data.frame(barcode = NICHESv2_example$cell.list, val = 1L)
  expect_error(
    add_cell_meta(NICHESv2_example, bad.meta),
    regexp = "'cell' column"
  )
})

test_that("add_cell_meta: errors on non-NICHESObject input", {
  expect_error(
    add_cell_meta(list(), data.frame(cell = "A")),
    regexp = "NICHESObject"
  )
})