# ==============================================================================
# FILE:    tests/testthat/test-add_cell_meta.R
# PURPOSE: Unit tests for add_cell_meta().
#          Uses data(NICHESv2_example) throughout.
#          Full replacement semantics: $cell.meta becomes exactly new.meta;
#          columns absent from new.meta are dropped.
# ==============================================================================

data(NICHESv2_example)


# -- Input validation ----------------------------------------------------------

test_that("add_cell_meta: errors if niches.obj is not a NICHESObject", {
  expect_error(
    add_cell_meta(list(), data.frame(cell = "A")),
    regexp = "NICHESObject"
  )
})

test_that("add_cell_meta: errors if new.meta is not a data.frame or data.table", {
  expect_error(
    add_cell_meta(NICHESv2_example, new.meta = list(cell = "A")),
    regexp = "data.frame or data.table"
  )
})

test_that("add_cell_meta: errors if no cell column and no real rownames", {
  bad.meta <- data.frame(barcode = NICHESv2_example$cell.list, val = 1L)
  expect_error(
    add_cell_meta(NICHESv2_example, bad.meta, verbose = FALSE),
    regexp = "'cell' column"
  )
})

test_that("add_cell_meta: errors if cells are missing from new.meta", {
  # Supply only the first 10 cells -- should hard-error
  partial.meta <- data.frame(
    cell = NICHESv2_example$cell.list[1:10],
    x    = NICHESv2_example$cell.meta$x[1:10],
    y    = NICHESv2_example$cell.meta$y[1:10]
  )
  expect_error(
    add_cell_meta(NICHESv2_example, partial.meta, verbose = FALSE),
    regexp = "missing"
  )
})


# -- Rowname auto-detection ----------------------------------------------------

test_that("add_cell_meta: accepts rowname-style input (no cell column)", {
  rn.meta <- as.data.frame(NICHESv2_example$cell.meta)
  rownames(rn.meta) <- rn.meta$cell
  rn.meta$cell <- NULL
  result <- add_cell_meta(NICHESv2_example, rn.meta, verbose = FALSE)
  expect_true(inherits(result, "NICHESObject"))
  expect_true("cell" %in% colnames(result$cell.meta))
})


# -- Extra cells ---------------------------------------------------------------

test_that("add_cell_meta: emits message and drops extra cells", {
  base.meta     <- as.data.frame(NICHESv2_example$cell.meta)
  fake.row      <- base.meta[1L, ]
  fake.row$cell <- "FAKE_CELL"
  extra.meta    <- rbind(base.meta, fake.row)
  expect_message(
    add_cell_meta(NICHESv2_example, extra.meta, verbose = TRUE),
    regexp = "will be dropped"
  )
  result <- add_cell_meta(NICHESv2_example, extra.meta, verbose = FALSE)
  expect_equal(nrow(result$cell.meta), length(NICHESv2_example$cell.list))
})


# -- Full replacement semantics ------------------------------------------------

test_that("add_cell_meta: columns absent from new.meta are dropped", {
  # Supply only cell + x + y; celltype should disappear
  stripped.meta <- data.frame(
    cell = NICHESv2_example$cell.list,
    x    = NICHESv2_example$cell.meta$x,
    y    = NICHESv2_example$cell.meta$y
  )
  result <- add_cell_meta(NICHESv2_example, stripped.meta, verbose = FALSE)
  expect_false("celltype" %in% colnames(result$cell.meta))
})

test_that("add_cell_meta: new columns in new.meta are present after replacement", {
  new.meta <- as.data.frame(NICHESv2_example$cell.meta)
  new.meta$region <- sample(c("core", "edge"), nrow(new.meta), replace = TRUE)
  result <- add_cell_meta(NICHESv2_example, new.meta, verbose = FALSE)
  expect_true("region" %in% colnames(result$cell.meta))
})

test_that("add_cell_meta: existing column values are overwritten by new.meta", {
  new.meta <- as.data.frame(NICHESv2_example$cell.meta)
  new.meta$celltype <- "overridden"
  result <- add_cell_meta(NICHESv2_example, new.meta, verbose = FALSE)
  expect_true(all(result$cell.meta$celltype == "overridden"))
})


# -- Return type and slot integrity --------------------------------------------

test_that("add_cell_meta: returns a NICHESObject", {
  new.meta <- as.data.frame(NICHESv2_example$cell.meta)
  result   <- add_cell_meta(NICHESv2_example, new.meta, verbose = FALSE)
  expect_true(inherits(result, "NICHESObject"))
})

test_that("add_cell_meta: $cell.meta is keyed on cell", {
  new.meta <- as.data.frame(NICHESv2_example$cell.meta)
  result   <- add_cell_meta(NICHESv2_example, new.meta, verbose = FALSE)
  expect_equal(data.table::key(result$cell.meta), "cell")
})

test_that("add_cell_meta: $cell.meta row count equals length of $cell.list", {
  new.meta <- as.data.frame(NICHESv2_example$cell.meta)
  result   <- add_cell_meta(NICHESv2_example, new.meta, verbose = FALSE)
  expect_equal(nrow(result$cell.meta), length(result$cell.list))
})


# -- $edge.meta rebuild --------------------------------------------------------

test_that("add_cell_meta: new categorical column propagates to $edge.meta as Joint.*", {
  new.meta         <- as.data.frame(NICHESv2_example$cell.meta)
  new.meta$cluster <- sample(paste0("C", 1:3), nrow(new.meta), replace = TRUE)
  result           <- add_cell_meta(NICHESv2_example, new.meta, verbose = FALSE)
  expect_true("Joint.cluster" %in% colnames(result$edge.meta))
})

test_that("add_cell_meta: rebuild.edge.meta = FALSE skips $edge.meta update", {
  new.meta         <- as.data.frame(NICHESv2_example$cell.meta)
  new.meta$cluster <- "X"
  result <- add_cell_meta(NICHESv2_example, new.meta,
                           rebuild.edge.meta = FALSE, verbose = FALSE)
  expect_false("Joint.cluster" %in% colnames(result$edge.meta))
})


# -- Staleness notification ----------------------------------------------------

test_that("add_cell_meta: emits staleness message when neighborhood.composition is populated", {
  obj.with.comp <- add_neighborhood_composition(
    NICHESv2_example,
    cell.type.col = "celltype",
    verbose       = FALSE
  )
  new.meta <- as.data.frame(NICHESv2_example$cell.meta)
  expect_message(
    add_cell_meta(obj.with.comp, new.meta, verbose = TRUE),
    regexp = "stale"
  )
})

test_that("add_cell_meta: no staleness message when neighborhood.composition is NULL", {
  new.meta <- as.data.frame(NICHESv2_example$cell.meta)
  expect_no_message(
    add_cell_meta(NICHESv2_example, new.meta, verbose = TRUE)
  )
})
