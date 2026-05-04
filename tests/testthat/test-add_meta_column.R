# ==============================================================================
# FILE:    tests/testthat/test-add_meta_column.R
# PURPOSE: Unit tests for add_meta_column().
#          Uses data(NICHESv2_example) throughout.
#          Single-column append/overwrite; all other $cell.meta columns
#          are preserved unchanged.
# ==============================================================================

data(NICHESv2_example)

# Helper: build a full named vector covering all 80 cells
.make_vec <- function(values, nm = NICHESv2_example$cell.list) {
  setNames(values, nm)
}


# -- Input validation ----------------------------------------------------------

test_that("add_meta_column: errors if niches.obj is not a NICHESObject", {
  expect_error(
    add_meta_column(list(), "col", setNames("x", "cell1")),
    regexp = "NICHESObject"
  )
})

test_that("add_meta_column: errors if col.name is not a single character string", {
  vec <- .make_vec(rep("A", 80L))
  expect_error(
    add_meta_column(NICHESv2_example, col.name = 123L, col.values = vec),
    regexp = "col.name"
  )
  expect_error(
    add_meta_column(NICHESv2_example, col.name = c("a", "b"), col.values = vec),
    regexp = "col.name"
  )
  expect_error(
    add_meta_column(NICHESv2_example, col.name = NA_character_, col.values = vec),
    regexp = "col.name"
  )
})

test_that("add_meta_column: errors if col.values is not a named vector", {
  expect_error(
    add_meta_column(NICHESv2_example, "region", col.values = rep("A", 80L)),
    regexp = "named atomic vector"
  )
  expect_error(
    add_meta_column(NICHESv2_example, "region", col.values = list(a = 1)),
    regexp = "named atomic vector"
  )
})

test_that("add_meta_column: errors if cells are missing from col.values", {
  partial.vec <- .make_vec(rep("A", 10L), NICHESv2_example$cell.list[1:10])
  expect_error(
    add_meta_column(NICHESv2_example, "region", partial.vec, verbose = FALSE),
    regexp = "missing"
  )
})


# -- Overwrite notification ----------------------------------------------------

test_that("add_meta_column: emits message when overwriting existing column", {
  vec <- .make_vec(rep("TypeA", 80L))
  expect_message(
    add_meta_column(NICHESv2_example, "celltype", vec, verbose = TRUE),
    regexp = "overwritten"
  )
})

test_that("add_meta_column: no message when adding a new column (verbose = TRUE)", {
  vec <- .make_vec(rep("core", 80L))
  expect_no_message(
    add_meta_column(NICHESv2_example, "region", vec, verbose = TRUE)
  )
})


# -- Column value correctness --------------------------------------------------

test_that("add_meta_column: new column values are correct and aligned to $cell.list", {
  vals <- setNames(seq_len(80L), NICHESv2_example$cell.list)
  result <- add_meta_column(NICHESv2_example, "rank", vals, verbose = FALSE)
  # Values should be retrievable in $cell.list order
  expected <- vals[result$cell.meta$cell]
  expect_equal(result$cell.meta$rank, unname(expected))
})

test_that("add_meta_column: overwrite replaces column values correctly", {
  vec    <- .make_vec(rep("overridden", 80L))
  result <- add_meta_column(NICHESv2_example, "celltype", vec, verbose = FALSE)
  expect_true(all(result$cell.meta$celltype == "overridden"))
})

test_that("add_meta_column: extra names in col.values are silently ignored", {
  extra.vec <- c(.make_vec(rep("A", 80L)), setNames("B", "FAKE_CELL"))
  result <- add_meta_column(NICHESv2_example, "region", extra.vec,
                             verbose = FALSE)
  expect_equal(nrow(result$cell.meta), 80L)
})


# -- Preservation of other columns ---------------------------------------------

test_that("add_meta_column: columns not in col.name are preserved", {
  vec    <- .make_vec(rep("core", 80L))
  result <- add_meta_column(NICHESv2_example, "region", vec, verbose = FALSE)
  expect_true("celltype" %in% colnames(result$cell.meta))
  expect_true("x"        %in% colnames(result$cell.meta))
  expect_true("y"        %in% colnames(result$cell.meta))
})


# -- Return type ---------------------------------------------------------------

test_that("add_meta_column: returns a NICHESObject", {
  vec    <- .make_vec(rep("A", 80L))
  result <- add_meta_column(NICHESv2_example, "region", vec, verbose = FALSE)
  expect_true(inherits(result, "NICHESObject"))
})


# -- $edge.meta rebuild --------------------------------------------------------

test_that("add_meta_column: new categorical column propagates to $edge.meta", {
  vec    <- .make_vec(sample(c("core", "edge"), 80L, replace = TRUE))
  result <- add_meta_column(NICHESv2_example, "region", vec, verbose = FALSE)
  expect_true("Joint.region" %in% colnames(result$edge.meta))
})

test_that("add_meta_column: rebuild.edge.meta = FALSE skips $edge.meta update", {
  vec    <- .make_vec(rep("core", 80L))
  result <- add_meta_column(NICHESv2_example, "region", vec,
                             rebuild.edge.meta = FALSE, verbose = FALSE)
  expect_false("Joint.region" %in% colnames(result$edge.meta))
})


# -- Staleness notification ----------------------------------------------------

test_that("add_meta_column: emits staleness message when neighborhood.composition is populated", {
  obj.with.comp <- add_neighborhood_composition(
    NICHESv2_example,
    cell.type.col = "celltype",
    verbose       = FALSE
  )
  vec <- .make_vec(rep("TypeA", 80L))
  expect_message(
    add_meta_column(obj.with.comp, "celltype", vec, verbose = TRUE),
    regexp = "stale"
  )
})

test_that("add_meta_column: no staleness message when neighborhood.composition is NULL", {
  vec <- .make_vec(rep("core", 80L))
  expect_no_message(
    add_meta_column(NICHESv2_example, "region", vec, verbose = TRUE)
  )
})
