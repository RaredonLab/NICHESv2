# ==============================================================================
# FILE:    tests/testthat/test-aggregate_NICHESObject.R
# PURPOSE: Unit tests for aggregate_NICHESObject().
#          Fixture: data(NICHESv2_example) -- neighborhood slots always
#          populated, no aggregation slots pre-populated.
# ==============================================================================

data(NICHESv2_example)


# -- Input validation ----------------------------------------------------------

test_that("aggregate_NICHESObject: errors if niches.obj is not a NICHESObject", {
  expect_error(
    aggregate_NICHESObject(list()),
    regexp = "NICHESObject"
  )
})


# -- cell.type.col = NULL: edge aggregation only -------------------------------

test_that("aggregate_NICHESObject: cell.type.col = NULL populates only neighborhood.edge.agg", {
  obj <- aggregate_NICHESObject(NICHESv2_example, cell.type.col = NULL,
                                verbose = FALSE)
  expect_true(inherits(obj, "NICHESObject"))
  expect_false(is.null(obj$aggregations$neighborhood.edge.agg))
  expect_true(is.null(obj$aggregations$neighborhood.composition))
})


# -- Full pipeline -------------------------------------------------------------

test_that("aggregate_NICHESObject: cell.type.col provided populates both slots", {
  obj <- aggregate_NICHESObject(NICHESv2_example,
                                cell.type.col = "celltype",
                                verbose       = FALSE)
  expect_false(is.null(obj$aggregations$neighborhood.edge.agg))
  expect_false(is.null(obj$aggregations$neighborhood.composition))
})


# -- Idempotency ---------------------------------------------------------------

test_that("aggregate_NICHESObject: calling twice replaces both slots and does not double rows", {
  obj1 <- aggregate_NICHESObject(NICHESv2_example,
                                 cell.type.col = "celltype",
                                 verbose       = FALSE)
  obj2 <- aggregate_NICHESObject(obj1,
                                 cell.type.col = "celltype",
                                 verbose       = FALSE)
  expect_equal(
    nrow(obj2$aggregations$neighborhood.edge.agg),
    nrow(obj1$aggregations$neighborhood.edge.agg)
  )
  expect_equal(
    nrow(obj2$aggregations$neighborhood.composition),
    nrow(obj1$aggregations$neighborhood.composition)
  )
})
