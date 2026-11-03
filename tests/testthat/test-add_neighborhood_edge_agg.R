# ==============================================================================
# FILE:    tests/testthat/test-add_neighborhood_edge_agg.R
# PURPOSE: Unit tests for add_neighborhood_edge_agg().
#          Fixture: data(NICHESv2_example) -- neighborhood slots always
#          populated, no aggregation slots pre-populated.
# ==============================================================================

data(NICHESv2_example)

# Run once at file level; reused across structural and arithmetic tests.
# Copy-on-modify semantics in R mean test_that blocks cannot corrupt this.
agg.obj <- add_neighborhood_edge_agg(NICHESv2_example, verbose = FALSE)
agg.dt  <- agg.obj$aggregations$neighborhood.edge.agg


# -- Input validation ----------------------------------------------------------

test_that("add_neighborhood_edge_agg: errors if niches.obj is not a NICHESObject", {
  expect_error(
    add_neighborhood_edge_agg(list()),
    regexp = "NICHESObject"
  )
})

test_that("add_neighborhood_edge_agg: errors if neighborhood slots are NULL", {
  obj.no.nbr <- NICHESv2_example
  obj.no.nbr$neighborhood.edge.list <- NULL
  expect_error(
    add_neighborhood_edge_agg(obj.no.nbr, verbose = FALSE),
    regexp = "Neighborhood slots are NULL"
  )
})


# -- Smoke test ----------------------------------------------------------------

test_that("add_neighborhood_edge_agg: returns a NICHESObject with slot populated", {
  expect_true(inherits(agg.obj, "NICHESObject"))
  expect_false(is.null(agg.dt))
})


# -- Slot type and key ---------------------------------------------------------

test_that("add_neighborhood_edge_agg: neighborhood.edge.agg is a data.table", {
  expect_true(data.table::is.data.table(agg.dt))
})

test_that("add_neighborhood_edge_agg: keyed on (cell, LRM)", {
  expect_equal(data.table::key(agg.dt), c("cell", "LRM"))
})


# -- Column structure ----------------------------------------------------------

test_that("add_neighborhood_edge_agg: all 27 expected columns are present", {
  expected.cols <- c(
    "cell", "LRM",
    "count.in", "count.out", "count.cross",
    "count.self.autocrine", "count.other.autocrine", "count.total",
    "score.in.sum",               "score.out.sum",               "score.cross.sum",
    "score.self.autocrine.sum",   "score.other.autocrine.sum",   "score.total.sum",
    "score.in.mean",              "score.out.mean",               "score.cross.mean",
    "score.self.autocrine.mean",  "score.other.autocrine.mean",  "score.total.mean",
    "score.in.cmean",             "score.out.cmean",              "score.cross.cmean",
    "score.self.autocrine.cmean", "score.other.autocrine.cmean", "score.total.cmean",
    "score.inout.ratio"
  )
  expect_true(all(expected.cols %in% colnames(agg.dt)))
  expect_equal(ncol(agg.dt), 27L)
})


# -- Count column sanity -------------------------------------------------------

test_that("add_neighborhood_edge_agg: count columns are non-negative integers", {
  count.cols <- c("count.in", "count.out", "count.cross",
                  "count.self.autocrine", "count.other.autocrine", "count.total")
  for (col in count.cols) {
    expect_true(is.integer(agg.dt[[col]]),
                label = paste(col, "is integer"))
    expect_true(all(agg.dt[[col]] >= 0L),
                label = paste(col, "is non-negative"))
  }
})

test_that("add_neighborhood_edge_agg: count.total equals sum of five category counts", {
  expect_true(all(
    agg.dt$count.total ==
      agg.dt$count.in + agg.dt$count.out + agg.dt$count.cross +
      agg.dt$count.self.autocrine + agg.dt$count.other.autocrine
  ))
})

test_that("add_neighborhood_edge_agg: count.self.autocrine is always 0 or 1", {
  expect_true(all(agg.dt$count.self.autocrine %in% c(0L, 1L)))
})


# -- score.inout.ratio correctness ---------------------------------------------

test_that("add_neighborhood_edge_agg: score.inout.ratio is NA iff score.out.sum == 0", {
  expect_equal(
    is.na(agg.dt$score.inout.ratio),
    agg.dt$score.out.sum == 0
  )
})

test_that("add_neighborhood_edge_agg: score.inout.ratio equals score.in.sum / score.out.sum where defined", {
  rows.defined <- !is.na(agg.dt$score.inout.ratio)
  expect_equal(
    agg.dt$score.inout.ratio[rows.defined],
    agg.dt$score.in.sum[rows.defined] / agg.dt$score.out.sum[rows.defined]
  )
})


# -- score.cat.cmean NA rule ---------------------------------------------------

test_that("add_neighborhood_edge_agg: score.cat.cmean is NA iff count.cat == 0", {
  cats <- c("in", "out", "cross", "self.autocrine", "other.autocrine", "total")
  for (cat in cats) {
    cmean.col <- paste0("score.", cat, ".cmean")
    count.col <- paste0("count.", cat)
    expect_equal(
      is.na(agg.dt[[cmean.col]]),
      agg.dt[[count.col]] == 0L,
      label = paste(cmean.col, "NA rule")
    )
  }
})


# -- Cell coverage -------------------------------------------------------------

# Slot is sparse: only (cell, LRM) pairs with at least one non-zero-score
# neighborhood edge are stored. Not all 80 cells are guaranteed to appear.
# Test that every cell that does appear is a valid cell from neighborhood.meta.
test_that("add_neighborhood_edge_agg: all output cells are valid (subset of neighborhood.meta)", {
  valid.cells <- NICHESv2_example$neighborhood.meta$cell
  expect_true(all(agg.dt$cell %in% valid.cells))
})


# -- Idempotency ---------------------------------------------------------------

test_that("add_neighborhood_edge_agg: calling twice replaces slot and does not double rows", {
  obj2 <- add_neighborhood_edge_agg(agg.obj, verbose = FALSE)
  expect_equal(
    nrow(obj2$aggregations$neighborhood.edge.agg),
    nrow(agg.dt)
  )
})


# -- Nonspatial mode -----------------------------------------------------------

# add_neighborhood_edge_agg has no mode guard; nonspatial is supported.
# Interpretation shifts (co-occurring cell pairs, not spatial microenvironment)
# but the computation runs without error.
test_that("add_neighborhood_edge_agg: runs without error in nonspatial mode", {
  obj.ns <- NICHESv2_example
  attr(obj.ns, "mode") <- "nonspatial"
  result <- add_neighborhood_edge_agg(obj.ns, verbose = FALSE)
  expect_true(data.table::is.data.table(
    result$aggregations$neighborhood.edge.agg
  ))
})
