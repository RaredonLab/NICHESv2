# ==============================================================================
# FILE:    tests/testthat/test-add_neighborhood_composition.R
# PURPOSE: Unit tests for add_neighborhood_composition().
#          Fixture: data(NICHESv2_example) -- spatial mode, cell types are
#          "TypeA", "TypeB", "TypeC" in $cell.meta$celltype.
#          gsub("[^A-Za-z0-9.]", ".", type) leaves these names unchanged.
# ==============================================================================

data(NICHESv2_example)

# Run once at file level with defaults (include.self = TRUE, normalize = TRUE).
comp.obj <- add_neighborhood_composition(
  NICHESv2_example,
  cell.type.col = "celltype",
  verbose       = FALSE
)
comp.dt <- comp.obj$aggregations$neighborhood.composition


# -- Input validation ----------------------------------------------------------

test_that("add_neighborhood_composition: errors if niches.obj is not a NICHESObject", {
  expect_error(
    add_neighborhood_composition(list()),
    regexp = "NICHESObject"
  )
})

test_that("add_neighborhood_composition: errors if mode is nonspatial", {
  obj.ns <- NICHESv2_example
  attr(obj.ns, "mode") <- "nonspatial"
  expect_error(
    add_neighborhood_composition(obj.ns, cell.type.col = "celltype",
                                 verbose = FALSE),
    regexp = "nonspatial"
  )
})

test_that("add_neighborhood_composition: errors if cell.type.col is NULL", {
  expect_error(
    add_neighborhood_composition(NICHESv2_example, cell.type.col = NULL,
                                 verbose = FALSE),
    regexp = "must be provided"
  )
})

test_that("add_neighborhood_composition: errors if cell.type.col is not in cell.meta", {
  expect_error(
    add_neighborhood_composition(NICHESv2_example, cell.type.col = "nosuchcol",
                                 verbose = FALSE),
    regexp = "not found"
  )
})


# -- Smoke test ----------------------------------------------------------------

test_that("add_neighborhood_composition: returns a NICHESObject with slot populated", {
  expect_true(inherits(comp.obj, "NICHESObject"))
  expect_false(is.null(comp.dt))
})


# -- Slot type and key ---------------------------------------------------------

test_that("add_neighborhood_composition: neighborhood.composition is a data.table", {
  expect_true(data.table::is.data.table(comp.dt))
})

test_that("add_neighborhood_composition: keyed on cell", {
  expect_equal(data.table::key(comp.dt), "cell")
})


# -- Column structure ----------------------------------------------------------

test_that("add_neighborhood_composition: expected count and prop columns present for all three types", {
  expected.count <- c("TypeA.neighbor.count", "TypeB.neighbor.count",
                      "TypeC.neighbor.count")
  expected.prop  <- c("TypeA.neighbor.prop",  "TypeB.neighbor.prop",
                      "TypeC.neighbor.prop")
  expect_true(all(expected.count %in% colnames(comp.dt)))
  expect_true(all(expected.prop  %in% colnames(comp.dt)))
})


# -- Cell coverage -------------------------------------------------------------

test_that("add_neighborhood_composition: all 80 cells present and zero-filled", {
  expect_equal(nrow(comp.dt), 80L)
  expect_true(all(NICHESv2_example$cell.list %in% comp.dt$cell))
})


# -- Proportion validity -------------------------------------------------------

test_that("add_neighborhood_composition: prop columns sum to 1 per cell with neighbors, 0 for isolated cells", {
  count.cols <- c("TypeA.neighbor.count", "TypeB.neighbor.count",
                  "TypeC.neighbor.count")
  prop.cols  <- c("TypeA.neighbor.prop",  "TypeB.neighbor.prop",
                  "TypeC.neighbor.prop")
  row.totals <- rowSums(comp.dt[, count.cols, with = FALSE])
  row.props  <- rowSums(comp.dt[, prop.cols,  with = FALSE])
  has.neighbor <- row.totals > 0
  expect_true(all(abs(row.props[has.neighbor] - 1) < 1e-10))
  expect_true(all(row.props[!has.neighbor] == 0))
})


# -- Count/prop arithmetic consistency -----------------------------------------

test_that("add_neighborhood_composition: prop equals count / row total per cell", {
  types      <- c("TypeA", "TypeB", "TypeC")
  count.cols <- paste0(types, ".neighbor.count")
  prop.cols  <- paste0(types, ".neighbor.prop")
  row.totals <- rowSums(comp.dt[, count.cols, with = FALSE])
  has.neighbor <- row.totals > 0
  for (i in seq_along(types)) {
    expected.prop <- ifelse(
      has.neighbor,
      comp.dt[[count.cols[i]]] / row.totals,
      0
    )
    expect_equal(comp.dt[[prop.cols[i]]], expected.prop,
                 label = paste(prop.cols[i], "arithmetic"))
  }
})


# -- Idempotency ---------------------------------------------------------------

test_that("add_neighborhood_composition: calling twice replaces slot and does not double rows", {
  obj2 <- add_neighborhood_composition(
    comp.obj,
    cell.type.col = "celltype",
    verbose       = FALSE
  )
  expect_equal(
    nrow(obj2$aggregations$neighborhood.composition),
    nrow(comp.dt)
  )
})


# -- include.self direction of effect ------------------------------------------

# Optional: tests direction of effect only, not exact values, because exact
# counts depend on autocrine self-loop structure in the synthetic data.
test_that("add_neighborhood_composition: include.self = FALSE gives row totals <= include.self = TRUE", {
  obj.no.self <- add_neighborhood_composition(
    NICHESv2_example,
    cell.type.col = "celltype",
    include.self  = FALSE,
    verbose       = FALSE
  )
  dt.no.self <- obj.no.self$aggregations$neighborhood.composition
  count.cols <- c("TypeA.neighbor.count", "TypeB.neighbor.count",
                  "TypeC.neighbor.count")
  totals.with    <- rowSums(comp.dt[,    count.cols, with = FALSE])
  totals.without <- rowSums(dt.no.self[, count.cols, with = FALSE])
  expect_true(all(totals.with >= totals.without))
})
