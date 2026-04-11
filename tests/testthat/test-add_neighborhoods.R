## Tests for add_neighborhoods()
##
## Assumes the standard test NICHESObject fixture is available via the
## helper-*.R convention (auto-sourced by testthat). The fixture is expected
## to be a spatial-mode NICHESObject with a non-empty $edge.list.

# ---------------------------------------------------------------------------
# Helper: build a minimal NICHESObject just for these tests if no fixture
# ---------------------------------------------------------------------------
# We try to use an existing helper if one is defined; otherwise fall back
# to a trivial hand-built object carrying only the slots add_neighborhoods()
# reads and writes.
make_minimal_niches_obj <- function() {
  el <- data.table::data.table(
    SendingCell   = c("A", "A", "B", "C", "D"),
    ReceivingCell = c("B", "C", "C", "D", "D")
  )
  el[, edge := paste0(SendingCell, "|", ReceivingCell)]
  data.table::setcolorder(el, c("edge", "SendingCell", "ReceivingCell"))
  data.table::setkey(el, edge)

  obj <- list(
    cell.list              = c("A", "B", "C", "D"),
    edge.list              = el,
    neighborhood.cell.list = NULL,
    neighborhood.edge.list = NULL,
    neighborhood.meta      = NULL
  )
  class(obj) <- "NICHESObject"
  attr(obj, "mode") <- "spatial"
  obj
}

# ---------------------------------------------------------------------------
# Input validation
# ---------------------------------------------------------------------------
test_that("errors on non-NICHESObject input", {
  expect_error(add_neighborhoods(list()), "must be a NICHESObject")
})

test_that("errors on empty edge.list", {
  obj <- make_minimal_niches_obj()
  obj$edge.list <- obj$edge.list[0]  # zero rows
  expect_error(add_neighborhoods(obj), "edge.list is empty")
})

# ---------------------------------------------------------------------------
# Slot population
# ---------------------------------------------------------------------------
test_that("populates the three neighborhood slots", {
  obj <- make_minimal_niches_obj()
  expect_null(obj$neighborhood.cell.list)
  expect_null(obj$neighborhood.edge.list)
  expect_null(obj$neighborhood.meta)

  obj2 <- add_neighborhoods(obj, verbose = FALSE)

  expect_type(obj2$neighborhood.cell.list, "list")
  expect_type(obj2$neighborhood.edge.list, "list")
  expect_s3_class(obj2$neighborhood.meta, "data.table")
  expect_identical(data.table::key(obj2$neighborhood.meta), "cell")

  # Other slots left untouched
  expect_identical(obj$cell.list, obj2$cell.list)
  expect_identical(obj$edge.list, obj2$edge.list)
})

test_that("returned object is still a NICHESObject", {
  obj <- make_minimal_niches_obj()
  obj2 <- add_neighborhoods(obj, verbose = FALSE)
  expect_s3_class(obj2, "NICHESObject")
})

test_that("sets edge.filter.mode attribute on the object", {
  obj  <- make_minimal_niches_obj()
  obj2 <- add_neighborhoods(obj, edge.filter.mode = "central_cell_only",
                            verbose = FALSE)
  expect_identical(attr(obj2, "edge.filter.mode"), "central_cell_only")
})

test_that("all cells from cell.list appear in meta", {
  obj  <- make_minimal_niches_obj()
  obj2 <- add_neighborhoods(obj, verbose = FALSE)
  # Every cell that participates in an edge must appear
  expect_true(all(obj$cell.list %in% obj2$neighborhood.meta$cell))
})

# ---------------------------------------------------------------------------
# Overwrite behavior
# ---------------------------------------------------------------------------
test_that("overwrites existing neighborhood slots with a message", {
  obj  <- make_minimal_niches_obj()
  obj2 <- add_neighborhoods(obj, edge.filter.mode = "both_in_neighborhood",
                            verbose = FALSE)

  # Now re-extract with a different mode; must emit an overwrite message
  expect_message(
    obj3 <- add_neighborhoods(obj2, edge.filter.mode = "central_cell_only",
                              verbose = TRUE),
    "Overwriting existing neighborhood slots"
  )
  expect_identical(attr(obj3, "edge.filter.mode"), "central_cell_only")

  # Results should differ between the two modes on this toy graph:
  # both_in_neighborhood for A -> {A|B, A|C, B|C}
  # central_cell_only    for A -> {A|B, A|C}
  expect_false(identical(
    sort(obj2$neighborhood.edge.list[["A"]]),
    sort(obj3$neighborhood.edge.list[["A"]])
  ))
})

# ---------------------------------------------------------------------------
# Parity with direct extract_neighborhoods() call
# ---------------------------------------------------------------------------
test_that("add_neighborhoods() is a pure wrapper over extract_neighborhoods()", {
  obj    <- make_minimal_niches_obj()
  direct <- extract_neighborhoods(obj$edge.list, verbose = FALSE)
  viaobj <- add_neighborhoods(obj, verbose = FALSE)

  # Sort within each cell to get a canonical form for comparison
  canon.list <- function(x) lapply(x[sort(names(x))], sort)
  expect_identical(canon.list(direct$neighborhood.cell.list),
                   canon.list(viaobj$neighborhood.cell.list))
  expect_identical(canon.list(direct$neighborhood.edge.list),
                   canon.list(viaobj$neighborhood.edge.list))

  # meta tables: sort rows by cell then compare
  d.meta <- direct$neighborhood.meta[order(cell)]
  v.meta <- viaobj$neighborhood.meta[order(cell)]
  expect_equal(d.meta$cell,        v.meta$cell)
  expect_equal(d.meta$n.neighbors, v.meta$n.neighbors)
  expect_equal(d.meta$n.edges,     v.meta$n.edges)
})
