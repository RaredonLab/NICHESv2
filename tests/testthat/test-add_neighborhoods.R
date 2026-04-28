# ==============================================================================
# FILE:    tests/testthat/test-add_neighborhoods.R
# PURPOSE: Unit tests for add_neighborhoods().
#          Builds a minimal NICHESObject inline — no shared fixtures needed.
# ==============================================================================

# File-level helper: builds the smallest valid NICHESObject that
# add_neighborhoods() can operate on (only reads $edge.list, writes 3 slots).
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
    neighborhood.meta      = NULL,
    aggregations           = list()
  )
  class(obj) <- "NICHESObject"
  attr(obj, "mode") <- "spatial"
  obj
}

# ── Input validation ───────────────────────────────────────────────────────────

test_that("add_neighborhoods: errors on non-NICHESObject input", {
  expect_error(add_neighborhoods(list()), "must be a NICHESObject")
})

test_that("add_neighborhoods: errors on empty edge.list", {
  obj <- make_minimal_niches_obj()
  obj$edge.list <- obj$edge.list[0L]
  expect_error(add_neighborhoods(obj), "edge.list is empty")
})

# ── Slot population ────────────────────────────────────────────────────────────

test_that("add_neighborhoods: populates the three neighborhood slots", {
  obj <- make_minimal_niches_obj()
  expect_null(obj$neighborhood.cell.list)
  expect_null(obj$neighborhood.edge.list)
  expect_null(obj$neighborhood.meta)
  
  obj2 <- add_neighborhoods(obj, n.cores = 1L, verbose = FALSE)
  
  expect_type(obj2$neighborhood.cell.list, "list")
  expect_type(obj2$neighborhood.edge.list, "list")
  expect_s3_class(obj2$neighborhood.meta, "data.table")
  expect_identical(data.table::key(obj2$neighborhood.meta), "cell")
  
  # Non-neighborhood slots must be untouched
  expect_identical(obj$cell.list,  obj2$cell.list)
  expect_identical(obj$edge.list,  obj2$edge.list)
})

test_that("add_neighborhoods: returned object is still a NICHESObject", {
  obj  <- make_minimal_niches_obj()
  obj2 <- add_neighborhoods(obj, n.cores = 1L, verbose = FALSE)
  expect_s3_class(obj2, "NICHESObject")
})

test_that("add_neighborhoods: sets edge.filter.mode attribute on the object", {
  obj  <- make_minimal_niches_obj()
  obj2 <- add_neighborhoods(obj, edge.filter.mode = "central_cell_only",
                            n.cores = 1L, verbose = FALSE)
  expect_identical(attr(obj2, "edge.filter.mode"), "central_cell_only")
})

test_that("add_neighborhoods: all cells appear in neighborhood.meta", {
  obj  <- make_minimal_niches_obj()
  obj2 <- add_neighborhoods(obj, n.cores = 1L, verbose = FALSE)
  expect_true(all(obj$cell.list %in% obj2$neighborhood.meta$cell))
})

# ── Overwrite behaviour ────────────────────────────────────────────────────────

test_that("add_neighborhoods: hard error when slots already populated and overwrite = FALSE", {
  obj  <- make_minimal_niches_obj()
  obj2 <- add_neighborhoods(obj, n.cores = 1L, verbose = FALSE)
  expect_error(
    add_neighborhoods(obj2, n.cores = 1L, verbose = FALSE),
    regexp = "overwrite = TRUE"
  )
})

test_that("add_neighborhoods: overwrite = TRUE recomputes slots and emits message", {
  obj  <- make_minimal_niches_obj()
  obj2 <- add_neighborhoods(obj,  edge.filter.mode = "both_in_neighborhood",
                            n.cores = 1L, verbose = FALSE)
  
  expect_message(
    obj3 <- add_neighborhoods(obj2, edge.filter.mode = "central_cell_only",
                              overwrite = TRUE,
                              n.cores   = 1L,
                              verbose   = TRUE),
    regexp = "Recomputing neighborhood slots"
  )
  expect_identical(attr(obj3, "edge.filter.mode"), "central_cell_only")
  
  # Results must differ between modes on this toy graph:
  # both_in_neighborhood for A → {A|B, A|C, B|C}
  # central_cell_only    for A → {A|B, A|C}
  expect_false(identical(
    sort(obj2$neighborhood.edge.list[["A"]]),
    sort(obj3$neighborhood.edge.list[["A"]])
  ))
})

# ── Parity with direct extract_neighborhoods() call ───────────────────────────

test_that("add_neighborhoods: results match direct extract_neighborhoods() call", {
  obj    <- make_minimal_niches_obj()
  direct <- extract_neighborhoods(obj$edge.list, n.cores = 1L, verbose = FALSE)
  viaobj <- add_neighborhoods(obj,               n.cores = 1L, verbose = FALSE)
  
  canon.list <- function(x) lapply(x[sort(names(x))], sort)
  
  expect_identical(canon.list(direct$neighborhood.cell.list),
                   canon.list(viaobj$neighborhood.cell.list))
  expect_identical(canon.list(direct$neighborhood.edge.list),
                   canon.list(viaobj$neighborhood.edge.list))
  
  d.meta <- direct$neighborhood.meta[order(cell)]
  v.meta <- viaobj$neighborhood.meta[order(cell)]
  expect_equal(d.meta$cell,        v.meta$cell)
  expect_equal(d.meta$n.neighbors, v.meta$n.neighbors)
  expect_equal(d.meta$n.edges,     v.meta$n.edges)
})