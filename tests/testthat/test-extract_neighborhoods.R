# ==============================================================================
# FILE:    tests/testthat/test-extract_neighborhoods.R
# PURPOSE: Unit tests for extract_neighborhoods().
#          Operates on a plain edge.list data.table — no NICHESObject needed.
#          Uses a small hand-built deterministic graph so expected outputs can
#          be computed by hand and verified exactly.
# ==============================================================================

# Top-level helper: reused by every block below.
# Graph (directed edges as stored): A->B, A->C, B->C, C->D, D->D (autocrine)
# Symmetric adjacency derived by extract_neighborhoods():
#   A: {B, C}      B: {A, C}      C: {A, B, D}      D: {C, D}
make_toy_edge_list <- function() {
  dt <- data.table::data.table(
    SendingCell   = c("A", "A", "B", "C", "D"),
    ReceivingCell = c("B", "C", "C", "D", "D")
  )
  dt[, edge := paste0(SendingCell, "|", ReceivingCell)]
  data.table::setcolorder(dt, c("edge", "SendingCell", "ReceivingCell"))
  data.table::setkey(dt, edge)
  dt
}

# ── Input validation ───────────────────────────────────────────────────────────

test_that("extract_neighborhoods: errors on non-data.table input", {
  expect_error(
    extract_neighborhoods(data.frame(edge = "A|B")),
    "must be a data.table"
  )
})

test_that("extract_neighborhoods: errors on missing required columns", {
  bad <- data.table::data.table(from = "A", to = "B")
  expect_error(extract_neighborhoods(bad), "missing required columns")
})

test_that("extract_neighborhoods: errors on empty edge list", {
  empty <- data.table::data.table(
    edge          = character(0),
    SendingCell   = character(0),
    ReceivingCell = character(0)
  )
  expect_error(extract_neighborhoods(empty), "zero rows")
})

test_that("extract_neighborhoods: errors on invalid edge.filter.mode", {
  el <- make_toy_edge_list()
  expect_error(
    extract_neighborhoods(el, edge.filter.mode = "bogus", n.cores = 1L),
    "edge.filter.mode must be one of"
  )
})

# ── Output structure and attributes ───────────────────────────────────────────

test_that("extract_neighborhoods: returns three named slots with correct types", {
  el  <- make_toy_edge_list()
  res <- extract_neighborhoods(el, n.cores = 1L, verbose = FALSE)
  
  expect_named(res, c("neighborhood.cell.list",
                      "neighborhood.edge.list",
                      "neighborhood.meta"))
  expect_type(res$neighborhood.cell.list, "list")
  expect_type(res$neighborhood.edge.list, "list")
  expect_s3_class(res$neighborhood.meta, "data.table")
  expect_named(res$neighborhood.meta, c("cell", "n.neighbors", "n.edges"))
  expect_identical(data.table::key(res$neighborhood.meta), "cell")
})

test_that("extract_neighborhoods: attributes edge.filter.mode and remove.autocrine are set", {
  el  <- make_toy_edge_list()
  res <- extract_neighborhoods(el, n.cores = 1L, verbose = FALSE)
  
  expect_identical(attr(res, "edge.filter.mode"), "both_in_neighborhood")
  expect_false(attr(res, "remove.autocrine"))
})

test_that("extract_neighborhoods: all cells appear as keys in all three outputs", {
  el        <- make_toy_edge_list()
  res       <- extract_neighborhoods(el, n.cores = 1L, verbose = FALSE)
  all.cells <- sort(unique(c(el$SendingCell, el$ReceivingCell)))
  
  expect_setequal(names(res$neighborhood.cell.list), all.cells)
  expect_setequal(names(res$neighborhood.edge.list), all.cells)
  expect_setequal(res$neighborhood.meta$cell,        all.cells)
})

# ── Adjacency correctness ──────────────────────────────────────────────────────

test_that("extract_neighborhoods: bidirectional adjacency is correct on toy graph", {
  el  <- make_toy_edge_list()
  res <- extract_neighborhoods(el, n.cores = 1L, verbose = FALSE)
  
  expect_setequal(res$neighborhood.cell.list[["A"]], c("B", "C"))
  expect_setequal(res$neighborhood.cell.list[["B"]], c("A", "C"))
  expect_setequal(res$neighborhood.cell.list[["C"]], c("A", "B", "D"))
  # D's autocrine edge makes D appear in its own neighbor set
  expect_setequal(res$neighborhood.cell.list[["D"]], c("C", "D"))
})

test_that("extract_neighborhoods: remove.autocrine=TRUE strips self from adjacency", {
  el  <- make_toy_edge_list()
  res <- extract_neighborhoods(el, remove.autocrine = TRUE,
                               n.cores = 1L, verbose = FALSE)
  
  expect_setequal(res$neighborhood.cell.list[["D"]], c("C"))
  expect_true(attr(res, "remove.autocrine"))
})

# ── Edge filter modes ──────────────────────────────────────────────────────────
# A's neighborhood = {A, B, C}
#   both_in_neighborhood : both endpoints in {A,B,C}  → A|B, A|C, B|C
#   any_in_neighborhood  : either endpoint in {A,B,C} → A|B, A|C, B|C, C|D
#   central_cell_only    : A is an endpoint            → A|B, A|C

test_that("extract_neighborhoods: edge.filter.mode = both_in_neighborhood", {
  el  <- make_toy_edge_list()
  res <- extract_neighborhoods(el, edge.filter.mode = "both_in_neighborhood",
                               n.cores = 1L, verbose = FALSE)
  expect_setequal(res$neighborhood.edge.list[["A"]], c("A|B", "A|C", "B|C"))
})

test_that("extract_neighborhoods: edge.filter.mode = any_in_neighborhood", {
  el  <- make_toy_edge_list()
  res <- extract_neighborhoods(el, edge.filter.mode = "any_in_neighborhood",
                               n.cores = 1L, verbose = FALSE)
  expect_setequal(res$neighborhood.edge.list[["A"]], c("A|B", "A|C", "B|C", "C|D"))
})

test_that("extract_neighborhoods: edge.filter.mode = central_cell_only", {
  el  <- make_toy_edge_list()
  res <- extract_neighborhoods(el, edge.filter.mode = "central_cell_only",
                               n.cores = 1L, verbose = FALSE)
  expect_setequal(res$neighborhood.edge.list[["A"]], c("A|B", "A|C"))
  # D has its own autocrine edge plus incoming C|D
  expect_setequal(res$neighborhood.edge.list[["D"]], c("C|D", "D|D"))
})

# ── Meta table count parity ────────────────────────────────────────────────────

test_that("extract_neighborhoods: neighborhood.meta counts match list lengths", {
  el   <- make_toy_edge_list()
  res  <- extract_neighborhoods(el, n.cores = 1L, verbose = FALSE)
  meta <- res$neighborhood.meta
  
  for (this.cell in meta$cell) {
    expect_equal(
      meta[cell == this.cell, n.neighbors],
      length(res$neighborhood.cell.list[[this.cell]])
    )
    expect_equal(
      meta[cell == this.cell, n.edges],
      length(res$neighborhood.edge.list[[this.cell]])
    )
  }
})

# ── Parallel == serial (regression guard for the unname() fix) ─────────────────

test_that("extract_neighborhoods: parallel result matches serial result exactly", {
  skip_on_cran()
  if (parallel::detectCores() < 2L) skip("only one core available")
  
  el       <- make_toy_edge_list()
  r.serial <- extract_neighborhoods(el, n.cores = 1L, verbose = FALSE)
  r.par    <- extract_neighborhoods(el, n.cores = 2L, verbose = FALSE)
  
  # Cell names must be clean — no "1." chunk prefixes from c() on named lists
  expect_false(any(grepl("^[0-9]+\\.", names(r.par$neighborhood.cell.list))))
  expect_false(any(grepl("^[0-9]+\\.", names(r.par$neighborhood.edge.list))))
  
  canon <- function(res) {
    list(
      cells = lapply(res$neighborhood.cell.list[sort(names(res$neighborhood.cell.list))], sort),
      edges = lapply(res$neighborhood.edge.list[sort(names(res$neighborhood.edge.list))], sort)
    )
  }
  expect_identical(canon(r.serial), canon(r.par))
})