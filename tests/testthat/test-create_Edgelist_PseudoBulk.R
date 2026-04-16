# ==============================================================================
# FILE:    tests/testthat/test-create_Edgelist_PseudoBulk.R
# PURPOSE: Unit tests for create_Edgelist_PseudoBulk().
#          All inputs synthesised inline — no shared fixtures needed.
# ==============================================================================

# File-level helper: 15 cells across 3 types — small, deterministic, reusable.
make_pseudobulk_meta <- function() {
  data.frame(
    row.names = paste0("cell_", seq_len(15L)),
    celltype  = rep(c("TypeA", "TypeB", "TypeC"), each = 5L),
    stringsAsFactors = FALSE
  )
}

# ── Structure and schema ───────────────────────────────────────────────────────

test_that("create_Edgelist_PseudoBulk: returns a data.table with from and to", {
  meta <- make_pseudobulk_meta()
  el   <- create_Edgelist_PseudoBulk(meta, cell.type.col = "celltype")
  
  expect_s3_class(el, "data.table")
  expect_true(all(c("from", "to") %in% colnames(el)))
})

test_that("create_Edgelist_PseudoBulk: from/to values are type labels not barcodes", {
  meta         <- make_pseudobulk_meta()
  el           <- create_Edgelist_PseudoBulk(meta, cell.type.col = "celltype")
  type.labels  <- unique(meta$celltype)
  all.barcodes <- rownames(meta)
  
  expect_true(all(el$from %in% type.labels))
  expect_true(all(el$to   %in% type.labels))
  expect_false(any(el$from %in% all.barcodes))
})

test_that("create_Edgelist_PseudoBulk: produces exactly T^2 edges for T types", {
  meta    <- make_pseudobulk_meta()
  el      <- create_Edgelist_PseudoBulk(meta, cell.type.col = "celltype")
  n.types <- length(unique(meta$celltype))
  
  expect_equal(nrow(el), n.types^2L)
})

test_that("create_Edgelist_PseudoBulk: all ordered type pairs are present", {
  meta        <- make_pseudobulk_meta()
  el          <- create_Edgelist_PseudoBulk(meta, cell.type.col = "celltype")
  type.labels <- sort(unique(meta$celltype))
  
  expected.from <- rep(type.labels, each  = length(type.labels))
  expected.to   <- rep(type.labels, times = length(type.labels))
  actual        <- el[order(el$from, el$to), ]
  
  expect_equal(actual$from, sort(expected.from))
  expect_equal(actual$to,   expected.to[order(expected.from)])
})

test_that("create_Edgelist_PseudoBulk: autocrine edges (from == to) are included", {
  meta      <- make_pseudobulk_meta()
  el        <- create_Edgelist_PseudoBulk(meta, cell.type.col = "celltype")
  autocrine <- el[el$from == el$to, ]
  n.types   <- length(unique(meta$celltype))
  
  expect_equal(nrow(autocrine), n.types)
})

test_that("create_Edgelist_PseudoBulk: output is deterministic", {
  meta <- make_pseudobulk_meta()
  el1  <- create_Edgelist_PseudoBulk(meta, cell.type.col = "celltype")
  el2  <- create_Edgelist_PseudoBulk(meta, cell.type.col = "celltype")
  
  expect_equal(el1$from, el2$from)
  expect_equal(el1$to,   el2$to)
})

# ── Attributes ────────────────────────────────────────────────────────────────

test_that("create_Edgelist_PseudoBulk: mode attribute is 'pseudobulk'", {
  meta <- make_pseudobulk_meta()
  el   <- create_Edgelist_PseudoBulk(meta, cell.type.col = "celltype")
  
  expect_equal(attr(el, "mode"), "pseudobulk")
})

test_that("create_Edgelist_PseudoBulk: n.cells attribute matches nrow(meta)", {
  meta <- make_pseudobulk_meta()
  el   <- create_Edgelist_PseudoBulk(meta, cell.type.col = "celltype")
  
  expect_equal(attr(el, "n.cells"), nrow(meta))
})

test_that("create_Edgelist_PseudoBulk: cell.names attribute is sorted unique types", {
  meta <- make_pseudobulk_meta()
  el   <- create_Edgelist_PseudoBulk(meta, cell.type.col = "celltype")
  
  expect_equal(attr(el, "cell.names"), sort(unique(meta$celltype)))
})

test_that("create_Edgelist_PseudoBulk: param attribute records cell.type.col", {
  meta <- make_pseudobulk_meta()
  el   <- create_Edgelist_PseudoBulk(meta, cell.type.col = "celltype")
  
  expect_equal(attr(el, "param"), "celltype")
})

# ── Input validation ───────────────────────────────────────────────────────────

test_that("create_Edgelist_PseudoBulk: errors on non-data.frame input", {
  expect_error(
    create_Edgelist_PseudoBulk(matrix(1:4, 2L), cell.type.col = "celltype"),
    regexp = "data.frame"
  )
})

test_that("create_Edgelist_PseudoBulk: errors when cell.type.col not in meta", {
  meta <- make_pseudobulk_meta()
  expect_error(
    create_Edgelist_PseudoBulk(meta, cell.type.col = "NOTACOL"),
    regexp = "not found"
  )
})

test_that("create_Edgelist_PseudoBulk: errors when fewer than 2 cell types", {
  mono.meta <- data.frame(
    celltype  = rep("OnlyType", 5L),
    row.names = paste0("cell_", 1:5)
  )
  expect_error(
    create_Edgelist_PseudoBulk(mono.meta, cell.type.col = "celltype"),
    regexp = "2 distinct"
  )
})