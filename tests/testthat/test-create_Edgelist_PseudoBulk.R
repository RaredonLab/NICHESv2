# ==============================================================================
# FILE:    tests/testthat/test-create_Edgelist_PseudoBulk.R
# ==============================================================================

# ── Structure and schema ───────────────────────────────────────────────────────

test_that("create_Edgelist_PseudoBulk: returns a data.table", {
  .setup_fixtures()
  el <- create_Edgelist_PseudoBulk(meta.data, cell.type.col = "celltype")
  expect_s3_class(el, "data.table")
  expect_true(all(c("from", "to") %in% colnames(el)))
})

test_that("create_Edgelist_PseudoBulk: from/to values are type labels not barcodes", {
  .setup_fixtures()
  el          <- create_Edgelist_PseudoBulk(meta.data, cell.type.col = "celltype")
  type.labels <- unique(meta.data$celltype)
  all.barcodes <- rownames(meta.data)
  expect_true(all(el$from %in% type.labels))
  expect_true(all(el$to   %in% type.labels))
  # Confirm they are NOT barcodes (no from value should be a raw barcode)
  expect_false(any(el$from %in% all.barcodes))
})

test_that("create_Edgelist_PseudoBulk: produces exactly T² edges", {
  .setup_fixtures()
  el     <- create_Edgelist_PseudoBulk(meta.data, cell.type.col = "celltype")
  n.types <- length(unique(meta.data$celltype))
  expect_equal(nrow(el), n.types^2)
})

test_that("create_Edgelist_PseudoBulk: all type pairs are present", {
  .setup_fixtures()
  el          <- create_Edgelist_PseudoBulk(meta.data, cell.type.col = "celltype")
  type.labels <- sort(unique(meta.data$celltype))
  expected.from <- rep(type.labels, each  = length(type.labels))
  expected.to   <- rep(type.labels, times = length(type.labels))
  # After sorting both, they should match exactly
  actual <- el[order(el$from, el$to), ]
  expect_equal(actual$from, sort(expected.from))
  expect_equal(actual$to,   expected.to[order(expected.from)])
})

test_that("create_Edgelist_PseudoBulk: autocrine edges (from == to) are included", {
  .setup_fixtures()
  el         <- create_Edgelist_PseudoBulk(meta.data, cell.type.col = "celltype")
  autocrine  <- el[el$from == el$to, ]
  n.types    <- length(unique(meta.data$celltype))
  expect_equal(nrow(autocrine), n.types)
})

test_that("create_Edgelist_PseudoBulk: output is deterministic (no randomness)", {
  .setup_fixtures()
  el1 <- create_Edgelist_PseudoBulk(meta.data, cell.type.col = "celltype")
  el2 <- create_Edgelist_PseudoBulk(meta.data, cell.type.col = "celltype")
  expect_equal(el1$from, el2$from)
  expect_equal(el1$to,   el2$to)
})

# ── Metadata attributes ────────────────────────────────────────────────────────

test_that("create_Edgelist_PseudoBulk: mode attribute is 'pseudobulk'", {
  .setup_fixtures()
  el <- create_Edgelist_PseudoBulk(meta.data, cell.type.col = "celltype")
  expect_equal(attr(el, "mode"), "pseudobulk")
})

test_that("create_Edgelist_PseudoBulk: n.cells attribute matches nrow(meta.data)", {
  .setup_fixtures()
  el <- create_Edgelist_PseudoBulk(meta.data, cell.type.col = "celltype")
  expect_equal(attr(el, "n.cells"), nrow(meta.data))
})

test_that("create_Edgelist_PseudoBulk: cell.names attribute is sorted unique types", {
  .setup_fixtures()
  el          <- create_Edgelist_PseudoBulk(meta.data, cell.type.col = "celltype")
  expect_equal(attr(el, "cell.names"), sort(unique(meta.data$celltype)))
})

test_that("create_Edgelist_PseudoBulk: param attribute records cell.type.col", {
  .setup_fixtures()
  el <- create_Edgelist_PseudoBulk(meta.data, cell.type.col = "celltype")
  expect_equal(attr(el, "param"), "celltype")
})

# ── Input validation ───────────────────────────────────────────────────────────

test_that("create_Edgelist_PseudoBulk: errors on non-data.frame input", {
  expect_error(
    create_Edgelist_PseudoBulk(matrix(1:4, 2), cell.type.col = "celltype"),
    regexp = "data.frame"
  )
})

test_that("create_Edgelist_PseudoBulk: errors when cell.type.col not in meta.data", {
  .setup_fixtures()
  expect_error(
    create_Edgelist_PseudoBulk(meta.data, cell.type.col = "NOTACOL"),
    regexp = "not found"
  )
})

test_that("create_Edgelist_PseudoBulk: errors when fewer than 2 cell types", {
  mono.meta <- data.frame(celltype = rep("OnlyType", 5),
                           row.names = paste0("cell_", 1:5))
  expect_error(
    create_Edgelist_PseudoBulk(mono.meta, cell.type.col = "celltype"),
    regexp = "2 distinct"
  )
})
