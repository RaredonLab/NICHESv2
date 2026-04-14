# ==============================================================================
# FILE:    tests/testthat/test-create_Edgelist_SamplingByType.R
# ==============================================================================

# ── Structure and schema ───────────────────────────────────────────────────────

test_that("create_Edgelist_SamplingByType: returns a data.table", {
  .setup_fixtures()
  el <- create_Edgelist_SamplingByType(meta.data, cell.type.col = "celltype", seed = 42L)
  expect_s3_class(el, "data.table")
  expect_true(all(c("from", "to") %in% colnames(el)))
})

test_that("create_Edgelist_SamplingByType: from/to values are cell barcodes", {
  .setup_fixtures()
  el           <- create_Edgelist_SamplingByType(meta.data, cell.type.col = "celltype", seed = 42L)
  all.barcodes <- rownames(meta.data)
  expect_true(all(el$from %in% all.barcodes))
  expect_true(all(el$to   %in% all.barcodes))
})

test_that("create_Edgelist_SamplingByType: from/to values are not type labels", {
  .setup_fixtures()
  el          <- create_Edgelist_SamplingByType(meta.data, cell.type.col = "celltype", seed = 42L)
  type.labels <- unique(meta.data$celltype)
  # Edge identifiers should be barcodes, not the type name strings
  expect_false(all(el$from %in% type.labels))
})

test_that("create_Edgelist_SamplingByType: all T² type crossings are represented", {
  .setup_fixtures()
  el          <- create_Edgelist_SamplingByType(meta.data, cell.type.col = "celltype", seed = 42L)
  type.labels <- unique(meta.data$celltype)
  # For each ordered type pair, at least one edge should exist
  for (s in type.labels) {
    s.barcodes <- rownames(meta.data)[meta.data$celltype == s]
    for (r in type.labels) {
      r.barcodes <- rownames(meta.data)[meta.data$celltype == r]
      pair.edges <- el[el$from %in% s.barcodes & el$to %in% r.barcodes, ]
      expect_gt(nrow(pair.edges), 0L)
    }
  }
})

test_that("create_Edgelist_SamplingByType: per-crossing edge count respects min(n_i, n_j)", {
  .setup_fixtures()
  el          <- create_Edgelist_SamplingByType(meta.data, cell.type.col = "celltype", seed = 42L)
  type.labels <- unique(meta.data$celltype)
  for (s in type.labels) {
    s.barcodes <- rownames(meta.data)[meta.data$celltype == s]
    for (r in type.labels) {
      r.barcodes <- rownames(meta.data)[meta.data$celltype == r]
      cap        <- min(length(s.barcodes), length(r.barcodes))
      pair.edges <- el[el$from %in% s.barcodes & el$to %in% r.barcodes, ]
      expect_lte(nrow(pair.edges), cap)
    }
  }
})

test_that("create_Edgelist_SamplingByType: max.cells caps edges per crossing", {
  .setup_fixtures()
  cap <- 2L
  el  <- create_Edgelist_SamplingByType(meta.data, cell.type.col = "celltype",
                                         seed = 42L, max.cells = cap)
  type.labels <- unique(meta.data$celltype)
  for (s in type.labels) {
    s.barcodes <- rownames(meta.data)[meta.data$celltype == s]
    for (r in type.labels) {
      r.barcodes <- rownames(meta.data)[meta.data$celltype == r]
      pair.edges <- el[el$from %in% s.barcodes & el$to %in% r.barcodes, ]
      expect_lte(nrow(pair.edges), cap)
    }
  }
})

# ── Reproducibility ────────────────────────────────────────────────────────────

test_that("create_Edgelist_SamplingByType: same seed gives identical output", {
  .setup_fixtures()
  el1 <- create_Edgelist_SamplingByType(meta.data, cell.type.col = "celltype", seed = 42L)
  el2 <- create_Edgelist_SamplingByType(meta.data, cell.type.col = "celltype", seed = 42L)
  expect_equal(el1$from, el2$from)
  expect_equal(el1$to,   el2$to)
})

test_that("create_Edgelist_SamplingByType: different seed gives different output", {
  .setup_fixtures()
  el1 <- create_Edgelist_SamplingByType(meta.data, cell.type.col = "celltype", seed = 42L)
  el2 <- create_Edgelist_SamplingByType(meta.data, cell.type.col = "celltype", seed = 99L)
  # It is astronomically unlikely these are identical
  expect_false(identical(el1$from, el2$from))
})

# ── Metadata attributes ────────────────────────────────────────────────────────

test_that("create_Edgelist_SamplingByType: mode attribute is 'sampled'", {
  .setup_fixtures()
  el <- create_Edgelist_SamplingByType(meta.data, cell.type.col = "celltype", seed = 42L)
  expect_equal(attr(el, "mode"), "sampled")
})

test_that("create_Edgelist_SamplingByType: n.cells attribute matches nrow(meta.data)", {
  .setup_fixtures()
  el <- create_Edgelist_SamplingByType(meta.data, cell.type.col = "celltype", seed = 42L)
  expect_equal(attr(el, "n.cells"), nrow(meta.data))
})

test_that("create_Edgelist_SamplingByType: param records seed and max.cells", {
  .setup_fixtures()
  el <- create_Edgelist_SamplingByType(meta.data, cell.type.col = "celltype",
                                        seed = 42L, max.cells = 10L)
  p  <- attr(el, "param")
  expect_equal(p$seed,          42L)
  expect_equal(p$max.cells,     10L)
  expect_equal(p$cell.type.col, "celltype")
})

# ── Input validation ───────────────────────────────────────────────────────────

test_that("create_Edgelist_SamplingByType: errors on non-data.frame input", {
  expect_error(
    create_Edgelist_SamplingByType(matrix(1:4, 2), cell.type.col = "celltype"),
    regexp = "data.frame"
  )
})

test_that("create_Edgelist_SamplingByType: errors when cell.type.col not in meta.data", {
  .setup_fixtures()
  expect_error(
    create_Edgelist_SamplingByType(meta.data, cell.type.col = "NOTACOL"),
    regexp = "not found"
  )
})

test_that("create_Edgelist_SamplingByType: errors when fewer than 2 cell types", {
  mono.meta <- data.frame(celltype = rep("OnlyType", 5),
                           row.names = paste0("cell_", 1:5))
  expect_error(
    create_Edgelist_SamplingByType(mono.meta, cell.type.col = "celltype"),
    regexp = "2 distinct"
  )
})

test_that("create_Edgelist_SamplingByType: errors on invalid max.cells", {
  .setup_fixtures()
  expect_error(
    create_Edgelist_SamplingByType(meta.data, cell.type.col = "celltype",
                                    seed = 42L, max.cells = 0L),
    regexp = "max.cells"
  )
  expect_error(
    create_Edgelist_SamplingByType(meta.data, cell.type.col = "celltype",
                                    seed = 42L, max.cells = -5L),
    regexp = "max.cells"
  )
})
