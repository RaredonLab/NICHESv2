# ==============================================================================
# FILE:    tests/testthat/test-create_Edgelist_SamplingByType.R
# PURPOSE: Unit tests for create_Edgelist_SamplingByType().
#          All inputs synthesised inline — no shared fixtures needed.
# ==============================================================================

# File-level helper: 12 cells across 2 types, 6 each — enough for all
# per-crossing cap and sampling assertions.
make_sampling_meta <- function() {
  data.frame(
    row.names = paste0("cell_", seq_len(12L)),
    celltype  = rep(c("TypeA", "TypeB"), each = 6L),
    stringsAsFactors = FALSE
  )
}

# ── Structure and schema ───────────────────────────────────────────────────────

test_that("create_Edgelist_SamplingByType: returns a data.table with from and to", {
  meta <- make_sampling_meta()
  el   <- create_Edgelist_SamplingByType(meta, cell.type.col = "celltype",
                                         seed = 42L)
  expect_s3_class(el, "data.table")
  expect_true(all(c("from", "to") %in% colnames(el)))
})

test_that("create_Edgelist_SamplingByType: from/to values are cell barcodes", {
  meta         <- make_sampling_meta()
  el           <- create_Edgelist_SamplingByType(meta, cell.type.col = "celltype",
                                                 seed = 42L)
  all.barcodes <- rownames(meta)
  
  expect_true(all(el$from %in% all.barcodes))
  expect_true(all(el$to   %in% all.barcodes))
})

test_that("create_Edgelist_SamplingByType: from/to values are not type labels", {
  meta        <- make_sampling_meta()
  el          <- create_Edgelist_SamplingByType(meta, cell.type.col = "celltype",
                                                seed = 42L)
  type.labels <- unique(meta$celltype)
  
  expect_false(all(el$from %in% type.labels))
})

test_that("create_Edgelist_SamplingByType: all T^2 type crossings are represented", {
  meta        <- make_sampling_meta()
  el          <- create_Edgelist_SamplingByType(meta, cell.type.col = "celltype",
                                                seed = 42L)
  type.labels <- unique(meta$celltype)
  
  for (s in type.labels) {
    s.barcodes <- rownames(meta)[meta$celltype == s]
    for (r in type.labels) {
      r.barcodes <- rownames(meta)[meta$celltype == r]
      pair.edges <- el[el$from %in% s.barcodes & el$to %in% r.barcodes, ]
      expect_gt(nrow(pair.edges), 0L)
    }
  }
})

test_that("create_Edgelist_SamplingByType: per-crossing edge count respects min(n_i, n_j)", {
  meta        <- make_sampling_meta()
  el          <- create_Edgelist_SamplingByType(meta, cell.type.col = "celltype",
                                                seed = 42L)
  type.labels <- unique(meta$celltype)
  
  for (s in type.labels) {
    s.barcodes <- rownames(meta)[meta$celltype == s]
    for (r in type.labels) {
      r.barcodes <- rownames(meta)[meta$celltype == r]
      cap        <- min(length(s.barcodes), length(r.barcodes))
      pair.edges <- el[el$from %in% s.barcodes & el$to %in% r.barcodes, ]
      expect_lte(nrow(pair.edges), cap)
    }
  }
})

test_that("create_Edgelist_SamplingByType: max.cells caps edges per crossing", {
  meta        <- make_sampling_meta()
  cap         <- 2L
  el          <- create_Edgelist_SamplingByType(meta, cell.type.col = "celltype",
                                                seed = 42L, max.cells = cap)
  type.labels <- unique(meta$celltype)
  
  for (s in type.labels) {
    s.barcodes <- rownames(meta)[meta$celltype == s]
    for (r in type.labels) {
      r.barcodes <- rownames(meta)[meta$celltype == r]
      pair.edges <- el[el$from %in% s.barcodes & el$to %in% r.barcodes, ]
      expect_lte(nrow(pair.edges), cap)
    }
  }
})

# ── Reproducibility ────────────────────────────────────────────────────────────

test_that("create_Edgelist_SamplingByType: same seed gives identical output", {
  meta <- make_sampling_meta()
  el1  <- create_Edgelist_SamplingByType(meta, cell.type.col = "celltype",
                                         seed = 42L)
  el2  <- create_Edgelist_SamplingByType(meta, cell.type.col = "celltype",
                                         seed = 42L)
  expect_equal(el1$from, el2$from)
  expect_equal(el1$to,   el2$to)
})

test_that("create_Edgelist_SamplingByType: different seed gives different output", {
  meta <- make_sampling_meta()
  el1  <- create_Edgelist_SamplingByType(meta, cell.type.col = "celltype",
                                         seed = 42L, max.cells = 3L)
  el2  <- create_Edgelist_SamplingByType(meta, cell.type.col = "celltype",
                                         seed = 99L, max.cells = 3L)
  expect_false(identical(el1$from, el2$from))
})

# ── Attributes ────────────────────────────────────────────────────────────────

test_that("create_Edgelist_SamplingByType: mode attribute is 'nonspatial'", {
  meta <- make_sampling_meta()
  el   <- create_Edgelist_SamplingByType(meta, cell.type.col = "celltype",
                                         seed = 42L)
  expect_equal(attr(el, "mode"), "nonspatial")
})

test_that("create_Edgelist_SamplingByType: n.cells attribute matches nrow(meta)", {
  meta <- make_sampling_meta()
  el   <- create_Edgelist_SamplingByType(meta, cell.type.col = "celltype",
                                         seed = 42L)
  expect_equal(attr(el, "n.cells"), nrow(meta))
})

test_that("create_Edgelist_SamplingByType: param attribute records seed, max.cells, cell.type.col", {
  meta <- make_sampling_meta()
  el   <- create_Edgelist_SamplingByType(meta, cell.type.col = "celltype",
                                         seed = 42L, max.cells = 10L)
  p    <- attr(el, "param")
  
  expect_equal(p$seed,          42L)
  expect_equal(p$max.cells,     10L)
  expect_equal(p$cell.type.col, "celltype")
})

# ── Input validation ───────────────────────────────────────────────────────────

test_that("create_Edgelist_SamplingByType: errors on non-data.frame input", {
  expect_error(
    create_Edgelist_SamplingByType(matrix(1:4, 2L), cell.type.col = "celltype"),
    regexp = "data.frame"
  )
})

test_that("create_Edgelist_SamplingByType: errors when cell.type.col not in meta", {
  meta <- make_sampling_meta()
  expect_error(
    create_Edgelist_SamplingByType(meta, cell.type.col = "NOTACOL"),
    regexp = "not found"
  )
})

test_that("create_Edgelist_SamplingByType: errors when fewer than 2 cell types", {
  mono.meta <- data.frame(
    celltype  = rep("OnlyType", 5L),
    row.names = paste0("cell_", 1:5)
  )
  expect_error(
    create_Edgelist_SamplingByType(mono.meta, cell.type.col = "celltype"),
    regexp = "2 distinct"
  )
})

test_that("create_Edgelist_SamplingByType: errors on invalid max.cells values", {
  meta <- make_sampling_meta()
  expect_error(
    create_Edgelist_SamplingByType(meta, cell.type.col = "celltype",
                                   seed = 42L, max.cells = 0L),
    regexp = "max.cells"
  )
  expect_error(
    create_Edgelist_SamplingByType(meta, cell.type.col = "celltype",
                                   seed = 42L, max.cells = -5L),
    regexp = "max.cells"
  )
})

test_that("create_Edgelist_SamplingByType restores the caller's random state", {
  meta <- make_sampling_meta()
  set.seed(99L)
  before <- get(".Random.seed", envir = globalenv())
  create_Edgelist_SamplingByType(meta, cell.type.col = "celltype", seed = 42L)
  expect_identical(get(".Random.seed", envir = globalenv()), before)
})