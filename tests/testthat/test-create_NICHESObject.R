# ==============================================================================
# FILE:    tests/testthat/test-create_NICHESObject.R
# PURPOSE: Unit tests for create_NICHESObject().
#          All inputs synthesised inline — no shared fixtures needed.
#          The one fixture-dependent test (checking attr(niches.obj,"mode"))
#          is dropped: the mode attribute is covered by the sample.col
#          functional tests in Section 4.
# ==============================================================================

# ── File-level helpers ─────────────────────────────────────────────────────────

.make_count_mtx <- function(n.cells, n.genes = 20L, seed = 1L) {
  set.seed(seed)
  m <- Matrix::rsparsematrix(n.genes, n.cells,
                             density = 0.3,
                             rand.x  = function(n) rpois(n, 2))
  rownames(m) <- paste0("gene", seq_len(n.genes))
  colnames(m) <- paste0("cell", seq_len(n.cells))
  m
}

.make_spatial_meta <- function(n.cells, seed = 1L) {
  set.seed(seed)
  data.frame(
    row.names = paste0("cell", seq_len(n.cells)),
    x         = runif(n.cells, 0, 100),
    y         = runif(n.cells, 0, 100),
    cell.type = sample(c("TypeA", "TypeB"), n.cells, replace = TRUE),
    stringsAsFactors = FALSE
  )
}

.make_sampled_meta <- function(n.cells, seed = 1L) {
  set.seed(seed)
  data.frame(
    row.names = paste0("cell", seq_len(n.cells)),
    cell.type = sample(c("TypeA", "TypeB"), n.cells, replace = TRUE),
    stringsAsFactors = FALSE
  )
}

.make_tiny_LRM_db <- function() {
  data.frame(
    ligand   = c("gene1", "gene3"),
    receptor = c("gene2", "gene4"),
    stringsAsFactors = FALSE
  )
}

# Inputs for sample.col tests: 40 cells split across 2 samples of 20 each.
# Returns a list(count.mtx, meta.data, lrm.db) — the NICHESObject itself is
# built inside each test_that block to preserve block independence.
.make_sample_col_inputs <- function() {
  n.per.sample <- 20L
  n.total      <- n.per.sample * 2L
  n.genes      <- 10L
  
  set.seed(99L)
  count.mtx <- Matrix::rsparsematrix(
    n.genes, n.total, density = 0.4,
    rand.x = function(n) rpois(n, 3)
  )
  rownames(count.mtx) <- paste0("gene", seq_len(n.genes))
  colnames(count.mtx) <- paste0("cell", seq_len(n.total))
  
  meta.data <- data.frame(
    row.names = paste0("cell", seq_len(n.total)),
    x         = c(runif(n.per.sample, 0, 50), runif(n.per.sample, 0, 50)),
    y         = c(runif(n.per.sample, 0, 50), runif(n.per.sample, 0, 50)),
    cell.type = sample(c("TypeA", "TypeB"), n.total, replace = TRUE),
    sample.id = rep(c("S1", "S2"), each = n.per.sample),
    stringsAsFactors = FALSE
  )
  
  list(count.mtx = count.mtx,
       meta.data  = meta.data,
       lrm.db     = .make_tiny_LRM_db())
}

# ── Section 1: spatial mode validation ────────────────────────────────────────

test_that("create_NICHESObject: errors on unknown mode", {
  expect_error(
    create_NICHESObject(
      count.mtx = .make_count_mtx(10L),
      meta.data = .make_spatial_meta(10L),
      LRM.db    = .make_tiny_LRM_db(),
      mode      = "banana"
    ),
    regexp = "arg"
  )
})

test_that("create_NICHESObject: spatial mode errors when x/y absent", {
  expect_error(
    create_NICHESObject(
      count.mtx = .make_count_mtx(10L),
      meta.data = .make_sampled_meta(10L),
      LRM.db    = .make_tiny_LRM_db(),
      mode      = "spatial",
      k         = 3L
    ),
    regexp = "x.*y|y.*x"
  )
})

test_that("create_NICHESObject: spatial mode errors when both k and rad provided", {
  expect_error(
    create_NICHESObject(
      count.mtx = .make_count_mtx(10L),
      meta.data = .make_spatial_meta(10L),
      LRM.db    = .make_tiny_LRM_db(),
      mode      = "spatial",
      k         = 3L,
      rad       = 50
    ),
    regexp = "k.*rad|rad.*k"
  )
})

test_that("create_NICHESObject: spatial mode errors when neither k nor rad provided", {
  expect_error(
    create_NICHESObject(
      count.mtx = .make_count_mtx(10L),
      meta.data = .make_spatial_meta(10L),
      LRM.db    = .make_tiny_LRM_db(),
      mode      = "spatial"
    ),
    regexp = "k.*rad|rad.*k"
  )
})

# ── Section 2: sampled mode validation ────────────────────────────────────────

test_that("create_NICHESObject: sampled mode errors when cell.type.col is NULL", {
  expect_error(
    create_NICHESObject(
      count.mtx     = .make_count_mtx(20L),
      meta.data     = .make_sampled_meta(20L),
      LRM.db        = .make_tiny_LRM_db(),
      mode          = "sampled",
      cell.type.col = NULL
    ),
    regexp = "cell.type.col"
  )
})

test_that("create_NICHESObject: sampled mode errors when cell.type.col not in meta.data", {
  expect_error(
    create_NICHESObject(
      count.mtx     = .make_count_mtx(20L),
      meta.data     = .make_sampled_meta(20L),
      LRM.db        = .make_tiny_LRM_db(),
      mode          = "sampled",
      cell.type.col = "nonexistent"
    ),
    regexp = "cell.type.col"
  )
})

# ── Section 3: sample.col argument validation ──────────────────────────────────

test_that("create_NICHESObject: sample.col must be a single character string", {
  expect_error(
    create_NICHESObject(
      count.mtx  = .make_count_mtx(10L),
      meta.data  = .make_spatial_meta(10L),
      LRM.db     = .make_tiny_LRM_db(),
      mode       = "spatial",
      k          = 3L,
      sample.col = c("col1", "col2")
    ),
    regexp = "sample.col"
  )
})

test_that("create_NICHESObject: sample.col must exist in meta.data", {
  expect_error(
    create_NICHESObject(
      count.mtx  = .make_count_mtx(10L),
      meta.data  = .make_spatial_meta(10L),
      LRM.db     = .make_tiny_LRM_db(),
      mode       = "spatial",
      k          = 3L,
      sample.col = "nonexistent.column"
    ),
    regexp = "sample.col"
  )
})

# ── Section 4: sample.col functional tests ────────────────────────────────────

test_that("create_NICHESObject: sample.col path returns a NICHESObject", {
  inp <- .make_sample_col_inputs()
  obj <- suppressWarnings(create_NICHESObject(
    count.mtx  = inp$count.mtx,
    meta.data  = inp$meta.data,
    LRM.db     = inp$lrm.db,
    mode       = "spatial",
    sample.col = "sample.id",
    k          = 3L,
    n.cores    = 1L,
    verbose    = FALSE
  ))
  expect_true(is_NICHESObject(obj))
})

test_that("create_NICHESObject: sample.col attribute is set on merged object", {
  inp <- .make_sample_col_inputs()
  obj <- suppressWarnings(create_NICHESObject(
    count.mtx  = inp$count.mtx,
    meta.data  = inp$meta.data,
    LRM.db     = inp$lrm.db,
    mode       = "spatial",
    sample.col = "sample.id",
    k          = 3L,
    n.cores    = 1L,
    verbose    = FALSE
  ))
  expect_equal(attr(obj, "sample.col"), "sample.id")
})

test_that("create_NICHESObject: mode attribute preserved through sample.col path", {
  inp <- .make_sample_col_inputs()
  obj <- suppressWarnings(create_NICHESObject(
    count.mtx  = inp$count.mtx,
    meta.data  = inp$meta.data,
    LRM.db     = inp$lrm.db,
    mode       = "spatial",
    sample.col = "sample.id",
    k          = 3L,
    n.cores    = 1L,
    verbose    = FALSE
  ))
  expect_equal(attr(obj, "mode"), "spatial")
})

test_that("create_NICHESObject: sample.col path — all cells present in $cell.list", {
  inp <- .make_sample_col_inputs()
  obj <- suppressWarnings(create_NICHESObject(
    count.mtx  = inp$count.mtx,
    meta.data  = inp$meta.data,
    LRM.db     = inp$lrm.db,
    mode       = "spatial",
    sample.col = "sample.id",
    k          = 3L,
    n.cores    = 1L,
    verbose    = FALSE
  ))
  expect_setequal(obj$cell.list, rownames(inp$meta.data))
})

test_that("create_NICHESObject: no cross-sample edges in merged object", {
  inp <- .make_sample_col_inputs()
  obj <- suppressWarnings(create_NICHESObject(
    count.mtx  = inp$count.mtx,
    meta.data  = inp$meta.data,
    LRM.db     = inp$lrm.db,
    mode       = "spatial",
    sample.col = "sample.id",
    k          = 3L,
    n.cores    = 1L,
    verbose    = FALSE
  ))
  barcode.to.sample <- setNames(inp$meta.data$sample.id,
                                rownames(inp$meta.data))
  sending.sample   <- barcode.to.sample[obj$edge.list$SendingCell]
  receiving.sample <- barcode.to.sample[obj$edge.list$ReceivingCell]
  expect_true(all(sending.sample == receiving.sample))
})

test_that("create_NICHESObject: $cell.meta retains sample.id column", {
  inp <- .make_sample_col_inputs()
  obj <- suppressWarnings(create_NICHESObject(
    count.mtx  = inp$count.mtx,
    meta.data  = inp$meta.data,
    LRM.db     = inp$lrm.db,
    mode       = "spatial",
    sample.col = "sample.id",
    k          = 3L,
    n.cores    = 1L,
    verbose    = FALSE
  ))
  expect_true("sample.id" %in% colnames(obj$cell.meta))
})

test_that("create_NICHESObject: $edge.meta has Sending.sample.id and Receiving.sample.id", {
  inp <- .make_sample_col_inputs()
  obj <- suppressWarnings(create_NICHESObject(
    count.mtx  = inp$count.mtx,
    meta.data  = inp$meta.data,
    LRM.db     = inp$lrm.db,
    mode       = "spatial",
    sample.col = "sample.id",
    k          = 3L,
    n.cores    = 1L,
    verbose    = FALSE
  ))
  expect_true("Sending.sample.id"   %in% colnames(obj$edge.meta))
  expect_true("Receiving.sample.id" %in% colnames(obj$edge.meta))
})

test_that("create_NICHESObject: NULL sample.col produces no sample.col attribute", {
  inp <- .make_sample_col_inputs()
  obj <- suppressWarnings(create_NICHESObject(
    count.mtx  = inp$count.mtx,
    meta.data  = inp$meta.data,
    LRM.db     = inp$lrm.db,
    mode       = "spatial",
    sample.col = NULL,
    k          = 3L,
    n.cores    = 1L,
    verbose    = FALSE
  ))
  expect_null(attr(obj, "sample.col"))
})

# ── Section 5: sample.col edge cases ──────────────────────────────────────────

test_that("create_NICHESObject: single-cell sample emits warning; remaining samples merged", {
  n.genes <- 8L
  set.seed(7L)
  count.mtx <- Matrix::rsparsematrix(n.genes, 22L, density = 0.5,
                                     rand.x = function(n) rpois(n, 2))
  rownames(count.mtx) <- paste0("gene", seq_len(n.genes))
  colnames(count.mtx) <- paste0("cell", seq_len(22L))
  
  meta.data <- data.frame(
    row.names = paste0("cell", seq_len(22L)),
    x         = runif(22L, 0, 100),
    y         = runif(22L, 0, 100),
    sample.id = c("singleton", rep("normal", 21L)),
    stringsAsFactors = FALSE
  )
  
  expect_warning(
    obj <- create_NICHESObject(
      count.mtx  = count.mtx,
      meta.data  = meta.data,
      LRM.db     = .make_tiny_LRM_db(),
      mode       = "spatial",
      sample.col = "sample.id",
      k          = 3L,
      n.cores    = 1L,
      verbose    = FALSE
    ),
    regexp = "fewer than 2 cells"
  )
  expect_true(is_NICHESObject(obj))
  expect_equal(length(obj$cell.list), 21L)
})

# ── Section 6: neighborhood slots ─────────────────────────────────────────────

test_that("create_NICHESObject: neighborhood slots non-NULL by default", {
  inp <- .make_sample_col_inputs()
  obj <- suppressWarnings(create_NICHESObject(
    count.mtx  = inp$count.mtx,
    meta.data  = inp$meta.data,
    LRM.db     = inp$lrm.db,
    mode       = "spatial",
    k          = 3L,
    n.cores    = 1L,
    verbose    = FALSE
  ))
  expect_false(is.null(obj$neighborhood.cell.list))
  expect_false(is.null(obj$neighborhood.edge.list))
  expect_false(is.null(obj$neighborhood.meta))
})

test_that("create_NICHESObject: edge.filter.mode attribute set correctly", {
  inp <- .make_sample_col_inputs()
  obj <- suppressWarnings(create_NICHESObject(
    count.mtx        = inp$count.mtx,
    meta.data        = inp$meta.data,
    LRM.db           = inp$lrm.db,
    mode             = "spatial",
    k                = 3L,
    edge.filter.mode = "central_cell_only",
    n.cores          = 1L,
    verbose          = FALSE
  ))
  expect_equal(attr(obj, "edge.filter.mode"), "central_cell_only")
})

test_that("create_NICHESObject: remove.autocrine attribute set correctly", {
  inp <- .make_sample_col_inputs()
  obj <- suppressWarnings(create_NICHESObject(
    count.mtx        = inp$count.mtx,
    meta.data        = inp$meta.data,
    LRM.db           = inp$lrm.db,
    mode             = "spatial",
    k                = 3L,
    remove.autocrine = TRUE,
    n.cores          = 1L,
    verbose          = FALSE
  ))
  expect_true(attr(obj, "remove.autocrine"))
})

test_that("create_NICHESObject: edge.filter.mode = NULL leaves neighborhood slots NULL", {
  inp <- .make_sample_col_inputs()
  obj <- suppressWarnings(create_NICHESObject(
    count.mtx        = inp$count.mtx,
    meta.data        = inp$meta.data,
    LRM.db           = inp$lrm.db,
    mode             = "spatial",
    k                = 3L,
    edge.filter.mode = NULL,
    n.cores          = 1L,
    verbose          = FALSE
  ))
  expect_null(obj$neighborhood.cell.list)
  expect_null(obj$neighborhood.edge.list)
  expect_null(obj$neighborhood.meta)
})