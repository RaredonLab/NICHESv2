## test-create_NICHESObject.R
## Tests for create_NICHESObject() — including new sample.col path.
## Existing tests for spatial and sampled modes are preserved below.
## New sample.col tests are in the final section.

# ---------------------------------------------------------------------------
# Shared fixture (loaded from inst/testdata via helper-fixtures.R)
# ---------------------------------------------------------------------------
# niches.obj — spatial NICHESObject, attr "mode" == "spatial"
# Loaded automatically by testthat via helper-fixtures.R

# ---------------------------------------------------------------------------
# Helpers used in this file only
# ---------------------------------------------------------------------------

## Build a minimal count matrix for N cells and G genes
.make_count_mtx <- function(n.cells, n.genes = 20L, seed = 1L) {
  set.seed(seed)
  m <- Matrix::rsparsematrix(n.genes, n.cells,
                             density = 0.3, rand.x = function(n) rpois(n, 2))
  rownames(m) <- paste0("gene", seq_len(n.genes))
  colnames(m) <- paste0("cell", seq_len(n.cells))
  m
}

## Build minimal spatial meta.data (with x, y)
.make_spatial_meta <- function(n.cells, seed = 1L) {
  set.seed(seed)
  data.frame(
    row.names  = paste0("cell", seq_len(n.cells)),
    x          = runif(n.cells, 0, 100),
    y          = runif(n.cells, 0, 100),
    cell.type  = sample(c("TypeA", "TypeB"), n.cells, replace = TRUE),
    stringsAsFactors = FALSE
  )
}

## Build minimal sampled meta.data (no x, y)
.make_sampled_meta <- function(n.cells, seed = 1L) {
  set.seed(seed)
  data.frame(
    row.names  = paste0("cell", seq_len(n.cells)),
    cell.type  = sample(c("TypeA", "TypeB"), n.cells, replace = TRUE),
    stringsAsFactors = FALSE
  )
}

## Tiny two-gene, two-LRM database — lowercase to match default ligand.col/receptor.col
.make_tiny_LRM_db <- function() {
  data.frame(
    ligand   = c("gene1", "gene3"),
    receptor = c("gene2", "gene4"),
    stringsAsFactors = FALSE
  )
}


# ===========================================================================
# SECTION 1: Existing spatial mode tests (unchanged)
# ===========================================================================

test_that("fixture NICHESObject has correct mode attribute", {
  expect_equal(attr(niches.obj, "mode"), "spatial")
})

test_that("create_NICHESObject errors on unknown mode", {
  expect_error(
    create_NICHESObject(
      count.mtx = .make_count_mtx(10),
      meta.data = .make_spatial_meta(10),
      LRM.db    = .make_tiny_LRM_db(),
      mode      = "banana"
    ),
    regexp = "arg"   # match.arg error
  )
})

test_that("spatial mode errors when x/y absent", {
  meta.no.xy <- .make_sampled_meta(10)   # no x, y
  expect_error(
    create_NICHESObject(
      count.mtx = .make_count_mtx(10),
      meta.data = meta.no.xy,
      LRM.db    = .make_tiny_LRM_db(),
      mode      = "spatial",
      k         = 3L
    ),
    regexp = "x.*y|y.*x"
  )
})

test_that("spatial mode errors when both k and rad are provided", {
  expect_error(
    create_NICHESObject(
      count.mtx = .make_count_mtx(10),
      meta.data = .make_spatial_meta(10),
      LRM.db    = .make_tiny_LRM_db(),
      mode      = "spatial",
      k         = 3L,
      rad       = 50
    ),
    regexp = "k.*rad|rad.*k"
  )
})

test_that("spatial mode errors when neither k nor rad is provided", {
  expect_error(
    create_NICHESObject(
      count.mtx = .make_count_mtx(10),
      meta.data = .make_spatial_meta(10),
      LRM.db    = .make_tiny_LRM_db(),
      mode      = "spatial"
    ),
    regexp = "k.*rad|rad.*k"
  )
})


# ===========================================================================
# SECTION 2: Existing sampled mode tests (unchanged)
# ===========================================================================

test_that("sampled mode errors when cell.type.col is NULL", {
  expect_error(
    create_NICHESObject(
      count.mtx     = .make_count_mtx(20),
      meta.data     = .make_sampled_meta(20),
      LRM.db        = .make_tiny_LRM_db(),
      mode          = "sampled",
      cell.type.col = NULL
    ),
    regexp = "cell.type.col"
  )
})

test_that("sampled mode errors when cell.type.col not in meta.data", {
  expect_error(
    create_NICHESObject(
      count.mtx     = .make_count_mtx(20),
      meta.data     = .make_sampled_meta(20),
      LRM.db        = .make_tiny_LRM_db(),
      mode          = "sampled",
      cell.type.col = "nonexistent"
    ),
    regexp = "cell.type.col"
  )
})


# ===========================================================================
# SECTION 3: NEW — sample.col validation tests
# ===========================================================================

test_that("sample.col must be a single character string", {
  expect_error(
    create_NICHESObject(
      count.mtx  = .make_count_mtx(10),
      meta.data  = .make_spatial_meta(10),
      LRM.db     = .make_tiny_LRM_db(),
      mode       = "spatial",
      k          = 3L,
      sample.col = c("col1", "col2")   # vector, not scalar
    ),
    regexp = "sample.col"
  )
})

test_that("sample.col must exist in meta.data", {
  expect_error(
    create_NICHESObject(
      count.mtx  = .make_count_mtx(10),
      meta.data  = .make_spatial_meta(10),
      LRM.db     = .make_tiny_LRM_db(),
      mode       = "spatial",
      k          = 3L,
      sample.col = "nonexistent.column"
    ),
    regexp = "sample.col"
  )
})


# ===========================================================================
# SECTION 4: NEW — sample.col functional tests
# ===========================================================================

## Shared setup: 40 cells split across 2 samples (20 each)
## Uses spatial mode with k = 3 and a tiny LRM database

local({

  n.per.sample <- 20L
  n.total      <- n.per.sample * 2L
  n.genes      <- 10L
  lrm.db       <- .make_tiny_LRM_db()

  set.seed(99L)
  count.mtx <- Matrix::rsparsematrix(
    n.genes, n.total,
    density = 0.4,
    rand.x  = function(n) rpois(n, 3)
  )
  rownames(count.mtx) <- paste0("gene", seq_len(n.genes))
  colnames(count.mtx) <- paste0("cell", seq_len(n.total))

  meta.data <- data.frame(
    row.names  = paste0("cell", seq_len(n.total)),
    x          = c(runif(n.per.sample, 0, 50),   runif(n.per.sample, 0, 50)),
    y          = c(runif(n.per.sample, 0, 50),   runif(n.per.sample, 0, 50)),
    cell.type  = sample(c("TypeA", "TypeB"), n.total, replace = TRUE),
    sample.id  = rep(c("S1", "S2"), each = n.per.sample),
    stringsAsFactors = FALSE
  )

  test_that("sample.col path returns a NICHESObject", {
    obj <- suppressWarnings(
      create_NICHESObject(
        count.mtx  = count.mtx,
        meta.data  = meta.data,
        LRM.db     = lrm.db,
        mode       = "spatial",
        sample.col = "sample.id",
        k          = 3L,
        n.cores    = 1L,
        verbose    = FALSE
      )
    )
    expect_true(is_NICHESObject(obj))
  })

  test_that("sample.col attribute is set on the merged object", {
    obj <- suppressWarnings(
      create_NICHESObject(
        count.mtx  = count.mtx,
        meta.data  = meta.data,
        LRM.db     = lrm.db,
        mode       = "spatial",
        sample.col = "sample.id",
        k          = 3L,
        n.cores    = 1L,
        verbose    = FALSE
      )
    )
    expect_equal(attr(obj, "sample.col"), "sample.id")
  })

  test_that("mode attribute is preserved through sample.col path", {
    obj <- suppressWarnings(
      create_NICHESObject(
        count.mtx  = count.mtx,
        meta.data  = meta.data,
        LRM.db     = lrm.db,
        mode       = "spatial",
        sample.col = "sample.id",
        k          = 3L,
        n.cores    = 1L,
        verbose    = FALSE
      )
    )
    expect_equal(attr(obj, "mode"), "spatial")
  })

  test_that("sample.col path: all cells present in merged $cell.list", {
    obj <- suppressWarnings(
      create_NICHESObject(
        count.mtx  = count.mtx,
        meta.data  = meta.data,
        LRM.db     = lrm.db,
        mode       = "spatial",
        sample.col = "sample.id",
        k          = 3L,
        n.cores    = 1L,
        verbose    = FALSE
      )
    )
    expect_setequal(obj$cell.list, rownames(meta.data))
  })

  test_that("no cross-sample edges exist in the merged object", {
    obj <- suppressWarnings(
      create_NICHESObject(
        count.mtx  = count.mtx,
        meta.data  = meta.data,
        LRM.db     = lrm.db,
        mode       = "spatial",
        sample.col = "sample.id",
        k          = 3L,
        n.cores    = 1L,
        verbose    = FALSE
      )
    )

    # Determine which sample each cell belongs to
    barcode.to.sample <- setNames(meta.data$sample.id, rownames(meta.data))

    # Every edge's sending and receiving cell must be in the same sample
    edge.dt <- obj$edge.list
    sending.sample   <- barcode.to.sample[edge.dt$SendingCell]
    receiving.sample <- barcode.to.sample[edge.dt$ReceivingCell]
    expect_true(all(sending.sample == receiving.sample))
  })

  test_that("sample.col path: $cell.meta retains sample.id column", {
    obj <- suppressWarnings(
      create_NICHESObject(
        count.mtx  = count.mtx,
        meta.data  = meta.data,
        LRM.db     = lrm.db,
        mode       = "spatial",
        sample.col = "sample.id",
        k          = 3L,
        n.cores    = 1L,
        verbose    = FALSE
      )
    )
    expect_true("sample.id" %in% colnames(obj$cell.meta))
  })

  test_that("sample.col path: $edge.meta has Sending.sample.id column", {
    obj <- suppressWarnings(
      create_NICHESObject(
        count.mtx  = count.mtx,
        meta.data  = meta.data,
        LRM.db     = lrm.db,
        mode       = "spatial",
        sample.col = "sample.id",
        k          = 3L,
        n.cores    = 1L,
        verbose    = FALSE
      )
    )
    expect_true("Sending.sample.id" %in% colnames(obj$edge.meta))
    expect_true("Receiving.sample.id" %in% colnames(obj$edge.meta))
  })

  test_that("sample.col path: object without sample.col has no sample.col attr", {
    obj.no.sample <- suppressWarnings(
      create_NICHESObject(
        count.mtx  = count.mtx,
        meta.data  = meta.data,
        LRM.db     = lrm.db,
        mode       = "spatial",
        sample.col = NULL,
        k          = 3L,
        n.cores    = 1L,
        verbose    = FALSE
      )
    )
    expect_null(attr(obj.no.sample, "sample.col"))
  })

})   # end local()


# ===========================================================================
# SECTION 5: NEW — sample.col edge cases
# ===========================================================================

test_that("sample.col: single-cell sample emits warning, remaining samples merged", {
  n.genes <- 8L
  set.seed(7L)
  count.mtx <- Matrix::rsparsematrix(n.genes, 22L, density = 0.5,
                                     rand.x = function(n) rpois(n, 2))
  rownames(count.mtx) <- paste0("gene", seq_len(n.genes))
  colnames(count.mtx) <- paste0("cell", seq_len(22L))

  meta.data <- data.frame(
    row.names  = paste0("cell", seq_len(22L)),
    x          = runif(22L, 0, 100),
    y          = runif(22L, 0, 100),
    sample.id  = c("singleton", rep("normal", 21L)),   # 1 cell in sample "singleton"
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

  # The merged object should only contain cells from "normal"
  expect_true(is_NICHESObject(obj))
  expect_equal(length(obj$cell.list), 21L)
})
