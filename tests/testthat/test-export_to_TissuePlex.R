# ==============================================================================
# FILE:    tests/testthat/test-export_to_TissuePlex.R
# PURPOSE: Unit tests for export_to_TissuePlex().
#          Uses data(NICHESv2_example) for the full NICHESObject.
#          All blocks skip when arrow is not installed.
# ==============================================================================

# TissuePlex schema: the 16 required output columns in order.
.tissueplex.cols <- c(
  "edge", "sending_cell", "receiving_cell", "is_autocrine",
  "lrm", "lrm_id", "ligand", "receptor",
  "score", "score_norm",
  "x1", "y1", "x2", "y2",
  "sending_type", "receiving_type"
)

# -- Input validation ----------------------------------------------------------

test_that("export_to_TissuePlex: errors on non-NICHESObject input", {
  skip_if_not_installed("arrow")
  tmp <- tempfile(fileext = ".parquet")
  on.exit(unlink(tmp), add = TRUE)
  expect_error(
    export_to_TissuePlex(list(), tmp),
    regexp = "NICHESObject"
  )
})

test_that("export_to_TissuePlex: errors when arrow is not available", {
  data(NICHESv2_example)
  tmp <- tempfile(fileext = ".parquet")
  on.exit(unlink(tmp), add = TRUE)
  mockr <- function(...) FALSE
  local_mocked_bindings(requireNamespace = mockr, .package = "base")
  expect_error(
    export_to_TissuePlex(NICHESv2_example, tmp),
    regexp = "arrow"
  )
})

test_that("export_to_TissuePlex: errors on empty edge.list", {
  skip_if_not_installed("arrow")
  data(NICHESv2_example)
  obj <- NICHESv2_example
  obj$edge.list <- data.table::data.table(
    edge = character(0), SendingCell = character(0), ReceivingCell = character(0)
  )
  tmp <- tempfile(fileext = ".parquet")
  on.exit(unlink(tmp), add = TRUE)
  expect_error(
    export_to_TissuePlex(obj, tmp),
    regexp = "edge.list is empty"
  )
})

test_that("export_to_TissuePlex: errors on empty edge.meta", {
  skip_if_not_installed("arrow")
  data(NICHESv2_example)
  obj <- NICHESv2_example
  obj$edge.meta <- data.table::data.table()
  tmp <- tempfile(fileext = ".parquet")
  on.exit(unlink(tmp), add = TRUE)
  expect_error(
    export_to_TissuePlex(obj, tmp),
    regexp = "edge.meta is empty"
  )
})

test_that("export_to_TissuePlex: errors when coordinate columns are missing", {
  skip_if_not_installed("arrow")
  data(NICHESv2_example)
  tmp <- tempfile(fileext = ".parquet")
  on.exit(unlink(tmp), add = TRUE)
  expect_error(
    export_to_TissuePlex(NICHESv2_example, tmp, x.col = "no_such_col"),
    regexp = "Coordinate columns not found"
  )
})

# -- Basic export with defaults ------------------------------------------------

test_that("export_to_TissuePlex: returns output.path invisibly", {
  skip_if_not_installed("arrow")
  data(NICHESv2_example)
  tmp <- tempfile(fileext = ".parquet")
  on.exit(unlink(tmp), add = TRUE)
  result <- expect_invisible(
    export_to_TissuePlex(NICHESv2_example, tmp)
  )
  expect_equal(result, tmp)
})

test_that("export_to_TissuePlex: creates a readable parquet file", {
  skip_if_not_installed("arrow")
  data(NICHESv2_example)
  tmp <- tempfile(fileext = ".parquet")
  on.exit(unlink(tmp), add = TRUE)
  export_to_TissuePlex(NICHESv2_example, tmp)
  expect_true(file.exists(tmp))
  out <- arrow::read_parquet(tmp)
  expect_true(is.data.frame(out))
  expect_gt(nrow(out), 0L)
})

test_that("export_to_TissuePlex: output has exactly the 16 TissuePlex columns", {
  skip_if_not_installed("arrow")
  data(NICHESv2_example)
  tmp <- tempfile(fileext = ".parquet")
  on.exit(unlink(tmp), add = TRUE)
  export_to_TissuePlex(NICHESv2_example, tmp)
  out <- arrow::read_parquet(tmp)
  expect_equal(names(out), .tissueplex.cols)
})

test_that("export_to_TissuePlex: every edge in edge.list is represented", {
  skip_if_not_installed("arrow")
  data(NICHESv2_example)
  tmp <- tempfile(fileext = ".parquet")
  on.exit(unlink(tmp), add = TRUE)
  export_to_TissuePlex(NICHESv2_example, tmp)
  out <- arrow::read_parquet(tmp)
  expect_true(
    all(NICHESv2_example$edge.list$edge %in% out$edge)
  )
})

# -- Column types --------------------------------------------------------------

test_that("export_to_TissuePlex: column types match TissuePlex spec", {
  skip_if_not_installed("arrow")
  data(NICHESv2_example)
  tmp <- tempfile(fileext = ".parquet")
  on.exit(unlink(tmp), add = TRUE)
  export_to_TissuePlex(NICHESv2_example, tmp)
  out <- data.table::as.data.table(arrow::read_parquet(tmp))
  expect_type(out$edge,           "character")
  expect_type(out$sending_cell,   "character")
  expect_type(out$receiving_cell, "character")
  expect_type(out$is_autocrine,   "logical")
  expect_type(out$lrm,            "character")
  expect_type(out$lrm_id,         "integer")
  expect_type(out$ligand,         "character")
  expect_type(out$receptor,       "character")
  expect_type(out$score,          "double")
  expect_type(out$score_norm,     "double")
  expect_type(out$x1,             "double")
  expect_type(out$y1,             "double")
  expect_type(out$x2,             "double")
  expect_type(out$y2,             "double")
  expect_type(out$sending_type,   "character")
  expect_type(out$receiving_type, "character")
})

# -- is_autocrine flag ---------------------------------------------------------

test_that("export_to_TissuePlex: is_autocrine is TRUE iff sender == receiver", {
  skip_if_not_installed("arrow")
  data(NICHESv2_example)
  tmp <- tempfile(fileext = ".parquet")
  on.exit(unlink(tmp), add = TRUE)
  export_to_TissuePlex(NICHESv2_example, tmp)
  out <- data.table::as.data.table(arrow::read_parquet(tmp))
  expect_equal(out$is_autocrine, out$sending_cell == out$receiving_cell)
})

# -- LRM splitting -------------------------------------------------------------

test_that("export_to_TissuePlex: ligand and receptor are correctly split from lrm", {
  skip_if_not_installed("arrow")
  data(NICHESv2_example)
  tmp <- tempfile(fileext = ".parquet")
  on.exit(unlink(tmp), add = TRUE)
  export_to_TissuePlex(NICHESv2_example, tmp)
  out <- data.table::as.data.table(arrow::read_parquet(tmp))
  scored <- out[!is.na(lrm)]
  expected <- data.table::tstrsplit(scored$lrm, "|", fixed = TRUE)
  expect_equal(scored$ligand,   expected[[1L]])
  expect_equal(scored$receptor, expected[[2L]])
})

test_that("export_to_TissuePlex: lrm_id is a dense 1..N integer index", {
  skip_if_not_installed("arrow")
  data(NICHESv2_example)
  tmp <- tempfile(fileext = ".parquet")
  on.exit(unlink(tmp), add = TRUE)
  export_to_TissuePlex(NICHESv2_example, tmp)
  out <- data.table::as.data.table(arrow::read_parquet(tmp))
  valid.ids <- out[!is.na(lrm_id), sort(unique(lrm_id))]
  expect_equal(valid.ids, seq_along(valid.ids))
})

# -- Unscored edges get placeholder rows ---------------------------------------

test_that("export_to_TissuePlex: unscored edges have NA for LRM fields", {
  skip_if_not_installed("arrow")
  data(NICHESv2_example)
  tmp <- tempfile(fileext = ".parquet")
  on.exit(unlink(tmp), add = TRUE)
  export_to_TissuePlex(NICHESv2_example, tmp)
  out <- data.table::as.data.table(arrow::read_parquet(tmp))
  scored.edges   <- unique(NICHESv2_example$edge.data$edge)
  unscored.rows  <- out[!edge %in% scored.edges]
  if (nrow(unscored.rows) > 0L) {
    expect_true(all(is.na(unscored.rows$lrm)))
    expect_true(all(is.na(unscored.rows$score)))
    expect_true(all(is.na(unscored.rows$score_norm)))
    expect_true(all(is.na(unscored.rows$lrm_id)))
    expect_true(all(is.na(unscored.rows$ligand)))
    expect_true(all(is.na(unscored.rows$receptor)))
  }
})

# -- celltype.col variations ---------------------------------------------------

test_that("export_to_TissuePlex: celltype columns populated when celltype.col matches", {
  skip_if_not_installed("arrow")
  data(NICHESv2_example)
  tmp <- tempfile(fileext = ".parquet")
  on.exit(unlink(tmp), add = TRUE)
  export_to_TissuePlex(NICHESv2_example, tmp, celltype.col = "celltype")
  out <- data.table::as.data.table(arrow::read_parquet(tmp))
  expect_false(all(is.na(out$sending_type)))
  expect_false(all(is.na(out$receiving_type)))
})

test_that("export_to_TissuePlex: celltype.col = NULL produces NA type columns", {
  skip_if_not_installed("arrow")
  data(NICHESv2_example)
  tmp <- tempfile(fileext = ".parquet")
  on.exit(unlink(tmp), add = TRUE)
  export_to_TissuePlex(NICHESv2_example, tmp, celltype.col = NULL)
  out <- data.table::as.data.table(arrow::read_parquet(tmp))
  expect_true(all(is.na(out$sending_type)))
  expect_true(all(is.na(out$receiving_type)))
})

test_that("export_to_TissuePlex: nonexistent celltype.col emits message and yields NA", {
  skip_if_not_installed("arrow")
  data(NICHESv2_example)
  tmp <- tempfile(fileext = ".parquet")
  on.exit(unlink(tmp), add = TRUE)
  expect_message(
    export_to_TissuePlex(NICHESv2_example, tmp, celltype.col = "no_such_col"),
    regexp = "not found"
  )
  out <- data.table::as.data.table(arrow::read_parquet(tmp))
  expect_true(all(is.na(out$sending_type)))
  expect_true(all(is.na(out$receiving_type)))
})

# -- Coordinates ---------------------------------------------------------------

test_that("export_to_TissuePlex: x1/y1/x2/y2 are non-NA for all rows", {
  skip_if_not_installed("arrow")
  data(NICHESv2_example)
  tmp <- tempfile(fileext = ".parquet")
  on.exit(unlink(tmp), add = TRUE)
  export_to_TissuePlex(NICHESv2_example, tmp)
  out <- data.table::as.data.table(arrow::read_parquet(tmp))
  expect_false(anyNA(out$x1))
  expect_false(anyNA(out$y1))
  expect_false(anyNA(out$x2))
  expect_false(anyNA(out$y2))
})

# -- n.threads parameter -------------------------------------------------------

test_that("export_to_TissuePlex: n.threads restores original thread count", {
  skip_if_not_installed("arrow")
  data(NICHESv2_example)
  tmp <- tempfile(fileext = ".parquet")
  on.exit(unlink(tmp), add = TRUE)
  before <- data.table::getDTthreads()
  export_to_TissuePlex(NICHESv2_example, tmp, n.threads = 1L)
  after <- data.table::getDTthreads()
  expect_equal(before, after)
})

# -- Output directory creation -------------------------------------------------

test_that("export_to_TissuePlex: creates parent directories automatically", {
  skip_if_not_installed("arrow")
  data(NICHESv2_example)
  nested <- file.path(tempdir(), "tp_test", "sub1", "sub2", "edges.parquet")
  on.exit(unlink(file.path(tempdir(), "tp_test"), recursive = TRUE), add = TRUE)
  export_to_TissuePlex(NICHESv2_example, nested)
  expect_true(file.exists(nested))
})

# -- Row count identity --------------------------------------------------------

test_that("export_to_TissuePlex: scored rows + unscored rows = total output", {
  skip_if_not_installed("arrow")
  data(NICHESv2_example)
  tmp <- tempfile(fileext = ".parquet")
  on.exit(unlink(tmp), add = TRUE)
  export_to_TissuePlex(NICHESv2_example, tmp)
  out <- data.table::as.data.table(arrow::read_parquet(tmp))
  scored.edges   <- unique(NICHESv2_example$edge.data$edge)
  n.scored.rows  <- nrow(out[edge %in% scored.edges])
  n.unscored.rows <- nrow(out[!edge %in% scored.edges])
  expect_equal(n.scored.rows + n.unscored.rows, nrow(out))
})
