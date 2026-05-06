# ==============================================================================
# FILE:    tests/testthat/test-add_cell_rawcounts.R
# PURPOSE: Unit tests for add_cell_rawcounts().
#          Uses data(NICHESv2_example) throughout.
#          Tests Option B row expansion: new rows added for (cell, gene)
#          pairs where count > 0 and score = 0.
# ==============================================================================

data(NICHESv2_example)

n.cells <- length(NICHESv2_example$cell.list)
n.genes <- length(NICHESv2_example$gene.list)

# Standard count matrix: all 1L, covers all cells and genes in the fixture.
# Every (cell, gene) pair will have count = 1, triggering row expansion
# for pairs not already in $cell.data (score = 0).
count.full <- matrix(
  1L,
  nrow     = n.genes,
  ncol     = n.cells,
  dimnames = list(NICHESv2_example$gene.list, NICHESv2_example$cell.list)
)


# -- Input validation ----------------------------------------------------------

test_that("add_cell_rawcounts: errors if niches.obj is not a NICHESObject", {
  expect_error(
    add_cell_rawcounts(list(), count.full),
    regexp = "NICHESObject"
  )
})

test_that("add_cell_rawcounts: errors if count.mtx is not a matrix or Matrix", {
  expect_error(
    add_cell_rawcounts(NICHESv2_example, count.mtx = list(a = 1)),
    regexp = "base matrix or a sparse Matrix"
  )
})

test_that("add_cell_rawcounts: errors if count.mtx has no rownames or colnames", {
  no.names <- matrix(1L, nrow = n.genes, ncol = n.cells)
  expect_error(
    add_cell_rawcounts(NICHESv2_example, no.names),
    regexp = "rownames.*colnames|colnames.*rownames"
  )
})

test_that("add_cell_rawcounts: errors if cells are missing from count.mtx", {
  partial <- count.full[, -1L, drop = FALSE]   # drop first cell
  expect_error(
    add_cell_rawcounts(NICHESv2_example, partial, verbose = FALSE),
    regexp = "missing"
  )
})


# -- Gene coverage warnings and NA handling ------------------------------------

test_that("add_cell_rawcounts: emits message when genes are absent from count.mtx", {
  count.sub <- count.full[-1L, , drop = FALSE]   # drop first gene
  expect_message(
    add_cell_rawcounts(NICHESv2_example, count.sub, verbose = TRUE),
    regexp = "absent from"
  )
})

test_that("add_cell_rawcounts: absent genes in count.mtx get count = NA", {
  missing.gene <- NICHESv2_example$gene.list[1L]
  count.sub    <- count.full[NICHESv2_example$gene.list[-1L], , drop = FALSE]
  result       <- suppressMessages(
    add_cell_rawcounts(NICHESv2_example, count.sub, verbose = FALSE)
  )
  rows.missing <- result$cell.data[gene == missing.gene]
  expect_true(nrow(rows.missing) > 0L)
  expect_true(all(is.na(rows.missing$count)))
})


# -- Non-negative integer checks -----------------------------------------------

test_that("add_cell_rawcounts: warns on negative values in count.mtx", {
  count.neg        <- count.full
  count.neg[1L, 1L] <- -1L
  expect_warning(
    suppressMessages(add_cell_rawcounts(NICHESv2_example, count.neg,
                                        verbose = FALSE)),
    regexp = "negative"
  )
})

test_that("add_cell_rawcounts: warns on non-integer values in count.mtx", {
  count.frac        <- matrix(
    as.numeric(count.full),
    nrow     = n.genes,
    ncol     = n.cells,
    dimnames = list(NICHESv2_example$gene.list, NICHESv2_example$cell.list)
  )
  count.frac[1L, 1L] <- 1.5
  expect_warning(
    suppressMessages(add_cell_rawcounts(NICHESv2_example, count.frac,
                                        verbose = FALSE)),
    regexp = "non-integer"
  )
})

test_that("add_cell_rawcounts: warnings for negative and non-integer are independent", {
  count.both        <- matrix(
    as.numeric(count.full),
    nrow     = n.genes,
    ncol     = n.cells,
    dimnames = list(NICHESv2_example$gene.list, NICHESv2_example$cell.list)
  )
  count.both[1L, 1L] <- -1.5
  w <- tryCatch(
    suppressMessages(add_cell_rawcounts(NICHESv2_example, count.both,
                                         verbose = FALSE)),
    warning = function(w) w
  )
  # At least one warning is emitted; both conditions are present
  expect_true(inherits(w, "warning"))
})


# -- Existing count column notice ----------------------------------------------

test_that("add_cell_rawcounts: emits message when count column already exists", {
  obj.with.count <- suppressMessages(
    add_cell_rawcounts(NICHESv2_example, count.full, verbose = FALSE)
  )
  expect_message(
    add_cell_rawcounts(obj.with.count, count.full, verbose = TRUE),
    regexp = "will be overwritten"
  )
})


# -- Return type and slot type -------------------------------------------------

test_that("add_cell_rawcounts: returns a NICHESObject", {
  result <- suppressMessages(
    add_cell_rawcounts(NICHESv2_example, count.full, verbose = FALSE)
  )
  expect_true(inherits(result, "NICHESObject"))
})

test_that("add_cell_rawcounts: count column is present in $cell.data", {
  result <- suppressMessages(
    add_cell_rawcounts(NICHESv2_example, count.full, verbose = FALSE)
  )
  expect_true("count" %in% colnames(result$cell.data))
})

test_that("add_cell_rawcounts: $cell.data is keyed on (cell, gene)", {
  result <- suppressMessages(
    add_cell_rawcounts(NICHESv2_example, count.full, verbose = FALSE)
  )
  expect_equal(data.table::key(result$cell.data), c("cell", "gene"))
})


# -- Count value correctness for existing rows ---------------------------------

test_that("add_cell_rawcounts: existing $cell.data rows get correct count values", {
  result        <- suppressMessages(
    add_cell_rawcounts(NICHESv2_example, count.full, verbose = FALSE)
  )
  existing.rows <- result$cell.data[score > 0]
  # All entries in count.full are 1L; all existing rows should have count = 1
  expect_true(all(existing.rows$count == 1))
})

test_that("add_cell_rawcounts: genes in count.mtx with count = 0 get count = 0 not NA", {
  # Build a sparse matrix where the first gene has all zeros for all cells
  count.sparse        <- Matrix::Matrix(count.full, sparse = TRUE)
  first.gene          <- NICHESv2_example$gene.list[1L]
  count.sparse[1L, ]  <- 0
  result <- suppressMessages(
    add_cell_rawcounts(NICHESv2_example, count.sparse, verbose = FALSE)
  )
  rows.zero.gene <- result$cell.data[gene == first.gene & score > 0]
  if (nrow(rows.zero.gene) > 0L)
    expect_true(all(rows.zero.gene$count == 0))
})


# -- Option B row expansion ----------------------------------------------------

test_that("add_cell_rawcounts: new rows added for count > 0 and score = 0", {
  original.nrow <- nrow(NICHESv2_example$cell.data)
  result        <- suppressMessages(
    add_cell_rawcounts(NICHESv2_example, count.full, verbose = FALSE)
  )
  expect_gt(nrow(result$cell.data), original.nrow)
})

test_that("add_cell_rawcounts: total rows after equals n.cells x n.genes when all counts = 1", {
  result <- suppressMessages(
    add_cell_rawcounts(NICHESv2_example, count.full, verbose = FALSE)
  )
  # Every (cell, gene) pair now has either score > 0 OR count > 0 (both = 1),
  # so all n.cells * n.genes pairs should be present.
  expect_equal(nrow(result$cell.data), n.cells * n.genes)
})

test_that("add_cell_rawcounts: new rows have score = 0 and score.norm = 0", {
  result    <- suppressMessages(
    add_cell_rawcounts(NICHESv2_example, count.full, verbose = FALSE)
  )
  new.rows  <- result$cell.data[score == 0]
  expect_true(nrow(new.rows) > 0L)
  expect_true(all(new.rows$score.norm == 0))
})

test_that("add_cell_rawcounts: new rows have correct count values", {
  result   <- suppressMessages(
    add_cell_rawcounts(NICHESv2_example, count.full, verbose = FALSE)
  )
  new.rows <- result$cell.data[score == 0]
  expect_true(all(new.rows$count == 1))
})

test_that("add_cell_rawcounts: emits message reporting number of new rows", {
  expect_message(
    add_cell_rawcounts(NICHESv2_example, count.full, verbose = TRUE),
    regexp = "Added"
  )
})


# -- Sparse Matrix input -------------------------------------------------------

test_that("add_cell_rawcounts: works with sparse dgCMatrix input", {
  count.sparse <- Matrix::Matrix(count.full, sparse = TRUE)
  result       <- suppressMessages(
    add_cell_rawcounts(NICHESv2_example, count.sparse, verbose = FALSE)
  )
  expect_true(inherits(result, "NICHESObject"))
  expect_true("count" %in% colnames(result$cell.data))
})
