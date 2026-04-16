# ==============================================================================
# FILE:    tests/testthat/test-NICHESObject_methods.R
# PURPOSE: Unit tests for S3 methods on NICHESObject:
#            print, summary, dim, length, [
#          Uses data(NICHESv2_example) — loaded inline per block.
# ==============================================================================

# ── print ──────────────────────────────────────────────────────────────────────

test_that("print.NICHESObject: returns the object invisibly", {
  data(NICHESv2_example)
  out <- withVisible(print(NICHESv2_example))
  expect_false(out$visible)
  expect_true(is_NICHESObject(out$value))
})

test_that("print.NICHESObject: emits output to console", {
  data(NICHESv2_example)
  expect_output(print(NICHESv2_example))
})

# ── summary ────────────────────────────────────────────────────────────────────

test_that("summary.NICHESObject: returns the object invisibly", {
  data(NICHESv2_example)
  out <- withVisible(summary(NICHESv2_example))
  expect_false(out$visible)
  expect_true(is_NICHESObject(out$value))
})

test_that("summary.NICHESObject: emits output to console", {
  data(NICHESv2_example)
  expect_output(summary(NICHESv2_example))
})

# ── dim ───────────────────────────────────────────────────────────────────────

test_that("dim.NICHESObject: returns a named integer vector of length 4", {
  data(NICHESv2_example)
  d <- dim(NICHESv2_example)
  expect_type(d, "integer")
  expect_length(d, 4L)
  expect_named(d, c("cells", "genes", "edges", "LRMs"))
})

test_that("dim.NICHESObject: cells dimension equals length of $cell.list", {
  data(NICHESv2_example)
  expect_equal(dim(NICHESv2_example)[["cells"]],
               length(NICHESv2_example$cell.list))
})

test_that("dim.NICHESObject: LRMs dimension equals length of $LRM.list", {
  data(NICHESv2_example)
  expect_equal(dim(NICHESv2_example)[["LRMs"]],
               length(NICHESv2_example$LRM.list))
})

# ── length ─────────────────────────────────────────────────────────────────────

test_that("length.NICHESObject: returns a single integer", {
  data(NICHESv2_example)
  n <- length(NICHESv2_example)
  expect_type(n, "integer")
  expect_length(n, 1L)
})

test_that("length.NICHESObject: equals number of cells", {
  data(NICHESv2_example)
  expect_equal(length(NICHESv2_example),
               length(NICHESv2_example$cell.list))
})

# ── [ subsetting ───────────────────────────────────────────────────────────────

test_that("[.NICHESObject: returns a NICHESObject", {
  data(NICHESv2_example)
  keep    <- NICHESv2_example$cell.list[1:10]
  sub.obj <- NICHESv2_example[keep]
  expect_true(is_NICHESObject(sub.obj))
})

test_that("[.NICHESObject: $cell.list contains exactly the requested cells", {
  data(NICHESv2_example)
  keep    <- NICHESv2_example$cell.list[1:10]
  sub.obj <- NICHESv2_example[keep]
  expect_setequal(sub.obj$cell.list, keep)
})

test_that("[.NICHESObject: $cell.meta has one row per retained cell", {
  data(NICHESv2_example)
  keep    <- NICHESv2_example$cell.list[1:10]
  sub.obj <- NICHESv2_example[keep]
  expect_equal(nrow(sub.obj$cell.meta), 10L)
  expect_true(all(sub.obj$cell.meta$cell %in% keep))
})

test_that("[.NICHESObject: both endpoints of every retained edge are in kept cells", {
  data(NICHESv2_example)
  keep    <- NICHESv2_example$cell.list[1:10]
  sub.obj <- NICHESv2_example[keep]
  expect_true(all(sub.obj$edge.list$SendingCell   %in% keep))
  expect_true(all(sub.obj$edge.list$ReceivingCell %in% keep))
})

test_that("[.NICHESObject: $edge.data references only retained edges", {
  data(NICHESv2_example)
  keep    <- NICHESv2_example$cell.list[1:10]
  sub.obj <- NICHESv2_example[keep]
  expect_true(all(sub.obj$edge.data$edge %in% sub.obj$edge.list$edge))
})

test_that("[.NICHESObject: $edge.meta references only retained edges", {
  data(NICHESv2_example)
  keep    <- NICHESv2_example$cell.list[1:10]
  sub.obj <- NICHESv2_example[keep]
  expect_true(all(sub.obj$edge.meta$edge %in% sub.obj$edge.list$edge))
})

test_that("[.NICHESObject: $cell.data references only retained cells", {
  data(NICHESv2_example)
  keep    <- NICHESv2_example$cell.list[1:10]
  sub.obj <- NICHESv2_example[keep]
  expect_true(all(sub.obj$cell.data$cell %in% keep))
})

test_that("[.NICHESObject: $gene.list and $LRM.list are unchanged by subsetting", {
  data(NICHESv2_example)
  keep    <- NICHESv2_example$cell.list[1:10]
  sub.obj <- NICHESv2_example[keep]
  expect_equal(sub.obj$gene.list, NICHESv2_example$gene.list)
  expect_equal(sub.obj$LRM.list,  NICHESv2_example$LRM.list)
})

test_that("[.NICHESObject: mode attribute is preserved after subsetting", {
  data(NICHESv2_example)
  keep    <- NICHESv2_example$cell.list[1:10]
  sub.obj <- NICHESv2_example[keep]
  expect_equal(attr(sub.obj, "mode"), attr(NICHESv2_example, "mode"))
})

test_that("[.NICHESObject: subsetting with all cells returns equivalent object", {
  data(NICHESv2_example)
  sub.obj <- NICHESv2_example[NICHESv2_example$cell.list]
  expect_equal(length(sub.obj), length(NICHESv2_example))
  expect_setequal(sub.obj$cell.list, NICHESv2_example$cell.list)
})