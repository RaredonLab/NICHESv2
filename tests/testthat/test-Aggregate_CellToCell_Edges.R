# ==============================================================================
# FILE:    tests/testthat/test-Aggregate_CellToCell_Edges.R
# PURPOSE: Unit tests for Aggregate_CellToCell_Edges().
#          Uses data(NICHESv2_example) for the full NICHESObject — loaded
#          inline inside each block to preserve block independence.
# ==============================================================================

test_that("Aggregate_CellToCell_Edges: returns list with $sending and $receiving", {
  data(NICHESv2_example)
  agg <- Aggregate_CellToCell_Edges(NICHESv2_example)
  expect_type(agg, "list")
  expect_true(all(c("sending", "receiving") %in% names(agg)))
})

test_that("Aggregate_CellToCell_Edges: full grid dimensions are n.cells x n.LRM", {
  data(NICHESv2_example)
  agg     <- Aggregate_CellToCell_Edges(NICHESv2_example)
  n.cells <- length(NICHESv2_example$cell.list)
  n.lrm   <- length(NICHESv2_example$LRM.list)
  expect_equal(nrow(agg$sending),   n.cells * n.lrm)
  expect_equal(nrow(agg$receiving), n.cells * n.lrm)
})

test_that("Aggregate_CellToCell_Edges: required columns present in both tables", {
  data(NICHESv2_example)
  agg  <- Aggregate_CellToCell_Edges(NICHESv2_example)
  cols <- c("cell", "LRM", "ligand", "receptor",
            "n.edges", "score.sum", "score.mean", "score.max", "score.stat")
  expect_true(all(cols %in% colnames(agg$sending)))
  expect_true(all(cols %in% colnames(agg$receiving)))
})

test_that("Aggregate_CellToCell_Edges: silent combinations are zero-filled not NA", {
  data(NICHESv2_example)
  agg <- Aggregate_CellToCell_Edges(NICHESv2_example)
  expect_false(anyNA(agg$sending$score.sum))
  expect_false(anyNA(agg$receiving$score.sum))
})

test_that("Aggregate_CellToCell_Edges: stat='mean' sets score.stat = score.mean", {
  data(NICHESv2_example)
  agg <- Aggregate_CellToCell_Edges(NICHESv2_example, stat = "mean")
  expect_equal(agg$sending$score.stat, agg$sending$score.mean)
})

test_that("Aggregate_CellToCell_Edges: stat='max' sets score.stat = score.max", {
  data(NICHESv2_example)
  agg <- Aggregate_CellToCell_Edges(NICHESv2_example, stat = "max")
  expect_equal(agg$sending$score.stat, agg$sending$score.max)
})

test_that("Aggregate_CellToCell_Edges: score.sum >= score.max for all rows", {
  data(NICHESv2_example)
  agg <- Aggregate_CellToCell_Edges(NICHESv2_example)
  expect_true(all(agg$sending$score.sum >= agg$sending$score.max))
})

test_that("Aggregate_CellToCell_Edges: silent rows have n.edges == 0", {
  data(NICHESv2_example)
  agg         <- Aggregate_CellToCell_Edges(NICHESv2_example, stat = "sum")
  silent.rows <- agg$sending[agg$sending$score.sum == 0, ]
  expect_true(all(silent.rows$n.edges == 0))
})

test_that("Aggregate_CellToCell_Edges: cell.meta columns are joined onto output", {
  data(NICHESv2_example)
  agg <- Aggregate_CellToCell_Edges(NICHESv2_example)
  expect_true("celltype" %in% colnames(agg$sending))
  expect_true("x"        %in% colnames(agg$sending))
  expect_true("y"        %in% colnames(agg$sending))
})

test_that("Aggregate_CellToCell_Edges: errors on non-NICHESObject input", {
  expect_error(
    Aggregate_CellToCell_Edges(list(a = 1)),
    regexp = "NICHESObject"
  )
})