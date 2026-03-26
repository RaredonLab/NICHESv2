# ==============================================================================
# FILE:    tests/testthat/test-Aggregate_CellToCell_Edges.R
# ==============================================================================

test_that("Aggregate_CellToCell_Edges: returns list with $sending and $receiving", {
  .setup_fixtures()
  agg <- Aggregate_CellToCell_Edges(niches.obj)
  expect_type(agg, "list")
  expect_true(all(c("sending", "receiving") %in% names(agg)))
})

test_that("Aggregate_CellToCell_Edges: full grid dimensions are n.cells x n.LRM", {
  .setup_fixtures()
  agg     <- Aggregate_CellToCell_Edges(niches.obj)
  n.cells <- length(niches.obj$cell.list)
  n.lrm   <- length(niches.obj$LRM.list)
  expect_equal(nrow(agg$sending),   n.cells * n.lrm)
  expect_equal(nrow(agg$receiving), n.cells * n.lrm)
})

test_that("Aggregate_CellToCell_Edges: required columns present", {
  .setup_fixtures()
  agg  <- Aggregate_CellToCell_Edges(niches.obj)
  cols <- c("cell", "LRM", "ligand", "receptor",
            "n.edges", "score.sum", "score.mean", "score.max", "score.stat")
  expect_true(all(cols %in% colnames(agg$sending)))
  expect_true(all(cols %in% colnames(agg$receiving)))
})

test_that("Aggregate_CellToCell_Edges: silent combinations are zero-filled (not NA)", {
  .setup_fixtures()
  agg <- Aggregate_CellToCell_Edges(niches.obj)
  expect_false(anyNA(agg$sending$score.sum))
  expect_false(anyNA(agg$receiving$score.sum))
})

test_that("Aggregate_CellToCell_Edges: stat='mean' sets score.stat = score.mean", {
  .setup_fixtures()
  agg <- Aggregate_CellToCell_Edges(niches.obj, stat = "mean")
  expect_equal(agg$sending$score.stat, agg$sending$score.mean)
})

test_that("Aggregate_CellToCell_Edges: stat='max' sets score.stat = score.max", {
  .setup_fixtures()
  agg <- Aggregate_CellToCell_Edges(niches.obj, stat = "max")
  expect_equal(agg$sending$score.stat, agg$sending$score.max)
})

test_that("Aggregate_CellToCell_Edges: score.sum >= score.max (non-negative scores)", {
  .setup_fixtures()
  agg <- Aggregate_CellToCell_Edges(niches.obj)
  expect_true(all(agg$sending$score.sum >= agg$sending$score.max))
})

test_that("Aggregate_CellToCell_Edges: sending n.edges equals number of active neighbours", {
  .setup_fixtures()
  agg  <- Aggregate_CellToCell_Edges(niches.obj, stat = "sum")
  # For any cell with no signal edges for a given LRM, n.edges should be 0
  silent.rows <- agg$sending[agg$sending$score.sum == 0, ]
  expect_true(all(silent.rows$n.edges == 0))
})

test_that("Aggregate_CellToCell_Edges: cell.meta columns are joined", {
  .setup_fixtures()
  agg <- Aggregate_CellToCell_Edges(niches.obj)
  expect_true("celltype" %in% colnames(agg$sending))
  expect_true("x"        %in% colnames(agg$sending))
  expect_true("y"        %in% colnames(agg$sending))
})

test_that("Aggregate_CellToCell_Edges: errors on non-NICHESObject input", {
  expect_error(Aggregate_CellToCell_Edges(list(a = 1)), regexp = "NICHESObject")
})
