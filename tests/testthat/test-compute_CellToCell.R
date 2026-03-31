# ==============================================================================
# FILE:    tests/testthat/test-compute_CellToCell.R
# ==============================================================================

test_that("compute_CellToCell: returns EdgeLRResult with correct slots", {
  .setup_fixtures()
  elr <- compute_CellToCell(
    count.mtx = count.mtx,
    edge.list = edge.list,
    LRM.db    = LRM.db
  )
  expect_s3_class(elr, "EdgeLRResult")
  expect_true(all(c("edges", "cell.meta", "lr.meta", "params") %in% names(elr)))
})

test_that("compute_CellToCell: $edges has required columns", {
  .setup_fixtures()
  elr <- compute_CellToCell(count.mtx, edge.list, LRM.db)
  expect_true(all(c("from", "to", "ligand", "receptor", "LRM", "score") %in% colnames(elr$edges)))
})

test_that("compute_CellToCell: LRM column uses '|' separator", {
  .setup_fixtures()
  elr <- compute_CellToCell(count.mtx, edge.list, LRM.db)
  expect_true(all(grepl("|", elr$edges$LRM, fixed = TRUE)))
  expect_false(any(grepl("\\.", elr$edges$LRM)))
})

test_that("compute_CellToCell: all scores are strictly positive", {
  .setup_fixtures()
  elr <- compute_CellToCell(count.mtx, edge.list, LRM.db)
  expect_true(all(elr$edges$score > 0))
})

test_that("compute_CellToCell: product method gives correct score", {
  mx <- matrix(c(2, 0, 1, 0,
                  0, 3, 0, 2,
                  1, 1, 0, 3), nrow = 3, byrow = TRUE,
               dimnames = list(c("GeneA", "GeneB", "GeneC"), paste0("C", 1:4)))
  mx <- as(mx, "CsparseMatrix")
  lrm.df <- data.frame(ligand = "GeneA", receptor = "GeneC")
  el.df  <- data.frame(from = c("C1", "C2"), to = c("C3", "C4"))
  # C1→C3: GeneA[C1]=2, GeneC[C3]=0 → product = 0 (filtered out)
  # C2→C4: GeneA[C2]=0, GeneC[C4]=3 → product = 0 (filtered out)
  elr <- compute_CellToCell(mx, el.df, lrm.df, method = "product")
  expect_equal(nrow(elr$edges), 0L)

  # Now a non-zero edge
  el2 <- data.frame(from = "C1", to = "C2")
  # C1→C2: GeneA[C1]=2, GeneC[C2]=1 → product = 2
  elr2 <- compute_CellToCell(mx, el2, lrm.df, method = "product")
  expect_equal(nrow(elr2$edges), 1L)
  expect_equal(elr2$edges$score, 2)
})

test_that("compute_CellToCell: minimum method gives correct score", {
  mx <- matrix(c(4, 3), nrow = 2, dimnames = list(c("GeneA", "GeneB"), "C1"))
  mx <- as(mx, "CsparseMatrix")
  lrm.df <- data.frame(ligand = "GeneA", receptor = "GeneB")
  el.df  <- data.frame(from = "C1", to = "C1")
  elr <- compute_CellToCell(mx, el.df, lrm.df, method = "minimum")
  expect_equal(elr$edges$score, 3)  # min(4, 3)
})

test_that("compute_CellToCell: complex ligand (product) scored correctly", {
  # Complex ligand: GeneA_GeneB → GeneC
  # score = expr(GeneA)[from] * expr(GeneB)[from] * expr(GeneC)[to]
  mx <- matrix(c(2, 3, 5), nrow = 3,
               dimnames = list(c("GeneA", "GeneB", "GeneC"), "C1"))
  mx <- as(mx, "CsparseMatrix")
  lrm.df <- data.frame(ligand = "GeneA_GeneB", receptor = "GeneC")
  el.df  <- data.frame(from = "C1", to = "C1")
  elr <- compute_CellToCell(mx, el.df, lrm.df, method = "product")
  expect_equal(elr$edges$score, 2 * 3 * 5)  # 30
})

test_that("compute_CellToCell: complex ligand (minimum) scored correctly", {
  mx <- matrix(c(2, 3, 5), nrow = 3,
               dimnames = list(c("GeneA", "GeneB", "GeneC"), "C1"))
  mx <- as(mx, "CsparseMatrix")
  lrm.df <- data.frame(ligand = "GeneA_GeneB", receptor = "GeneC")
  el.df  <- data.frame(from = "C1", to = "C1")
  elr <- compute_CellToCell(mx, el.df, lrm.df, method = "minimum")
  expect_equal(elr$edges$score, 2)  # min(2, 3, 5)
})

test_that("compute_CellToCell: drops LRM pairs with missing genes", {
  .setup_fixtures()
  bad.lrm <- data.frame(ligand = c("Gene1", "NOTEXIST"),
                         receptor = c("Gene2", "Gene3"))
  expect_message(
    elr <- compute_CellToCell(count.mtx, edge.list, bad.lrm),
    regexp = "Dropping 1"
  )
})

test_that("compute_CellToCell: lr.meta carries extra columns", {
  .setup_fixtures()
  elr <- compute_CellToCell(count.mtx, edge.list, LRM.db)
  expect_true("pathway" %in% colnames(elr$lr.meta))
})

test_that("compute_CellToCell: no coordinates message absent for spatial edgelist", {
  # Spatial edgelists carry coordinates — no message should fire
  .setup_fixtures()
  expect_no_message(
    compute_CellToCell(count.mtx, edge.list, LRM.db),
    message = "No coordinates"
  )
})

test_that("compute_CellToCell: accepts edgelist without coordinates (non-spatial)", {
  # A plain data.frame with from/to only (no x.from etc.) should work fine
  .setup_fixtures()
  bare.el <- data.frame(from = edge.list$from, to = edge.list$to)
  elr <- compute_CellToCell(count.mtx, bare.el, LRM.db)
  expect_s3_class(elr, "EdgeLRResult")
  # coord columns should be absent from $edges
  expect_false("x.from" %in% colnames(elr$edges))
})
