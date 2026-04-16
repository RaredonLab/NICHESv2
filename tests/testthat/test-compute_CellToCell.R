# ==============================================================================
# FILE:    tests/testthat/test-compute_CellToCell.R
# PURPOSE: Unit tests for compute_CellToCell().
#          All inputs synthesised inline — no shared fixtures needed.
# ==============================================================================

# ── File-level helpers ─────────────────────────────────────────────────────────

# 20 genes x 10 cells, all non-zero (density + rpois offset ensures no zeros
# after abs), so product scores are guaranteed to be positive.
make_ctc_count_mtx <- function() {
  set.seed(1L)
  m <- Matrix::rsparsematrix(20L, 10L, density = 0.6,
                             rand.x = function(n) rpois(n, 3) + 1)
  rownames(m) <- paste0("gene", seq_len(20L))
  colnames(m) <- paste0("cell", seq_len(10L))
  m
}

# 3 LRM pairs; pathway column exercises the "extra columns preserved" path.
make_ctc_LRM_db <- function() {
  data.frame(
    ligand   = c("gene1", "gene3", "gene5"),
    receptor = c("gene2", "gene4", "gene6"),
    pathway  = c("PathA", "PathB", "PathC"),
    stringsAsFactors = FALSE
  )
}

# Spatial edge list (carries x.from / y.from / x.to / y.to).
make_ctc_spatial_el <- function() {
  data.frame(
    from   = paste0("cell", 1:5),
    to     = paste0("cell", 6:10),
    x.from = seq_len(5L) * 1.0,
    y.from = seq_len(5L) * 1.0,
    x.to   = seq_len(5L) * 2.0,
    y.to   = seq_len(5L) * 2.0,
    stringsAsFactors = FALSE
  )
}

# Bare edge list (from / to only — non-spatial path).
make_ctc_bare_el <- function() {
  data.frame(
    from = paste0("cell", 1:5),
    to   = paste0("cell", 6:10),
    stringsAsFactors = FALSE
  )
}

# ── Output structure ───────────────────────────────────────────────────────────

test_that("compute_CellToCell: returns plain named list with edges and lr.meta", {
  elr <- compute_CellToCell(
    count.mtx = make_ctc_count_mtx(),
    edge.list = make_ctc_spatial_el(),
    LRM.db    = make_ctc_LRM_db(),
    n.cores   = 1L
  )
  expect_type(elr, "list")
  expect_identical(sort(names(elr)), c("edges", "lr.meta"))
})

test_that("compute_CellToCell: $edges has required columns", {
  elr <- compute_CellToCell(
    count.mtx = make_ctc_count_mtx(),
    edge.list = make_ctc_spatial_el(),
    LRM.db    = make_ctc_LRM_db(),
    n.cores   = 1L
  )
  expect_true(all(c("from", "to", "ligand", "receptor", "LRM", "score") %in%
                    colnames(elr$edges)))
})

test_that("compute_CellToCell: LRM column uses '|' separator", {
  elr <- compute_CellToCell(
    count.mtx = make_ctc_count_mtx(),
    edge.list = make_ctc_spatial_el(),
    LRM.db    = make_ctc_LRM_db(),
    n.cores   = 1L
  )
  expect_true(all(grepl("|", elr$edges$LRM, fixed = TRUE)))
  expect_false(any(grepl("\\.", elr$edges$LRM)))
})

test_that("compute_CellToCell: all returned scores are strictly positive", {
  # Zero-score edges are filtered out by the function; everything remaining
  # must be > 0.
  elr <- compute_CellToCell(
    count.mtx = make_ctc_count_mtx(),
    edge.list = make_ctc_spatial_el(),
    LRM.db    = make_ctc_LRM_db(),
    n.cores   = 1L
  )
  expect_true(all(elr$edges$score > 0))
})

# ── Exact score arithmetic ─────────────────────────────────────────────────────

test_that("compute_CellToCell: product method gives correct score", {
  mx <- matrix(c(2, 0, 1, 0,
                 0, 3, 0, 2,
                 1, 1, 0, 3), nrow = 3, byrow = TRUE,
               dimnames = list(c("GeneA", "GeneB", "GeneC"), paste0("C", 1:4)))
  mx     <- as(mx, "CsparseMatrix")
  lrm.df <- data.frame(ligand = "GeneA", receptor = "GeneC")
  el.df  <- data.frame(from = c("C1", "C2"), to = c("C3", "C4"))
  
  # C1->C3: GeneA[C1]=2, GeneC[C3]=0 -> 0 (filtered out)
  # C2->C4: GeneA[C2]=0, GeneC[C4]=3 -> 0 (filtered out)
  elr <- compute_CellToCell(mx, el.df, lrm.df, method = "product", n.cores = 1L)
  expect_equal(nrow(elr$edges), 0L)
  
  # C1->C2: GeneA[C1]=2, GeneC[C2]=1 -> product = 2
  elr2 <- compute_CellToCell(mx, data.frame(from = "C1", to = "C2"),
                             lrm.df, method = "product", n.cores = 1L)
  expect_equal(nrow(elr2$edges), 1L)
  expect_equal(elr2$edges$score, 2)
})

test_that("compute_CellToCell: minimum method gives correct score", {
  mx <- matrix(c(4, 3), nrow = 2,
               dimnames = list(c("GeneA", "GeneB"), "C1"))
  mx     <- as(mx, "CsparseMatrix")
  lrm.df <- data.frame(ligand = "GeneA", receptor = "GeneB")
  el.df  <- data.frame(from = "C1", to = "C1")
  
  elr <- compute_CellToCell(mx, el.df, lrm.df, method = "minimum", n.cores = 1L)
  expect_equal(elr$edges$score, 3)   # min(4, 3)
})

test_that("compute_CellToCell: complex ligand (product) scored correctly", {
  # score = gene1[from] * gene2[from] * gene3[to]
  mx <- matrix(c(2, 3, 5), nrow = 3,
               dimnames = list(c("GeneA", "GeneB", "GeneC"), "C1"))
  mx     <- as(mx, "CsparseMatrix")
  lrm.df <- data.frame(ligand = "GeneA_GeneB", receptor = "GeneC")
  el.df  <- data.frame(from = "C1", to = "C1")
  
  elr <- compute_CellToCell(mx, el.df, lrm.df, method = "product", n.cores = 1L)
  expect_equal(elr$edges$score, 2 * 3 * 5)   # 30
})

test_that("compute_CellToCell: complex ligand (minimum) scored correctly", {
  mx <- matrix(c(2, 3, 5), nrow = 3,
               dimnames = list(c("GeneA", "GeneB", "GeneC"), "C1"))
  mx     <- as(mx, "CsparseMatrix")
  lrm.df <- data.frame(ligand = "GeneA_GeneB", receptor = "GeneC")
  el.df  <- data.frame(from = "C1", to = "C1")
  
  elr <- compute_CellToCell(mx, el.df, lrm.df, method = "minimum", n.cores = 1L)
  expect_equal(elr$edges$score, 2)   # min(2, 3, 5)
})

# ── LRM filtering and extra columns ───────────────────────────────────────────

test_that("compute_CellToCell: drops LRM pairs with missing genes and emits message", {
  bad.lrm <- data.frame(
    ligand   = c("gene1", "NOTEXIST"),
    receptor = c("gene2", "gene3")
  )
  expect_message(
    compute_CellToCell(make_ctc_count_mtx(), make_ctc_spatial_el(),
                       bad.lrm, n.cores = 1L),
    regexp = "Dropping 1"
  )
})

test_that("compute_CellToCell: lr.meta preserves extra columns from LRM.db", {
  elr <- compute_CellToCell(
    count.mtx = make_ctc_count_mtx(),
    edge.list = make_ctc_spatial_el(),
    LRM.db    = make_ctc_LRM_db(),
    n.cores   = 1L
  )
  expect_true("pathway" %in% colnames(elr$lr.meta))
})

# ── Spatial vs non-spatial edge list ──────────────────────────────────────────

test_that("compute_CellToCell: spatial edgelist passes coordinates through to $edges", {
  elr <- compute_CellToCell(
    count.mtx = make_ctc_count_mtx(),
    edge.list = make_ctc_spatial_el(),
    LRM.db    = make_ctc_LRM_db(),
    n.cores   = 1L
  )
  expect_true("x.from" %in% colnames(elr$edges))
  expect_true("y.to"   %in% colnames(elr$edges))
})

test_that("compute_CellToCell: bare from/to edgelist accepted; no coord columns in $edges", {
  elr <- compute_CellToCell(
    count.mtx = make_ctc_count_mtx(),
    edge.list = make_ctc_bare_el(),
    LRM.db    = make_ctc_LRM_db(),
    n.cores   = 1L
  )
  expect_type(elr, "list")
  expect_identical(sort(names(elr)), c("edges", "lr.meta"))
  expect_false("x.from" %in% colnames(elr$edges))
})