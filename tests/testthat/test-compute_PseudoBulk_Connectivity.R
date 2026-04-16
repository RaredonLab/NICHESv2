# ==============================================================================
# FILE:    tests/testthat/test-compute_PseudoBulk_Connectivity.R
# PURPOSE: Unit tests for compute_PseudoBulk_Connectivity().
#          All inputs synthesised inline — no shared fixtures needed.
# ==============================================================================

# ── File-level helpers ─────────────────────────────────────────────────────────

# 15 genes x 12 cells, high density + rpois offset guarantees non-zero
# pseudobulk averages for both cell types.
make_pb_count_mtx <- function() {
  set.seed(1L)
  m <- Matrix::rsparsematrix(15L, 12L, density = 0.8,
                             rand.x = function(n) rpois(n, 5) + 1)
  rownames(m) <- paste0("gene", seq_len(15L))
  colnames(m) <- paste0("cell", seq_len(12L))
  m
}

# 12 cells split evenly across 2 types; rownames match count matrix colnames.
make_pb_meta <- function() {
  data.frame(
    row.names = paste0("cell", seq_len(12L)),
    celltype  = rep(c("TypeA", "TypeB"), each = 6L),
    stringsAsFactors = FALSE
  )
}

# 3 LRM pairs drawn from the first 6 genes; pathway column for extra-col test.
make_pb_LRM_db <- function() {
  data.frame(
    ligand   = c("gene1", "gene3", "gene5"),
    receptor = c("gene2", "gene4", "gene6"),
    pathway  = c("PathA", "PathB", "PathC"),
    stringsAsFactors = FALSE
  )
}

# ── Output structure ───────────────────────────────────────────────────────────

test_that("compute_PseudoBulk_Connectivity: returns a data.table", {
  pb <- compute_PseudoBulk_Connectivity(
    make_pb_count_mtx(), make_pb_meta(), make_pb_LRM_db(),
    cell.type.col = "celltype", n.cores = 1L
  )
  expect_s3_class(pb, "data.table")
})

test_that("compute_PseudoBulk_Connectivity: NOT a NICHESObject", {
  pb <- compute_PseudoBulk_Connectivity(
    make_pb_count_mtx(), make_pb_meta(), make_pb_LRM_db(),
    cell.type.col = "celltype", n.cores = 1L
  )
  expect_false(is_NICHESObject(pb))
})

test_that("compute_PseudoBulk_Connectivity: has all required columns", {
  pb <- compute_PseudoBulk_Connectivity(
    make_pb_count_mtx(), make_pb_meta(), make_pb_LRM_db(),
    cell.type.col = "celltype", n.cores = 1L
  )
  expect_true(all(c("SendingType", "ReceivingType", "LRM",
                    "ligand", "receptor", "score", "score.norm") %in%
                    colnames(pb)))
})

test_that("compute_PseudoBulk_Connectivity: LRM column uses '|' separator", {
  pb <- compute_PseudoBulk_Connectivity(
    make_pb_count_mtx(), make_pb_meta(), make_pb_LRM_db(),
    cell.type.col = "celltype", n.cores = 1L
  )
  expect_true(all(grepl("|", pb$LRM, fixed = TRUE)))
})

test_that("compute_PseudoBulk_Connectivity: SendingType and ReceivingType are type labels", {
  meta        <- make_pb_meta()
  pb          <- compute_PseudoBulk_Connectivity(
    make_pb_count_mtx(), meta, make_pb_LRM_db(),
    cell.type.col = "celltype", n.cores = 1L
  )
  type.labels <- unique(meta$celltype)
  expect_true(all(pb$SendingType   %in% type.labels))
  expect_true(all(pb$ReceivingType %in% type.labels))
})

test_that("compute_PseudoBulk_Connectivity: all scores are strictly positive", {
  pb <- compute_PseudoBulk_Connectivity(
    make_pb_count_mtx(), make_pb_meta(), make_pb_LRM_db(),
    cell.type.col = "celltype", n.cores = 1L
  )
  expect_true(all(pb$score > 0))
})

# ── Normalisation ──────────────────────────────────────────────────────────────

test_that("compute_PseudoBulk_Connectivity: score.norm sums to 1 per type-pair (prop mode)", {
  pb  <- compute_PseudoBulk_Connectivity(
    make_pb_count_mtx(), make_pb_meta(), make_pb_LRM_db(),
    cell.type.col    = "celltype",
    normalize.method = "prop",
    n.cores          = 1L
  )
  pair.sums <- pb[, .(s = sum(score.norm)), by = .(SendingType, ReceivingType)]
  expect_true(all(abs(pair.sums$s - 1) < 1e-9))
})

test_that("compute_PseudoBulk_Connectivity: score.norm equals score when normalize.method='none'", {
  pb <- compute_PseudoBulk_Connectivity(
    make_pb_count_mtx(), make_pb_meta(), make_pb_LRM_db(),
    cell.type.col    = "celltype",
    normalize.method = "none",
    n.cores          = 1L
  )
  expect_equal(pb$score, pb$score.norm)
})

# ── Determinism and method contrast ───────────────────────────────────────────

test_that("compute_PseudoBulk_Connectivity: output is deterministic", {
  pb1 <- compute_PseudoBulk_Connectivity(
    make_pb_count_mtx(), make_pb_meta(), make_pb_LRM_db(),
    cell.type.col = "celltype", verbose = FALSE, n.cores = 1L
  )
  pb2 <- compute_PseudoBulk_Connectivity(
    make_pb_count_mtx(), make_pb_meta(), make_pb_LRM_db(),
    cell.type.col = "celltype", verbose = FALSE, n.cores = 1L
  )
  expect_equal(pb1$score, pb2$score)
})

test_that("compute_PseudoBulk_Connectivity: product vs minimum give different scores", {
  pb.prod <- compute_PseudoBulk_Connectivity(
    make_pb_count_mtx(), make_pb_meta(), make_pb_LRM_db(),
    cell.type.col = "celltype", method = "product",
    verbose = FALSE, n.cores = 1L
  )
  pb.min <- compute_PseudoBulk_Connectivity(
    make_pb_count_mtx(), make_pb_meta(), make_pb_LRM_db(),
    cell.type.col = "celltype", method = "minimum",
    verbose = FALSE, n.cores = 1L
  )
  expect_false(identical(pb.prod$score, pb.min$score))
})

test_that("compute_PseudoBulk_Connectivity: output contains valid type labels and is non-empty", {
  meta        <- make_pb_meta()
  pb          <- compute_PseudoBulk_Connectivity(
    make_pb_count_mtx(), meta, make_pb_LRM_db(),
    cell.type.col = "celltype", n.cores = 1L
  )
  type.labels <- sort(unique(meta$celltype))
  expect_gt(nrow(pb), 0L)
  expect_true(all(pb$SendingType   %in% type.labels))
  expect_true(all(pb$ReceivingType %in% type.labels))
  expect_gt(nrow(unique(pb[, .(SendingType, ReceivingType)])), 0L)
})

# ── Input validation ───────────────────────────────────────────────────────────

test_that("compute_PseudoBulk_Connectivity: errors when cell.type.col missing", {
  expect_error(
    compute_PseudoBulk_Connectivity(
      make_pb_count_mtx(), make_pb_meta(), make_pb_LRM_db(),
      cell.type.col = "NOTACOL", n.cores = 1L
    ),
    regexp = "not found"
  )
})

test_that("compute_PseudoBulk_Connectivity: errors on non-matrix count.mtx", {
  expect_error(
    compute_PseudoBulk_Connectivity(
      as.data.frame(as.matrix(make_pb_count_mtx())),
      make_pb_meta(), make_pb_LRM_db(),
      cell.type.col = "celltype", n.cores = 1L
    ),
    regexp = "matrix or sparseMatrix"
  )
})

test_that("compute_PseudoBulk_Connectivity: errors when no shared cells between matrix and meta", {
  bad.meta <- make_pb_meta()
  rownames(bad.meta) <- paste0("WRONG_", rownames(bad.meta))
  expect_error(
    compute_PseudoBulk_Connectivity(
      make_pb_count_mtx(), bad.meta, make_pb_LRM_db(),
      cell.type.col = "celltype", n.cores = 1L
    ),
    regexp = "No shared cells"
  )
})