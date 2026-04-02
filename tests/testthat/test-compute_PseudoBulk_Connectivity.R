# ==============================================================================
# FILE:    tests/testthat/test-compute_PseudoBulk_Connectivity.R
# ==============================================================================

# ── Output structure ───────────────────────────────────────────────────────────

test_that("compute_PseudoBulk_Connectivity: returns a data.table", {
  .setup_fixtures()
  pb <- compute_PseudoBulk_Connectivity(count.mtx, meta.data, LRM.db,
                                         cell.type.col = "celltype")
  expect_s3_class(pb, "data.table")
})

test_that("compute_PseudoBulk_Connectivity: NOT a NICHESObject", {
  .setup_fixtures()
  pb <- compute_PseudoBulk_Connectivity(count.mtx, meta.data, LRM.db,
                                         cell.type.col = "celltype")
  expect_false(is_NICHESObject(pb))
})

test_that("compute_PseudoBulk_Connectivity: has required columns", {
  .setup_fixtures()
  pb <- compute_PseudoBulk_Connectivity(count.mtx, meta.data, LRM.db,
                                         cell.type.col = "celltype")
  expect_true(all(c("SendingType", "ReceivingType", "LRM",
                    "ligand", "receptor", "score", "score.norm") %in% colnames(pb)))
})

test_that("compute_PseudoBulk_Connectivity: LRM uses '|' separator", {
  .setup_fixtures()
  pb <- compute_PseudoBulk_Connectivity(count.mtx, meta.data, LRM.db,
                                         cell.type.col = "celltype")
  expect_true(all(grepl("|", pb$LRM, fixed = TRUE)))
})

test_that("compute_PseudoBulk_Connectivity: SendingType and ReceivingType are type labels", {
  .setup_fixtures()
  pb          <- compute_PseudoBulk_Connectivity(count.mtx, meta.data, LRM.db,
                                                  cell.type.col = "celltype")
  type.labels <- unique(meta.data$celltype)
  expect_true(all(pb$SendingType   %in% type.labels))
  expect_true(all(pb$ReceivingType %in% type.labels))
})

test_that("compute_PseudoBulk_Connectivity: all scores are strictly positive", {
  .setup_fixtures()
  pb <- compute_PseudoBulk_Connectivity(count.mtx, meta.data, LRM.db,
                                         cell.type.col = "celltype")
  expect_true(all(pb$score > 0))
})

test_that("compute_PseudoBulk_Connectivity: score.norm sums to 1 per type-pair (prop mode)", {
  .setup_fixtures()
  pb  <- compute_PseudoBulk_Connectivity(count.mtx, meta.data, LRM.db,
                                          cell.type.col = "celltype",
                                          normalize.method = "prop")
  tol <- 1e-9
  pair.sums <- pb[, .(s = sum(score.norm)), by = .(SendingType, ReceivingType)]
  expect_true(all(abs(pair.sums$s - 1) < tol))
})

test_that("compute_PseudoBulk_Connectivity: score.norm equals score when normalize.method='none'", {
  .setup_fixtures()
  pb <- compute_PseudoBulk_Connectivity(count.mtx, meta.data, LRM.db,
                                         cell.type.col = "celltype",
                                         normalize.method = "none")
  expect_equal(pb$score, pb$score.norm)
})

test_that("compute_PseudoBulk_Connectivity: output contains valid type labels", {
  .setup_fixtures()
  pb          <- compute_PseudoBulk_Connectivity(count.mtx, meta.data, LRM.db,
                                                  cell.type.col = "celltype")
  type.labels <- sort(unique(meta.data$celltype))
  # Output must be non-empty and contain only valid type labels
  expect_gt(nrow(pb), 0L)
  expect_true(all(pb$SendingType   %in% type.labels))
  expect_true(all(pb$ReceivingType %in% type.labels))
  # At least some type pairs must be represented
  present.pairs <- unique(pb[, .(SendingType, ReceivingType)])
  expect_gt(nrow(present.pairs), 0L)
})

test_that("compute_PseudoBulk_Connectivity: is deterministic (no randomness)", {
  .setup_fixtures()
  pb1 <- compute_PseudoBulk_Connectivity(count.mtx, meta.data, LRM.db,
                                          cell.type.col = "celltype", verbose = FALSE)
  pb2 <- compute_PseudoBulk_Connectivity(count.mtx, meta.data, LRM.db,
                                          cell.type.col = "celltype", verbose = FALSE)
  expect_equal(pb1$score, pb2$score)
})

test_that("compute_PseudoBulk_Connectivity: product vs minimum give different scores", {
  .setup_fixtures()
  pb.prod <- compute_PseudoBulk_Connectivity(count.mtx, meta.data, LRM.db,
                                              cell.type.col = "celltype",
                                              method = "product", verbose = FALSE)
  pb.min  <- compute_PseudoBulk_Connectivity(count.mtx, meta.data, LRM.db,
                                              cell.type.col = "celltype",
                                              method = "minimum", verbose = FALSE)
  expect_false(identical(pb.prod$score, pb.min$score))
})

# ── Input validation ───────────────────────────────────────────────────────────

test_that("compute_PseudoBulk_Connectivity: errors when cell.type.col missing", {
  .setup_fixtures()
  expect_error(
    compute_PseudoBulk_Connectivity(count.mtx, meta.data, LRM.db,
                                     cell.type.col = "NOTACOL"),
    regexp = "not found"
  )
})

test_that("compute_PseudoBulk_Connectivity: errors on non-matrix count.mtx", {
  .setup_fixtures()
  expect_error(
    compute_PseudoBulk_Connectivity(as.data.frame(as.matrix(count.mtx)),
                                     meta.data, LRM.db,
                                     cell.type.col = "celltype"),
    regexp = "matrix or sparseMatrix"
  )
})

test_that("compute_PseudoBulk_Connectivity: errors when no shared cells", {
  .setup_fixtures()
  bad.meta <- meta.data
  rownames(bad.meta) <- paste0("WRONG_", rownames(bad.meta))
  expect_error(
    compute_PseudoBulk_Connectivity(count.mtx, bad.meta, LRM.db,
                                     cell.type.col = "celltype"),
    regexp = "No shared cells"
  )
})
