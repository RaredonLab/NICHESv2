# ==============================================================================
# FILE:    tests/testthat/test-create_NICHESObject.R
# ==============================================================================

# ── Spatial mode (fixture-based) ──────────────────────────────────────────────

test_that("create_NICHESObject: returns a NICHESObject", {
  .setup_fixtures()
  expect_s3_class(niches.obj, "NICHESObject")
  expect_true(is_NICHESObject(niches.obj))
})

test_that("create_NICHESObject: all 13 slots present", {
  .setup_fixtures()
  required.slots <- c(
    "cell.list", "cell.data", "cell.meta",
    "edge.list", "edge.data", "edge.meta",
    "neighborhood.cell.list", "neighborhood.edge.list", "neighborhood.meta",
    "gene.list", "gene.meta", "LRM.list", "LRM.meta"
  )
  expect_true(all(required.slots %in% names(niches.obj)))
})

test_that("create_NICHESObject: neighborhood slots are NULL", {
  .setup_fixtures()
  expect_null(niches.obj$neighborhood.cell.list)
  expect_null(niches.obj$neighborhood.edge.list)
  expect_null(niches.obj$neighborhood.meta)
})

test_that("create_NICHESObject: cell.list matches count.mtx columns", {
  .setup_fixtures()
  expect_equal(sort(niches.obj$cell.list), sort(colnames(count.mtx)))
})

test_that("create_NICHESObject: gene.list matches count.mtx rows", {
  .setup_fixtures()
  expect_equal(niches.obj$gene.list, rownames(count.mtx))
})

test_that("create_NICHESObject: cell.data column schema is correct", {
  .setup_fixtures()
  expect_true(all(c("cell", "gene", "score", "score.norm") %in% colnames(niches.obj$cell.data)))
})

test_that("create_NICHESObject: cell.data is sparse (no zero scores)", {
  .setup_fixtures()
  expect_true(all(niches.obj$cell.data$score > 0))
})

test_that("create_NICHESObject: cell.data score.norm sums to 1 per cell", {
  .setup_fixtures()
  tol <- 1e-9
  dt  <- niches.obj$cell.data
  cell.sums <- dt[, .(s = sum(score.norm)), by = cell]
  expect_true(all(abs(cell.sums$s - 1) < tol))
})

test_that("create_NICHESObject: edge.list column schema is correct", {
  .setup_fixtures()
  expect_true(all(c("edge", "SendingCell", "ReceivingCell") %in% colnames(niches.obj$edge.list)))
})

test_that("create_NICHESObject: edge column is 'SendingCell|ReceivingCell'", {
  .setup_fixtures()
  el <- niches.obj$edge.list
  reconstructed <- paste(el$SendingCell, el$ReceivingCell, sep = "|")
  expect_equal(el$edge, reconstructed)
})

test_that("create_NICHESObject: edge.data column schema is correct", {
  .setup_fixtures()
  expect_true(all(c("edge", "LRM", "score", "score.norm") %in% colnames(niches.obj$edge.data)))
})

test_that("create_NICHESObject: edge.data score.norm sums to 1 per edge (prop mode)", {
  .setup_fixtures()
  tol  <- 1e-9
  ed   <- niches.obj$edge.data
  edge.sums <- ed[, .(s = sum(score.norm)), by = edge]
  expect_true(all(abs(edge.sums$s - 1) < tol))
})

test_that("create_NICHESObject: LRM.list uses '|' separator", {
  .setup_fixtures()
  expect_true(all(grepl("|", niches.obj$LRM.list, fixed = TRUE)))
})

test_that("create_NICHESObject: LRM.meta has required columns", {
  .setup_fixtures()
  expect_true(all(c("LRM", "ligand", "receptor") %in% colnames(niches.obj$LRM.meta)))
})

test_that("create_NICHESObject: edge.meta has Sending/Receiving/Joint cols for celltype", {
  .setup_fixtures()
  em.cols <- colnames(niches.obj$edge.meta)
  expect_true("Sending.celltype"  %in% em.cols)
  expect_true("Receiving.celltype" %in% em.cols)
  expect_true("Joint.celltype"    %in% em.cols)
})

test_that("create_NICHESObject: Joint.celltype is 'SendingType|ReceivingType'", {
  .setup_fixtures()
  em <- niches.obj$edge.meta
  expected <- paste(em$Sending.celltype, em$Receiving.celltype, sep = "|")
  expect_equal(em$Joint.celltype, expected)
})

test_that("create_NICHESObject: mode attribute is 'spatial' on spatial fixture", {
  .setup_fixtures()
  expect_equal(attr(niches.obj, "mode"), "spatial")
})

test_that("create_NICHESObject: errors when x/y missing in spatial mode", {
  .setup_fixtures()
  bad.meta <- meta.data
  bad.meta$x <- NULL
  expect_error(
    create_NICHESObject(count.mtx, bad.meta, LRM.db),
    regexp = "requires columns 'x' and 'y'"
  )
})

test_that("create_NICHESObject: [.NICHESObject subsetting preserves class", {
  .setup_fixtures()
  sub.cells <- niches.obj$cell.list[1:20]
  sub.obj   <- niches.obj[sub.cells]
  expect_s3_class(sub.obj, "NICHESObject")
  expect_equal(length(sub.obj$cell.list), 20L)
})

test_that("create_NICHESObject: dim() returns named vector", {
  .setup_fixtures()
  d <- dim(niches.obj)
  expect_equal(names(d), c("cells", "genes", "edges", "LRMs"))
  expect_equal(d["cells"], c(cells = length(niches.obj$cell.list)))
})

# ── Pseudobulk mode ───────────────────────────────────────────────────────────

test_that("create_NICHESObject pseudobulk: returns a NICHESObject", {
  .setup_fixtures()
  pb.obj <- create_NICHESObject(
    count.mtx     = count.mtx,
    meta.data     = meta.data,
    LRM.db        = LRM.db,
    mode          = "pseudobulk",
    cell.type.col = "celltype"
  )
  expect_s3_class(pb.obj, "NICHESObject")
  expect_true(is_NICHESObject(pb.obj))
})

test_that("create_NICHESObject pseudobulk: mode attribute is 'pseudobulk'", {
  .setup_fixtures()
  pb.obj <- create_NICHESObject(count.mtx, meta.data, LRM.db,
                                 mode = "pseudobulk", cell.type.col = "celltype")
  expect_equal(attr(pb.obj, "mode"), "pseudobulk")
})

test_that("create_NICHESObject pseudobulk: $cell.list holds original barcodes", {
  .setup_fixtures()
  pb.obj <- create_NICHESObject(count.mtx, meta.data, LRM.db,
                                 mode = "pseudobulk", cell.type.col = "celltype")
  expect_equal(sort(pb.obj$cell.list), sort(colnames(count.mtx)))
})

test_that("create_NICHESObject pseudobulk: edge nodes are type labels not barcodes", {
  .setup_fixtures()
  pb.obj    <- create_NICHESObject(count.mtx, meta.data, LRM.db,
                                    mode = "pseudobulk", cell.type.col = "celltype")
  type.labels <- sort(unique(meta.data$celltype))
  # SendingCell and ReceivingCell values must all be type labels
  expect_true(all(pb.obj$edge.list$SendingCell   %in% type.labels))
  expect_true(all(pb.obj$edge.list$ReceivingCell %in% type.labels))
})

test_that("create_NICHESObject pseudobulk: T² edges are produced", {
  .setup_fixtures()
  pb.obj    <- create_NICHESObject(count.mtx, meta.data, LRM.db,
                                    mode = "pseudobulk", cell.type.col = "celltype")
  n.types   <- length(unique(meta.data$celltype))
  expect_equal(nrow(pb.obj$edge.list), n.types^2)
})

test_that("create_NICHESObject pseudobulk: autocrine type edges present", {
  .setup_fixtures()
  pb.obj <- create_NICHESObject(count.mtx, meta.data, LRM.db,
                                 mode = "pseudobulk", cell.type.col = "celltype")
  # Every type should have a self-loop
  autocrine <- pb.obj$edge.list[pb.obj$edge.list$SendingCell ==
                                  pb.obj$edge.list$ReceivingCell, ]
  n.types   <- length(unique(meta.data$celltype))
  expect_equal(nrow(autocrine), n.types)
})

test_that("create_NICHESObject pseudobulk: errors when cell.type.col missing", {
  .setup_fixtures()
  expect_error(
    create_NICHESObject(count.mtx, meta.data, LRM.db, mode = "pseudobulk"),
    regexp = "cell.type.col"
  )
})

# ── Sampled mode ──────────────────────────────────────────────────────────────

test_that("create_NICHESObject sampled: returns a NICHESObject", {
  .setup_fixtures()
  samp.obj <- create_NICHESObject(
    count.mtx     = count.mtx,
    meta.data     = meta.data,
    LRM.db        = LRM.db,
    mode          = "sampled",
    cell.type.col = "celltype",
    seed          = 42L
  )
  expect_s3_class(samp.obj, "NICHESObject")
  expect_true(is_NICHESObject(samp.obj))
})

test_that("create_NICHESObject sampled: mode attribute is 'sampled'", {
  .setup_fixtures()
  samp.obj <- create_NICHESObject(count.mtx, meta.data, LRM.db,
                                   mode = "sampled", cell.type.col = "celltype",
                                   seed = 42L)
  expect_equal(attr(samp.obj, "mode"), "sampled")
})

test_that("create_NICHESObject sampled: $cell.list holds original barcodes", {
  .setup_fixtures()
  samp.obj <- create_NICHESObject(count.mtx, meta.data, LRM.db,
                                   mode = "sampled", cell.type.col = "celltype",
                                   seed = 42L)
  expect_equal(sort(samp.obj$cell.list), sort(colnames(count.mtx)))
})

test_that("create_NICHESObject sampled: edge nodes are cell barcodes", {
  .setup_fixtures()
  samp.obj <- create_NICHESObject(count.mtx, meta.data, LRM.db,
                                   mode = "sampled", cell.type.col = "celltype",
                                   seed = 42L)
  all.barcodes <- colnames(count.mtx)
  expect_true(all(samp.obj$edge.list$SendingCell   %in% all.barcodes))
  expect_true(all(samp.obj$edge.list$ReceivingCell %in% all.barcodes))
})

test_that("create_NICHESObject sampled: errors when cell.type.col missing", {
  .setup_fixtures()
  expect_error(
    create_NICHESObject(count.mtx, meta.data, LRM.db, mode = "sampled"),
    regexp = "cell.type.col"
  )
})
