# ==============================================================================
# FILE:    tests/testthat/test-merge_NICHESObjects.R
# PURPOSE: Unit tests for merge_NICHESObjects().
#          All NICHESObjects are synthesised inline via .make_minimal_NICHESObject()
#          — no shared fixtures needed.
#          obj.A/B/C are defined at file level as deterministic read-only
#          objects; no block mutates them, so this does not violate the
#          block-isolation requirement.
# ==============================================================================

# ── File-level helper ──────────────────────────────────────────────────────────

.make_minimal_NICHESObject <- function(barcodes,
                                       mode       = "spatial",
                                       gene.names = c("gene1", "gene2",
                                                      "gene3", "gene4"),
                                       sample.id  = NULL) {
  n.cells <- length(barcodes)
  n.genes <- length(gene.names)
  
  cell.meta.cols <- list(
    cell = barcodes,
    x    = seq_len(n.cells) * 10.0,
    y    = seq_len(n.cells) *  5.0
  )
  if (!is.null(sample.id))
    cell.meta.cols[["sample.id"]] <- rep(sample.id, n.cells)
  
  cell.meta <- do.call(data.table::data.table, cell.meta.cols)
  data.table::setkeyv(cell.meta, "cell")
  
  cell.data <- data.table::data.table(
    cell       = barcodes,
    gene       = gene.names[1L],
    score      = rpois(n.cells, lambda = 3),
    score.norm = 1.0
  )
  cell.data <- cell.data[score > 0]
  data.table::setkeyv(cell.data, c("cell", "gene"))
  
  edges <- paste(barcodes, barcodes, sep = "|")
  
  edge.list <- data.table::data.table(
    edge          = edges,
    SendingCell   = barcodes,
    ReceivingCell = barcodes
  )
  data.table::setkeyv(edge.list, "edge")
  
  edge.data <- data.table::data.table(
    edge       = edges,
    LRM        = "gene1|gene2",
    score      = rpois(n.cells, lambda = 2),
    score.norm = 1.0
  )
  edge.data <- edge.data[score > 0]
  data.table::setkeyv(edge.data, c("edge", "LRM"))
  
  edge.meta <- data.table::data.table(
    edge          = edges,
    SendingCell   = barcodes,
    ReceivingCell = barcodes,
    Sending.x     = seq_len(n.cells) * 10.0,
    Receiving.x   = seq_len(n.cells) * 10.0,
    Joint.x       = NA_real_
  )
  data.table::setkeyv(edge.meta, "edge")
  
  gene.meta <- data.table::data.table(gene = gene.names)
  data.table::setkeyv(gene.meta, "gene")
  
  LRM.meta <- data.table::data.table(
    LRM      = "gene1|gene2",
    ligand   = "gene1",
    receptor = "gene2"
  )
  data.table::setkeyv(LRM.meta, "LRM")
  
  obj <- NICHESv2:::.new_NICHESObject(
    cell.list              = barcodes,
    cell.data              = cell.data,
    cell.meta              = cell.meta,
    edge.list              = edge.list,
    edge.data              = edge.data,
    edge.meta              = edge.meta,
    neighborhood.cell.list = NULL,
    neighborhood.edge.list = NULL,
    neighborhood.meta      = NULL,
    gene.list              = gene.names,
    gene.meta              = gene.meta,
    LRM.list               = "gene1|gene2",
    LRM.meta               = LRM.meta,
    segmentation           = NULL
  )
  attr(obj, "mode") <- mode
  obj
}

# File-level objects: deterministic, read-only, never mutated by any block.
set.seed(42L)
obj.A <- .make_minimal_NICHESObject(paste0("s1_cell", 1:5), sample.id = "S1")
obj.B <- .make_minimal_NICHESObject(paste0("s2_cell", 1:4), sample.id = "S2")
obj.C <- .make_minimal_NICHESObject(paste0("s3_cell", 1:6), sample.id = "S3")
attr(obj.A, "sample.col") <- "sample.id"
attr(obj.B, "sample.col") <- "sample.id"
attr(obj.C, "sample.col") <- "sample.id"

# ── Input validation ───────────────────────────────────────────────────────────

test_that("merge_NICHESObjects: empty list is a hard error", {
  expect_error(merge_NICHESObjects(list()), regexp = "empty")
})

test_that("merge_NICHESObjects: non-NICHESObject in list is a hard error", {
  expect_error(
    merge_NICHESObjects(obj.A, list(1, 2, 3)),
    regexp = "NICHESObject"
  )
})

test_that("merge_NICHESObjects: mode mismatch is a hard error", {
  set.seed(1L)
  obj.sampled <- .make_minimal_NICHESObject(paste0("sc_cell", 1:3),
                                            mode = "nonspatial")
  expect_error(
    merge_NICHESObjects(obj.A, obj.sampled),
    regexp = "mode"
  )
})

test_that("merge_NICHESObjects: duplicate barcodes is a hard error", {
  set.seed(2L)
  obj.dup <- .make_minimal_NICHESObject(paste0("s1_cell", 1:5))
  attr(obj.dup, "mode") <- "spatial"
  expect_error(
    merge_NICHESObjects(obj.A, obj.dup),
    regexp = "[Dd]uplicate"
  )
})

# ── Single-object passthrough ──────────────────────────────────────────────────

test_that("merge_NICHESObjects: single-object list returns that object unchanged", {
  result <- suppressMessages(merge_NICHESObjects(list(obj.A)))
  expect_true(is_NICHESObject(result))
  expect_equal(result$cell.list, obj.A$cell.list)
})

# ── Calling conventions ────────────────────────────────────────────────────────

test_that("merge_NICHESObjects: positional args — merge_NICHESObjects(a, b, c)", {
  result <- suppressMessages(suppressWarnings(
    merge_NICHESObjects(obj.A, obj.B, obj.C)
  ))
  expect_true(is_NICHESObject(result))
})

test_that("merge_NICHESObjects: pre-built list — merge_NICHESObjects(list(a, b, c))", {
  result <- suppressMessages(suppressWarnings(
    merge_NICHESObjects(list(obj.A, obj.B, obj.C))
  ))
  expect_true(is_NICHESObject(result))
})

test_that("merge_NICHESObjects: named variable list — merge_NICHESObjects(obj.list)", {
  obj.list <- list(obj.A, obj.B, obj.C)
  result   <- suppressMessages(suppressWarnings(
    merge_NICHESObjects(obj.list)
  ))
  expect_true(is_NICHESObject(result))
})

test_that("merge_NICHESObjects: two-object list works", {
  result <- suppressMessages(suppressWarnings(
    merge_NICHESObjects(list(obj.A, obj.B))
  ))
  expect_true(is_NICHESObject(result))
})

# ── Cell slot correctness ──────────────────────────────────────────────────────

test_that("merge_NICHESObjects: $cell.list contains all barcodes from all objects", {
  result <- suppressMessages(suppressWarnings(
    merge_NICHESObjects(obj.A, obj.B, obj.C)
  ))
  expected <- c(obj.A$cell.list, obj.B$cell.list, obj.C$cell.list)
  expect_setequal(result$cell.list, expected)
})

test_that("merge_NICHESObjects: $cell.list length equals sum of individual lengths", {
  result <- suppressMessages(suppressWarnings(
    merge_NICHESObjects(obj.A, obj.B, obj.C)
  ))
  expect_equal(
    length(result$cell.list),
    length(obj.A$cell.list) + length(obj.B$cell.list) + length(obj.C$cell.list)
  )
})

test_that("merge_NICHESObjects: $cell.meta has one row per cell", {
  result <- suppressMessages(suppressWarnings(
    merge_NICHESObjects(obj.A, obj.B)
  ))
  expect_equal(
    nrow(result$cell.meta),
    length(obj.A$cell.list) + length(obj.B$cell.list)
  )
})

test_that("merge_NICHESObjects: $cell.meta is keyed on 'cell'", {
  result <- suppressMessages(suppressWarnings(
    merge_NICHESObjects(obj.A, obj.B)
  ))
  expect_equal(data.table::key(result$cell.meta), "cell")
})

# ── Edge slot correctness ──────────────────────────────────────────────────────

test_that("merge_NICHESObjects: all edge endpoints are in $cell.list", {
  result <- suppressMessages(suppressWarnings(
    merge_NICHESObjects(obj.A, obj.B)
  ))
  expect_true(all(result$edge.list$SendingCell   %in% result$cell.list))
  expect_true(all(result$edge.list$ReceivingCell %in% result$cell.list))
})

test_that("merge_NICHESObjects: $edge.data is keyed on (edge, LRM)", {
  result <- suppressMessages(suppressWarnings(
    merge_NICHESObjects(obj.A, obj.B)
  ))
  expect_equal(data.table::key(result$edge.data), c("edge", "LRM"))
})

# ── Gene / LRM union ──────────────────────────────────────────────────────────

test_that("merge_NICHESObjects: $gene.list is the union of individual gene lists", {
  set.seed(3L)
  obj.B.extra <- .make_minimal_NICHESObject(
    paste0("extra_cell", 1:3),
    gene.names = c("gene3", "gene4", "gene5", "gene6")
  )
  result <- suppressMessages(suppressWarnings(
    merge_NICHESObjects(obj.A, obj.B.extra)
  ))
  expect_setequal(result$gene.list, union(obj.A$gene.list, obj.B.extra$gene.list))
})

test_that("merge_NICHESObjects: $LRM.list is the union of individual LRM lists", {
  result <- suppressMessages(suppressWarnings(
    merge_NICHESObjects(obj.A, obj.B)
  ))
  expect_setequal(result$LRM.list, union(obj.A$LRM.list, obj.B$LRM.list))
})

test_that("merge_NICHESObjects: $gene.meta has no duplicate gene rows", {
  result <- suppressMessages(suppressWarnings(
    merge_NICHESObjects(obj.A, obj.B)
  ))
  expect_equal(nrow(result$gene.meta), data.table::uniqueN(result$gene.meta$gene))
})

test_that("merge_NICHESObjects: $LRM.meta has no duplicate LRM rows", {
  result <- suppressMessages(suppressWarnings(
    merge_NICHESObjects(obj.A, obj.B)
  ))
  expect_equal(nrow(result$LRM.meta), data.table::uniqueN(result$LRM.meta$LRM))
})

# ── Attribute propagation ──────────────────────────────────────────────────────

test_that("merge_NICHESObjects: mode attribute is preserved", {
  result <- suppressMessages(suppressWarnings(
    merge_NICHESObjects(obj.A, obj.B)
  ))
  expect_equal(attr(result, "mode"), "spatial")
})

test_that("merge_NICHESObjects: sample.col attribute propagated when consistent", {
  result <- suppressMessages(suppressWarnings(
    merge_NICHESObjects(obj.A, obj.B, obj.C)
  ))
  expect_equal(attr(result, "sample.col"), "sample.id")
})

test_that("merge_NICHESObjects: inconsistent sample.col produces warning and no attr", {
  obj.diff <- obj.B
  attr(obj.diff, "sample.col") <- "different.col"
  expect_warning(
    result <- suppressMessages(merge_NICHESObjects(obj.A, obj.diff)),
    regexp = "sample.col"
  )
  expect_null(attr(result, "sample.col"))
})

test_that("merge_NICHESObjects: no sample.col attr when none of the objects have it", {
  set.seed(4L)
  obj.plain1 <- .make_minimal_NICHESObject(paste0("plain_cell",  1:3))
  obj.plain2 <- .make_minimal_NICHESObject(paste0("plain2_cell", 1:3))
  attr(obj.plain1, "mode") <- "spatial"
  attr(obj.plain2, "mode") <- "spatial"
  
  result <- suppressMessages(merge_NICHESObjects(obj.plain1, obj.plain2))
  expect_null(attr(result, "sample.col"))
})

# ── Segmentation handling ──────────────────────────────────────────────────────

test_that("merge_NICHESObjects: all-NULL segmentation stays NULL in merged object", {
  result <- suppressMessages(suppressWarnings(
    merge_NICHESObjects(obj.A, obj.B)
  ))
  expect_null(result$segmentation)
})

test_that("merge_NICHESObjects: partial segmentation warns and returns NULL", {
  skip_if_not_installed("sf")
  
  polys <- lapply(seq_along(obj.A$cell.list), function(i) {
    sf::st_polygon(list(matrix(
      c(i, i, i+1, i+1, i, 0, 1, 1, 0, 0),
      ncol = 2
    )))
  })
  seg.sf <- sf::st_sf(
    cell     = obj.A$cell.list,
    geometry = sf::st_sfc(polys)
  )
  obj.A.seg <- obj.A
  obj.A.seg$segmentation <- seg.sf
  
  expect_warning(
    result <- suppressMessages(merge_NICHESObjects(obj.A.seg, obj.B)),
    regexp = "segmentation"
  )
  expect_null(result$segmentation)
})