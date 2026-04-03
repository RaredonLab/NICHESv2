## test-merge_NICHESObjects.R
## Tests for merge_NICHESObjects() and the underlying .merge_NICHESObject_list()

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

## Build a minimal NICHESObject directly via .new_NICHESObject()
## so these tests are independent of create_NICHESObject() correctness.
.make_minimal_NICHESObject <- function(barcodes,
                                       mode       = "spatial",
                                       gene.names = c("gene1", "gene2",
                                                      "gene3", "gene4"),
                                       sample.id  = NULL) {

  n.cells <- length(barcodes)
  n.genes <- length(gene.names)

  # $cell.list
  cell.list <- barcodes

  # $cell.meta — minimal: cell, x, y (spatial); optionally a sample column
  cell.meta.cols <- list(
    cell = barcodes,
    x    = seq_len(n.cells) * 10.0,
    y    = seq_len(n.cells) * 5.0
  )
  if (!is.null(sample.id)) {
    cell.meta.cols[["sample.id"]] <- rep(sample.id, n.cells)
  }
  cell.meta <- do.call(data.table, cell.meta.cols)
  setkeyv(cell.meta, "cell")

  # $cell.data — sparse: just first gene has expression for all cells
  cell.data <- data.table(
    cell       = barcodes,
    gene       = gene.names[1L],
    score      = rpois(n.cells, lambda = 3),
    score.norm = 1.0   # simplification
  )
  cell.data <- cell.data[score > 0]
  setkeyv(cell.data, c("cell", "gene"))

  # $edge.list — a few self-loops (autocrine) for simplicity
  edge.list <- data.table(
    edge          = paste(barcodes, barcodes, sep = "|"),
    SendingCell   = barcodes,
    ReceivingCell = barcodes
  )
  setkeyv(edge.list, "edge")

  # $edge.data
  edges <- paste(barcodes, barcodes, sep = "|")
  edge.data <- data.table(
    edge       = edges,
    LRM        = "gene1|gene2",
    score      = rpois(n.cells, lambda = 2),
    score.norm = 1.0
  )
  edge.data <- edge.data[score > 0]
  setkeyv(edge.data, c("edge", "LRM"))

  # $edge.meta
  edge.meta <- data.table(
    edge                 = edges,
    SendingCell          = barcodes,
    ReceivingCell        = barcodes,
    Sending.x            = seq_len(n.cells) * 10.0,
    Receiving.x          = seq_len(n.cells) * 10.0,
    Joint.x              = NA_real_
  )
  setkeyv(edge.meta, "edge")

  # $gene.list / $gene.meta
  gene.list <- gene.names
  gene.meta <- data.table(gene = gene.names)
  setkeyv(gene.meta, "gene")

  # $LRM.list / $LRM.meta
  LRM.list <- "gene1|gene2"
  LRM.meta <- data.table(LRM = "gene1|gene2", ligand = "gene1",
                          receptor = "gene2")
  setkeyv(LRM.meta, "LRM")

  obj <- .new_NICHESObject(
    cell.list              = cell.list,
    cell.data              = cell.data,
    cell.meta              = cell.meta,
    edge.list              = edge.list,
    edge.data              = edge.data,
    edge.meta              = edge.meta,
    neighborhood.cell.list = NULL,
    neighborhood.edge.list = NULL,
    neighborhood.meta      = NULL,
    gene.list              = gene.list,
    gene.meta              = gene.meta,
    LRM.list               = LRM.list,
    LRM.meta               = LRM.meta,
    segmentation           = NULL
  )

  attr(obj, "mode") <- mode
  obj
}

set.seed(42L)
obj.A <- .make_minimal_NICHESObject(paste0("s1_cell", 1:5),
                                    sample.id = "S1")
obj.B <- .make_minimal_NICHESObject(paste0("s2_cell", 1:4),
                                    sample.id = "S2")
obj.C <- .make_minimal_NICHESObject(paste0("s3_cell", 1:6),
                                    sample.id = "S3")
attr(obj.A, "sample.col") <- "sample.id"
attr(obj.B, "sample.col") <- "sample.id"
attr(obj.C, "sample.col") <- "sample.id"


# ===========================================================================
# SECTION 1: Input validation
# ===========================================================================

test_that("empty obj.list is a hard error", {
  expect_error(merge_NICHESObjects(list()), regexp = "empty")
})

test_that("non-NICHESObject in list is a hard error", {
  expect_error(
    merge_NICHESObjects(obj.A, list(1, 2, 3)),
    regexp = "NICHESObject"
  )
})

test_that("mode mismatch is a hard error", {
  obj.sampled <- .make_minimal_NICHESObject(paste0("sc_cell", 1:3),
                                            mode = "sampled")
  expect_error(
    merge_NICHESObjects(obj.A, obj.sampled),
    regexp = "mode"
  )
})

test_that("duplicate barcodes is a hard error", {
  obj.dup <- .make_minimal_NICHESObject(paste0("s1_cell", 1:5))  # same as obj.A
  attr(obj.dup, "mode") <- "spatial"
  expect_error(
    merge_NICHESObjects(obj.A, obj.dup),
    regexp = "[Dd]uplicate"
  )
})


# ===========================================================================
# SECTION 2: Single-object passthrough
# ===========================================================================

test_that("single-object list returns that object unchanged", {
  result <- suppressMessages(merge_NICHESObjects(list(obj.A)))
  expect_true(is_NICHESObject(result))
  expect_equal(result$cell.list, obj.A$cell.list)
})


# ===========================================================================
# SECTION 3: Calling conventions — all three styles work
# ===========================================================================

test_that("direct positional arguments: merge_NICHESObjects(a, b, c)", {
  result <- suppressMessages(
    suppressWarnings(merge_NICHESObjects(obj.A, obj.B, obj.C))
  )
  expect_true(is_NICHESObject(result))
})

test_that("pre-built list: merge_NICHESObjects(list(a, b, c))", {
  result <- suppressMessages(
    suppressWarnings(merge_NICHESObjects(list(obj.A, obj.B, obj.C)))
  )
  expect_true(is_NICHESObject(result))
})

test_that("explicit list() call: merge_NICHESObjects(obj.list)", {
  obj.list <- list(obj.A, obj.B, obj.C)
  result   <- suppressMessages(
    suppressWarnings(merge_NICHESObjects(obj.list))
  )
  expect_true(is_NICHESObject(result))
})

test_that("two-object list: merge_NICHESObjects(list(a, b))", {
  result <- suppressMessages(
    suppressWarnings(merge_NICHESObjects(list(obj.A, obj.B)))
  )
  expect_true(is_NICHESObject(result))
})


# ===========================================================================
# SECTION 4: Cell slot correctness
# ===========================================================================

test_that("$cell.list contains all barcodes from all objects", {
  result <- suppressMessages(
    suppressWarnings(merge_NICHESObjects(obj.A, obj.B, obj.C))
  )
  expected <- c(obj.A$cell.list, obj.B$cell.list, obj.C$cell.list)
  expect_setequal(result$cell.list, expected)
})

test_that("$cell.list length equals sum of individual lengths", {
  result <- suppressMessages(
    suppressWarnings(merge_NICHESObjects(obj.A, obj.B, obj.C))
  )
  expected.n <- length(obj.A$cell.list) +
                length(obj.B$cell.list) +
                length(obj.C$cell.list)
  expect_equal(length(result$cell.list), expected.n)
})

test_that("$cell.meta has one row per cell", {
  result <- suppressMessages(
    suppressWarnings(merge_NICHESObjects(obj.A, obj.B))
  )
  expect_equal(nrow(result$cell.meta),
               length(obj.A$cell.list) + length(obj.B$cell.list))
})

test_that("$cell.meta is keyed on 'cell'", {
  result <- suppressMessages(
    suppressWarnings(merge_NICHESObjects(obj.A, obj.B))
  )
  expect_equal(key(result$cell.meta), "cell")
})


# ===========================================================================
# SECTION 5: Edge slot correctness
# ===========================================================================

test_that("$edge.list has no cross-object edges (edges are self-contained)", {
  result <- suppressMessages(
    suppressWarnings(merge_NICHESObjects(obj.A, obj.B))
  )
  # All sending cells must be in $cell.list
  expect_true(all(result$edge.list$SendingCell %in% result$cell.list))
  expect_true(all(result$edge.list$ReceivingCell %in% result$cell.list))
})

test_that("$edge.data is keyed on (edge, LRM)", {
  result <- suppressMessages(
    suppressWarnings(merge_NICHESObjects(obj.A, obj.B))
  )
  expect_equal(key(result$edge.data), c("edge", "LRM"))
})


# ===========================================================================
# SECTION 6: Gene / LRM union
# ===========================================================================

test_that("merged $gene.list is the union of individual gene lists", {
  # Give obj.B a different gene set (partial overlap)
  obj.B.extra <- .make_minimal_NICHESObject(
    paste0("extra_cell", 1:3),
    gene.names = c("gene3", "gene4", "gene5", "gene6")  # 2 new genes
  )
  result <- suppressMessages(
    suppressWarnings(merge_NICHESObjects(obj.A, obj.B.extra))
  )
  expected.genes <- union(obj.A$gene.list, obj.B.extra$gene.list)
  expect_setequal(result$gene.list, expected.genes)
})

test_that("merged $LRM.list is the union of individual LRM lists", {
  result <- suppressMessages(
    suppressWarnings(merge_NICHESObjects(obj.A, obj.B))
  )
  expected <- union(obj.A$LRM.list, obj.B$LRM.list)
  expect_setequal(result$LRM.list, expected)
})

test_that("$gene.meta has no duplicate gene rows", {
  result <- suppressMessages(
    suppressWarnings(merge_NICHESObjects(obj.A, obj.B))
  )
  expect_equal(nrow(result$gene.meta), uniqueN(result$gene.meta$gene))
})

test_that("$LRM.meta has no duplicate LRM rows", {
  result <- suppressMessages(
    suppressWarnings(merge_NICHESObjects(obj.A, obj.B))
  )
  expect_equal(nrow(result$LRM.meta), uniqueN(result$LRM.meta$LRM))
})


# ===========================================================================
# SECTION 7: Attribute propagation
# ===========================================================================

test_that("mode attribute is preserved in merged object", {
  result <- suppressMessages(
    suppressWarnings(merge_NICHESObjects(obj.A, obj.B))
  )
  expect_equal(attr(result, "mode"), "spatial")
})

test_that("sample.col attribute propagated when consistent across objects", {
  result <- suppressMessages(
    suppressWarnings(merge_NICHESObjects(obj.A, obj.B, obj.C))
  )
  expect_equal(attr(result, "sample.col"), "sample.id")
})

test_that("inconsistent sample.col attributes produce a warning and no attr", {
  obj.diff <- obj.B
  attr(obj.diff, "sample.col") <- "different.col"
  expect_warning(
    result <- suppressMessages(merge_NICHESObjects(obj.A, obj.diff)),
    regexp = "sample.col"
  )
  expect_null(attr(result, "sample.col"))
})

test_that("no sample.col attr when none of the objects have it", {
  obj.no.attr <- .make_minimal_NICHESObject(paste0("plain_cell", 1:3))
  attr(obj.no.attr, "mode") <- "spatial"
  obj.no.attr2 <- .make_minimal_NICHESObject(paste0("plain2_cell", 1:3))
  attr(obj.no.attr2, "mode") <- "spatial"

  result <- suppressMessages(merge_NICHESObjects(obj.no.attr, obj.no.attr2))
  expect_null(attr(result, "sample.col"))
})


# ===========================================================================
# SECTION 8: Segmentation handling
# ===========================================================================

test_that("segmentation NULL in all objects → NULL in merged", {
  result <- suppressMessages(
    suppressWarnings(merge_NICHESObjects(obj.A, obj.B))
  )
  expect_null(result$segmentation)
})

test_that("partial segmentation (some NULL) warns and returns NULL", {
  skip_if_not_installed("sf")

  # Create a minimal sf polygon for obj.A cells
  polys <- lapply(seq_along(obj.A$cell.list), function(i) {
    sf::st_polygon(list(matrix(c(i, i, i+1, i+1, i,
                                 0,  1, 1,   0,   0),
                               ncol = 2)))
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
