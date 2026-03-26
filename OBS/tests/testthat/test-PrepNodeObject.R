# Tests for PrepNodeObject function
# tests/testthat/test-PrepNodeObject.R

test_that("PrepNodeObject filters cells by min.cells.per.ident", {
  skip_if_no_seurat()

  # Create object with unbalanced cell types
  test_obj <- create_test_seurat(n_cells = 60, n_celltypes = 3)

  # Manually create imbalance
  # Keep only 2 cells of CellType3
  cells_to_keep <- c(
    which(test_obj$cell_type == "CellType1"),
    which(test_obj$cell_type == "CellType2"),
    head(which(test_obj$cell_type == "CellType3"), 2)
  )
  test_obj <- subset(test_obj, cells = cells_to_keep)

  # Filter with min.cells.per.ident = 5
  result <- PrepNodeObject(
    node.object = test_obj,
    assay = "RNA",
    min.cells.per.ident = 5,
    min.cells.per.gene = NULL
  )

  # CellType3 should be removed (only 2 cells)
  remaining_types <- unique(result$cell_type)
  expect_false("CellType3" %in% remaining_types)
  expect_true("CellType1" %in% remaining_types)
  expect_true("CellType2" %in% remaining_types)
})


test_that("PrepNodeObject filters genes by min.cells.per.gene", {
  skip_if_no_seurat()

  test_obj <- create_test_seurat(n_cells = 50, n_genes = 100)

  # Set some genes to zero expression in most cells
  counts <- Seurat::GetAssayData(test_obj, layer = "counts")
  counts[1:10, ] <- 0  # First 10 genes have no expression
  counts[11:20, 1:2] <- 1  # Next 10 genes expressed in only 2 cells

  test_obj <- Seurat::CreateSeuratObject(
    counts = counts,
    meta.data = test_obj@meta.data
  )
  Seurat::Idents(test_obj) <- test_obj$cell_type

  # Filter genes expressed in < 5 cells
  result <- PrepNodeObject(
    node.object = test_obj,
    assay = "RNA",
    min.cells.per.ident = NULL,
    min.cells.per.gene = 5
  )

  # First 20 genes should be removed
  remaining_genes <- rownames(result)
  expect_false(any(paste0("GENE", 1:20) %in% remaining_genes))
  expect_true(any(paste0("GENE", 21:100) %in% remaining_genes))
})


test_that("PrepNodeObject handles both filters simultaneously", {
  skip_if_no_seurat()

  test_obj <- create_test_seurat(n_cells = 60, n_genes = 100, n_celltypes = 3)

  # Create imbalanced cell types
  cells_to_keep <- c(
    which(test_obj$cell_type == "CellType1"),
    which(test_obj$cell_type == "CellType2"),
    head(which(test_obj$cell_type == "CellType3"), 2)
  )
  test_obj <- subset(test_obj, cells = cells_to_keep)

  # Set some genes to low expression
  counts <- Seurat::GetAssayData(test_obj, layer = "counts")
  counts[1:10, ] <- 0
  test_obj <- Seurat::CreateSeuratObject(
    counts = counts,
    meta.data = test_obj@meta.data
  )
  Seurat::Idents(test_obj) <- test_obj$cell_type

  result <- PrepNodeObject(
    node.object = test_obj,
    assay = "RNA",
    min.cells.per.ident = 5,
    min.cells.per.gene = 3
  )

  # Check both filters applied
  expect_false("CellType3" %in% unique(result$cell_type))
  expect_false(any(paste0("GENE", 1:10) %in% rownames(result)))
})


test_that("PrepNodeObject preserves object when no filtering needed", {
  skip_if_no_seurat()

  test_obj <- create_test_seurat(n_cells = 60, n_genes = 100)
  original_cells <- ncol(test_obj)
  original_genes <- nrow(test_obj)

  # No filtering (very permissive thresholds)
  result <- PrepNodeObject(
    node.object = test_obj,
    assay = "RNA",
    min.cells.per.ident = 1,
    min.cells.per.gene = 1
  )

  # Should have same dimensions
  expect_equal(ncol(result), original_cells)
  expect_equal(nrow(result), original_genes)
})


test_that("PrepNodeObject validates input", {
  skip_if_no_seurat()

  test_obj <- create_test_seurat()

  # Test with invalid assay
  expect_error(
    PrepNodeObject(
      node.object = test_obj,
      assay = "INVALID_ASSAY",
      min.cells.per.ident = NULL,
      min.cells.per.gene = NULL
    )
  )
})


test_that("PrepNodeObject handles NULL thresholds", {
  skip_if_no_seurat()

  test_obj <- create_test_seurat()

  # Both NULL should return object as-is
  result <- PrepNodeObject(
    node.object = test_obj,
    assay = "RNA",
    min.cells.per.ident = NULL,
    min.cells.per.gene = NULL
  )

  expect_equal(ncol(result), ncol(test_obj))
  expect_equal(nrow(result), nrow(test_obj))
})
