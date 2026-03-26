# Tests for RunCellToCell function
# tests/testthat/test-RunCellToCell.R

test_that("RunCellToCell creates valid Seurat output", {
  skip_if_no_seurat()

  # Arrange
  test_obj <- create_test_seurat()
  test_gt <- create_test_groundtruth()

  # Act
  result <- RunCellToCell(
    node.object = test_obj,
    ground.truth = test_gt,
    assay = "RNA",
    meta.data.to.map = c("cell_type"),
    output_format = "seurat"
  )

  # Assert structure
  expect_s4_class(result, "Seurat")

  # Assert dimensions
  expect_equal(nrow(result), 5)  # 5 mechanisms in test ground truth
  expect_gt(ncol(result), 0)     # At least some cell pairs

  # Assert metadata columns exist
  expect_true("SendingType" %in% colnames(result@meta.data))
  expect_true("ReceivingType" %in% colnames(result@meta.data))
  expect_true("VectorType" %in% colnames(result@meta.data))
  expect_true("cell_type.Sending" %in% colnames(result@meta.data))
  expect_true("cell_type.Receiving" %in% colnames(result@meta.data))
  expect_true("cell_type.Joint" %in% colnames(result@meta.data))
})


test_that("RunCellToCell handles raw output format", {
  skip_if_no_seurat()

  test_obj <- create_test_seurat()
  test_gt <- create_test_groundtruth()

  result <- RunCellToCell(
    node.object = test_obj,
    ground.truth = test_gt,
    assay = "RNA",
    meta.data.to.map = c("cell_type"),
    output_format = "raw"
  )

  # Assert it's a list with expected components
  expect_type(result, "list")
  expect_true("CellToCellMatrix" %in% names(result))
  expect_true("metadata" %in% names(result))

  # Assert matrix properties
  expect_true(is.matrix(result$CellToCellMatrix) |
              inherits(result$CellToCellMatrix, "Matrix"))
  expect_equal(nrow(result$CellToCellMatrix), 5)

  # Assert metadata is dataframe
  expect_s3_class(result$metadata, "data.frame")
  expect_true("SendingType" %in% colnames(result$metadata))
})


test_that("RunCellToCell validates required parameters", {
  skip_if_no_seurat()

  test_obj <- create_test_seurat()
  test_gt <- create_test_groundtruth()

  # Test missing required parameters
  expect_error(
    RunCellToCell(
      node.object = test_obj,
      # missing ground.truth
      assay = "RNA"
    )
  )

  expect_error(
    RunCellToCell(
      # missing node.object
      ground.truth = test_gt,
      assay = "RNA"
    )
  )
})


test_that("RunCellToCell handles single cell type", {
  skip_if_no_seurat()

  test_obj <- create_test_seurat()
  test_gt <- create_test_groundtruth()

  # Subset to single cell type
  single_type_obj <- subset(test_obj, idents = "CellType1")

  result <- RunCellToCell(
    node.object = single_type_obj,
    ground.truth = test_gt,
    assay = "RNA",
    meta.data.to.map = c("cell_type"),
    output_format = "seurat"
  )

  expect_s4_class(result, "Seurat")
  expect_gt(ncol(result), 0)

  # All pairs should be CellType1—CellType1
  expect_true(all(result$VectorType == "CellType1—CellType1"))
})


test_that("RunCellToCell handles multiple metadata columns", {
  skip_if_no_seurat()

  test_obj <- create_test_seurat()
  test_obj$condition <- rep(c("Control", "Treatment"), length.out = ncol(test_obj))
  test_gt <- create_test_groundtruth()

  result <- RunCellToCell(
    node.object = test_obj,
    ground.truth = test_gt,
    assay = "RNA",
    meta.data.to.map = c("cell_type", "condition"),
    output_format = "seurat"
  )

  # Check both metadata columns are propagated
  expect_true("cell_type.Sending" %in% colnames(result@meta.data))
  expect_true("cell_type.Receiving" %in% colnames(result@meta.data))
  expect_true("condition.Sending" %in% colnames(result@meta.data))
  expect_true("condition.Receiving" %in% colnames(result@meta.data))
})


test_that("RunCellToCell creates proper interaction matrix", {
  skip_if_no_seurat()

  test_obj <- create_test_seurat(n_cells = 20, n_genes = 50)
  test_gt <- create_test_groundtruth(n_mechanisms = 3)

  result <- RunCellToCell(
    node.object = test_obj,
    ground.truth = test_gt,
    assay = "RNA",
    meta.data.to.map = c("cell_type"),
    output_format = "raw"
  )

  # Check matrix dimensions
  mat <- result$CellToCellMatrix
  expect_equal(nrow(mat), 3)  # 3 mechanisms

  # Check that values are non-negative
  expect_true(all(mat >= 0))

  # Check rownames are mechanism names
  expect_true(all(rownames(mat) %in% paste0("Mechanism_", 1:3)))
})
