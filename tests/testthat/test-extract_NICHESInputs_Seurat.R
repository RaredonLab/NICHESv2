# ==============================================================================
# FILE:    tests/testthat/test-extract_NICHESInputs_Seurat.R
# PURPOSE: Unit tests for extract_NICHESInputs_Seurat().
#          All Seurat objects are synthesised inline — no shared fixtures needed.
#          All blocks skip automatically when SeuratObject is not installed.
# ==============================================================================

test_that("extract_NICHESInputs_Seurat: returns list with count.mtx and meta.data", {
  skip_if_not_installed("SeuratObject")
  
  set.seed(42L)
  counts <- abs(Matrix::rsparsematrix(20L, 10L, density = 0.6))
  rownames(counts) <- paste0("gene", seq_len(20L))
  colnames(counts) <- paste0("cell", seq_len(10L))
  meta <- data.frame(
    cell.type = rep(c("A", "B"), 5L),
    row.names = colnames(counts)
  )
  seurat.obj <- SeuratObject::CreateSeuratObject(counts = counts,
                                                 meta.data = meta)
  
  result <- extract_NICHESInputs_Seurat(seurat.obj, verbose = FALSE)
  
  expect_type(result, "list")
  expect_named(result, c("count.mtx", "meta.data"))
})

test_that("extract_NICHESInputs_Seurat: count.mtx has correct dimensions and names", {
  skip_if_not_installed("SeuratObject")
  
  set.seed(42L)
  counts <- abs(Matrix::rsparsematrix(20L, 10L, density = 0.6))
  rownames(counts) <- paste0("gene", seq_len(20L))
  colnames(counts) <- paste0("cell", seq_len(10L))
  seurat.obj <- SeuratObject::CreateSeuratObject(counts = counts)
  
  result <- extract_NICHESInputs_Seurat(seurat.obj, verbose = FALSE)
  
  expect_true(inherits(result$count.mtx, "dgCMatrix"))
  expect_equal(nrow(result$count.mtx), 20L)
  expect_equal(ncol(result$count.mtx), 10L)
  expect_equal(rownames(result$count.mtx), paste0("gene", seq_len(20L)))
  expect_equal(colnames(result$count.mtx), paste0("cell", seq_len(10L)))
})

test_that("extract_NICHESInputs_Seurat: meta.data preserves user-supplied columns", {
  skip_if_not_installed("SeuratObject")
  
  set.seed(42L)
  counts <- abs(Matrix::rsparsematrix(5L, 5L, density = 0.8))
  rownames(counts) <- paste0("gene", seq_len(5L))
  colnames(counts) <- paste0("cell", seq_len(5L))
  meta <- data.frame(
    cell.type = c("A", "B", "A", "B", "A"),
    x         = c(1.1, 2.2, 3.3, 4.4, 5.5),
    y         = c(5.5, 4.4, 3.3, 2.2, 1.1),
    row.names = paste0("cell", seq_len(5L))
  )
  seurat.obj <- SeuratObject::CreateSeuratObject(counts = counts,
                                                 meta.data = meta)
  
  result <- extract_NICHESInputs_Seurat(seurat.obj, verbose = FALSE)
  
  expect_true(is.data.frame(result$meta.data))
  expect_equal(rownames(result$meta.data), paste0("cell", seq_len(5L)))
  expect_true("cell.type" %in% colnames(result$meta.data))
  expect_true("x"         %in% colnames(result$meta.data))
  expect_true("y"         %in% colnames(result$meta.data))
})

test_that("extract_NICHESInputs_Seurat: errors on missing assay with available names listed", {
  skip_if_not_installed("SeuratObject")
  
  set.seed(42L)
  counts <- abs(Matrix::rsparsematrix(5L, 5L, density = 0.8))
  rownames(counts) <- paste0("gene", seq_len(5L))
  colnames(counts) <- paste0("cell", seq_len(5L))
  seurat.obj <- SeuratObject::CreateSeuratObject(counts = counts)
  
  expect_error(
    extract_NICHESInputs_Seurat(seurat.obj, assay = "NOTEXIST", verbose = FALSE),
    "Available assays"
  )
})

test_that("extract_NICHESInputs_Seurat: errors on missing layer with available names listed", {
  skip_if_not_installed("SeuratObject")
  
  set.seed(42L)
  counts <- abs(Matrix::rsparsematrix(5L, 5L, density = 0.8))
  rownames(counts) <- paste0("gene", seq_len(5L))
  colnames(counts) <- paste0("cell", seq_len(5L))
  seurat.obj <- SeuratObject::CreateSeuratObject(counts = counts)
  
  expect_error(
    extract_NICHESInputs_Seurat(seurat.obj, layer = "NOTLAYER", verbose = FALSE),
    regexp = "not found"
  )
})

test_that("extract_NICHESInputs_Seurat: verbose=TRUE emits spatial hint message", {
  skip_if_not_installed("SeuratObject")
  
  set.seed(42L)
  counts <- abs(Matrix::rsparsematrix(5L, 5L, density = 0.8))
  rownames(counts) <- paste0("gene", seq_len(5L))
  colnames(counts) <- paste0("cell", seq_len(5L))
  seurat.obj <- SeuratObject::CreateSeuratObject(counts = counts)
  
  expect_message(
    extract_NICHESInputs_Seurat(seurat.obj, verbose = TRUE),
    regexp = "spatial"
  )
})

test_that("extract_NICHESInputs_Seurat: verbose=FALSE suppresses all messages", {
  skip_if_not_installed("SeuratObject")
  
  set.seed(42L)
  counts <- abs(Matrix::rsparsematrix(5L, 5L, density = 0.8))
  rownames(counts) <- paste0("gene", seq_len(5L))
  colnames(counts) <- paste0("cell", seq_len(5L))
  seurat.obj <- SeuratObject::CreateSeuratObject(counts = counts)
  
  expect_no_message(
    extract_NICHESInputs_Seurat(seurat.obj, verbose = FALSE)
  )
})

test_that("extract_NICHESInputs_Seurat: class-based V4/V5 dispatch produces correct shape", {
  skip_if_not_installed("SeuratObject")
  
  # The key invariant: extract_NICHESInputs_Seurat() detects assay class via
  # inherits(assay.obj, "Assay5"), NOT by checking the SeuratObject package
  # version. Whichever version is installed, CreateSeuratObject() produces the
  # correct class for that version, and the function must handle it correctly.
  set.seed(42L)
  counts <- abs(Matrix::rsparsematrix(10L, 8L, density = 0.7))
  rownames(counts) <- paste0("gene", seq_len(10L))
  colnames(counts) <- paste0("cell", seq_len(8L))
  seurat.obj <- SeuratObject::CreateSeuratObject(counts = counts)
  
  assay.class <- class(seurat.obj@assays[["RNA"]])
  result      <- extract_NICHESInputs_Seurat(seurat.obj, verbose = FALSE)
  
  expect_equal(dim(result$count.mtx), c(10L, 8L))
  expect_true(assay.class %in% c("Assay", "Assay5"))
})