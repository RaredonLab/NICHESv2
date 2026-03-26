# NICHES Package Testing Guide

## Overview

This guide explains how to organize, write, and run tests for the NICHES package using the `testthat` framework—the standard for R package testing.

---

## Directory Structure

```
NICHES/
├── tests/
│   ├── testthat.R                    # Entry point (auto-runs all tests)
│   ├── testthat/
│   │   ├── fixtures/                 # Small test data (modern approach)
│   │   │   ├── test_seurat.rds       # Minimal Seurat object
│   │   │   ├── test_matrix.rds       # Small expression matrix
│   │   │   ├── test_metadata.rds     # Metadata dataframe
│   │   │   └── test_groundtruth.rds  # Minimal L-R database
│   │   │
│   │   ├── helper-data.R             # Functions to load/generate test data
│   │   ├── setup.R                   # Runs once before all tests
│   │   │
│   │   ├── test-RunNICHES.R          # Tests for RunNICHES.R
│   │   ├── test-RunCellToCell.R      # Tests for RunCellToCell.R
│   │   ├── test-LoadGroundTruth.R    # Tests for LoadGroundTruth.R
│   │   ├── test-PrepNodeObject.R     # Tests for PrepNodeObject.R
│   │   ├── test-ComputeEdgelist.R    # Tests for ComputeEdgelist.R
│   │   ├── test-visualization.R      # Tests for CircuitPlot, NetworkPlot
│   │   └── ...                       # One test file per R source file
│   │
│   └── manual/                       # Manual testing scripts (not automated)
│       └── interactive_examples.R
│
└── inst/
    └── testdata/                     # Alternative: larger test datasets
        └── example_dataset.rds       # (if fixtures are too small)
```

---

## Key Principles

### 1. Test File Naming Convention

- **Pattern**: `test-<source-file>.R` corresponds to `R/<source-file>.R`
- **Examples**:
  - `R/RunNICHES.R` → `tests/testthat/test-RunNICHES.R`
  - `R/LoadFantom5.R` → `tests/testthat/test-LoadFantom5.R`
  - `R/CircuitPlot.R` → `tests/testthat/test-CircuitPlot.R`

### 2. Test Organization Within Files

```r
# tests/testthat/test-RunCellToCell.R

test_that("RunCellToCell handles valid Seurat input", {
  # Arrange: Set up test data
  test_obj <- load_test_seurat()
  test_gt <- load_test_groundtruth()

  # Act: Run the function
  result <- RunCellToCell(
    node.object = test_obj,
    ground.truth = test_gt,
    assay = "RNA",
    output_format = "seurat"
  )

  # Assert: Check expectations
  expect_s4_class(result, "Seurat")
  expect_true("CellToCell" %in% names(result@assays))
  expect_gt(ncol(result), 0)
})

test_that("RunCellToCell validates input parameters", {
  test_obj <- load_test_seurat()

  expect_error(
    RunCellToCell(node.object = NULL, ground.truth = list()),
    "node.object cannot be NULL"
  )
})

test_that("RunCellToCell returns raw format correctly", {
  # ... test raw output format
})
```

### 3. Test Data Size

Keep test data **MINIMAL**:
- **Seurat objects**: 50-200 cells, 100-500 genes
- **Expression matrices**: Small enough to run tests in < 1 second
- **Ground truth**: 5-10 L-R mechanisms
- **Goal**: Fast tests that run in milliseconds

---

## Test Data Storage Options

### Option 1: `tests/testthat/fixtures/` (RECOMMENDED)

**Best for**: Small, static test data files

**Structure**:
```
tests/testthat/fixtures/
├── minimal_seurat.rds      # 100 cells × 200 genes
├── minimal_matrix.rds       # Expression matrix
├── minimal_metadata.rds     # Metadata dataframe
├── minimal_groundtruth.rds  # 10 L-R pairs
└── spatial_coords.rds       # X/Y coordinates
```

**Usage**:
```r
# In test files, load with:
test_seurat <- readRDS(test_path("fixtures/minimal_seurat.rds"))
```

**Advantages**:
- Self-contained (data lives with tests)
- Easy to create/modify
- Standard testthat convention

---

### Option 2: `inst/testdata/`

**Best for**: Larger datasets that might also be useful for examples

**Structure**:
```
inst/testdata/
├── example_dataset.rds      # ~1000 cells for integration tests
└── example_spatial.rds      # Spatial data example
```

**Usage**:
```r
# In test files:
data_path <- system.file("testdata", "example_dataset.rds", package = "NICHES")
test_data <- readRDS(data_path)
```

**Advantages**:
- Available after package installation
- Can be used in vignettes/examples
- Good for "real-world" integration tests

---

### Option 3: Programmatic Generation in `helper-data.R`

**Best for**: Dynamic data that's simple to generate

**Structure**:
```r
# tests/testthat/helper-data.R

#' Create minimal Seurat object for testing
#' @return Seurat object with 50 cells, 100 genes, 3 cell types
create_test_seurat <- function() {
  library(Seurat)

  # Generate fake expression data
  counts <- matrix(
    rpois(5000, lambda = 2),  # 100 genes × 50 cells
    nrow = 100,
    dimnames = list(
      paste0("Gene_", 1:100),
      paste0("Cell_", 1:50)
    )
  )

  # Create metadata
  metadata <- data.frame(
    cell_type = rep(c("TypeA", "TypeB", "TypeC"), length.out = 50),
    sample = "test_sample",
    row.names = colnames(counts)
  )

  # Create Seurat object
  obj <- CreateSeuratObject(
    counts = counts,
    meta.data = metadata
  )

  Idents(obj) <- obj$cell_type
  return(obj)
}

#' Create minimal ground truth for testing
create_test_groundtruth <- function() {
  list(
    source.subunits = matrix(
      c("Gene_1", "Gene_3", "Gene_5", "Gene_7", "Gene_9"),
      ncol = 1,
      dimnames = list(
        paste0("Mechanism_", 1:5),
        "subunit_1"
      )
    ),
    target.subunits = matrix(
      c("Gene_2", "Gene_4", "Gene_6", "Gene_8", "Gene_10"),
      ncol = 1,
      dimnames = list(
        paste0("Mechanism_", 1:5),
        "subunit_1"
      )
    )
  )
}
```

**Usage**:
```r
# In test files (helper-data.R is auto-loaded by testthat):
test_obj <- create_test_seurat()
test_gt <- create_test_groundtruth()
```

**Advantages**:
- No files to maintain
- Easy to customize per test
- Self-documenting (code shows exactly what's in the data)

---

## Recommended Approach for NICHES

**Use a hybrid approach**:

1. **`tests/testthat/fixtures/`**: Store a few pre-built objects
   - Real Seurat object from actual data (heavily downsampled)
   - Real FANTOM5 subset (5-10 mechanisms)
   - Pre-computed edgelist for spatial tests

2. **`tests/testthat/helper-data.R`**: Generate simple test objects
   - Quick mock objects for unit tests
   - Edge cases (empty objects, single cell type, etc.)

3. **`inst/testdata/`**: Optional larger datasets
   - For integration tests
   - For vignette examples

---

## Setting Up testthat

### Step 1: Add testthat to DESCRIPTION

```r
# In DESCRIPTION file, add to Suggests:
Suggests:
  knitr,
  rmarkdown,
  testthat (>= 3.0.0)
```

### Step 2: Initialize testthat structure

```r
# From R console in package root:
usethis::use_testthat()
```

This creates:
- `tests/testthat.R` (entry point)
- `tests/testthat/` directory

### Step 3: Create your first test

```r
usethis::use_test("RunNICHES")
```

This creates `tests/testthat/test-RunNICHES.R`

---

## Example Test Files

### Example 1: Testing Core Transformation Function

```r
# tests/testthat/test-RunCellToCell.R

test_that("RunCellToCell creates valid Seurat output", {
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
  expect_true("CellToCell" %in% names(result@assays))

  # Assert dimensions
  expect_equal(nrow(result), 5)  # 5 mechanisms in test ground truth
  expect_gt(ncol(result), 0)     # At least some cell pairs

  # Assert metadata
  expect_true("SendingType" %in% colnames(result@meta.data))
  expect_true("ReceivingType" %in% colnames(result@meta.data))
  expect_true("VectorType" %in% colnames(result@meta.data))
  expect_true("cell_type.Sending" %in% colnames(result@meta.data))
})

test_that("RunCellToCell handles raw output format", {
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
})

test_that("RunCellToCell validates inputs", {
  test_obj <- create_test_seurat()
  test_gt <- create_test_groundtruth()

  # Test NULL object
  expect_error(
    RunCellToCell(node.object = NULL, ground.truth = test_gt),
    "node.object"
  )

  # Test NULL ground truth
  expect_error(
    RunCellToCell(node.object = test_obj, ground.truth = NULL),
    "ground.truth"
  )

  # Test invalid output format
  expect_error(
    RunCellToCell(
      node.object = test_obj,
      ground.truth = test_gt,
      output_format = "invalid"
    ),
    "output_format"
  )
})

test_that("RunCellToCell handles edge cases", {
  test_obj <- create_test_seurat()
  test_gt <- create_test_groundtruth()

  # Single cell type
  single_type_obj <- subset(test_obj, idents = "TypeA")
  result <- RunCellToCell(
    node.object = single_type_obj,
    ground.truth = test_gt,
    assay = "RNA",
    output_format = "seurat"
  )

  expect_s4_class(result, "Seurat")
  # All pairs should be TypeA-TypeA
  expect_true(all(result$VectorType == "TypeA—TypeA"))
})
```

---

### Example 2: Testing S3 Dispatch

```r
# tests/testthat/test-RunNICHES.R

test_that("RunNICHES dispatches to correct method for Seurat", {
  test_obj <- create_test_seurat()

  result <- RunNICHES(
    object = test_obj,
    species = "human",
    LR.database = "fantom5",
    CellToCell = TRUE,
    CellToSystem = FALSE,
    SystemToCell = FALSE
  )

  expect_type(result, "list")
  expect_true("CellToCell" %in% names(result))
  expect_s4_class(result$CellToCell, "Seurat")
})

test_that("RunNICHES dispatches to correct method for matrix", {
  test_obj <- create_test_seurat()

  # Extract components
  expr_matrix <- GetAssayData(test_obj, layer = "counts")
  metadata <- test_obj@meta.data

  result <- RunNICHES(
    object = expr_matrix,
    meta.data.df = metadata,
    cell_types = "cell_type",
    species = "human",
    LR.database = "fantom5",
    CellToCell = TRUE
  )

  expect_type(result, "list")
  expect_true("CellToCell" %in% names(result))
})

test_that("RunNICHES enables multiple organization modes", {
  test_obj <- create_test_seurat()

  result <- RunNICHES(
    object = test_obj,
    species = "human",
    LR.database = "fantom5",
    CellToCell = TRUE,
    CellToSystem = TRUE,
    SystemToCell = TRUE
  )

  expect_length(result, 3)
  expect_true(all(c("CellToCell", "CellToSystem", "SystemToCell") %in% names(result)))
})
```

---

### Example 3: Testing Ground Truth Loading

```r
# tests/testthat/test-LoadGroundTruth.R

test_that("LoadGroundTruth loads FANTOM5 correctly", {
  skip_if_not_installed("NICHES")  # Skip if package not installed

  gt <- LoadGroundTruth(
    LR.database = "fantom5",
    species = "human",
    input_rownames = paste0("Gene_", 1:100)
  )

  expect_type(gt, "list")
  expect_true("source.subunits" %in% names(gt))
  expect_true("target.subunits" %in% names(gt))
  expect_true(is.matrix(gt$source.subunits))
  expect_true(is.matrix(gt$target.subunits))
})

test_that("LoadGroundTruth filters to input genes", {
  # Create test data with known genes
  input_genes <- c("TGFB1", "TGFBR1", "IL6", "IL6R")

  gt <- LoadGroundTruth(
    LR.database = "fantom5",
    species = "human",
    input_rownames = input_genes
  )

  # All genes in ground truth should be in input_genes
  all_genes <- c(
    gt$source.subunits[!is.na(gt$source.subunits)],
    gt$target.subunits[!is.na(gt$target.subunits)]
  )

  expect_true(all(all_genes %in% input_genes))
})

test_that("LoadGroundTruth handles custom database", {
  custom_db <- data.frame(
    ligand = c("Gene_1", "Gene_3"),
    receptor = c("Gene_2", "Gene_4")
  )

  gt <- LoadGroundTruth(
    LR.database = "custom",
    custom_LR_database = custom_db,
    input_rownames = paste0("Gene_", 1:10)
  )

  expect_equal(nrow(gt$source.subunits), 2)
  expect_equal(nrow(gt$target.subunits), 2)
})
```

---

### Example 4: Testing Spatial Functions

```r
# tests/testthat/test-ComputeEdgelist.R

test_that("ComputeEdgelist computes k-NN correctly", {
  # Create test object with spatial coordinates
  test_obj <- create_test_seurat()
  test_obj$x <- runif(ncol(test_obj), 0, 100)
  test_obj$y <- runif(ncol(test_obj), 0, 100)

  edgelist <- ComputeEdgelist(
    node.object = test_obj,
    position.x = "x",
    position.y = "y",
    k = 5,
    nn.method = "aoz"
  )

  expect_s3_class(edgelist, "data.frame")
  expect_true("from" %in% colnames(edgelist))
  expect_true("to" %in% colnames(edgelist))
  expect_gt(nrow(edgelist), 0)

  # Each cell should have approximately k neighbors
  neighbors_per_cell <- table(edgelist$from)
  expect_true(all(neighbors_per_cell <= 5))
})

test_that("ComputeEdgelist validates spatial coordinates", {
  test_obj <- create_test_seurat()

  expect_error(
    ComputeEdgelist(
      node.object = test_obj,
      position.x = "nonexistent",
      position.y = "y",
      k = 5
    ),
    "position"
  )
})
```

---

### Example 5: Testing Visualization Functions

```r
# tests/testthat/test-CircuitPlot.R

test_that("CircuitPlot creates ggplot object", {
  # Arrange: Create test data
  transcr_obj <- create_test_seurat()

  # Create minimal CellToCell object
  ctc_obj <- RunCellToCell(
    node.object = transcr_obj,
    ground.truth = create_test_groundtruth(),
    assay = "RNA",
    output_format = "seurat"
  )

  # Act: Create plot
  p <- CircuitPlot(
    transcr.obj = transcr_obj,
    connect.obj = ctc_obj,
    feature = "Mechanism_1",
    group.by = "cell_type"
  )

  # Assert
  expect_s3_class(p, "ggplot")
})

test_that("CircuitPlot validates inputs", {
  transcr_obj <- create_test_seurat()
  ctc_obj <- create_test_seurat()  # Wrong type

  expect_error(
    CircuitPlot(
      transcr.obj = transcr_obj,
      connect.obj = ctc_obj,
      feature = "nonexistent_mechanism",
      group.by = "cell_type"
    )
  )
})
```

---

## Running Tests

### Run all tests

```r
# From R console:
devtools::test()

# Or from command line:
R CMD check .
```

### Run specific test file

```r
testthat::test_file("tests/testthat/test-RunCellToCell.R")
```

### Run tests with coverage

```r
covr::package_coverage()
```

### Run tests on package check

```r
devtools::check()  # Runs tests as part of R CMD check
```

---

## Test Coverage Goals

Aim for **80%+ code coverage** with priority on:

1. **Core transformation functions** (100% coverage)
   - RunCellToCell, RunCellToSystem, RunSystemToCell
   - RunCellToCellSpatial, RunCellToNeighborhood, RunNeighborhoodToCell

2. **S3 dispatch** (100% coverage)
   - All RunNICHES methods

3. **Ground truth loading** (90%+ coverage)
   - LoadGroundTruth, LoadFantom5, LoadOmniPath, LoadCustom

4. **Helper functions** (80%+ coverage)
   - PrepNodeObject, ComputeEdgelist, GetSeuratAssay

5. **Visualization** (50%+ coverage)
   - Focus on error handling, not rendering

---

## Best Practices

### 1. Keep Tests Fast

- Use minimal data (50-200 cells)
- Mock external dependencies
- Skip slow tests in CI with `skip_on_ci()`

### 2. Test Behavior, Not Implementation

```r
# GOOD: Test what the function does
expect_equal(ncol(result), 50)

# BAD: Test how it does it
expect_true(exists("temp_var_inside_function"))
```

### 3. Use Descriptive Test Names

```r
# GOOD
test_that("RunCellToCell returns Seurat object with correct metadata", {})

# BAD
test_that("test 1", {})
```

### 4. Test Edge Cases

- Empty inputs
- Single cell type
- Single mechanism
- NA values
- Very large inputs (performance)

### 5. Use Fixtures for Complex Setup

```r
# tests/testthat/fixtures/complex_setup.R
# This runs once and results are cached
saveRDS(expensive_computation(), "complex_object.rds")
```

---

## Continuous Integration

Add GitHub Actions for automated testing:

```yaml
# .github/workflows/R-CMD-check.yaml
on: [push, pull_request]

jobs:
  R-CMD-check:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v2
      - uses: r-lib/actions/setup-r@v2
      - name: Install dependencies
        run: |
          install.packages(c("devtools", "testthat"))
          devtools::install_deps(dependencies = TRUE)
      - name: Check
        run: devtools::check()
```

---

## Summary

### Quick Start Checklist

- [ ] Add `testthat (>= 3.0.0)` to DESCRIPTION Suggests
- [ ] Run `usethis::use_testthat()` to create infrastructure
- [ ] Create `tests/testthat/helper-data.R` with test data generators
- [ ] Create `tests/testthat/fixtures/` for static test data
- [ ] Write test files following `test-<source>.R` convention
- [ ] Aim for 80%+ coverage on core functions
- [ ] Run `devtools::test()` frequently during development
- [ ] Set up CI/CD with GitHub Actions

### File Organization Summary

```
tests/
├── testthat.R                    # Auto-generated entry point
├── testthat/
│   ├── fixtures/                 # Small test data files
│   ├── helper-data.R             # Data generation functions
│   ├── test-RunNICHES.R          # Test S3 dispatch
│   ├── test-RunCellToCell.R      # Test transformations
│   ├── test-LoadGroundTruth.R    # Test database loading
│   └── test-visualization.R      # Test plots
```

---

**Now you're ready to write comprehensive tests for NICHES!** 🧪
