# NICHES Package Tests

This directory contains the test suite for the NICHES package using the `testthat` framework.

## Quick Start

### Run all tests

```r
# From R console in package root:
devtools::test()

# Or use testthat directly:
testthat::test_dir("tests/testthat")
```

### Run specific test file

```r
testthat::test_file("tests/testthat/test-RunCellToCell.R")
```

### Check test coverage

```r
covr::package_coverage()
```

## Directory Structure

```
tests/
├── README.md                     # This file
├── testthat.R                    # Entry point (auto-runs all tests)
├── testthat/
│   ├── fixtures/                 # Small test data files (empty for now)
│   ├── helper-data.R             # Test data generation functions
│   ├── test-RunCellToCell.R      # Tests for RunCellToCell()
│   ├── test-LoadGroundTruth.R    # Tests for ground truth loading
│   ├── test-PrepNodeObject.R     # Tests for PrepNodeObject()
│   └── ...                       # Add more test files here
```

## Test Data Functions

The `helper-data.R` file provides functions for generating test data:

- `create_test_seurat()` - Minimal Seurat object (50 cells, 100 genes, 3 cell types)
- `create_test_spatial_seurat()` - Seurat with spatial coordinates
- `create_test_groundtruth()` - Minimal L-R database (5 mechanisms)
- `create_test_matrix()` - Expression matrix without Seurat
- `create_test_metadata()` - Metadata dataframe
- `create_test_edgelist()` - Spatial edgelist for neighbor analysis
- `load_real_fantom5_subset()` - Load real FANTOM5 data (if available)

### Example Usage

```r
# In your test file:
test_that("my function works", {
  # Create test data
  test_obj <- create_test_seurat()
  test_gt <- create_test_groundtruth()

  # Run your function
  result <- my_function(test_obj, test_gt)

  # Make assertions
  expect_s4_class(result, "Seurat")
})
```

## Writing New Tests

### 1. Create a new test file

Follow the naming convention: `test-<source-file>.R`

```r
# For R/RunSystemToCell.R, create:
# tests/testthat/test-RunSystemToCell.R
```

### 2. Structure your tests

```r
test_that("descriptive test name", {
  # Arrange: Set up test data
  test_obj <- create_test_seurat()

  # Act: Run the function
  result <- my_function(test_obj)

  # Assert: Check expectations
  expect_s4_class(result, "Seurat")
  expect_equal(ncol(result), 50)
})
```

### 3. Test multiple scenarios

- Valid inputs
- Invalid inputs (should error)
- Edge cases (empty, single cell type, etc.)
- Different parameter combinations

## Existing Test Files

### test-RunCellToCell.R

Tests for the core `RunCellToCell()` transformation function:
- ✓ Creates valid Seurat output
- ✓ Handles raw output format
- ✓ Validates required parameters
- ✓ Handles single cell type
- ✓ Handles multiple metadata columns
- ✓ Creates proper interaction matrix

### test-LoadGroundTruth.R

Tests for ground truth loading system:
- ✓ Loads FANTOM5 correctly
- ✓ Handles custom database
- ✓ Filters mechanisms correctly
- ✓ Handles multi-subunit complexes
- ✓ Removes mechanisms with missing genes
- ✓ Works for different species

### test-PrepNodeObject.R

Tests for preprocessing function:
- ✓ Filters cells by min.cells.per.ident
- ✓ Filters genes by min.cells.per.gene
- ✓ Handles both filters simultaneously
- ✓ Preserves object when no filtering needed
- ✓ Validates input
- ✓ Handles NULL thresholds

## Test Coverage Goals

| Function Category | Target Coverage | Priority |
|-------------------|----------------|----------|
| Core transformations (Run*) | 100% | High |
| S3 dispatch (RunNICHES) | 100% | High |
| Ground truth loading | 90% | High |
| Helper functions | 80% | Medium |
| Visualization | 50% | Low |

## Adding Test Data

### Option 1: Add to helper-data.R

For programmatically generated data:

```r
# In tests/testthat/helper-data.R
create_my_test_data <- function() {
  # Generate test data
  return(data)
}
```

### Option 2: Add to fixtures/

For static data files:

```r
# Save data to fixtures/
saveRDS(my_test_data, "tests/testthat/fixtures/my_test_data.rds")

# Load in tests:
test_data <- readRDS(test_path("fixtures/my_test_data.rds"))
```

## Running Tests in Development

```r
# Load package for interactive testing
devtools::load_all()

# Run tests
devtools::test()

# Run specific test
testthat::test_file("tests/testthat/test-RunCellToCell.R")

# Run tests with reporter that shows details
devtools::test(reporter = "progress")
```

## Continuous Integration

Tests run automatically on:
- Every commit (if CI is set up)
- Every pull request
- Before CRAN submission

## Common Test Patterns

### Testing S3 dispatch

```r
test_that("function dispatches correctly", {
  seurat_obj <- create_test_seurat()
  matrix_obj <- create_test_matrix()

  # Both should work
  result1 <- my_function(seurat_obj)
  result2 <- my_function(matrix_obj)

  expect_s4_class(result1, "Seurat")
  expect_type(result2, "list")
})
```

### Testing error handling

```r
test_that("function validates inputs", {
  expect_error(
    my_function(NULL),
    "cannot be NULL"
  )

  expect_error(
    my_function(invalid_input),
    "must be a Seurat object"
  )
})
```

### Testing edge cases

```r
test_that("function handles edge cases", {
  # Empty input
  empty_obj <- create_test_seurat(n_cells = 0)
  expect_error(my_function(empty_obj))

  # Single cell
  single_cell <- create_test_seurat(n_cells = 1)
  result <- my_function(single_cell)
  expect_equal(ncol(result), 1)
})
```

## Debugging Tests

### Run test with browser

```r
# In your test file, add:
test_that("my test", {
  browser()  # Execution will pause here
  # ... rest of test
})
```

### Get detailed output

```r
devtools::test(reporter = "location")
```

### Check specific test

```r
# Run just one test_that() block
testthat::test_file(
  "tests/testthat/test-RunCellToCell.R",
  filter = "creates valid Seurat output"
)
```

## Resources

- [R Packages - Testing chapter](https://r-pkgs.org/testing-basics.html)
- [testthat documentation](https://testthat.r-lib.org/)
- [Writing tests for Bioconductor](https://bioconductor.org/developers/how-to/unitTesting-guidelines/)
- See `TESTING_GUIDE.md` in package root for comprehensive guide

## Next Steps

To continue building the test suite:

1. Write tests for remaining Run* functions:
   - test-RunCellToSystem.R
   - test-RunSystemToCell.R
   - test-RunCellToCellSpatial.R
   - test-RunCellToNeighborhood.R
   - test-RunNeighborhoodToCell.R

2. Write tests for S3 dispatch:
   - test-RunNICHES.R

3. Write tests for spatial functions:
   - test-ComputeEdgelist.R

4. Write tests for utilities:
   - test-GetSeuratAssay.R
   - test-ReturnCelltypes.R
   - test-FilterGroundTruth.R

5. Write tests for visualization:
   - test-CircuitPlot.R
   - test-NetworkPlot.R
   - test-DotPlotSCC.R

Aim for 80%+ overall test coverage!
