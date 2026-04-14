# ==============================================================================
# FILE:    tests/testthat/test-utilities.R
# PURPOSE: Tests for add_cell_meta(), load_LRM_database() (custom path)
# ==============================================================================

# ── add_cell_meta ──────────────────────────────────────────────────────────────

test_that("add_cell_meta: adds a new column to $cell.meta", {
  .setup_fixtures()
  new.col <- data.frame(
    cell    = niches.obj$cell.list,
    cluster = sample(1:5, length(niches.obj$cell.list), replace = TRUE),
    stringsAsFactors = FALSE
  )
  updated <- add_cell_meta(niches.obj, new.col)
  expect_true("cluster" %in% colnames(updated$cell.meta))
})

test_that("add_cell_meta: overwrites an existing column", {
  .setup_fixtures()
  override <- data.frame(
    cell     = niches.obj$cell.list,
    celltype = "overridden",
    stringsAsFactors = FALSE
  )
  updated <- add_cell_meta(niches.obj, override)
  expect_true(all(updated$cell.meta$celltype == "overridden"))
})

test_that("add_cell_meta: rebuilds $edge.meta with new Joint column", {
  .setup_fixtures()
  new.col <- data.frame(
    cell    = niches.obj$cell.list,
    cluster = paste0("C", sample(1:3, length(niches.obj$cell.list), replace = TRUE)),
    stringsAsFactors = FALSE
  )
  updated <- add_cell_meta(niches.obj, new.col)
  expect_true("Joint.cluster" %in% colnames(updated$edge.meta))
})

test_that("add_cell_meta: returns NICHESObject", {
  .setup_fixtures()
  dummy <- data.frame(cell = niches.obj$cell.list[1], extra = "x")
  updated <- add_cell_meta(niches.obj, dummy)
  expect_s3_class(updated, "NICHESObject")
})

test_that("add_cell_meta: errors when 'cell' column missing from new.meta", {
  .setup_fixtures()
  bad.meta <- data.frame(barcode = niches.obj$cell.list, val = 1)
  expect_error(add_cell_meta(niches.obj, bad.meta), regexp = "'cell' column")
})

test_that("add_cell_meta: errors on non-NICHESObject input", {
  expect_error(add_cell_meta(list(), data.frame(cell = "A")), regexp = "NICHESObject")
})

# ── load_LRM_database (custom data.frame path) ────────────────────────────────

test_that("load_LRM_database: custom data.frame standardises column names", {
  .setup_fixtures()
  custom.db <- data.frame(
    MyLigand   = c("Gene1", "Gene3"),
    MyReceptor = c("Gene2", "Gene4"),
    pathway    = c("PathA", "PathB"),
    stringsAsFactors = FALSE
  )
  result <- load_LRM_database(
    custom.db,
    ligand.col   = "MyLigand",
    receptor.col = "MyReceptor",
    verbose      = FALSE
  )
  expect_true(all(c("ligand", "receptor") %in% colnames(result)))
  expect_equal(nrow(result), 2L)
})

test_that("load_LRM_database: extra columns are preserved", {
  custom.db <- data.frame(
    ligand   = "Gene1",
    receptor = "Gene2",
    score    = 0.95,
    stringsAsFactors = FALSE
  )
  result <- load_LRM_database(custom.db, verbose = FALSE)
  expect_true("score" %in% colnames(result))
})

test_that("load_LRM_database: NA rows are silently removed", {
  custom.db <- data.frame(
    ligand   = c("Gene1", NA,      "Gene3"),
    receptor = c("Gene2", "Gene4", NA),
    stringsAsFactors = FALSE
  )
  result <- load_LRM_database(custom.db, verbose = FALSE)
  expect_equal(nrow(result), 1L)
})

test_that("load_LRM_database: errors on wrong column names", {
  custom.db <- data.frame(L = "Gene1", R = "Gene2")
  expect_error(
    load_LRM_database(custom.db, ligand.col = "ligand", receptor.col = "receptor",
                      verbose = FALSE),
    regexp = "not found"
  )
})

test_that("load_LRM_database: errors on unknown string shorthand", {
  expect_error(load_LRM_database("connectomeDB"), regexp = "Unknown built-in")
})

test_that("load_LRM_database: errors on invalid species for fantom5", {
  expect_error(
    load_LRM_database("fantom5", species = "zebrafish", verbose = FALSE),
    regexp = "species"
  )
})

test_that("load_LRM_database: errors on invalid species for omnipath", {
  expect_error(
    load_LRM_database("omnipath", species = "pig", verbose = FALSE),
    regexp = "species"
  )
})

