# ==============================================================================
# FILE:    tests/testthat/test-load_LRM_database.R
# PURPOSE: Unit tests for load_LRM_database().
#          All inputs synthesised inline — no shared fixtures needed.
# ==============================================================================

# ── Custom data.frame path ─────────────────────────────────────────────────────

test_that("load_LRM_database: custom data.frame standardises column names", {
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

# ── FANTOM5 path ───────────────────────────────────────────────────────────────

test_that("load_LRM_database: fantom5 human returns data.frame with ligand and receptor", {
  result <- load_LRM_database("fantom5", species = "human", verbose = FALSE)
  expect_true(is.data.frame(result))
  expect_true(all(c("ligand", "receptor") %in% colnames(result)))
  expect_gt(nrow(result), 0L)
})

# ── ConnectomeDB2025 path (default database) ───────────────────────────────────

test_that("load_LRM_database: default call returns valid data.frame", {
  result <- load_LRM_database(verbose = FALSE)
  expect_true(is.data.frame(result))
  expect_true(all(c("ligand", "receptor") %in% colnames(result)))
  expect_gt(nrow(result), 0L)
})

test_that("load_LRM_database: connectomedb2025 human returns data.frame with ligand and receptor", {
  result <- load_LRM_database("connectomedb2025", species = "human", verbose = FALSE)
  expect_true(is.data.frame(result))
  expect_true(all(c("ligand", "receptor") %in% colnames(result)))
  expect_gt(nrow(result), 0L)
})

test_that("load_LRM_database: connectomedb2025 pig is supported", {
  result <- load_LRM_database("connectomedb2025", species = "pig", verbose = FALSE)
  expect_true(is.data.frame(result))
  expect_gt(nrow(result), 0L)
})

test_that("load_LRM_database: errors on invalid species for connectomedb2025", {
  expect_error(
    load_LRM_database("connectomedb2025", species = "zebrafish", verbose = FALSE),
    regexp = "species"
  )
})

# ── Input validation ───────────────────────────────────────────────────────────

test_that("load_LRM_database: errors on wrong column names in custom data.frame", {
  custom.db <- data.frame(L = "Gene1", R = "Gene2")
  expect_error(
    load_LRM_database(custom.db, ligand.col = "ligand", receptor.col = "receptor",
                      verbose = FALSE),
    regexp = "not found"
  )
})

test_that("load_LRM_database: errors on unknown string shorthand", {
  expect_error(
    load_LRM_database("connectomeDB"),
    regexp = "Unknown built-in"
  )
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