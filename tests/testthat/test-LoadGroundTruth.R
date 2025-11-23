# Tests for LoadGroundTruth and related functions
# tests/testthat/test-LoadGroundTruth.R

test_that("LoadGroundTruth loads FANTOM5 correctly", {
  skip_if_not_installed("NICHES")

  # Create gene list that matches some real genes
  input_genes <- c("TGFB1", "TGFBR1", "TGFBR2", "IL6", "IL6R", "IL6ST",
                   "TNF", "TNFRSF1A", "TNFRSF1B", "GENE1", "GENE2")

  gt <- LoadGroundTruth(
    LR.database = "fantom5",
    species = "human",
    input_rownames = input_genes
  )

  # Assert structure
  expect_type(gt, "list")
  expect_true("source.subunits" %in% names(gt))
  expect_true("target.subunits" %in% names(gt))
  expect_true(is.matrix(gt$source.subunits))
  expect_true(is.matrix(gt$target.subunits))

  # Assert same number of mechanisms for ligands and receptors
  expect_equal(nrow(gt$source.subunits), nrow(gt$target.subunits))

  # Assert that only input genes are present
  all_ligands <- gt$source.subunits[!is.na(gt$source.subunits)]
  all_receptors <- gt$target.subunits[!is.na(gt$target.subunits)]
  all_genes <- unique(c(all_ligands, all_receptors))

  expect_true(all(all_genes %in% input_genes))
})


test_that("LoadGroundTruth handles custom database", {
  # Create custom L-R database
  custom_db <- data.frame(
    ligand = c("GENE1", "GENE3", "GENE5"),
    receptor = c("GENE2", "GENE4", "GENE6"),
    stringsAsFactors = FALSE
  )

  input_genes <- paste0("GENE", 1:10)

  gt <- LoadGroundTruth(
    LR.database = "custom",
    custom_LR_database = custom_db,
    input_rownames = input_genes
  )

  # Assert correct number of mechanisms
  expect_equal(nrow(gt$source.subunits), 3)
  expect_equal(nrow(gt$target.subunits), 3)

  # Assert correct ligands and receptors
  expect_true(all(gt$source.subunits[, 1] %in% c("GENE1", "GENE3", "GENE5")))
  expect_true(all(gt$target.subunits[, 1] %in% c("GENE2", "GENE4", "GENE6")))
})


test_that("LoadGroundTruth filters mechanisms correctly", {
  custom_db <- data.frame(
    ligand = c("GENE1", "GENE3", "GENE99"),  # GENE99 not in input
    receptor = c("GENE2", "GENE98", "GENE6"), # GENE98 not in input
    stringsAsFactors = FALSE
  )

  input_genes <- paste0("GENE", 1:10)

  gt <- LoadGroundTruth(
    LR.database = "custom",
    custom_LR_database = custom_db,
    input_rownames = input_genes
  )

  # Only first mechanism should remain (GENE1—GENE2)
  # GENE3—GENE98 filtered (receptor missing)
  # GENE99—GENE6 filtered (ligand missing)
  expect_equal(nrow(gt$source.subunits), 1)
  expect_equal(gt$source.subunits[1, 1], "GENE1")
  expect_equal(gt$target.subunits[1, 1], "GENE2")
})


test_that("LoadCustom handles multi-subunit complexes", {
  # Custom database with multi-subunit complexes (underscore-separated)
  custom_db <- data.frame(
    ligand = c("GENE1", "GENE3_GENE4"),  # Second is 2-subunit
    receptor = c("GENE2", "GENE5"),
    stringsAsFactors = FALSE
  )

  gt <- LoadCustom(custom_LR_database = custom_db)

  # Assert structure
  expect_equal(nrow(gt$source.subunits), 2)
  expect_equal(ncol(gt$source.subunits), 2)  # Max 2 subunits

  # Check first mechanism (single subunit)
  expect_equal(gt$source.subunits[1, 1], "GENE1")
  expect_true(is.na(gt$source.subunits[1, 2]))

  # Check second mechanism (two subunits)
  expect_equal(gt$source.subunits[2, 1], "GENE3")
  expect_equal(gt$source.subunits[2, 2], "GENE4")
})


test_that("FilterGroundTruth removes mechanisms with missing genes", {
  # Create ground truth with some genes not in input
  gt <- list(
    source.subunits = matrix(
      c("GENE1", "GENE3", "GENE99"),
      ncol = 1,
      dimnames = list(paste0("Mech_", 1:3), "subunit_1")
    ),
    target.subunits = matrix(
      c("GENE2", "GENE4", "GENE6"),
      ncol = 1,
      dimnames = list(paste0("Mech_", 1:3), "subunit_1")
    )
  )

  input_genes <- paste0("GENE", 1:10)

  filtered_gt <- FilterGroundTruth(
    ground.truth = gt,
    input_rownames = input_genes
  )

  # GENE99 is not in input, so Mech_3 should be removed
  expect_equal(nrow(filtered_gt$source.subunits), 2)
  expect_true("Mech_1" %in% rownames(filtered_gt$source.subunits))
  expect_true("Mech_2" %in% rownames(filtered_gt$source.subunits))
  expect_false("Mech_3" %in% rownames(filtered_gt$source.subunits))
})


test_that("LoadFantom5 works for different species", {
  skip_if_not_installed("NICHES")

  # Test each species
  for (species in c("human", "mouse", "rat", "pig")) {
    gt <- tryCatch({
      LoadFantom5(species = species)
    }, error = function(e) NULL)

    # Skip if data not available
    if (is.null(gt)) {
      skip(paste("FANTOM5 data not available for", species))
    }

    expect_type(gt, "list")
    expect_true("source.subunits" %in% names(gt))
    expect_true("target.subunits" %in% names(gt))
    expect_gt(nrow(gt$source.subunits), 0)
  }
})
