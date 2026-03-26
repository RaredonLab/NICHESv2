# ==============================================================================
# FILE:    tests/testthat/test-CellToCell_InSitu.R
# ==============================================================================

test_that("CellToCell_InSitu_Single_LRM: returns a ggplot", {
  .setup_fixtures()
  lrm <- niches.obj$LRM.list[1]
  p   <- CellToCell_InSitu_Single_LRM(niches.obj, LRM = lrm)
  expect_s3_class(p, "ggplot")
})

test_that("CellToCell_InSitu_Single_LRM: plot title defaults to 'LRM Connectivity'", {
  .setup_fixtures()
  lrm <- niches.obj$LRM.list[1]
  p   <- CellToCell_InSitu_Single_LRM(niches.obj, LRM = lrm)
  expect_equal(p$labels$title, paste(lrm, "Connectivity"))
})

test_that("CellToCell_InSitu_Single_LRM: custom title is applied", {
  .setup_fixtures()
  lrm   <- niches.obj$LRM.list[1]
  p     <- CellToCell_InSitu_Single_LRM(niches.obj, LRM = lrm, title = "My Title")
  expect_equal(p$labels$title, "My Title")
})

test_that("CellToCell_InSitu_Single_LRM: errors on wrong LRM format (no '|')", {
  .setup_fixtures()
  expect_error(
    CellToCell_InSitu_Single_LRM(niches.obj, LRM = "MIF.CD74"),
    regexp = "ligand\\|receptor"
  )
})

test_that("CellToCell_InSitu_Single_LRM: errors on LRM not in $LRM.list", {
  .setup_fixtures()
  expect_error(
    CellToCell_InSitu_Single_LRM(niches.obj, LRM = "NOTREAL|GENE"),
    regexp = "not found"
  )
})

test_that("CellToCell_InSitu_Single_LRM: errors on non-NICHESObject input", {
  expect_error(
    CellToCell_InSitu_Single_LRM(list(), LRM = "A|B"),
    regexp = "NICHESObject"
  )
})

test_that("CellToCell_InSitu: saves PNGs and returns summary data.frame", {
  .setup_fixtures()
  tmp.dir <- withr::local_tempdir()

  # Plot first 2 LRMs only
  result <- CellToCell_InSitu(
    niches.obj = niches.obj,
    LRM.list   = niches.obj$LRM.list[1:2],
    folder     = tmp.dir,
    dpi        = 72   # low for speed in tests
  )

  expect_s3_class(result, "data.frame")
  expect_true(all(c("LRM", "status", "path") %in% colnames(result)))
  expect_equal(nrow(result), 2L)
  expect_true(all(result$status == "saved"))
  expect_true(all(file.exists(result$path)))
})

test_that("CellToCell_InSitu: selected mode saves to 'selected_LRMs' subfolder", {
  .setup_fixtures()
  tmp.dir <- withr::local_tempdir()
  result  <- CellToCell_InSitu(
    niches.obj = niches.obj,
    LRM.list   = niches.obj$LRM.list[1],
    folder     = tmp.dir,
    dpi        = 72
  )
  # Path should be inside selected_LRMs/
  expect_true(all(grepl("selected_LRMs", result$path, fixed = TRUE)))
})

test_that("CellToCell_InSitu: NULL LRM.list plots all pairs", {
  .setup_fixtures()
  # Only run if there are ≤5 LRM pairs to keep test fast
  n.lrm <- length(niches.obj$LRM.list)
  skip_if(n.lrm > 5, "Too many LRM pairs for a full-run InSitu test.")

  tmp.dir <- withr::local_tempdir()
  result  <- CellToCell_InSitu(niches.obj, LRM.list = NULL, folder = tmp.dir, dpi = 72)
  expect_equal(nrow(result), n.lrm)
})
