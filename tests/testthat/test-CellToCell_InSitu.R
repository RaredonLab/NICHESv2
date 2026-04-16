# ==============================================================================
# FILE:    tests/testthat/test-CellToCell_InSitu.R
# PURPOSE: Unit tests for CellToCell_InSitu() and CellToCell_InSitu_Single_LRM().
#          Uses data(NICHESv2_example) for the full NICHESObject — loaded
#          inline inside each block to preserve block independence.
#          Blocks that write to disk use withr::local_tempdir() for cleanup.
# ==============================================================================

# ── CellToCell_InSitu_Single_LRM: basic output ────────────────────────────────

test_that("CellToCell_InSitu_Single_LRM: returns a ggplot", {
  data(NICHESv2_example)
  lrm <- NICHESv2_example$LRM.list[1L]
  p   <- CellToCell_InSitu_Single_LRM(NICHESv2_example, LRM = lrm)
  expect_s3_class(p, "ggplot")
})

test_that("CellToCell_InSitu_Single_LRM: plot title defaults to '<LRM> Connectivity'", {
  data(NICHESv2_example)
  lrm <- NICHESv2_example$LRM.list[1L]
  p   <- CellToCell_InSitu_Single_LRM(NICHESv2_example, LRM = lrm)
  expect_equal(p$labels$title, paste(lrm, "Connectivity"))
})

test_that("CellToCell_InSitu_Single_LRM: custom title is applied", {
  data(NICHESv2_example)
  lrm <- NICHESv2_example$LRM.list[1L]
  p   <- CellToCell_InSitu_Single_LRM(NICHESv2_example, LRM = lrm,
                                      title = "My Title")
  expect_equal(p$labels$title, "My Title")
})

# ── CellToCell_InSitu_Single_LRM: input validation ────────────────────────────

test_that("CellToCell_InSitu_Single_LRM: errors on wrong LRM format (no '|')", {
  data(NICHESv2_example)
  expect_error(
    CellToCell_InSitu_Single_LRM(NICHESv2_example, LRM = "MIF.CD74"),
    regexp = "ligand\\|receptor"
  )
})

test_that("CellToCell_InSitu_Single_LRM: errors on LRM not in $LRM.list", {
  data(NICHESv2_example)
  expect_error(
    CellToCell_InSitu_Single_LRM(NICHESv2_example, LRM = "NOTREAL|GENE"),
    regexp = "not found"
  )
})

test_that("CellToCell_InSitu_Single_LRM: errors on non-NICHESObject input", {
  expect_error(
    CellToCell_InSitu_Single_LRM(list(), LRM = "A|B"),
    regexp = "NICHESObject"
  )
})

# ── CellToCell_InSitu_Single_LRM: value.col ───────────────────────────────────

test_that("CellToCell_InSitu_Single_LRM: value.col = 'score.norm' returns a ggplot", {
  data(NICHESv2_example)
  lrm <- NICHESv2_example$LRM.list[1L]
  p   <- CellToCell_InSitu_Single_LRM(NICHESv2_example, LRM = lrm,
                                      value.col = "score.norm")
  expect_s3_class(p, "ggplot")
})

test_that("CellToCell_InSitu_Single_LRM: errors on invalid value.col", {
  data(NICHESv2_example)
  lrm <- NICHESv2_example$LRM.list[1L]
  expect_error(
    CellToCell_InSitu_Single_LRM(NICHESv2_example, LRM = lrm,
                                 value.col = "raw.count"),
    regexp = "value.col"
  )
})

# ── CellToCell_InSitu_Single_LRM: scale.upper ─────────────────────────────────

test_that("CellToCell_InSitu_Single_LRM: user-supplied scale.upper returns a ggplot", {
  data(NICHESv2_example)
  lrm <- NICHESv2_example$LRM.list[1L]
  p   <- CellToCell_InSitu_Single_LRM(NICHESv2_example, LRM = lrm,
                                      scale.upper = 5)
  expect_s3_class(p, "ggplot")
})

test_that("CellToCell_InSitu_Single_LRM: errors on non-positive scale.upper", {
  data(NICHESv2_example)
  lrm <- NICHESv2_example$LRM.list[1L]
  expect_error(
    CellToCell_InSitu_Single_LRM(NICHESv2_example, LRM = lrm,
                                 scale.upper = -1),
    regexp = "scale.upper"
  )
})

test_that("CellToCell_InSitu_Single_LRM: errors on non-numeric scale.upper", {
  data(NICHESv2_example)
  lrm <- NICHESv2_example$LRM.list[1L]
  expect_error(
    CellToCell_InSitu_Single_LRM(NICHESv2_example, LRM = lrm,
                                 scale.upper = "high"),
    regexp = "scale.upper"
  )
})

test_that("CellToCell_InSitu_Single_LRM: NULL scale.upper auto-computes without error", {
  data(NICHESv2_example)
  lrm <- NICHESv2_example$LRM.list[1L]
  p   <- CellToCell_InSitu_Single_LRM(NICHESv2_example, LRM = lrm,
                                      scale.upper = NULL)
  expect_s3_class(p, "ggplot")
})

# ── CellToCell_InSitu: batch function ─────────────────────────────────────────

test_that("CellToCell_InSitu: saves PNGs and returns summary data.frame", {
  data(NICHESv2_example)
  tmp.dir <- withr::local_tempdir()
  result  <- CellToCell_InSitu(
    niches.obj = NICHESv2_example,
    LRM.list   = NICHESv2_example$LRM.list[1:2],
    folder     = tmp.dir,
    dpi        = 72L
  )
  expect_s3_class(result, "data.frame")
  expect_true(all(c("LRM", "status", "path") %in% colnames(result)))
  expect_equal(nrow(result), 2L)
  expect_true(all(result$status == "saved"))
  expect_true(all(file.exists(result$path)))
})

test_that("CellToCell_InSitu: saves to 'selected_LRMs' subfolder", {
  data(NICHESv2_example)
  tmp.dir <- withr::local_tempdir()
  result  <- CellToCell_InSitu(
    niches.obj = NICHESv2_example,
    LRM.list   = NICHESv2_example$LRM.list[1L],
    folder     = tmp.dir,
    dpi        = 72L
  )
  expect_true(all(grepl("selected_LRMs", result$path, fixed = TRUE)))
})

test_that("CellToCell_InSitu: NULL LRM.list plots all pairs", {
  data(NICHESv2_example)
  n.lrm <- length(NICHESv2_example$LRM.list)
  skip_if(n.lrm > 5L, "Too many LRM pairs for a full-run InSitu test.")
  
  tmp.dir <- withr::local_tempdir()
  result  <- CellToCell_InSitu(NICHESv2_example, LRM.list = NULL,
                               folder = tmp.dir, dpi = 72L)
  expect_equal(nrow(result), n.lrm)
})

test_that("CellToCell_InSitu: value.col = 'score.norm' saves PNG without error", {
  data(NICHESv2_example)
  tmp.dir <- withr::local_tempdir()
  result  <- CellToCell_InSitu(
    niches.obj = NICHESv2_example,
    LRM.list   = NICHESv2_example$LRM.list[1L],
    folder     = tmp.dir,
    value.col  = "score.norm",
    dpi        = 72L
  )
  expect_true(result$status[1L] == "saved")
  expect_true(file.exists(result$path[1L]))
})

test_that("CellToCell_InSitu: user-supplied scale.upper saves PNG without error", {
  data(NICHESv2_example)
  tmp.dir <- withr::local_tempdir()
  result  <- CellToCell_InSitu(
    niches.obj  = NICHESv2_example,
    LRM.list    = NICHESv2_example$LRM.list[1L],
    folder      = tmp.dir,
    scale.upper = 5,
    dpi         = 72L
  )
  expect_true(result$status[1L] == "saved")
  expect_true(file.exists(result$path[1L]))
})