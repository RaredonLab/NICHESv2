# ==============================================================================
# FILE:    tests/testthat.R
# PURPOSE: Standard testthat entry point. Run via devtools::test() or
#          R CMD check. Generates test fixtures on first run if absent.
# ==============================================================================

library(testthat)
library(NICHESv2)

# Auto-generate fixtures if not yet present (convenient for first-time runs)
fixture.path <- system.file("testdata", "niches.obj.rds",
                             package = "NICHESv2", mustWork = FALSE)
if (!nchar(fixture.path) || !file.exists(fixture.path)) {
  gen.script <- system.file("testdata", "generate_testdata.R",
                             package = "NICHESv2", mustWork = FALSE)
  if (nchar(gen.script) && file.exists(gen.script)) {
    message("NICHESv2 tests: generating test fixtures ...")
    source(gen.script)
  }
}

test_check("NICHESv2")
