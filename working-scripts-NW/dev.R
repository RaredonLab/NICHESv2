# dev.R — developer workflow
#
# HOW TO USE:
#   1. Open NICHESv2.Rproj in RStudio
#   2. Run blocks individually with Ctrl+Enter
#
# FIRST TIME ONLY:  source fixtures line, then the daily loop
# DAILY LOOP:       load_all()  →  test()
# BEFORE COMMIT:    document()  →  check()

library(devtools)

# =============================================================================
# 0. Generate test fixtures — run ONCE before the first test()
# =============================================================================

source("inst/testdata/generate_testdata.R")

# =============================================================================
# 1. Generate NAMESPACE + man/ pages from roxygen2 comments
# =============================================================================

devtools::document()

# =============================================================================
# 2. Load the package into your session (no install needed)
# =============================================================================

devtools::load_all(".")

# =============================================================================
# 3. Run the test suite
# =============================================================================
options(rlang_backtrace_on_error = "none")
devtools::test()

# =============================================================================
# 4. Full R CMD CHECK — run before every commit
# =============================================================================

devtools::check()


niches.obj <- readRDS("~/Desktop/for_Docker/RQ32878-002_Slide_2_PPLR_726-2/Rat-PPLR-Run1/NICHESv2_PPLR_rad25.rds")

cell.type.col <- 'Type.6'

niches.obj <- aggregate_NICHESObject(
  niches.obj,
  cell.type.col = cell.type.col,
  n.cores       = 4L
)
