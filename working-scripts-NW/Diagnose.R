devtools::document()    # regenerate man/ pages
?add_cell_meta          # Examples section should now render
?add_neighborhoods      # Examples section should now render
?Explore_Radius_Neighborhood   # Examples section should now render

# Confirm examples actually run (this is what R CMD check --as-cran does for \donttest)
devtools::run_examples()


devtools::document()
?Aggregate_CellToCell_Edges        # Examples section now executable
?compute_PseudoBulk_Connectivity   # Examples section now executable

# Confirm both examples actually run end-to-end
devtools::run_examples(run_donttest = TRUE)

devtools::document()
?add_segmentation       # @seealso section now appears, example is more pedagogical
?load_LRM_database      # two example blocks visible (one runs, one skipped)

# Confirm \donttest examples still run end-to-end
devtools::run_examples(run_donttest = TRUE)

devtools::document()
?create_NICHESObject     # confirm sections render, examples appear, @return enumerates slots
devtools::run_examples(run_donttest = TRUE)  # all 3 example variants must run clean
devtools::test()         # full suite must still pass (function body unchanged, so should)

devtools::document()
?compute_CellToCell                          # confirm 3 sections render, example appears, FANTOM5 line correct
devtools::run_examples(run_donttest = TRUE)  # the new example must run clean
devtools::test()                             # full suite must still pass

devtools::document()
?CellToCell_InSitu
devtools::run_examples(run_donttest = TRUE)
devtools::test()

devtools::document()
?CellToCell_InSitu_Single_LRM
devtools::run_examples(run_donttest = TRUE)
devtools::test()


devtools::document()
devtools::test()                         # all tests must still pass
devtools::run_examples(run_donttest = TRUE)  # examples for these 4 must still run


# source("inst/testdata/generate_testdata.R")
# source("data-raw/build_NICHESv2_example.R")

# 1. What's in the test fixture?
count.mtx <- readRDS("inst/testdata/count.mtx.rds")
LRM.db    <- readRDS("inst/testdata/LRM.db.rds")
edge.list <- readRDS("inst/testdata/edge.list.rds")

cat("count.mtx dimensions:", dim(count.mtx), "\n")
cat("Sample of rownames(count.mtx):\n"); print(head(rownames(count.mtx), 10))

cat("\nLRM.db dimensions:", dim(LRM.db), "\n")
cat("LRM.db column names:", colnames(LRM.db), "\n")
cat("Sample of LRM.db:\n"); print(head(LRM.db))

# 2. Which test exactly is failing — show the call?

# In package root
devtools::load_all(".")    # Re-load the package fresh
exists(".genes_available", where = asNamespace("NICHESv2"))
# expect: TRUE

NICHESv2:::.genes_available("Gene1", c("Gene1", "Gene2"))
# expect: TRUE

devtools::test(filter = "compute_CellToCell")

devtools::check() 
