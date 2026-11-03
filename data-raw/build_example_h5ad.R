# build_example_h5ad.R
#
# PURPOSE: Writes inst/extdata/NICHESv2_example.h5ad, a small AnnData file built
#          from the synthetic NICHESv2_inputs dataset. It backs the example in
#          ?extract_NICHESInputs_AnnData. Reuses the spec-compliant writer from
#          the AnnData test file so the two cannot drift apart.
#
# USAGE:   Run from the package root:
#          Rscript data-raw/build_example_h5ad.R

library(Matrix)  # needed for dim() on the dgCMatrix inside the writer

load("data/NICHESv2_inputs.rda")

test.exprs <- parse("tests/testthat/test-extract_NICHESInputs_AnnData.R")
for (e in test.exprs) {
  if (is.call(e) && identical(as.character(e[[2L]]), ".write_test_h5ad")) eval(e)
}

out <- "inst/extdata/NICHESv2_example.h5ad"
dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
unlink(out)

.write_test_h5ad(out, NICHESv2_inputs$count.mtx, NICHESv2_inputs$meta.data)
rhdf5::h5closeAll()

cat(out, ":", file.size(out), "bytes\n")
