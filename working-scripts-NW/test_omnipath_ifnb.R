# test_omnipath_ifnb.R ---------------------------------------------------------
# End-to-end test of the OmniPath resilient-loading fix (NICHESv1 PR #72 /
# issue #70, ported into NICHESv2's load_LRM_database.R).
#
# Goal: confirm that create_NICHESObject() runs cleanly with LRM.db = "omnipath"
# on a real, public Seurat dataset (ifnb interferon-beta-stimulated human PBMCs),
# exercising the full pipeline: OmniPath fetch -> edge list -> LR scoring ->
# NICHESObject assembly -> neighborhood extraction -> aggregation.
#
# The dataset is non-spatial, so we use mode = "nonspatial" (v1-style stratified
# random sampling), which requires a cell.type.col.
#
# Run from the package root:
#   Rscript working-scripts-NW/test_omnipath_ifnb.R
# ------------------------------------------------------------------------------

# -- 0. Dependencies -----------------------------------------------------------
# Seurat / SeuratData provide the public ifnb dataset; OmnipathR is needed for
# the *primary* OmniPath path (the REST fallback is base-R only).
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratData)
})

# Load the package under test (dev version, straight from source).
devtools::load_all(quiet = TRUE)

set.seed(42)

# -- 1. Get the ifnb dataset ---------------------------------------------------
if (!"ifnb" %in% SeuratData::AvailableData()[, "Dataset"] ||
    is.na(SeuratData::AvailableData()["ifnb.SeuratData", "Installed"]) ||
    !isTRUE(SeuratData::AvailableData()["ifnb.SeuratData", "Installed"])) {
  message("Installing ifnb dataset via SeuratData ...")
  SeuratData::InstallData("ifnb")
}

ifnb <- SeuratData::LoadData("ifnb")
ifnb <- UpdateSeuratObject(ifnb)   # harmless if already current; normalizes V4/V5

# ifnb metadata carries:
#   - seurat_annotations : cell-type labels (used as cell.type.col)
#   - stim               : CTRL vs STIM condition
cat("\n--- ifnb meta.data columns ---\n")
print(colnames(ifnb@meta.data))
cat("\n--- cell-type counts (seurat_annotations) ---\n")
print(table(ifnb$seurat_annotations, useNA = "ifany"))

# Optional: subset to keep the test fast. Comment out to run on the full object.
n.sub <- 2000L
if (ncol(ifnb) > n.sub) {
  ifnb <- ifnb[, sample(colnames(ifnb), n.sub)]
  message(sprintf("Subsetted to %d cells for a fast smoke test.", ncol(ifnb)))
}

# Normalize so LR scoring runs on expression-scaled values (v1 convention).
ifnb <- NormalizeData(ifnb, verbose = FALSE)

# -- 2. Extract NICHES inputs --------------------------------------------------
# Pull the normalized "data" layer; extract does NOT add coordinates (fine for
# non-spatial mode).
inputs <- extract_NICHESInputs_Seurat(ifnb, assay = "RNA", layer = "data")

# Drop cells with an NA cell-type label (nonspatial mode stratifies on it).
keep <- !is.na(inputs$meta.data$seurat_annotations)
inputs$count.mtx <- inputs$count.mtx[, keep, drop = FALSE]
inputs$meta.data <- inputs$meta.data[keep, , drop = FALSE]

# -- 3. Load OmniPath directly (isolates the fix before the full pipeline) -----
cat("\n========== STEP A: load_LRM_database('omnipath') ==========\n")
op.db <- load_LRM_database("omnipath", species = "human", verbose = TRUE)
cat(sprintf("OmniPath LR pairs loaded: %d\n", nrow(op.db)))
stopifnot(nrow(op.db) > 0,
          all(c("ligand", "receptor") %in% colnames(op.db)))

# -- 4. Full pipeline with LRM.db = "omnipath" ---------------------------------
cat("\n========== STEP B: create_NICHESObject(LRM.db = 'omnipath') ==========\n")
niches <- create_NICHESObject(
  count.mtx     = inputs$count.mtx,
  meta.data     = inputs$meta.data,
  LRM.db        = "omnipath",
  mode          = "nonspatial",
  species       = "human",
  cell.type.col = "seurat_annotations",
  max.cells     = 100L,      # cap edges per cell-type crossing (nonspatial)
  method        = "product",
  n.cores       = 1L,
  verbose       = TRUE
)

# -- 5. Sanity checks on the resulting object ----------------------------------
cat("\n========== STEP C: sanity checks ==========\n")
cat("class            : ", paste(class(niches), collapse = ", "), "\n")
cat("mode attr        : ", attr(niches, "mode"), "\n")
cat("n cells          : ", length(niches$cell.list), "\n")
cat("n edges          : ", length(niches$edge.list), "\n")
cat("n LRMs           : ", length(niches$LRM.list), "\n")
cat("edge.data rows   : ", nrow(niches$edge.data), "\n")

stopifnot(
  inherits(niches, "NICHESObject"),
  length(niches$edge.list) > 0,
  length(niches$LRM.list)  > 0,
  nrow(niches$edge.data)   > 0
)

# -- 6. Aggregation (optional; exercises the rest of the pipeline) -------------
# Non-fatal: the OmniPath fix is already proven by STEP B/C. This just confirms
# downstream aggregation runs on an OmniPath-scored object.
cat("\n========== STEP D: aggregate_NICHESObject() (optional) ==========\n")
tryCatch({
  niches <- aggregate_NICHESObject(niches,
                                   cell.type.col = "seurat_annotations",
                                   verbose = TRUE)
  cat("aggregation slots: ", paste(names(niches$aggregations), collapse = ", "), "\n")
}, error = function(e) {
  cat("STEP D skipped (non-fatal):", conditionMessage(e), "\n")
})

cat("\n\n*** OmniPath end-to-end test PASSED (STEPS A-C) ***\n")
