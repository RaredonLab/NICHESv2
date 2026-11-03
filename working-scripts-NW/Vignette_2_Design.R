library(SeuratData)
InstallData("ifnb")
data("ifnb")

sink("ifnb.obs.txt")
    
    # Object overview
    ifnb
    
    # How many cells, what conditions?
    table(ifnb@meta.data$stim)             # should show CTRL and STIM
    
    # What cell type column is available?
    head(ifnb@meta.data, 3)      # look for the annotation column name
    
    # All metadata column names
    colnames(ifnb@meta.data)
    
    # Seurat version (V4 vs V5 matters for extract_NICHESInputs_Seurat)
    packageVersion("SeuratObject")
    inherits(ifnb@assays[["RNA"]], "Assay5")   # TRUE = V5, FALSE = V4

sink()

devtools::load_all()

inputs <- extract_NICHESInputs_Seurat(ifnb, assay = 'RNA', layer = 'counts')
dim(inputs$count.mtx)
colnames(inputs$meta.data)

# PLACEHOLDER — swap <cell.type.col> with the real column name after Step 2
niches.ifnb <- create_NICHESObject(
  count.mtx     = inputs$count.mtx,
  meta.data     = inputs$meta.data,
  LRM.db        = load_LRM_database("fantom5", species = "human"),
  mode          = "sampled",
  species       = "human",
  sample.col    = "stim",
  cell.type.col = "seurat_annotations",   # ← fill in after Step 2
  method        = "product",
  n.cores       = 4L,
  verbose       = TRUE
)

sink("ifnb.obs.txt")
print(niches.ifnb)
head(niches.ifnb$edge.data, 3)
head(niches.ifnb$edge.meta, 3)
sink()

pb <- compute_PseudoBulk_Connectivity(
  count.mtx     = inputs$count.mtx,
  meta.data     = inputs$meta.data,
  LRM.db        = "fantom5",
  cell.type.col = "seurat_annotations",
  method        = "product",
  n.cores       = 4L,
  verbose       = TRUE
)

# Share these:
colnames(pb)
head(pb, 5)
dim(pb)




devtools::build_vignettes()