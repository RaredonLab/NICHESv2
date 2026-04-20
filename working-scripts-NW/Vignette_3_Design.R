SeuratData::InstallData("stxBrain")   # ~150 MB, one-time
library(SeuratData)
library(Seurat)
data("anterior1")                      # the anterior section object

sink('./working-scripts-NW/check_stxBrain.txt')
anterior1
colnames(anterior1@meta.data)
head(anterior1@meta.data, 3)
inherits(anterior1@assays[["Spatial"]], "Assay5")
sink()

sink('./working-scripts-NW/check_stxBrain.txt')
tmp <- anterior1@images[[1]]@boundaries$centroids@coords
head(tmp, 3)
colnames(tmp)

devtools::load_all()
inputs.stx <- extract_NICHESInputs_Seurat(
  anterior1, assay = "Spatial", layer = "counts", verbose = TRUE
)
colnames(inputs.stx$meta.data)
head(inputs.stx$meta.data, 3)

sink()


# Retrieve barcodes and attach to coordinate matrix
coords <- anterior1@images$anterior1@boundaries$centroids@coords
barcodes <- anterior1@images$anterior1@boundaries$centroids@cells
rownames(coords) <- barcodes

# Verify alignment with meta.data
all(barcodes %in% rownames(inputs.stx$meta.data))   # should be TRUE
sum(is.na(coords[rownames(inputs.stx$meta.data), "x"]))  # should be 0

# Safe named join
inputs.stx$meta.data$x <- coords[rownames(inputs.stx$meta.data), "x"]
inputs.stx$meta.data$y <- coords[rownames(inputs.stx$meta.data), "y"]

# Spot-check
head(inputs.stx$meta.data[, c("x", "y")], 3)



# Confirm no NAs
sum(is.na(inputs.stx$meta.data$x))

# Then create the object
niches.stx <- create_NICHESObject(
  count.mtx  = inputs.stx$count.mtx,
  meta.data  = inputs.stx$meta.data,
  LRM.db     = load_LRM_database("fantom5", species = "mouse"),
  mode       = "spatial",
  species    = "mouse",
  sample.col = NULL,
  k          = 6L,          # start with kNN=6, we will explore radius after
  method     = "product",
  n.cores    = 4L,
  verbose    = TRUE
)

sink('./working-scripts-NW/check_stxBrain_niches.txt')
print(niches.stx)
head(niches.stx$edge.data, 3)
sink()

# Top 20 LRM pairs by total score
top.lrm <- niches.stx$edge.data[, .(total.score = sum(score)), by = LRM
][order(-total.score)][1:20]
print(top.lrm)


CellToCell_InSitu_Single_LRM(niches.obj = niches.stx, LRM = top.lrm$LRM[1])
CellToCell_InSitu_Single_LRM(niches.obj = niches.stx, LRM = top.lrm$LRM[2])
CellToCell_InSitu_Single_LRM(niches.obj = niches.stx, LRM = top.lrm$LRM[3])
CellToCell_InSitu_Single_LRM(niches.obj = niches.stx, LRM = top.lrm$LRM[4])

CellToCell_InSitu_Single_LRM(
  niches.obj    = niches.stx,
  LRM           = "Gnas|Adcy1",
  bg.alpha      = 0.3,
  bg.linewidth  = 0.2,
  sig.linewidth = 0.3,
  sig.alpha = 0.7
)


coords.df <- data.frame(
  x = inputs.stx$meta.data$x,
  y = inputs.stx$meta.data$y
)

Explore_Radius_Neighborhood(coords.df, r.max = 800)
