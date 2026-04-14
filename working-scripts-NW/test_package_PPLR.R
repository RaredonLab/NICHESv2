library(Seurat)
PPLR.all <- readRDS("~/Desktop/TEMP/PPLR.all.integrated.ann.updated.rds")
sample.data <- subset(PPLR.all, subset = Sample.global == 'PNX7_M1_Sp')

lrm.df <- read.csv("~/Desktop/LR_DB_Rat.csv")
LRM.db <- load_LRM_database(lrm.df, ligand.col = "Ligand", receptor.col = "Receptor")

count.mtx <- sample.data@assays$Xenium$counts
meta.data <- sample.data@meta.data

niches.obj <- create_NICHESObject(
  count.mtx        = count.mtx,
  meta.data        = meta.data,
  LRM.db           = "fantom5",       # "omnipath", "fantom5", or a data.frame
  k                = NULL,
  rad              = 25,
  method           = NULL,
  normalize.method = "prop",
  n.cores          = 4,
  verbose          = TRUE
)



