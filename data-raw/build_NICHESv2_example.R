# ==============================================================================
# FILE:    data-raw/build_NICHESv2_example.R
# PURPOSE: Builds the user-facing example datasets bundled with NICHESv2:
#            data/NICHESv2_inputs.rda   -- raw inputs (list of three)
#            data/NICHESv2_example.rda  -- prebuilt NICHESObject
#
# USAGE:   Run from package root after any change to the bundled data spec:
#            source("data-raw/build_NICHESv2_example.R")
#          Then commit the resulting data/*.rda files. End users access via:
#            data(NICHESv2_inputs)
#            data(NICHESv2_example)
#
# This file is excluded from the package build via .Rbuildignore.
# Deterministic seed ensures bundled data are reproducible.
# ==============================================================================

suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
  library(NICHESv2)
})

set.seed(42)

# -- Parameters (matches inst/testdata/generate_testdata.R) -------------------
N.CELLS <- 80L
N.GENES <- 50L
DENSITY <- 0.25

# -- Genes ---------------------------------------------------------------------
gene.names <- paste0("Gene", seq_len(N.GENES))
lrm.genes  <- gene.names[1:30]

# -- Count matrix --------------------------------------------------------------
n.nonzero <- round(N.GENES * N.CELLS * DENSITY)
nonzero.i <- sample(N.GENES, n.nonzero, replace = TRUE)
nonzero.j <- sample(N.CELLS, n.nonzero, replace = TRUE)
nonzero.x <- rpois(n.nonzero, lambda = 3) + 1L

count.mtx <- sparseMatrix(
  i        = nonzero.i,
  j        = nonzero.j,
  x        = as.numeric(nonzero.x),
  dims     = c(N.GENES, N.CELLS),
  dimnames = list(gene.names, paste0("Cell_", seq_len(N.CELLS)))
)
count.mtx <- as(count.mtx, "CsparseMatrix")

# -- Cell metadata: 3 spatial clusters ----------------------------------------
cluster.centers <- list(
  TypeA = c(0.25, 0.25),
  TypeB = c(0.75, 0.25),
  TypeC = c(0.50, 0.80)
)
celltypes <- rep(names(cluster.centers), times = c(30L, 25L, 25L))
x.coords  <- numeric(N.CELLS)
y.coords  <- numeric(N.CELLS)
for (ct in names(cluster.centers)) {
  idx <- which(celltypes == ct)
  cx  <- cluster.centers[[ct]]
  x.coords[idx] <- rnorm(length(idx), mean = cx[1], sd = 0.12)
  y.coords[idx] <- rnorm(length(idx), mean = cx[2], sd = 0.12)
}

meta.data <- data.frame(
  x          = x.coords,
  y          = y.coords,
  nCount_RNA = round(colSums(count.mtx)),
  celltype   = celltypes,
  row.names  = paste0("Cell_", seq_len(N.CELLS)),
  stringsAsFactors = FALSE
)

# -- LRM database: 15 pairs (12 simple + 3 with complexes) --------------------
single.ligands   <- lrm.genes[c(1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23)]
single.receptors <- lrm.genes[c(2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24)]

complex.rows <- data.frame(
  ligand   = c(paste(lrm.genes[25], lrm.genes[26], sep = "_"),
               lrm.genes[27],
               paste(lrm.genes[29], lrm.genes[30], sep = "_")),
  receptor = c(lrm.genes[27],
               paste(lrm.genes[28], lrm.genes[29], sep = "_"),
               paste(lrm.genes[28], lrm.genes[30], sep = "_")),
  pathway  = c("PathwayX", "PathwayY", "PathwayZ"),
  stringsAsFactors = FALSE
)

LRM.db <- rbind(
  data.frame(
    ligand   = single.ligands,
    receptor = single.receptors,
    pathway  = paste0("Pathway", LETTERS[1:12]),
    stringsAsFactors = FALSE
  ),
  complex.rows
)

# -- NICHESv2_inputs: the raw triple ------------------------------------------
NICHESv2_inputs <- list(
  count.mtx = count.mtx,
  meta.data = meta.data,
  LRM.db    = LRM.db
)

# -- NICHESv2_example: prebuilt NICHESObject ----------------------------------
# Neighborhood slots intentionally left NULL so add_neighborhoods()
# examples remain meaningful. Segmentation likewise NULL so
# add_segmentation() examples have something to attach.
NICHESv2_example <- create_NICHESObject(
  count.mtx        = count.mtx,
  meta.data        = meta.data,
  LRM.db           = LRM.db,
  mode             = "spatial",
  ligand.col       = "ligand",
  receptor.col     = "receptor",
  k                = 5L,
  method           = "product",
  normalize.method = "prop",
  n.cores          = 1L,
  verbose          = FALSE
)

# -- Save (xz compression keeps payload small) --------------------------------
if (!dir.exists("data")) dir.create("data")
save(NICHESv2_inputs,  file = "data/NICHESv2_inputs.rda",  compress = "xz")
save(NICHESv2_example, file = "data/NICHESv2_example.rda", compress = "xz")

message("[build_NICHESv2_example] Wrote:")
message("  data/NICHESv2_inputs.rda   -- list(count.mtx, meta.data, LRM.db)")
message("  data/NICHESv2_example.rda  -- NICHESObject (",
        length(NICHESv2_example$cell.list), " cells, ",
        length(NICHESv2_example$LRM.list),  " LRMs, ",
        nrow(NICHESv2_example$edge.list),   " edges)")
