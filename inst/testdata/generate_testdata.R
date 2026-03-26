# ==============================================================================
# FILE:    inst/testdata/generate_testdata.R
# PURPOSE: Generates all .rds test fixtures used by tests/testthat/*.
#          Run this script once before running devtools::test():
#
#            source("inst/testdata/generate_testdata.R")
#
# OUTPUTS (written to inst/testdata/):
#   count.mtx.rds    — 50-gene × 80-cell sparse dgCMatrix
#   meta.data.rds    — 80-row data.frame: x, y, nCount_RNA, celltype
#   LRM.db.rds       — 15-pair LRM data.frame (3 complexes)
#   edge.list.rds    — create_Edgelist_Spatial() output (k = 5)
#   niches.obj.rds   — NICHESObject (full pipeline)
# ==============================================================================

suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
})

# Ensure NICHESv2 functions are available (works both from source and installed)
if (!exists("create_NICHESObject", mode = "function")) {
  pkg.root <- normalizePath(file.path(dirname(sys.frame(1)$ofile), "../.."),
                            mustWork = FALSE)
  invisible(lapply(list.files(file.path(pkg.root, "R"), full.names = TRUE), source))
}

set.seed(42)
out.dir <- dirname(sys.frame(1)$ofile)
if (!dir.exists(out.dir)) out.dir <- "inst/testdata"

# ── Parameters ─────────────────────────────────────────────────────────────────
N.CELLS   <- 80L
N.GENES   <- 50L
N.LRM     <- 15L     # LRM pairs in the database
DENSITY   <- 0.25    # fraction of non-zero entries in count matrix

# ── Gene names (include enough for complex pairs) ──────────────────────────────
gene.names <- paste0("Gene", seq_len(N.GENES))

# Designate some genes for LRM pairs (first 30 reserved for LRM use)
lrm.genes  <- gene.names[1:30]

# ── Count matrix: sparse dgCMatrix ─────────────────────────────────────────────
n.nonzero <- round(N.GENES * N.CELLS * DENSITY)
nonzero.i <- sample(N.GENES,   n.nonzero, replace = TRUE)
nonzero.j <- sample(N.CELLS,   n.nonzero, replace = TRUE)
nonzero.x <- rpois(n.nonzero, lambda = 3) + 1L

count.mtx <- sparseMatrix(
  i        = nonzero.i,
  j        = nonzero.j,
  x        = as.numeric(nonzero.x),
  dims     = c(N.GENES, N.CELLS),
  dimnames = list(
    gene.names,
    paste0("Cell_", seq_len(N.CELLS))
  )
)
count.mtx <- as(count.mtx, "CsparseMatrix")

# ── Cell metadata ───────────────────────────────────────────────────────────────
# Simulate 3 spatial clusters
cluster.centers <- list(
  TypeA = c(0.25, 0.25),
  TypeB = c(0.75, 0.25),
  TypeC = c(0.50, 0.80)
)
celltypes <- rep(names(cluster.centers), times = c(30L, 25L, 25L))
x.coords  <- numeric(N.CELLS)
y.coords  <- numeric(N.CELLS)
for (ct in names(cluster.centers)) {
  idx        <- which(celltypes == ct)
  cx         <- cluster.centers[[ct]]
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

# ── LRM database (15 pairs, 3 complexes) ───────────────────────────────────────
# Single-gene pairs (12)
single.ligands   <- lrm.genes[c(1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23)]
single.receptors <- lrm.genes[c(2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24)]

# Complex pairs (3): ligand or receptor is a complex
complex.rows <- data.frame(
  ligand   = c(paste(lrm.genes[25], lrm.genes[26], sep = "_"),  # complex ligand
               lrm.genes[27],                                     # single ligand
               paste(lrm.genes[29], lrm.genes[30], sep = "_")),  # complex ligand
  receptor = c(lrm.genes[27],                                    # single receptor
               paste(lrm.genes[28], lrm.genes[29], sep = "_"),   # complex receptor
               paste(lrm.genes[28], lrm.genes[30], sep = "_")),  # complex receptor
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

# ── Edge list ───────────────────────────────────────────────────────────────────
coord.df  <- data.frame(x = meta.data$x, y = meta.data$y,
                         row.names = rownames(meta.data))
edge.list <- create_Edgelist_Spatial(coord = coord.df, k = 5L)

# ── Full NICHESObject ───────────────────────────────────────────────────────────
niches.obj <- create_NICHESObject(
  count.mtx        = count.mtx,
  meta.data        = meta.data,
  LRM.db           = LRM.db,
  ligand.col       = "ligand",
  receptor.col     = "receptor",
  k                = 5L,
  method           = "product",
  normalize.method = "prop",
  n.cores          = 1L,
  verbose          = TRUE
)

# ── Write fixtures ──────────────────────────────────────────────────────────────
saveRDS(count.mtx,  file.path(out.dir, "count.mtx.rds"))
saveRDS(meta.data,  file.path(out.dir, "meta.data.rds"))
saveRDS(LRM.db,     file.path(out.dir, "LRM.db.rds"))
saveRDS(edge.list,  file.path(out.dir, "edge.list.rds"))
saveRDS(niches.obj, file.path(out.dir, "niches.obj.rds"))

message("\n[generate_testdata] All fixtures saved to: ", out.dir)
message("  count.mtx.rds  — ", nrow(count.mtx),  " genes × ", ncol(count.mtx), " cells")
message("  meta.data.rds  — ", nrow(meta.data),  " cells")
message("  LRM.db.rds     — ", nrow(LRM.db),     " LRM pairs (", sum(grepl("_", LRM.db$ligand) | grepl("_", LRM.db$receptor)), " complexes)")
message("  edge.list.rds  — ", nrow(edge.list),  " edges")
message("  niches.obj.rds — NICHESObject with ",
        length(niches.obj$cell.list), " cells, ",
        length(niches.obj$LRM.list),  " LRMs")
