# NICHESv2 <img src="man/figures/logo.png" align="right" height="120" alt="NICHESv2 logo placeholder"/>

> **Spatial Cell-to-Cell Ligand-Receptor Communication Analysis**

NICHESv2 computes ligand-receptor (LR) signalling scores between spatially
neighbouring cells in single-cell spatial transcriptomics data.  It returns a
**`NICHESObject`** — a 13-slot S3 class — with sparse per-cell and per-edge
expression tables, flexible aggregation, and batch in situ visualisation.

---

## Key features

| Feature | Details |
|---------|---------|
| **Multi-gene complexes** | `"COL1A1_COL1A2"` style ligands / receptors handled automatically |
| **Built-in LR databases** | `"omnipath"` and `"fantom5"` loaded by name |
| **Sparse output** | Only non-zero signals stored; scales to 100k+ cells |
| **Per-edge normalisation** | `score.norm` normalised within each directed edge |
| **Parallel computation** | LR chunks dispatched via `mclapply` on Unix/macOS |
| **Flexible aggregation** | `Aggregate_CellToCell_Edges()` produces full cell × LRM grids |
| **In situ visualisation** | Per-LRM segment plots, dark theme, viridis colour scale |

---

## Installation

```r
# Install from local source
devtools::install("path/to/NICHESv2")

# Dependencies (installed automatically)
# data.table, Matrix, RANN, dbscan, ggplot2

# Optional: built-in LR databases (requires Bioconductor)
BiocManager::install("OmnipathR")
```

---

## Quick start

```r
library(NICHESv2)

# ── Inputs ──────────────────────────────────────────────────────────────────
# count.mtx : gene × cell sparse matrix (dgCMatrix or matrix)
#             rownames = gene symbols, colnames = cell barcodes
# meta.data : data.frame with rownames = cell barcodes
#             MUST contain columns "x" and "y" (spatial coordinates)

# ── 1. Build the NICHESObject ────────────────────────────────────────────────
niches.obj <- create_NICHESObject(
  count.mtx = count.mtx,
  meta.data = meta.data,     # must have x, y columns
  LRM.db    = "omnipath",    # or "fantom5" or your own data.frame
  k         = 6L,            # 6-nearest-neighbour spatial graph
  method    = "product"      # or "minimum"
)

# ── 2. Inspect ───────────────────────────────────────────────────────────────
print(niches.obj)
summary(niches.obj)
dim(niches.obj)       # cells / genes / edges / LRMs

# ── 3. Aggregate per-cell sending / receiving scores ─────────────────────────
agg <- Aggregate_CellToCell_Edges(niches.obj, stat = "mean")
head(agg$sending)
head(agg$receiving)

# ── 4. Visualise a single LRM pair in situ ────────────────────────────────────
p <- CellToCell_InSitu_Single_LRM(niches.obj, LRM = "MIF|CD74")
print(p)

# ── 5. Batch-export PNGs for all LRM pairs ────────────────────────────────────
CellToCell_InSitu(niches.obj, folder = "niches_plots", n.cores = 4L)

# ── 6. Subset to a cell population ────────────────────────────────────────────
b.cells <- niches.obj$cell.meta[celltype == "B_cell", cell]
sub.obj  <- niches.obj[b.cells]

# ── 7. Add cell type annotations after the fact ───────────────────────────────
new.annot <- data.frame(cell = my.cells, cluster = my.clusters)
niches.obj <- add_cell_meta(niches.obj, new.annot)
```

---

## Working with a Seurat object

```r
library(Seurat)
library(NICHESv2)

# Extract count matrix
count.mtx <- GetAssayData(seu, assay = "Spatial", slot = "counts")

# Extract metadata — add spatial coordinates if stored separately
meta.data        <- seu@meta.data
meta.data$x      <- GetTissueCoordinates(seu)[, "imagerow"]
meta.data$y      <- GetTissueCoordinates(seu)[, "imagecol"]

niches.obj <- create_NICHESObject(
  count.mtx = count.mtx,
  meta.data = meta.data,
  LRM.db    = "omnipath",
  k         = 6L
)
```

---

## NICHESObject slot reference

| # | Slot | Type | Contents |
|---|------|------|----------|
| 1 | `$cell.list` | `character` | Cell barcodes |
| 2 | `$cell.data` | `data.table` | `(cell, gene, score, score.norm)` — sparse |
| 3 | `$cell.meta` | `data.table` | Per-cell metadata: `x`, `y`, `celltype`, … |
| 4 | `$edge.list` | `data.table` | `(edge, SendingCell, ReceivingCell)` |
| 5 | `$edge.data` | `data.table` | `(edge, LRM, score, score.norm)` — sparse |
| 6 | `$edge.meta` | `data.table` | Per-edge metadata with `Sending.*`, `Receiving.*`, `Joint.*` cols |
| 7–9 | `$neighborhood.*` | reserved | `NULL` — populated externally |
| 10 | `$gene.list` | `character` | Gene names |
| 11 | `$gene.meta` | `data.table` | Per-gene annotation |
| 12 | `$LRM.list` | `character` | `"ligand\|receptor"` strings |
| 13 | `$LRM.meta` | `data.table` | `LRM`, `ligand`, `receptor`, pathway columns, … |

---

## Key conventions

| Convention | Rule |
|-----------|------|
| Function names | Underscore separator: `create_NICHESObject()` |
| Variable / argument names | Dot separator: `n.cores`, `bg.color` |
| Edge IDs | `"SendingCell\|ReceivingCell"` |
| LRM IDs | `"ligand\|receptor"` |
| Gene complexes | `"gene1_gene2"` (underscore in database) |
| Score columns | `score` (raw), `score.norm` (normalised) |

---

## Running tests

```r
# 1. Generate fixtures (once)
source("inst/testdata/generate_testdata.R")

# 2. Run the full test suite
devtools::test()

# 3. Full R CMD check
devtools::check()
```

---

## Citation

If you use NICHESv2 in published work, please cite the original NICHES paper
and the LR database you chose (FANTOM5 or OmniPath):

> Ramilowski J.A. et al. (2015). A draft network of ligand-receptor-mediated
> multicellular signalling in human. *Nature Communications*, 6, 7866.

> Türei D. et al. (2021). Integrated intra- and intercellular signaling
> knowledge for multicellular omics analysis. *Molecular Systems Biology*,
> 17(3), e9923.

---

## License

MIT © 2026 NICHESv2 Authors
