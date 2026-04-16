# NICHESv2

[![R-CMD-check](https://github.com/RaredonLab/NICHESv2/actions/workflows/R-CMD-check.yaml/badge.svg?branch=dev)](https://github.com/RaredonLab/NICHESv2/actions/workflows/R-CMD-check.yaml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)

> **Development version (0.0.9000).** This branch is shared with lab members
> and collaborators for pre-release testing. A formal 0.1.0 release will
> accompany the methods paper. Install instructions are below.

Spatial cell-to-cell ligand-receptor communication scoring for single-cell
transcriptomics data.

---

## Overview

NICHESv2 computes ligand-receptor (LR) signalling scores between pairs of
cells and organises results in a memory-efficient, long-format S3 object
(`NICHESObject`). Three edge-definition modes are supported — spatial
nearest-neighbour, type-stratified sampling, and pseudo-bulk — so the same
downstream workflow applies to both spatially resolved and conventional
scRNA-seq datasets.

Expression data can be supplied as a raw sparse matrix, a Seurat V4/V5
object, or an AnnData (`.h5ad`) file. Built-in LR database loaders cover
FANTOM5 (bundled, four species, no internet required) and OmniPath, and
any user-supplied data frame is also accepted. Multi-gene receptor complexes
are supported via underscore delimiting (`gene1_gene2`).

---

## Installation

```r
# 1. Install CRAN dependencies first (avoids timeouts on GitHub install)
install.packages(c(
  "data.table", "Matrix", "RANN", "dbscan", "irlba", "ggplot2", "scales"
))

# 2. Install NICHESv2 from the dev branch
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
remotes::install_github("RaredonLab/NICHESv2", ref = "dev")
```

**Optional — install only what you need:**

```r
# Seurat V4 / V5 input support
install.packages("SeuratObject")

# AnnData (.h5ad) input support — no Python required
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("rhdf5")

# OmniPath LR database
BiocManager::install("OmnipathR")

# Segmentation polygon support (sf geometries)
install.packages("sf")
```

---

## Quick start

### Explore the bundled example

```r
library(NICHESv2)

# Pre-built NICHESObject: 80 synthetic cells, spatial mode, 15 LR pairs
data(NICHESv2_example)

print(NICHESv2_example)          # slot summary and object attributes
dim(NICHESv2_example)            # cells x LRMs

# Aggregate per-edge scores to the cell level
agg <- Aggregate_CellToCell_Edges(NICHESv2_example, stat = "mean")
head(agg$sending)                # cell x LRM matrix, sending perspective
head(agg$receiving)              # cell x LRM matrix, receiving perspective

# Visualise one LR pair in situ
CellToCell_InSitu_Single_LRM(
  niches.obj = NICHESv2_example,
  LRM        = NICHESv2_example$LRM.list[1]
)
```

### Build a NICHESObject from raw inputs

```r
# The same inputs used to construct NICHESv2_example
data(NICHESv2_inputs)   # named list: $count.mtx, $meta.data, $LRM.db

obj <- create_NICHESObject(
  count.mtx        = NICHESv2_inputs$count.mtx,
  meta.data        = NICHESv2_inputs$meta.data,
  LRM.db           = NICHESv2_inputs$LRM.db,
  mode             = "spatial",
  ligand.col       = "ligand",
  receptor.col     = "receptor",
  k                = 5L,           # 5 nearest spatial neighbours per cell
  method           = "product",    # score = ligand expression x receptor expression
  normalize.method = "prop",
  n.cores          = 2L
)

print(obj)
```

For Seurat and AnnData input workflows, vignettes are in preparation.

---

## Features

- **Three edge modes:** spatial k-NN or radius-based neighbours (including
  autocrine self-loops), type-stratified random sampling, and pseudo-bulk
  all-type-pairs
- **LR databases:** FANTOM5 bundled for human, mouse, rat, and pig (no
  internet required); OmniPath via `OmnipathR`; user-supplied data frames
  accepted; multi-gene complexes supported via underscore delimiting
- **Memory-efficient storage:** scores held in long-format sparse
  `data.table`s inside the 14-slot `NICHESObject` S3 class — footprint
  scales with dataset density, not total dimensions
- **Multi-sample support:** `sample.col` argument in `create_NICHESObject()`
  processes samples independently and merges results; barcode uniqueness
  enforced as a hard error
- **Flexible data ingest:** raw sparse matrix, Seurat V4/V5
  (`extract_NICHESInputs_Seurat()`), or AnnData via `rhdf5` only —
  no Python, conda, or basilisk required
  (`extract_NICHESInputs_AnnData()`)
- **Cross-platform parallelism:** socket clusters only
  (`parallel::makeCluster` + `parLapply`); never `mclapply`
- **In-situ visualisation:** `CellToCell_InSitu()` exports batch PNGs for
  any list of LR pairs; `CellToCell_InSitu_Single_LRM()` for interactive
  single-pair inspection
- **Neighborhood analysis:** `add_neighborhoods()` populates per-cell edge
  and cell-neighbour lists; three `edge.filter.mode` options control which
  edges are included in each neighborhood
- **Aggregation:** `Aggregate_CellToCell_Edges()` returns zero-filled full
  cell × LRM grids from both sending and receiving perspectives

---

## Citation

NICHESv2 is described in a manuscript currently in preparation. In the
meantime, please cite the original NICHES method. BibTeX will be added here
upon publication.

```bibtex
@article{NICHESv2_preprint,
  title   = {{NICHESv2}: memory-efficient spatial cell-to-cell
             ligand-receptor communication scoring for
             single-cell transcriptomics},
  author  = {Wang, Nuoya and Raredon, {Micha Sam Brickman}},
  journal = {},
  year    = {2026},
  doi     = {}
}
```

---

## License

MIT © 2026 Nuoya Wang, Micha Sam Brickman Raredon
