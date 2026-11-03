# NICHESv2

[![R-CMD-check](https://github.com/RaredonLab/NICHESv2/actions/workflows/R-CMD-check.yaml/badge.svg?branch=main)](https://github.com/RaredonLab/NICHESv2/actions/workflows/R-CMD-check.yaml)
[![License: AGPL
v3](https://img.shields.io/badge/License-AGPL_v3-blue.svg)](https://raredonlab.github.io/NICHESv2/LICENSE.md)

NICHESv2 scores ligand-receptor (LR) signaling between individual pairs
of cells in spatial and non-spatial single-cell transcriptomics data. It
keeps every non-zero score, so you can trace a signal back to the cells
that produced it, and summarize signals by neighborhood when you need a
wider view.

## Installation

NICHESv2 is being prepared for CRAN. Until it is released there, install
it from GitHub:

``` r

install.packages("remotes")
remotes::install_github("RaredonLab/NICHESv2", ref = "main")
```

To read the vignette locally, add `build_vignettes = TRUE`, which needs
knitr, rmarkdown and pandoc.

These packages are optional. Install only the ones you need:

``` r

# Seurat V4 and V5 objects
install.packages("SeuratObject")

# AnnData (.h5ad) files, no Python required
install.packages("BiocManager")
BiocManager::install("rhdf5")

# OmniPath LR database
BiocManager::install("OmnipathR")

# TissuePlex export
install.packages("arrow")

# Cell segmentation polygons
install.packages("sf")
```

## Quick start

This example uses the small synthetic dataset that ships with the
package.

``` r

library(NICHESv2)

data(NICHESv2_inputs)

# Score LR pairs between mutual spatial neighbors
obj <- create_NICHESObject(
  count.mtx = NICHESv2_inputs$count.mtx,
  meta.data = NICHESv2_inputs$meta.data,
  LRM.db    = NICHESv2_inputs$LRM.db,
  mode      = "spatial",
  k         = 5L
)

# One row per non-zero (edge, LR pair) score
head(obj$edge.data)

# Summarize each cell's neighborhood
obj <- aggregate_NICHESObject(obj, cell.type.col = "celltype")
head(obj$aggregations$neighborhood.edge.agg)
head(obj$aggregations$neighborhood.composition)
```

For your own data, pass a database name to `LRM.db` and set `species`:

``` r

obj <- create_NICHESObject(
  count.mtx = count.mtx,
  meta.data = meta.data,
  LRM.db    = "connectomedb2025",
  species   = "mouse",
  mode      = "spatial",
  k         = 6L,
  n.cores   = 4L
)
```

The vignette walks through the full workflow:

``` r

vignette("NICHESv2")
```

The [website](https://raredonlab.github.io/NICHESv2/) has tutorials on
public spatial and non-spatial datasets, plus a page for reviewers.

## Features

- **Two analysis modes:** `"spatial"` connects cells by k-nearest
  neighbors or a fixed radius, and adds autocrine self-loops.
  `"nonspatial"` pairs cells by random sampling within each pair of cell
  types.
- **LR databases:** ConnectomeDB2025 and FANTOM5 are bundled for human,
  mouse, rat and pig, with no internet needed. OmniPath covers human,
  mouse and rat. You can also supply your own data frame, and multi-gene
  complexes are supported.
- **Memory-efficient storage:** the 15-slot `NICHESObject` stores scores
  in sparse long-format `data.table`s, so memory use follows how much
  signal exists rather than every possible combination.
- **Neighborhood aggregation:**
  [`aggregate_NICHESObject()`](https://raredonlab.github.io/NICHESv2/reference/aggregate_NICHESObject.md)
  gives directional scores per cell and LR pair (in, out, cross,
  self-autocrine and other-autocrine) and the cell-type makeup of each
  neighborhood.
- **Multiple samples:** `sample.col` processes each sample separately
  and merges the results, so no edges cross samples.
- **Flexible input:** use a sparse matrix, a Seurat V4 or V5 object, or
  an AnnData file (no Python required).
- **Works on every platform:** parallel scoring uses socket clusters
  through `n.cores`.
- **Interactive viewing:**
  [`export_to_TissuePlex()`](https://raredonlab.github.io/NICHESv2/reference/export_to_TissuePlex.md)
  writes edges to a Parquet file for viewing in
  [TissuePlex](https://raredonlab.github.io/TissuePlex/).

## Citation

A paper describing NICHESv2 is in preparation. Until it is published,
please cite the original NICHES method:

Raredon MSB, Yang J, Kothapalli N, et al. Comprehensive visualization of
cell-cell interactions in single-cell and spatial transcriptomics with
NICHES. *Bioinformatics*. 2023;39(1):btac775.
[doi:10.1093/bioinformatics/btac775](https://doi.org/10.1093/bioinformatics/btac775)

## License

NICHESv2 is dual-licensed.

- **Open source** under the [GNU Affero General Public License
  v3.0](https://raredonlab.github.io/NICHESv2/LICENSE.md). Use, modify
  and share it freely. If you distribute a modified version or offer it
  over a network, the AGPL requires that you make your modified source
  available under the same terms. Academic and non-profit use is
  unaffected in practice.
- **Commercial license** for use in products or services without the
  AGPL’s source-sharing obligations: contact the Raredon Lab.

© 2026 Nuoya Wang, Micha Sam Brickman Raredon
