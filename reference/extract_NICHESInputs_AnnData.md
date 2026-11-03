# Extract inputs for create_NICHESObject from an AnnData (.h5ad) file

Reads a `.h5ad` file using rhdf5 (no Python, no conda, no basilisk –
pure C HDF5 bindings, compatible with R \\\geq\\ 4.2) and returns the
count matrix and cell metadata as a named list ready to pass directly to
[`create_NICHESObject`](https://raredonlab.github.io/NICHESv2/reference/create_NICHESObject.md).

## Usage

``` r
extract_NICHESInputs_AnnData(
  file,
  layer = "counts",
  obs.cols = NULL,
  verbose = TRUE
)
```

## Arguments

- file:

  Character. Full path to the `.h5ad` file, including the file name and
  extension (e.g. `"/data/samples/sample1.h5ad"`).

- layer:

  Character. Which data to extract. Default `"counts"` (raw UMI counts
  stored in `adata.layers[["counts"]]`), a starting point rather than a
  recommendation: this choice affects which edges score highest. Use
  `"X"` to extract `adata.X` directly. A hard error with available layer
  names is thrown if the layer is absent.

- obs.cols:

  Character vector or `NULL`. Subset of `adata.obs` columns to keep.
  `NULL` (default) keeps all.

- verbose:

  Logical. If `TRUE` (default), prints a summary and a spatial
  coordinate reminder.

## Value

A named list:

- `count.mtx`:

  Sparse `dgCMatrix`, genes \\\times\\ cells. Row names are gene names;
  column names are cell barcodes.

- `meta.data`:

  Data frame from `adata.obs` (optionally column-subsetted). Row names
  are cell barcodes.

## Details

The AnnData layer structure maps to the `layer` argument as follows:

- `adata.X` \\\rightarrow\\ `layer = "X"`

- `adata.layers[["counts"]]` \\\rightarrow\\ `layer = "counts"`

- any other named layer \\\rightarrow\\ `layer = "<name>"`

Install rhdf5 with:


      BiocManager::install("rhdf5")

The function does **not** fetch spatial coordinates from `adata.obsm`.
For spatial mode, add `x` and `y` columns to `meta.data` yourself:


      library(rhdf5)
      coords <- t(rhdf5::h5read(file, "obsm/spatial"))
      inputs$meta.data$x <- coords[, 1]
      inputs$meta.data$y <- coords[, 2]

## See also

[`create_NICHESObject`](https://raredonlab.github.io/NICHESv2/reference/create_NICHESObject.md),
[`extract_NICHESInputs_Seurat`](https://raredonlab.github.io/NICHESv2/reference/extract_NICHESInputs_Seurat.md)

## Examples

``` r
if (requireNamespace("rhdf5", quietly = TRUE)) {
  h5ad.file <- system.file("extdata", "NICHESv2_example.h5ad",
                           package = "NICHESv2")
  inputs <- extract_NICHESInputs_AnnData(
    file    = h5ad.file,
    layer   = "counts",
    verbose = FALSE
  )
  data(NICHESv2_inputs)
  niches.obj <- create_NICHESObject(
    count.mtx = inputs$count.mtx,
    meta.data = inputs$meta.data,
    LRM.db    = NICHESv2_inputs$LRM.db,
    mode      = "spatial",
    k         = 5L,
    verbose   = FALSE
  )
  niches.obj
}
#> [create_Edgelist_Spatial] Mode: kNN    | Param: 5 | Cells: 80 | Inter-cell: 274 | Autocrine: 80
#> [compute_CellToCell] 354 edges x 15 LR pairs | method = 'product' | cores = 1
#> [compute_CellToCell] Done. 190 non-zero edge signals.
#> -- NICHESObject --------------------------------------------
#>   Cells        : 80  (80 non-empty)
#>   Genes        : 50  (50 valid)
#>   Edges        : 354  (141 non-empty)
#>   LRM pairs    : 15  (15 valid)
#>   Segmentation : none
#>   Neighborhoods: 80 cells  |  1364 edges  (mode: both_in_neighborhood, autocrine: included)
#>   Aggregations : none
#> 
#>   Slot summary:
#>     $cell.list     -- character vector of cell barcodes
#>     $cell.data     -- (cell, gene, score, score.norm) [sparse]
#>     $cell.meta     -- per-cell metadata
#>     $edge.list     -- (edge, SendingCell, ReceivingCell)
#>     $edge.data     -- (edge, LRM, score, score.norm) [sparse]
#>     $edge.meta     -- per-edge metadata with Sending/Receiving/Joint cols
#>     $neighborhood.cell.list -- named list: cell -> neighbor barcodes
#>     $neighborhood.edge.list -- named list: cell -> edge IDs
#>     $neighborhood.meta      -- (cell, n.neighbors, n.edges) [keyed]
#>     $gene.list     -- character vector of gene names
#>     $gene.meta     -- per-gene metadata
#>     $LRM.list      -- character vector of 'ligand|receptor' pairs
#>     $LRM.meta      -- per-LRM metadata
#>     $segmentation  -- sf object of cell polygons, or NULL
#>     $aggregations  -- named list of aggregation data.tables
#> ------------------------------------------------------------
```
