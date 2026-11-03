# Extract inputs for create_NICHESObject from a Seurat object

Extracts the count matrix and cell metadata from a Seurat object (V4 or
V5 internal structure) and returns them as a named list ready to pass
directly to
[`create_NICHESObject`](https://raredonlab.github.io/NICHESv2/reference/create_NICHESObject.md).
The object structure is detected from the assay class, not the installed
package version, so V4 objects loaded in a V5 environment are handled
correctly.

## Usage

``` r
extract_NICHESInputs_Seurat(
  seurat.obj,
  assay = "RNA",
  layer = "counts",
  verbose = TRUE
)
```

## Arguments

- seurat.obj:

  A Seurat object. Both V4 (`Assay`) and V5 (`Assay5`) internal
  structures are supported.

- assay:

  Character. The assay to extract from. Default `"RNA"`. Common
  alternatives: `"SCT"`, `"Spatial"`, `"Xenium"`.

- layer:

  Character. The data layer (V5) or slot (V4) to extract. Default
  `"counts"` (raw UMI counts), a starting point rather than a
  recommendation: this choice affects which edges score highest. Other
  options: `"data"` (log-normalized), `"scale.data"` (scaled; V4 only).

- verbose:

  Logical. If `TRUE` (default), prints a summary of what was extracted
  and a reminder about spatial coordinate requirements.

## Value

A named list with two elements:

- `count.mtx`:

  Sparse matrix (genes \\\times\\ cells), class `dgCMatrix`. Row names
  are gene names; column names are cell barcodes.

- `meta.data`:

  Data frame of cell metadata extracted directly from
  `seurat.obj@meta.data`, unmodified. Row names are cell barcodes.

## Details

The function does **not** fetch spatial coordinates from
`seurat.obj@images`. For spatial mode in
[`create_NICHESObject()`](https://raredonlab.github.io/NICHESv2/reference/create_NICHESObject.md),
columns named `x` and `y` must be present in `meta.data`. If your
coordinates are in `seurat.obj@images`, add them to the returned
`meta.data` yourself before calling
[`create_NICHESObject()`](https://raredonlab.github.io/NICHESv2/reference/create_NICHESObject.md):


      coords <- SeuratObject::GetTissueCoordinates(seurat.obj)
      inputs$meta.data$x <- coords$imagerow[match(rownames(inputs$meta.data),
                                                   rownames(coords))]
      inputs$meta.data$y <- coords$imagecol[match(rownames(inputs$meta.data),
                                                   rownames(coords))]

## See also

[`create_NICHESObject`](https://raredonlab.github.io/NICHESv2/reference/create_NICHESObject.md),
[`extract_NICHESInputs_AnnData`](https://raredonlab.github.io/NICHESv2/reference/extract_NICHESInputs_AnnData.md)

## Examples

``` r
if (requireNamespace("SeuratObject", quietly = TRUE)) {
  data(NICHESv2_inputs)
  seurat.obj <- SeuratObject::CreateSeuratObject(
    counts    = NICHESv2_inputs$count.mtx,
    meta.data = NICHESv2_inputs$meta.data
  )
  inputs <- extract_NICHESInputs_Seurat(seurat.obj, assay = "RNA",
                                        layer = "counts", verbose = FALSE)
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
