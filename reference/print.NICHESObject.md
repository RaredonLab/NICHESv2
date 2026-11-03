# Print a NICHESObject

Compact one-screen summary of a `NICHESObject`: cell, gene, edge, and
LRM-pair counts (with non-empty / valid counts where applicable),
segmentation status, and a slot legend. Use
[`summary.NICHESObject`](https://raredonlab.github.io/NICHESv2/reference/summary.NICHESObject.md)
for a longer breakdown including cell-type composition and signal
density.

## Usage

``` r
# S3 method for class 'NICHESObject'
print(x, ...)
```

## Arguments

- x:

  A `NICHESObject`.

- ...:

  Unused. Present for S3 method signature conformance.

## Value

The input `x`, returned invisibly. Called for the side effect of
printing to the console.

## See also

[`summary.NICHESObject`](https://raredonlab.github.io/NICHESv2/reference/summary.NICHESObject.md),
[`dim.NICHESObject`](https://raredonlab.github.io/NICHESv2/reference/dim.NICHESObject.md),
[`is_NICHESObject`](https://raredonlab.github.io/NICHESv2/reference/is_NICHESObject.md).

## Examples

``` r
# \donttest{
# Build a tiny NICHESObject for demonstration
set.seed(1L)
count.mtx <- Matrix::sparseMatrix(
  i = sample(6L, 30L, replace = TRUE),
  j = sample(8L, 30L, replace = TRUE),
  x = as.numeric(rpois(30L, lambda = 3) + 1L),
  dims = c(6L, 8L),
  dimnames = list(paste0("Gene", 1:6), paste0("Cell", 1:8))
)
meta.data <- data.frame(
  x = runif(8L), y = runif(8L),
  row.names = paste0("Cell", 1:8)
)
LRM.db <- data.frame(
  ligand   = c("Gene1", "Gene3"),
  receptor = c("Gene2", "Gene4")
)
obj <- create_NICHESObject(
  count.mtx = count.mtx, meta.data = meta.data, LRM.db = LRM.db,
  mode = "spatial", ligand.col = "ligand", receptor.col = "receptor",
  k = 3L, n.cores = 1L, verbose = FALSE
)
#> [create_Edgelist_Spatial] Mode: kNN    | Param: 3 | Cells: 8 | Inter-cell: 18 | Autocrine: 8
#> [compute_CellToCell] 26 edges x 2 LR pairs | method = 'product' | cores = 1
#> [compute_CellToCell] Done. 9 non-zero edge signals.

print(obj)
#> -- NICHESObject --------------------------------------------
#>   Cells        : 8  (7 non-empty)
#>   Genes        : 6  (6 valid)
#>   Edges        : 26  (9 non-empty)
#>   LRM pairs    : 2  (2 valid)
#>   Segmentation : none
#>   Neighborhoods: 8 cells  |  68 edges  (mode: both_in_neighborhood, autocrine: included)
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
obj            # implicit print
#> -- NICHESObject --------------------------------------------
#>   Cells        : 8  (7 non-empty)
#>   Genes        : 6  (6 valid)
#>   Edges        : 26  (9 non-empty)
#>   LRM pairs    : 2  (2 valid)
#>   Segmentation : none
#>   Neighborhoods: 8 cells  |  68 edges  (mode: both_in_neighborhood, autocrine: included)
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
# }
```
