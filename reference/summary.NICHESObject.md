# Summarize a NICHESObject

Long-form summary of a `NICHESObject`. Reports cell counts, cell-type
composition (when a `celltype` column is present in `$cell.meta`), gene
count, edge totals split into autocrine (self-loop) versus inter-cell,
LRM-pair count, edge-signal density (filled / possible cells of the edge
x LRM grid), how many of the three neighborhood slots are populated, and
segmentation status.

## Usage

``` r
# S3 method for class 'NICHESObject'
summary(object, ...)
```

## Arguments

- object:

  A `NICHESObject`.

- ...:

  Unused. Present for S3 method signature conformance.

## Value

The input `object`, returned invisibly. Called for the side effect of
printing to the console.

## Details

Use
[`print.NICHESObject`](https://raredonlab.github.io/NICHESv2/reference/print.NICHESObject.md)
for a shorter one-screen overview.

## See also

[`print.NICHESObject`](https://raredonlab.github.io/NICHESv2/reference/print.NICHESObject.md),
[`dim.NICHESObject`](https://raredonlab.github.io/NICHESv2/reference/dim.NICHESObject.md).

## Examples

``` r
# \donttest{
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
  celltype = rep(c("A", "B"), each = 4L),
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

summary(obj)
#> -- NICHESObject Summary ------------------------------------
#>   Cells: 8
#>   Cell types:
#>     A                    : 4
#>     B                    : 4
#> 
#>   Genes: 6
#> 
#>   Edges: 26
#>     Autocrine (self-loop) : 8
#>     Inter-cell            : 18
#> 
#>   LRM pairs: 2
#> 
#>   Edge signal density: 9 / 52 (17.31%)
#> 
#>   Neighborhood slots populated: 3 / 3
#> 
#>   Segmentation: none
#> 
#>   Aggregations: none
#> ------------------------------------------------------------
# }
```
