# Number of Cells in a NICHESObject

Returns the number of cells in a `NICHESObject`, identical to
`length(x$cell.list)`. Provided so that `length(obj)` works the way most
users expect for a per-cell collection.

## Usage

``` r
# S3 method for class 'NICHESObject'
length(x)
```

## Arguments

- x:

  A `NICHESObject`.

## Value

Integer scalar: number of cells.

## See also

[`dim.NICHESObject`](https://raredonlab.github.io/NICHESv2/reference/dim.NICHESObject.md)
for cell, gene, edge, and LRM counts together.

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

length(obj)   # 8
#> [1] 8
# }
```
