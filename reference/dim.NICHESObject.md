# Dimensions of a NICHESObject

Returns the four headline counts that characterize a `NICHESObject`: the
number of cells, genes, edges, and ligand–receptor (LRM) pairs. Provided
so that `dim(obj)` works on a `NICHESObject` the same way it does on a
matrix.

## Usage

``` r
# S3 method for class 'NICHESObject'
dim(x)
```

## Arguments

- x:

  A `NICHESObject`.

## Value

Named integer vector of length four with names `cells`, `genes`,
`edges`, `LRMs`.

## See also

[`length.NICHESObject`](https://raredonlab.github.io/NICHESv2/reference/length.NICHESObject.md)
for cell count alone;
[`summary.NICHESObject`](https://raredonlab.github.io/NICHESv2/reference/summary.NICHESObject.md)
for a fuller breakdown.

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

dim(obj)
#> cells genes edges  LRMs 
#>     8     6    26     2 
# }
```
