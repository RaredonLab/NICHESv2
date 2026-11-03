# Subset a NICHESObject to a Set of Cells

Filters all cell-indexed slots (`$cell.list`, `$cell.data`,
`$cell.meta`, the three neighborhood slots and every table in
`$aggregations`) to the cells in `i`, and all edge-indexed slots
(`$edge.list`, `$edge.data`, `$edge.meta`) to keep only edges where
**both** sending and receiving cells are in `i`. If `$segmentation` is
populated, it is filtered to retain only rows whose `cell` is in the
retained set. Slot-level invariants are re-validated by the internal
constructor before the subset is returned.

## Usage

``` r
# S3 method for class 'NICHESObject'
x[i, ...]
```

## Arguments

- x:

  A `NICHESObject`.

- i:

  Character vector of cell barcodes to retain.

- ...:

  Unused. Present for S3 method signature conformance.

## Value

A new `NICHESObject` restricted to the retained cells and their
fully-internal edges, with their neighborhood entries and aggregation
rows unchanged. Gene and LRM slots are unchanged.

## Details

Cell barcodes in `i` that are not present in `x$cell.list` are silently
dropped. If the intersection is empty the function errors rather than
returning an unusable empty object.

## Neighborhoods and aggregations after subsetting

A retained cell keeps its neighborhood entries and aggregation rows
unchanged, so they still describe its whole neighborhood in the full
object, including neighbors outside `i`. Aggregate the full object
first, then subset. Summing a cell's edges in `$edge.data` can give less
than its aggregation row, because the edge slots hold only edges between
retained cells; a
[`export_to_TissuePlex`](https://raredonlab.github.io/NICHESv2/reference/export_to_TissuePlex.md)
export of the subset shows those edges only.

Aggregations cannot be recomputed on a subset whose neighborhoods
include removed edges:
[`add_neighborhood_edge_agg`](https://raredonlab.github.io/NICHESv2/reference/add_neighborhood_edge_agg.md)
and
[`add_neighborhood_composition`](https://raredonlab.github.io/NICHESv2/reference/add_neighborhood_composition.md)
stop with an error. To analyze the subset as a separate tissue instead,
call `add_neighborhoods(obj, overwrite = TRUE)` first. This rebuilds the
neighborhoods from the retained edges and clears the carried-over
aggregation tables.

## See also

[`is_NICHESObject`](https://raredonlab.github.io/NICHESv2/reference/is_NICHESObject.md)
to validate the result.

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

# Subset to the first four cells; only fully-internal edges are kept
sub <- obj[paste0("Cell", 1:4)]
dim(sub)
#> cells genes edges  LRMs 
#>     4     6     6     2 
# }
```
