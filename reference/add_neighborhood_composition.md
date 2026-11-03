# Add Neighborhood Cell-Type Composition to a NICHESObject

For each cell, counts how many neighbors of each cell type fall within
its local neighborhood. Neighborhood membership is derived directly from
`$edge.list`: cell B is considered a neighbor of cell A if any directed
edge connects them (regardless of direction). Counts are zero-filled
against the full dataset's cell type universe so that every cell has a
value for every type.

## Usage

``` r
add_neighborhood_composition(
  niches.obj,
  cell.type.col = NULL,
  include.self = TRUE,
  normalize = TRUE,
  verbose = TRUE
)
```

## Arguments

- niches.obj:

  A `NICHESObject` with `mode = "spatial"`.

- cell.type.col:

  Character. Name of the column in `$cell.meta` holding cell type
  annotations. Hard error if `NULL` or not found.

- include.self:

  Logical. Whether to count the focal cell itself as one of its own
  neighbors. Default `TRUE`. When `TRUE`, cell A contributes 1 to its
  own type's count regardless of whether an autocrine self-loop edge
  exists.

- normalize:

  Logical. When `TRUE` (default), also produces
  `<PolishedType>.neighbor.prop` columns: each count divided by the
  total number of neighbors (including self when `include.self = TRUE`).
  Proportions sum to 1 per cell; cells with zero neighbors receive 0.

- verbose:

  Logical. Emit progress messages. Default `TRUE`.

## Value

The input `niches.obj` with `$aggregations$neighborhood.composition`
populated or replaced.

## Details

The result is stored as a single data.table at
`$aggregations$neighborhood.composition`, replacing any previously
computed version unconditionally.

## Column naming

Cell type labels are sanitized for use as column names via
`gsub("[^A-Za-z0-9.]", ".", type)`, which replaces spaces, slashes,
parentheses, and other special characters with dots. The resulting
columns follow the pattern `<PolishedType>.neighbor.count` (and
optionally `<PolishedType>.neighbor.prop`). Columns are ordered
alphabetically by polished type name. For example, cell type
`"CD8+ T cell"` becomes `"CD8..T.cell.neighbor.count"`.

## Spatial mode only

This function raises a hard error when
`attr(niches.obj, "mode") == "nonspatial"`. In nonspatial mode, edges
are statistically paired rather than spatially proximate, so
neighborhood composition would reflect sampling design rather than local
tissue structure. This is biologically meaningless and is therefore not
supported.

## Double-counting prevention

Because spatial graphs are typically bidirectional (edges A-\>B and
B-\>A both exist for proximate cells), the neighbor set for each focal
cell is formed as the *union* of all cells appearing on the other end of
any edge, deduplicated before counting. Each physical neighbor is
counted exactly once regardless of how many directed edges connect them.

## See also

[`add_neighborhood_edge_agg`](https://raredonlab.github.io/NICHESv2/reference/add_neighborhood_edge_agg.md),
[`add_neighborhoods`](https://raredonlab.github.io/NICHESv2/reference/add_neighborhoods.md),
[`aggregate_NICHESObject`](https://raredonlab.github.io/NICHESv2/reference/aggregate_NICHESObject.md)

## Examples

``` r
# \donttest{
data(NICHESv2_example)

obj <- add_neighborhood_composition(
  NICHESv2_example,
  cell.type.col = "celltype"
)
#> [add_neighborhood_composition] Computing composition across 80 cells | 3 cell types ...
#> [add_neighborhood_composition] Done. 80 cells | 3 cell types | columns: count + prop.
head(obj$aggregations$neighborhood.composition)
#> Key: <cell>
#>       cell TypeA.neighbor.count TypeB.neighbor.count TypeC.neighbor.count
#>     <char>                <int>                <int>                <int>
#> 1:  Cell_1                    5                    0                    0
#> 2: Cell_10                    6                    0                    0
#> 3: Cell_11                    6                    0                    0
#> 4: Cell_12                    4                    0                    0
#> 5: Cell_13                    2                    0                    0
#> 6: Cell_14                    6                    0                    0
#>    TypeA.neighbor.prop TypeB.neighbor.prop TypeC.neighbor.prop
#>                  <num>               <num>               <num>
#> 1:                   1                   0                   0
#> 2:                   1                   0                   0
#> 3:                   1                   0                   0
#> 4:                   1                   0                   0
#> 5:                   1                   0                   0
#> 6:                   1                   0                   0

# Without self-counting and without proportions
obj2 <- add_neighborhood_composition(
  NICHESv2_example,
  cell.type.col = "celltype",
  include.self  = FALSE,
  normalize     = FALSE
)
#> [add_neighborhood_composition] Computing composition across 80 cells | 3 cell types ...
#> [add_neighborhood_composition] Done. 80 cells | 3 cell types | columns: count.
# }
```
