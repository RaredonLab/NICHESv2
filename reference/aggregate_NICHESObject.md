# Run the Full Aggregation Pipeline on a NICHESObject

A convenience wrapper that calls all aggregation functions in the
correct dependency order, populating the `$aggregations` slot in a
single call. Both steps always fully recompute and replace their slots –
there is no incremental update mode. Individual functions can still be
called directly for finer control.

## Usage

``` r
aggregate_NICHESObject(
  niches.obj,
  cell.type.col = NULL,
  n.cores = 1L,
  chunk.size = 50L,
  verbose = TRUE
)
```

## Arguments

- niches.obj:

  A `NICHESObject` with neighborhood slots populated.

- cell.type.col:

  Character or `NULL`. Name of the column in `$cell.meta` holding cell
  type annotations, passed to
  [`add_neighborhood_composition`](https://raredonlab.github.io/NICHESv2/reference/add_neighborhood_composition.md).
  When `NULL` (default), step 2 is skipped.

- n.cores:

  Positive integer. Retained for API compatibility; passed to
  [`add_neighborhood_edge_agg`](https://raredonlab.github.io/NICHESv2/reference/add_neighborhood_edge_agg.md)
  but not currently used by the sparse matrix implementation. Default
  `1L`.

- chunk.size:

  Positive integer. Number of LRM pairs processed per sparse matrix
  batch in
  [`add_neighborhood_edge_agg`](https://raredonlab.github.io/NICHESv2/reference/add_neighborhood_edge_agg.md).
  Smaller values reduce peak memory; larger values are marginally
  faster. Default `50L`. Reduce to between `20L` and `30L` on
  memory-constrained machines (e.g. laptops); increase on HPC nodes with
  abundant RAM.

- verbose:

  Logical. Emit progress messages. Default `TRUE`.

## Value

The input `niches.obj` with `$aggregations` populated:

- `$neighborhood.edge.agg`:

  Always populated (five-category edge aggregation keyed on
  `(cell, LRM)`, with all three score methods and in/out ratio).

- `$neighborhood.composition`:

  Populated when `cell.type.col` is provided and `mode = "spatial"`.

## Pipeline steps

1.  [`add_neighborhood_edge_agg`](https://raredonlab.github.io/NICHESv2/reference/add_neighborhood_edge_agg.md):
    for each focal cell and each LRM pair, aggregates scores across all
    neighborhood edges using a five-way classification (in, out, cross,
    self.autocrine, other.autocrine). All three score methods (sum,
    mean, cmean) and the in/out ratio are always computed. LRM pairs are
    processed in batches of `chunk.size` to bound peak memory.

2.  [`add_neighborhood_composition`](https://raredonlab.github.io/NICHESv2/reference/add_neighborhood_composition.md):
    computes cell-type composition of each cell's local neighborhood
    from `$edge.list` and `$cell.meta`. *Skipped if
    `cell.type.col = NULL` or if the object mode is `"nonspatial"`.*

## Neighborhood slots requirement

Step 1 requires the three neighborhood slots to be populated.
[`create_NICHESObject`](https://raredonlab.github.io/NICHESv2/reference/create_NICHESObject.md)
populates them automatically using
`edge.filter.mode = "both_in_neighborhood"` and
`remove.autocrine = FALSE`, so this requirement is met in the standard
workflow.

## See also

[`add_neighborhood_edge_agg`](https://raredonlab.github.io/NICHESv2/reference/add_neighborhood_edge_agg.md),
[`add_neighborhood_composition`](https://raredonlab.github.io/NICHESv2/reference/add_neighborhood_composition.md)

## Examples

``` r
# \donttest{
data(NICHESv2_example)

# Edge aggregation only
obj <- aggregate_NICHESObject(NICHESv2_example, n.cores = 1L)
#> [aggregate_NICHESObject] Starting aggregation pipeline ...
#> [aggregate_NICHESObject] Step 1/2: neighborhood edge aggregation.
#> [add_neighborhood_edge_agg] Aggregating 80 cells x 15 LRMs x 1364 neighborhood edges (chunk.size = 50) ...
#>   [add_neighborhood_edge_agg] chunk 1 / 1 done.
#> [add_neighborhood_edge_agg] Done. 399 rows | 74 cells | 15 LRM pairs.
#> [aggregate_NICHESObject] Step 2/2: SKIPPED -- 'cell.type.col' not provided. Supply cell.type.col to compute neighborhood composition.
#> [aggregate_NICHESObject] Pipeline complete. Slots populated: neighborhood.edge.agg.
names(obj$aggregations)
#> [1] "neighborhood.edge.agg"

# With cell-type composition
obj <- aggregate_NICHESObject(
  NICHESv2_example,
  cell.type.col = "celltype",
  n.cores       = 1L
)
#> [aggregate_NICHESObject] Starting aggregation pipeline ...
#> [aggregate_NICHESObject] Step 1/2: neighborhood edge aggregation.
#> [add_neighborhood_edge_agg] Aggregating 80 cells x 15 LRMs x 1364 neighborhood edges (chunk.size = 50) ...
#>   [add_neighborhood_edge_agg] chunk 1 / 1 done.
#> [add_neighborhood_edge_agg] Done. 399 rows | 74 cells | 15 LRM pairs.
#> [aggregate_NICHESObject] Step 2/2: neighborhood composition.
#> [add_neighborhood_composition] Computing composition across 80 cells | 3 cell types ...
#> [add_neighborhood_composition] Done. 80 cells | 3 cell types | columns: count + prop.
#> [aggregate_NICHESObject] Pipeline complete. Slots populated: neighborhood.edge.agg, neighborhood.composition.
# }
```
