# Add Unified Neighborhood Edge Aggregation to a NICHESObject

For each focal cell and each LRM pair, aggregates LR communication
scores across all edges in the cell's neighborhood. Each neighborhood
edge is classified into one of five mutually exclusive categories
relative to the focal cell: `out` (focal cell sends), `in` (focal cell
receives), `cross` (neither endpoint is the focal cell),
`self.autocrine` (focal cell's own self-loop: SendingCell ==
ReceivingCell == focal cell), or `other.autocrine` (a neighbor cell's
self-loop: SendingCell == ReceivingCell != focal cell).

## Usage

``` r
add_neighborhood_edge_agg(
  niches.obj,
  n.cores = 1L,
  chunk.size = 50L,
  verbose = TRUE
)
```

## Arguments

- niches.obj:

  A `NICHESObject` with neighborhood slots populated.

- n.cores:

  Positive integer. Retained for API compatibility; not used in the
  current sparse matrix implementation. Default `1L`.

- chunk.size:

  Positive integer. Number of LRM pairs processed per sparse matrix
  batch. Smaller values reduce peak memory at the cost of slightly more
  overhead; larger values are faster but require more RAM. Default `50L`
  (about 12 batches for 558 LRMs). Reduce to between `20L` and `30L` on
  memory-constrained machines.

- verbose:

  Logical. Emit progress messages. Default `TRUE`.

## Value

The input `niches.obj` with `$aggregations$neighborhood.edge.agg`
populated (always replaced).

## Details

All three score methods (`sum`, `mean`, `cmean`) and the in/out ratio
are always computed. The slot is always fully replaced on each call –
there is no incremental update mode. This follows the same convention as
Seurat slot-populating functions: call the function when you want a
fresh result.

Results are stored as a sparse data.table keyed on `(cell, LRM)` at
`$aggregations$neighborhood.edge.agg`. Only `(cell, LRM)` pairs with at
least one non-zero-score neighborhood edge are included.

## Five-way edge classification

For a focal cell F and an edge (S -\> R) in F's neighborhood:

- `self.autocrine`:

  S == R == F. F's own self-loop edge. Checked first.
  `count.self.autocrine` is always 0 or 1 per (cell, LRM).

- `other.autocrine`:

  S == R != F. A neighbor cell's self-loop edge; F is not the
  self-looping cell.

- `out`:

  S == F and S != R. Focal cell is the sender.

- `in`:

  R == F and S != R. Focal cell is the receiver.

- `cross`:

  S != F and R != F and S != R. Edge between two neighbors; focal cell
  is not involved.

## Columns produced

Count columns (always produced):

- `count.in`, `count.out`, `count.cross`, `count.self.autocrine`,
  `count.other.autocrine`:

  Number of non-zero-score edges in each category for this LRM. Because
  `$edge.data` is sparse, every stored row represents an active edge.
  `count.self.autocrine` is always 0 or 1 per (cell, LRM).

- `count.total`:

  Sum of all five category counts.

Score columns (all three methods always produced, for each category):

- `score.{cat}.sum`:

  Total score across category edges.

- `score.{cat}.mean`:

  Sum divided by total neighborhood edges
  (`$neighborhood.meta$n.edges`), including zero-score edges.
  LRM-agnostic denominator.

- `score.{cat}.cmean`:

  Conditional mean: sum divided by the category count. `NA` when count
  is zero.

- `score.inout.ratio`:

  `score.in.sum / score.out.sum`. `NA` when `score.out.sum == 0`.

## See also

[`add_neighborhoods`](https://raredonlab.github.io/NICHESv2/reference/add_neighborhoods.md),
[`add_neighborhood_composition`](https://raredonlab.github.io/NICHESv2/reference/add_neighborhood_composition.md),
[`aggregate_NICHESObject`](https://raredonlab.github.io/NICHESv2/reference/aggregate_NICHESObject.md)

## Examples

``` r
# \donttest{
data(NICHESv2_example)

obj <- add_neighborhood_edge_agg(NICHESv2_example, n.cores = 1L)
#> [add_neighborhood_edge_agg] Aggregating 80 cells x 15 LRMs x 1364 neighborhood edges (chunk.size = 50) ...
#>   [add_neighborhood_edge_agg] chunk 1 / 1 done.
#> [add_neighborhood_edge_agg] Done. 399 rows | 74 cells | 15 LRM pairs.
head(obj$aggregations$neighborhood.edge.agg)
#> Key: <cell, LRM>
#>      cell           LRM count.in count.out count.cross count.self.autocrine
#>    <char>        <char>    <int>     <int>       <int>                <int>
#> 1: Cell_1 Gene11|Gene12        1         0           0                    1
#> 2: Cell_1 Gene13|Gene14        0         0           1                    0
#> 3: Cell_1 Gene15|Gene16        0         0           0                    0
#> 4: Cell_1 Gene17|Gene18        0         0           0                    0
#> 5: Cell_1   Gene3|Gene4        1         0           0                    0
#> 6: Cell_1  Gene9|Gene10        1         0           0                    1
#>    count.other.autocrine count.total score.in.sum score.out.sum score.cross.sum
#>                    <int>       <int>        <num>         <num>           <num>
#> 1:                     0           2           10             0               0
#> 2:                     0           1            0             0               8
#> 3:                     1           1            0             0               0
#> 4:                     1           1            0             0               0
#> 5:                     0           1           48             0               0
#> 6:                     0           2           21             0               0
#>    score.self.autocrine.sum score.other.autocrine.sum score.total.sum
#>                       <num>                     <num>           <num>
#> 1:                       12                         0              22
#> 2:                        0                         0               8
#> 3:                        0                        35              35
#> 4:                        0                         6               6
#> 5:                        0                         0              48
#> 6:                       14                         0              35
#>    score.in.mean score.out.mean score.cross.mean score.self.autocrine.mean
#>            <num>          <num>            <num>                     <num>
#> 1:     0.5882353              0        0.0000000                 0.7058824
#> 2:     0.0000000              0        0.4705882                 0.0000000
#> 3:     0.0000000              0        0.0000000                 0.0000000
#> 4:     0.0000000              0        0.0000000                 0.0000000
#> 5:     2.8235294              0        0.0000000                 0.0000000
#> 6:     1.2352941              0        0.0000000                 0.8235294
#>    score.other.autocrine.mean score.total.mean score.in.cmean score.out.cmean
#>                         <num>            <num>          <num>           <num>
#> 1:                  0.0000000        1.2941176             10              NA
#> 2:                  0.0000000        0.4705882             NA              NA
#> 3:                  2.0588235        2.0588235             NA              NA
#> 4:                  0.3529412        0.3529412             NA              NA
#> 5:                  0.0000000        2.8235294             48              NA
#> 6:                  0.0000000        2.0588235             21              NA
#>    score.cross.cmean score.self.autocrine.cmean score.other.autocrine.cmean
#>                <num>                      <num>                       <num>
#> 1:                NA                         12                          NA
#> 2:                 8                         NA                          NA
#> 3:                NA                         NA                          35
#> 4:                NA                         NA                           6
#> 5:                NA                         NA                          NA
#> 6:                NA                         14                          NA
#>    score.total.cmean score.inout.ratio
#>                <num>             <num>
#> 1:              11.0                NA
#> 2:               8.0                NA
#> 3:              35.0                NA
#> 4:               6.0                NA
#> 5:              48.0                NA
#> 6:              17.5                NA
# }
```
