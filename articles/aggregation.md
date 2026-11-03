# Neighborhood aggregation

Edges show individual cell pairs. Aggregation summarizes them for each
cell’s neighborhood, which gives you per-cell values to compare across a
tissue.

``` r

library(NICHESv2)
library(data.table)

data(NICHESv2_example)
obj <- aggregate_NICHESObject(
  NICHESv2_example,
  cell.type.col = "celltype",
  verbose = FALSE
)
names(obj$aggregations)
#> [1] "neighborhood.edge.agg"    "neighborhood.composition"
```

[`aggregate_NICHESObject()`](https://raredonlab.github.io/NICHESv2/reference/aggregate_NICHESObject.md)
runs
[`add_neighborhood_edge_agg()`](https://raredonlab.github.io/NICHESv2/reference/add_neighborhood_edge_agg.md),
then
[`add_neighborhood_composition()`](https://raredonlab.github.io/NICHESv2/reference/add_neighborhood_composition.md)
if you pass `cell.type.col`. Each run replaces the previous result.

## Edge aggregation

`$aggregations$neighborhood.edge.agg` has 27 columns, with one row for
each (cell, LRM) pair that has signal.

``` r

agg <- obj$aggregations$neighborhood.edge.agg
dim(agg)
#> [1] 399  27
```

### Five categories

Each edge in a neighborhood falls into one category relative to the
focal cell:

| Category | Meaning |
|----|----|
| `in` | A neighbor sends to the focal cell |
| `out` | The focal cell sends to a neighbor |
| `cross` | One neighbor sends to another; the focal cell isn’t involved |
| `self.autocrine` | The focal cell’s own self-loop |
| `other.autocrine` | A neighbor’s self-loop |

The `cross` category captures what a cell’s surroundings exchange, even
when the cell itself takes no part.

### Columns

`{category}` below stands for each of the five categories, or `total`
for all of them together.

| Column | Meaning |
|----|----|
| `count.{category}` | Number of edges with signal |
| `score.{category}.sum` | Total score |
| `score.{category}.mean` | Total score divided by all edges in the neighborhood, including those with no signal for this LRM |
| `score.{category}.cmean` | Total score divided by the count, or `NA` when the count is zero |
| `score.inout.ratio` | `score.in.sum` divided by `score.out.sum`, or `NA` when nothing is sent |

Use `mean` to compare the overall signal between neighborhoods of
different sizes. Use `cmean` to see how strong the signal is on the
edges that carry it.

``` r

agg[, .(cell, LRM, count.in, score.in.sum, score.in.mean, score.in.cmean)][1:3]
```

| cell   | LRM            | count.in | score.in.sum | score.in.mean | score.in.cmean |
|:-------|:---------------|---------:|-------------:|--------------:|---------------:|
| Cell_1 | Gene11\|Gene12 |        1 |           10 |         0.588 |             10 |
| Cell_1 | Gene13\|Gene14 |        0 |            0 |         0.000 |             NA |
| Cell_1 | Gene15\|Gene16 |        0 |            0 |         0.000 |             NA |

Edges with no signal for any LR pair still count toward the `mean`
denominator, which comes from `$neighborhood.meta$n.edges`.

### Memory

`chunk.size` sets how many LR pairs are processed at once, and defaults
to 50. Lower it if you run short of memory. Aggregation uses sparse
matrix multiplication, so `n.cores` has no effect on this step.

## Neighborhood composition

`$aggregations$neighborhood.composition` has one row per cell, with
counts and proportions of neighboring cell types. By default it counts
the focal cell too; pass `include.self = FALSE` to
[`add_neighborhood_composition()`](https://raredonlab.github.io/NICHESv2/reference/add_neighborhood_composition.md)
to leave it out.

``` r

comp <- obj$aggregations$neighborhood.composition
head(comp[, .(cell, TypeA.neighbor.prop, TypeB.neighbor.prop, TypeC.neighbor.prop)], 3)
```

| cell    | TypeA.neighbor.prop | TypeB.neighbor.prop | TypeC.neighbor.prop |
|:--------|--------------------:|--------------------:|--------------------:|
| Cell_1  |                   1 |                   0 |                   0 |
| Cell_10 |                   1 |                   0 |                   0 |
| Cell_11 |                   1 |                   0 |                   0 |

This table needs spatial data. In `"nonspatial"` mode it stops with an
error, because neighbor cell types would only reflect the sampling.

## Subsets

Aggregate the full object first, then subset. Neighborhood slots and
aggregation tables are indexed by cell, so each kept cell keeps its rows
unchanged. Those values still describe the cell’s whole neighborhood,
including neighbors you didn’t keep.

``` r

keep <- obj$cell.meta[celltype == "TypeA", cell]
sub <- obj[keep]

identical(
  sub$aggregations$neighborhood.edge.agg$score.in.sum,
  agg[cell %in% keep]$score.in.sum
)
#> [1] TRUE
```

The edge slots of a subset keep only edges between kept cells. Two
things follow:

- Summing a cell’s rows in `sub$edge.data` can give less than its
  aggregation row.
- A TissuePlex export of the subset shows only those kept edges.

### Recalculating on a subset

The removed edges are gone, so a subset can’t recompute full-tissue
values, and trying stops with an error:

``` r

add_neighborhood_edge_agg(sub, verbose = FALSE)
#> Error:
#> ! [add_neighborhood_edge_agg] This object is a subset whose neighborhoods include edges to removed cells, so aggregations cannot be recomputed on it. Aggregate the full object before subsetting (the results carry over), or run add_neighborhoods(obj, overwrite = TRUE) first to analyze the subset as a separate tissue.
```

Each kept cell already has its correct values, so you rarely need to
recalculate. To treat the subset as a separate tissue instead, rebuild
its neighborhoods first. This clears the carried-over tables, and the
new values use only the kept cells:

``` r

separate <- add_neighborhoods(sub, overwrite = TRUE, verbose = FALSE)
separate <- aggregate_NICHESObject(separate, verbose = FALSE)
```

## Changing neighborhoods

`add_neighborhoods(overwrite = TRUE)` can also apply another
`edge.filter.mode`. It clears existing aggregation tables, so run
aggregation again afterward.

``` r

narrow <- add_neighborhoods(
  NICHESv2_example,
  edge.filter.mode = "central_cell_only",
  overwrite = TRUE,
  verbose = FALSE
)
narrow <- add_neighborhood_edge_agg(narrow, verbose = FALSE)
```
