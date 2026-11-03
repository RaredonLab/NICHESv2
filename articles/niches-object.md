# The NICHESObject

A NICHESObject is an S3 list with 15 slots, which you read with `$`.
This guide uses `NICHESv2_example`, a small spatial object with 80
synthetic cells.

``` r

library(NICHESv2)
library(data.table)

data(NICHESv2_example)
obj <- NICHESv2_example
dim(obj)
#> cells genes edges  LRMs 
#>    80    50   354    15
```

The object is a plain list, so functions return a changed copy. Always
assign the result back: `obj <- add_meta_column(obj, ...)`.

## Why long format

Most cell pairs carry no signal for most LR pairs. NICHESv2 stores only
non-zero values, in keyed `data.table`s with one row per value. Memory
use then follows how much signal exists, not how many combinations are
possible.

## Cell slots

| Slot         | Contents                                      |
|--------------|-----------------------------------------------|
| `$cell.list` | Cell barcodes                                 |
| `$cell.data` | Non-zero expression, one row per (cell, gene) |
| `$cell.meta` | One row per cell                              |

``` r

head(obj$cell.data, 3)
```

| cell   | gene   | score | score.norm |
|:-------|:-------|------:|-----------:|
| Cell_1 | Gene10 |     7 |      0.106 |
| Cell_1 | Gene11 |     6 |      0.091 |
| Cell_1 | Gene12 |     2 |      0.030 |

In `$cell.data`, `score` is the value you supplied, and `score.norm` is
its share of that cell’s total.

## Edge slots

An edge is one directed cell pair, with an ID of the form
`"SendingCell|ReceivingCell"`.

| Slot         | Contents                                    |
|--------------|---------------------------------------------|
| `$edge.list` | One row per edge                            |
| `$edge.data` | Non-zero scores, one row per (edge, LRM)    |
| `$edge.meta` | One row per edge, with both cells’ metadata |

``` r

head(obj$edge.data, 3)
```

| edge             | LRM            | score | score.norm |
|:-----------------|:---------------|------:|-----------:|
| Cell_10\|Cell_1  | Gene9\|Gene10  |    21 |          1 |
| Cell_10\|Cell_10 | Gene15\|Gene16 |    35 |          1 |
| Cell_10\|Cell_4  | Gene9\|Gene10  |    18 |          1 |

An LRM ID has the form `"ligand|receptor"`. With the default
`normalize.method = "prop"`, `score.norm` is that LRM’s share of the
edge’s total score.

`$edge.meta` adds `Sending.` and `Receiving.` versions of each metadata
column, plus `Joint.` columns for text columns:

``` r

obj$edge.meta[, .(edges = .N), by = Joint.celltype][order(-edges)][1:3]
```

| Joint.celltype | edges |
|:---------------|------:|
| TypeA\|TypeA   |   132 |
| TypeB\|TypeB   |   115 |
| TypeC\|TypeC   |   101 |

## Neighborhood slots

These are filled in when the object is built, and all three are indexed
by cell.

| Slot                      | Contents                              |
|---------------------------|---------------------------------------|
| `$neighborhood.cell.list` | Each cell’s neighbors                 |
| `$neighborhood.edge.list` | The edges in each cell’s neighborhood |
| `$neighborhood.meta`      | Neighbor and edge counts per cell     |

``` r

head(obj$neighborhood.meta, 3)
```

| cell    | n.neighbors | n.edges |
|:--------|------------:|--------:|
| Cell_1  |           5 |      17 |
| Cell_10 |           6 |      30 |
| Cell_11 |           6 |      28 |

By default, a neighborhood holds every edge whose two cells are both in
that neighborhood. `add_neighborhoods(overwrite = TRUE)` changes this:

| `edge.filter.mode` | Edges kept for a focal cell |
|----|----|
| `"both_in_neighborhood"` (default) | Both cells are the focal cell or its neighbors |
| `"any_in_neighborhood"` | At least one cell is the focal cell or a neighbor |
| `"central_cell_only"` | One cell is the focal cell |

## Gene and LRM slots

| Slot                       | Contents                       |
|----------------------------|--------------------------------|
| `$gene.list`, `$gene.meta` | Genes and any annotations      |
| `$LRM.list`, `$LRM.meta`   | LR pairs and their annotations |

``` r

head(obj$LRM.meta, 3)
```

| LRM            | ligand | receptor | pathway  |
|:---------------|:-------|:---------|:---------|
| Gene11\|Gene12 | Gene11 | Gene12   | PathwayF |
| Gene13\|Gene14 | Gene13 | Gene14   | PathwayG |
| Gene15\|Gene16 | Gene15 | Gene16   | PathwayH |

## Extras

`$aggregations` starts empty and is filled by
[`aggregate_NICHESObject()`](https://raredonlab.github.io/NICHESv2/reference/aggregate_NICHESObject.md).
`$segmentation` holds optional cell polygons, added with
[`add_segmentation()`](https://raredonlab.github.io/NICHESv2/reference/add_segmentation.md).

## Attributes

Four attributes record how the object was built:

``` r

attr(obj, "mode")
#> [1] "spatial"
attr(obj, "sample.col")
#> NULL
attr(obj, "edge.filter.mode")
#> [1] "both_in_neighborhood"
attr(obj, "remove.autocrine")
#> [1] FALSE
```

## Changing metadata

To add or change one column, use
[`add_meta_column()`](https://raredonlab.github.io/NICHESv2/reference/add_meta_column.md)
with values named by barcode:

``` r

region <- setNames(
  ifelse(obj$cell.meta$x > 0.5, "right", "left"),
  obj$cell.meta$cell
)
obj <- add_meta_column(obj, "region", region, verbose = FALSE)
obj$edge.meta[, .(edges = .N), by = Sending.region]
```

| Sending.region | edges |
|:---------------|------:|
| left           |   188 |
| right          |   166 |

[`add_cell_meta()`](https://raredonlab.github.io/NICHESv2/reference/add_cell_meta.md)
replaces the whole table, so any column you leave out is lost, including
coordinates. Use it only when you mean to replace everything.

## Subsetting

`obj[cells]` takes a vector of barcodes. Each kept cell keeps its
neighborhood entries and aggregation rows unchanged, while the edge
slots keep only edges between kept cells.

``` r

type.a <- obj$cell.meta[celltype == "TypeA", cell]
sub <- obj[type.a]
dim(sub)
#> cells genes edges  LRMs 
#>    30    50   132    15
```

The [aggregation
guide](https://raredonlab.github.io/NICHESv2/articles/aggregation.md)
explains how to work with subsets.

## Merging

[`merge_NICHESObjects()`](https://raredonlab.github.io/NICHESv2/reference/merge_NICHESObjects.md)
combines objects from separate runs. They must share a mode and have no
barcodes in common. The merged object has no neighborhoods or
aggregations, so run
[`add_neighborhoods()`](https://raredonlab.github.io/NICHESv2/reference/add_neighborhoods.md)
afterward. To process several samples in one step, use `sample.col` in
[`create_NICHESObject()`](https://raredonlab.github.io/NICHESv2/reference/create_NICHESObject.md)
instead.
