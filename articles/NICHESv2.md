# Getting started with NICHESv2

NICHESv2 scores ligand-receptor (LR) signaling between individual pairs
of cells. For each directed cell pair and each LR pair, it asks whether
the sending cell expresses the ligand and the receiving cell expresses
the receptor. It keeps every non-zero score, so you can trace a signal
back to the exact cells that produced it.

This vignette runs the whole workflow on a small synthetic dataset that
ships with the package, so every step takes a moment and needs no
download.

``` r

library(NICHESv2)
library(data.table)
```

## Inputs

[`create_NICHESObject()`](https://raredonlab.github.io/NICHESv2/reference/create_NICHESObject.md)
needs three things: a gene-by-cell expression matrix, cell metadata, and
an LR database. The bundled `NICHESv2_inputs` list holds all three, for
50 genes and 80 cells.

``` r

data(NICHESv2_inputs)
names(NICHESv2_inputs)
#> [1] "count.mtx" "meta.data" "LRM.db"

dim(NICHESv2_inputs$count.mtx)
#> [1] 50 80
head(NICHESv2_inputs$meta.data, 3)
```

|        |     x |     y | nCount_RNA | celltype |
|:-------|------:|------:|-----------:|:---------|
| Cell_1 | 0.136 | 0.298 |         66 | TypeA    |
| Cell_2 | 0.238 | 0.172 |         53 | TypeA    |
| Cell_3 | 0.225 | 0.240 |         55 | TypeA    |

``` r

head(NICHESv2_inputs$LRM.db, 3)
```

| ligand | receptor | pathway  |
|:-------|:---------|:---------|
| Gene1  | Gene2    | PathwayA |
| Gene3  | Gene4    | PathwayB |
| Gene5  | Gene6    | PathwayC |

The column names of the matrix must match the row names of the metadata.
Spatial mode needs `x` and `y` columns in the metadata.

## Build a NICHESObject

Here each cell is connected to its spatial neighbors using `k = 5`. Only
mutual nearest neighbors become edges, and every cell also gets a
self-loop for autocrine signaling.

``` r

obj <- create_NICHESObject(
  count.mtx = NICHESv2_inputs$count.mtx,
  meta.data = NICHESv2_inputs$meta.data,
  LRM.db    = NICHESv2_inputs$LRM.db,
  mode      = "spatial",
  k         = 5L,
  method    = "product",
  verbose   = FALSE
)
dim(obj)
#> cells genes edges  LRMs 
#>    80    50   354    15
```

An edge is one directed cell pair. `method = "product"` multiplies the
sender’s ligand expression by the receiver’s receptor expression, and
`"minimum"` takes the smaller of the two. To use a fixed radius instead
of `k`, pass `rad`;
[`Explore_Radius_Neighborhood()`](https://raredonlab.github.io/NICHESv2/reference/Explore_Radius_Neighborhood.md)
helps you choose a value.

For larger datasets, set `n.cores` to score LR pairs in parallel. This
uses socket clusters, so it works on every platform.

## Inside the object

A NICHESObject is a list with 15 slots, which you read with `$`.

| Level | Slots |
|----|----|
| Cell | `$cell.list`, `$cell.data`, `$cell.meta` |
| Edge | `$edge.list`, `$edge.data`, `$edge.meta` |
| Neighborhood | `$neighborhood.cell.list`, `$neighborhood.edge.list`, `$neighborhood.meta` |
| Gene and LRM | `$gene.list`, `$gene.meta`, `$LRM.list`, `$LRM.meta` |
| Extras | `$segmentation`, `$aggregations` |

Scores are stored in long format, with one row for each non-zero (edge,
LRM) pair. This keeps memory use tied to how much signal exists, not to
the number of possible combinations.

``` r

head(obj$edge.data, 3)
```

| edge             | LRM            | score | score.norm |
|:-----------------|:---------------|------:|-----------:|
| Cell_10\|Cell_1  | Gene9\|Gene10  |    21 |          1 |
| Cell_10\|Cell_10 | Gene15\|Gene16 |    35 |          1 |
| Cell_10\|Cell_4  | Gene9\|Gene10  |    18 |          1 |

Edge IDs have the form `"SendingCell|ReceivingCell"` and LRM IDs have
the form `"ligand|receptor"`. `$edge.data` is a `data.table`, so you can
summarize it directly. These are the strongest LR pairs overall:

``` r

obj$edge.data[, .(total.score = sum(score)), by = LRM][order(-total.score)][1:3]
```

| LRM                          | total.score |
|:-----------------------------|------------:|
| Gene27\|Gene28_Gene29        |         632 |
| Gene29_Gene30\|Gene28_Gene30 |         600 |
| Gene25_Gene26\|Gene27        |         574 |

`$edge.meta` carries each cell’s metadata onto its edges, with
`Sending.` and `Receiving.` prefixes, which makes it easy to filter by
cell type:

``` r

obj$edge.meta[, .(edges = .N), by = .(Sending.celltype, Receiving.celltype)][1:3]
```

| Sending.celltype | Receiving.celltype | edges |
|:-----------------|:-------------------|------:|
| TypeA            | TypeA              |   132 |
| TypeB            | TypeB              |   115 |
| TypeB            | TypeC              |     2 |

## Aggregate by neighborhood

Each cell’s neighborhood is filled in when the object is built.
[`aggregate_NICHESObject()`](https://raredonlab.github.io/NICHESv2/reference/aggregate_NICHESObject.md)
summarizes each one per LRM, and with `cell.type.col` it also reports
which cell types surround each cell.

``` r

obj <- aggregate_NICHESObject(obj, cell.type.col = "celltype", verbose = FALSE)
names(obj$aggregations)
#> [1] "neighborhood.edge.agg"    "neighborhood.composition"
```

`$aggregations$neighborhood.edge.agg` has one row per (cell, LRM). It
sorts each neighborhood edge into five categories relative to the focal
cell: `in`, `out`, `cross` (between two neighbors), `self.autocrine` and
`other.autocrine`.

``` r

agg <- obj$aggregations$neighborhood.edge.agg
agg[, .(cell, LRM, score.in.sum, score.out.sum, score.cross.sum)][1:3]
```

| cell   | LRM            | score.in.sum | score.out.sum | score.cross.sum |
|:-------|:---------------|-------------:|--------------:|----------------:|
| Cell_1 | Gene11\|Gene12 |           10 |             0 |               0 |
| Cell_1 | Gene13\|Gene14 |            0 |             0 |               8 |
| Cell_1 | Gene15\|Gene16 |            0 |             0 |               0 |

`$aggregations$neighborhood.composition` has one row per cell:

``` r

comp <- obj$aggregations$neighborhood.composition
comp[, .(cell, TypeA.neighbor.prop, TypeB.neighbor.prop, TypeC.neighbor.prop)][1:3]
```

| cell    | TypeA.neighbor.prop | TypeB.neighbor.prop | TypeC.neighbor.prop |
|:--------|--------------------:|--------------------:|--------------------:|
| Cell_1  |                   1 |                   0 |                   0 |
| Cell_10 |                   1 |                   0 |                   0 |
| Cell_11 |                   1 |                   0 |                   0 |

## Subset to cells of interest

Aggregate first, then subset. Each retained cell keeps its neighborhood
entries and aggregation rows unchanged, so its values still describe its
whole neighborhood in the full tissue, including neighbors you didn’t
keep.

``` r

type.a <- obj$cell.meta[celltype == "TypeA", cell]
sub <- obj[type.a]
dim(sub)
#> cells genes edges  LRMs 
#>    30    50   132    15

identical(
  sub$aggregations$neighborhood.edge.agg[, score.in.sum],
  agg[cell %in% type.a, score.in.sum]
)
#> [1] TRUE
```

The edge slots of a subset hold only edges between retained cells.
Summing a cell’s rows in `sub$edge.data` can therefore give less than
its aggregation row.

Because those removed edges are gone, aggregations can’t be recomputed
on the subset, and trying stops with an error. To analyze the subset as
a separate tissue, rebuild its neighborhoods first:

``` r

separate <- add_neighborhoods(sub, overwrite = TRUE, verbose = FALSE)
separate <- aggregate_NICHESObject(separate, verbose = FALSE)
```

## Data without coordinates

For scRNA-seq data, `mode = "nonspatial"` pairs cells by random sampling
within each pair of cell types. It needs `cell.type.col`, and `seed`
makes the sampling repeatable.

``` r

ns <- create_NICHESObject(
  count.mtx     = NICHESv2_inputs$count.mtx,
  meta.data     = NICHESv2_inputs$meta.data,
  LRM.db        = NICHESv2_inputs$LRM.db,
  mode          = "nonspatial",
  cell.type.col = "celltype",
  verbose       = FALSE
)
dim(ns)
#> cells genes edges  LRMs 
#>    80    50   230    15
```

## Using your own data

### LR databases

Pass a database name to `LRM.db`, and set `species` to match your data.

``` r

db <- load_LRM_database("connectomedb2025", species = "human", verbose = FALSE)
nrow(db)
#> [1] 3550
```

`"connectomedb2025"` (the default in
[`load_LRM_database()`](https://raredonlab.github.io/NICHESv2/reference/load_LRM_database.md))
and `"fantom5"` are bundled for human, mouse, rat and pig. `"omnipath"`
needs internet access and the OmnipathR package, and covers human, mouse
and rat. You can also pass your own data frame. Multi-gene complexes are
written with underscores, such as `"COL1A1_COL1A2"`.

### Seurat objects

[`extract_NICHESInputs_Seurat()`](https://raredonlab.github.io/NICHESv2/reference/extract_NICHESInputs_Seurat.md)
works with Seurat V4 and V5 objects, so you don’t need to update an
older object first. It returns `count.mtx` and `meta.data`, ready for
[`create_NICHESObject()`](https://raredonlab.github.io/NICHESv2/reference/create_NICHESObject.md).
For spatial data, add the coordinates to `meta.data` yourself as `x` and
`y` columns;
[`SeuratObject::GetTissueCoordinates()`](https://satijalab.github.io/seurat-object/reference/GetTissueCoordinates.html)
returns them.

``` r

seu <- SeuratObject::CreateSeuratObject(
  counts    = NICHESv2_inputs$count.mtx,
  meta.data = NICHESv2_inputs$meta.data
)
inputs <- extract_NICHESInputs_Seurat(seu, layer = "counts", verbose = FALSE)
dim(inputs$count.mtx)
#> [1] 50 80
```

The `layer` you choose sets what the scores mean: raw counts give
absolute magnitudes, and normalized expression gives relative scores.

### AnnData files

[`extract_NICHESInputs_AnnData()`](https://raredonlab.github.io/NICHESv2/reference/extract_NICHESInputs_AnnData.md)
reads `.h5ad` files with rhdf5, so it doesn’t need Python.

``` r

h5ad <- system.file("extdata", "NICHESv2_example.h5ad", package = "NICHESv2")
inputs <- extract_NICHESInputs_AnnData(h5ad, verbose = FALSE)
dim(inputs$count.mtx)
#> [1] 50 80
```

### Exporting to TissuePlex

[`export_to_TissuePlex()`](https://raredonlab.github.io/NICHESv2/reference/export_to_TissuePlex.md)
writes one Parquet file of edges and needs the arrow package. To view
the file in tissue coordinates, see the [TissuePlex
manual](https://raredonlab.github.io/TissuePlex/).

``` r

out.file <- tempfile(fileext = ".parquet")
export_to_TissuePlex(obj, output.path = out.file)
unlink(out.file)
```

## Where to go next

See
[`?create_NICHESObject`](https://raredonlab.github.io/NICHESv2/reference/create_NICHESObject.md)
for every construction option,
[`?add_neighborhood_edge_agg`](https://raredonlab.github.io/NICHESv2/reference/add_neighborhood_edge_agg.md)
for all 27 aggregation columns, and
[`?"[.NICHESObject"`](https://raredonlab.github.io/NICHESv2/reference/sub-.NICHESObject.md)
for the subsetting rules.
[`?NICHESv2`](https://raredonlab.github.io/NICHESv2/reference/NICHESv2-package.md)
gives an overview.

The [website](https://raredonlab.github.io/NICHESv2/) has tutorials on
public spatial and non-spatial datasets.
