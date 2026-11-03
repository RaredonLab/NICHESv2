# Exporting to TissuePlex

[TissuePlex](https://github.com/RaredonLab/TissuePlex) is an interactive
viewer that draws cell-cell communication on the tissue image.
NICHESv2’s part ends with one file, `edges.parquet`, which TissuePlex
reads as an edge layer.

## Write the file

[`export_to_TissuePlex()`](https://raredonlab.github.io/NICHESv2/reference/export_to_TissuePlex.md)
needs the arrow package and a spatial object. Set `celltype.col` to your
cell type column, or `NULL` if you don’t have one.

``` r

library(NICHESv2)

data(NICHESv2_example)

out.file <- file.path(tempdir(), "edges.parquet")
export_to_TissuePlex(
  NICHESv2_example,
  output.path  = out.file,
  celltype.col = "celltype"
)
```

Every edge in `$edge.list` is written. Edges with signal get one row per
LR pair, and edges with no signal get one placeholder row, so TissuePlex
can still draw the full tissue graph.

If your coordinates aren’t in columns named `x` and `y`, set `x.col` and
`y.col`.

Exporting a subset writes only the edges between kept cells.

## View it in TissuePlex

Save the file as `edges.parquet` next to your dataset, then follow the
TissuePlex documentation from there:

- [User manual](https://raredonlab.github.io/TissuePlex/): installation
  and loading your data.
- [Worked
  examples](https://raredonlab.github.io/TissuePlex/vignettes.html):
  full walkthroughs.

TissuePlex needs your platform’s own output files as well, such as the
tissue image. These come from your data, not from NICHESv2.
