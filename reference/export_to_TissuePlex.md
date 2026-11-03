# Export a NICHESObject to TissuePlex edges.parquet

Writes a Parquet file that TissuePlex can load directly as its edge
layer. Every edge in `$edge.list` is exported. Edges present in
`$edge.data` produce one row per (edge, LRM) pair with real scores;
edges absent from `$edge.data` (no active LR signal) produce a single
placeholder row with `NA` for all LRM fields. This ensures TissuePlex
can render the full tissue graph even for edges with no scored LR pairs.

## Usage

``` r
export_to_TissuePlex(
  object,
  output.path,
  x.col = "x",
  y.col = "y",
  celltype.col = "celltype",
  n.threads = NULL
)
```

## Arguments

- object:

  A `NICHESObject`.

- output.path:

  Full path to the output `.parquet` file. Parent directories are
  created automatically.

- x.col:

  Name of the x-coordinate column in `$cell.meta` (default `"x"`).

- y.col:

  Name of the y-coordinate column in `$cell.meta` (default `"y"`).

- celltype.col:

  Name of the cell-type column in `$cell.meta`, or `NULL` to omit
  cell-type columns (default `"celltype"`).

- n.threads:

  Number of threads for `data.table` operations (`NULL` uses the current
  `getDTthreads()` setting). On a MacBook Pro, setting this to the
  number of performance cores (e.g. 8–12) significantly speeds up joins
  and string operations on large datasets.

## Value

Invisibly returns `output.path`.

## Details

NICHESv2 stops at this file. To view it, see the [TissuePlex
manual](https://raredonlab.github.io/TissuePlex/).

## Examples

``` r
if (requireNamespace("arrow", quietly = TRUE)) {
  data(NICHESv2_example)
  out.file <- tempfile(fileext = ".parquet")
  export_to_TissuePlex(NICHESv2_example, output.path = out.file)
  unlink(out.file)
}
#> Wrote 403 rows | 141 scored edges (15 LRMs) | 213 unscored edges -> /tmp/RtmpZzMH1L/file2141b7991a8.parquet
```
