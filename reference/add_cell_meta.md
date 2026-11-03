# Replace `$cell.meta` on a NICHESObject

Fully replaces `niches.obj$cell.meta` with `new.meta` and (by default)
rebuilds `$edge.meta` so the new columns propagate as `Sending.*` /
`Receiving.*` / `Joint.*` columns on every edge.

The replacement is complete: any column present in the existing
`$cell.meta` but absent from `new.meta` is dropped. The user is
responsible for supplying all columns they wish to retain, including
spatial coordinates (`x`, `y`).

`new.meta` must cover every cell in `niches.obj$cell.list`. Extra rows
(cells not in `$cell.list`) are silently dropped with a message when
`verbose = TRUE`.

## Usage

``` r
add_cell_meta(niches.obj, new.meta, rebuild.edge.meta = TRUE, verbose = TRUE)
```

## Arguments

- niches.obj:

  A `NICHESObject`.

- new.meta:

  A `data.frame` or `data.table` with either a `"cell"` column or
  rownames matching barcodes in `niches.obj$cell.list`.

- rebuild.edge.meta:

  Logical. If `TRUE` (default), regenerates `niches.obj$edge.meta` via
  `.build_edge_meta()` so new metadata columns appear on each edge. Set
  to `FALSE` when calling repeatedly in batch – call once manually at
  the end.

- verbose:

  Logical. If `TRUE` (default), emits messages for dropped extra cells
  and stale aggregation slots.

## Value

Updated `NICHESObject`.

## See also

[`add_meta_column`](https://raredonlab.github.io/NICHESv2/reference/add_meta_column.md)
to append or overwrite a single column without replacing the full table;
[`add_neighborhood_composition`](https://raredonlab.github.io/NICHESv2/reference/add_neighborhood_composition.md)
and
[`aggregate_NICHESObject`](https://raredonlab.github.io/NICHESv2/reference/aggregate_NICHESObject.md)
to recompute stale aggregations.

## Examples

``` r
# \donttest{
data(NICHESv2_example)

# Build a complete replacement metadata table
new.meta <- data.frame(
  cell     = NICHESv2_example$cell.list,
  x        = NICHESv2_example$cell.meta$x,
  y        = NICHESv2_example$cell.meta$y,
  celltype = NICHESv2_example$cell.meta$celltype,
  region   = sample(c("core", "edge"),
                    length(NICHESv2_example$cell.list),
                    replace = TRUE),
  stringsAsFactors = FALSE
)
obj.updated <- add_cell_meta(NICHESv2_example, new.meta)

# New column propagates to edges
head(obj.updated$edge.meta[, c("SendingCell", "ReceivingCell",
                               "Joint.region")])
#>    SendingCell ReceivingCell Joint.region
#>         <char>        <char>       <char>
#> 1:     Cell_10        Cell_1    edge|core
#> 2:     Cell_10       Cell_10    edge|edge
#> 3:     Cell_10       Cell_14    edge|core
#> 4:     Cell_10       Cell_25    edge|core
#> 5:     Cell_10        Cell_4    edge|core
#> 6:     Cell_10        Cell_6    edge|edge
# }
```
