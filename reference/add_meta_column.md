# Append or Overwrite a Single Column in `$cell.meta`

Adds one new column to `niches.obj$cell.meta` from a named vector, or
overwrites an existing column of the same name. All other columns in
`$cell.meta` are preserved unchanged.

Use
[`add_cell_meta`](https://raredonlab.github.io/NICHESv2/reference/add_cell_meta.md)
instead when you want to replace the entire metadata table.

`col.values` must be a named vector whose names are cell barcodes. It
must supply a value for every cell in `niches.obj$cell.list`. Extra
names (cells not in `$cell.list`) are silently ignored.

If `col.name` already exists in `$cell.meta`, a message is emitted and
the column is overwritten.

## Usage

``` r
add_meta_column(
  niches.obj,
  col.name,
  col.values,
  rebuild.edge.meta = TRUE,
  verbose = TRUE
)
```

## Arguments

- niches.obj:

  A `NICHESObject`.

- col.name:

  Character scalar. Name of the column to add or overwrite.

- col.values:

  Named vector. Names must be cell barcodes covering every cell in
  `niches.obj$cell.list`. The vector type determines the column type
  (character, numeric, integer, logical, factor).

- rebuild.edge.meta:

  Logical. If `TRUE` (default), regenerates `niches.obj$edge.meta` via
  `.build_edge_meta()` so the new column appears as `Sending.<col.name>`
  / `Receiving.<col.name>` (and `Joint.<col.name>` for categorical
  types) on every edge.

- verbose:

  Logical. If `TRUE` (default), emits a message when overwriting an
  existing column and when aggregation slots may be stale.

## Value

Updated `NICHESObject`.

## See also

[`add_cell_meta`](https://raredonlab.github.io/NICHESv2/reference/add_cell_meta.md)
to replace the full metadata table;
[`add_neighborhood_composition`](https://raredonlab.github.io/NICHESv2/reference/add_neighborhood_composition.md)
and
[`aggregate_NICHESObject`](https://raredonlab.github.io/NICHESv2/reference/aggregate_NICHESObject.md)
to recompute stale aggregations.

## Examples

``` r
# \donttest{
data(NICHESv2_example)

# Add a new numeric column
scores <- setNames(
  runif(length(NICHESv2_example$cell.list)),
  NICHESv2_example$cell.list
)
obj.updated <- add_meta_column(NICHESv2_example,
                               col.name   = "quality.score",
                               col.values = scores)
head(obj.updated$cell.meta[, c("cell", "quality.score")])
#> Key: <cell>
#>       cell quality.score
#>     <char>         <num>
#> 1:  Cell_1    0.52791704
#> 2: Cell_10    0.05861411
#> 3: Cell_11    0.99706914
#> 4: Cell_12    0.14903547
#> 5: Cell_13    0.51855664
#> 6: Cell_14    0.84612005

# Overwrite an existing categorical column
new.types <- setNames(
  sample(c("TypeA", "TypeB"), length(NICHESv2_example$cell.list),
         replace = TRUE),
  NICHESv2_example$cell.list
)
obj.updated <- add_meta_column(NICHESv2_example,
                               col.name   = "celltype",
                               col.values = new.types)
#> [add_meta_column] Column "celltype" already exists in $cell.meta and will be overwritten.
# }
```
