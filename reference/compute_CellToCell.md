# Compute Cell-to-Cell Ligand-Receptor Signals

Universal scoring engine for all NICHESv2 edge modes. For each edge and
each valid LR pair, computes a signal score. The edge list may come from
either edgelist constructor – spatial or type-sampled – and the count
matrix column names must be the cell barcodes it uses.

Supports **multi-gene complexes**: a ligand or receptor represented as
`"gene1_gene2"` is split on `"_"` and all component expressions are
combined via the chosen method (product or minimum) *across all ligand
and receptor components together*.

Built-in LRM databases can be loaded by name (see
[`load_LRM_database`](https://raredonlab.github.io/NICHESv2/reference/load_LRM_database.md)):

- `"fantom5"` – Ramilowski et al. 2015. Bundled with NICHESv2 (human,
  mouse, rat, pig). No internet or extra packages required.

- `"omnipath"` – OmniPath intercell network (Turei et al. 2021).
  Requires the OmnipathR package and internet access.

## Usage

``` r
compute_CellToCell(
  count.mtx,
  edge.list,
  LRM.db,
  ligand.col = "ligand",
  receptor.col = "receptor",
  method = c("product", "minimum"),
  n.cores = 1L
)
```

## Arguments

- count.mtx:

  Gene x cell expression matrix (sparse `dgCMatrix` or dense). Row names
  = gene symbols; column names = cell barcodes matching `edge.list$from`
  / `$to`.

- edge.list:

  `data.table` from
  [`create_Edgelist_Spatial`](https://raredonlab.github.io/NICHESv2/reference/create_Edgelist_Spatial.md)
  or
  [`create_Edgelist_SamplingByType`](https://raredonlab.github.io/NICHESv2/reference/create_Edgelist_SamplingByType.md);
  or any `data.frame` with `from` and `to` columns whose values match
  `colnames(count.mtx)`. Coordinate columns (`x.from`, `y.from`, `x.to`,
  `y.to`), when present, are propagated to the output.

- LRM.db:

  LR pair table (`data.frame` or `data.table`) **or** a character string
  (`"fantom5"` or `"omnipath"`) loaded via
  [`load_LRM_database`](https://raredonlab.github.io/NICHESv2/reference/load_LRM_database.md).
  When a string is passed, `ligand.col` and `receptor.col` are silently
  overridden to the canonical `"ligand"` and `"receptor"`.

- ligand.col:

  Character. Column name for ligand symbols in `LRM.db`. Default
  `"ligand"`. Ignored when `LRM.db` is a built-in database name.

- receptor.col:

  Character. Column name for receptor symbols in `LRM.db`. Default
  `"receptor"`. Ignored when `LRM.db` is a built-in database name.

- method:

  Per-edge LR scoring method. `"product"` (default) multiplies sender
  ligand expression by receiver receptor expression; `"minimum"` takes
  the elementwise minimum across all ligand and receptor components. See
  the Complex scoring section.

- n.cores:

  Integer. CPU cores for parallel LRM-chunk processing. Default `1L`.
  See the Parallelization section.

## Value

A plain named list with two elements:

- `$edges`:

  `data.table`: `from`, `to`, `ligand`, `receptor`, `LRM`
  (`"ligand|receptor"`), `score`. Additional columns `x.from`, `y.from`,
  `x.to`, `y.to` appear **only when** the input `edge.list` carries
  coordinates (i.e.\\ for spatial edge lists). One row per non-zero
  (edge x LRM) combination – silent pairs are omitted for sparsity.
  Keyed on `(from, to, LRM)`.

- `$lr.meta`:

  `data.table`: `ligand`, `receptor`, `LRM`, plus any extra annotation
  columns present in `LRM.db` (e.g.\\ `pathway`). One row per valid LRM
  pair after gene-availability filtering. Keyed on `(ligand, receptor)`.

Both elements are `data.table`s so downstream joins on `from`/`to`/`LRM`
are O(log n).

## Complex scoring example

LRM `"COL1A1_COL1A2 -> ITGB1"` with `method = "product"`: \$\$score =
\text{expr}\_{send}(COL1A1) \times \text{expr}\_{send}(COL1A2) \times
\text{expr}\_{recv}(ITGB1)\$\$ With `method = "minimum"`: \$\$score =
\min(\text{expr}\_{send}(COL1A1),\\ \text{expr}\_{send}(COL1A2),\\
\text{expr}\_{recv}(ITGB1))\$\$

## Direct vs wrapped use

Most users should call
[`create_NICHESObject`](https://raredonlab.github.io/NICHESv2/reference/create_NICHESObject.md),
which invokes `compute_CellToCell()` internally and assembles the result
into a full `NICHESObject`. Call this function directly when you only
need the raw `(edge, LRM, score)` table – for example, to score a custom
edge list (e.g.\\ a CCC graph from another method) without building the
surrounding object, or to integrate the score table into a non-NICHESv2
downstream pipeline.

## Parallelization

Work is partitioned across LRM pairs into `n.cores` chunks and processed
by a socket cluster
([`parallel::makeCluster`](https://rdrr.io/r/parallel/makeCluster.html)),
which works identically on Windows, macOS, and Linux. Each worker has
its `data.table` thread count throttled to 1 to prevent OpenMP
contention with the outer parallelism. The cluster is torn down via
`on.exit(stopCluster(cl))` so it is released even if scoring errors
midway.

## See also

[`create_NICHESObject`](https://raredonlab.github.io/NICHESv2/reference/create_NICHESObject.md)
for the wrapped pipeline that assembles a full `NICHESObject`;
[`create_Edgelist_Spatial`](https://raredonlab.github.io/NICHESv2/reference/create_Edgelist_Spatial.md),
[`create_Edgelist_SamplingByType`](https://raredonlab.github.io/NICHESv2/reference/create_Edgelist_SamplingByType.md)
for the edge-list constructors;
[`load_LRM_database`](https://raredonlab.github.io/NICHESv2/reference/load_LRM_database.md)
for the LR-database loader.

## Examples

``` r
# \donttest{
data(NICHESv2_inputs)

# Build a spatial kNN edge list (coord must have row names = barcodes)
coord <- data.frame(
  x = NICHESv2_inputs$meta.data$x,
  y = NICHESv2_inputs$meta.data$y,
  row.names = rownames(NICHESv2_inputs$meta.data)
)
edge.list <- create_Edgelist_Spatial(coord = coord, k = 5L, n.cores = 1L)
#> [create_Edgelist_Spatial] Mode: kNN    | Param: 5 | Cells: 80 | Inter-cell: 274 | Autocrine: 80

# Score every edge x LRM pair (returns a plain 2-element list)
result <- compute_CellToCell(
  count.mtx    = NICHESv2_inputs$count.mtx,
  edge.list    = edge.list,
  LRM.db       = NICHESv2_inputs$LRM.db,
  ligand.col   = "ligand",
  receptor.col = "receptor",
  method       = "product",
  n.cores      = 1L
)
#> [compute_CellToCell] 354 edges x 15 LR pairs | method = 'product' | cores = 1
#> [compute_CellToCell] Done. 190 non-zero edge signals.

str(result, max.level = 1)   # $edges and $lr.meta
#> List of 2
#>  $ edges  :Classes ‘data.table’ and 'data.frame':    190 obs. of  10 variables:
#>   ..- attr(*, ".internal.selfref")=<pointer: 0x55aced2a7a30> 
#>   ..- attr(*, "sorted")= chr [1:3] "from" "to" "LRM"
#>  $ lr.meta:Classes ‘data.table’ and 'data.frame':    15 obs. of  4 variables:
#>   ..- attr(*, ".internal.selfref")=<pointer: 0x55aced2a7a30> 
#>   ..- attr(*, "sorted")= chr [1:2] "ligand" "receptor"
head(result$edges)           # one row per non-zero (edge, LRM)
#> Key: <from, to, LRM>
#>       from      to     x.from    y.from       x.to      y.to ligand receptor
#>     <char>  <char>      <num>     <num>      <num>     <num> <char>   <char>
#> 1:  Cell_1  Cell_1 0.13612093 0.2979561 0.13612093 0.2979561 Gene11   Gene12
#> 2:  Cell_1  Cell_1 0.13612093 0.2979561 0.13612093 0.2979561  Gene9   Gene10
#> 3: Cell_10  Cell_1 0.08153932 0.2934977 0.13612093 0.2979561  Gene9   Gene10
#> 4: Cell_10 Cell_10 0.08153932 0.2934977 0.08153932 0.2934977 Gene15   Gene16
#> 5: Cell_10  Cell_4 0.08153932 0.2934977 0.03535238 0.2339655  Gene9   Gene10
#> 6: Cell_12 Cell_12 0.23674938 0.4357651 0.23674938 0.4357651 Gene11   Gene12
#>    score           LRM
#>    <num>        <char>
#> 1:    12 Gene11|Gene12
#> 2:    14  Gene9|Gene10
#> 3:    21  Gene9|Gene10
#> 4:    35 Gene15|Gene16
#> 5:    18  Gene9|Gene10
#> 6:    24 Gene11|Gene12
head(result$lr.meta)         # one row per valid LRM pair
#> Key: <ligand, receptor>
#>    ligand receptor  pathway           LRM
#>    <char>   <char>   <char>        <char>
#> 1:  Gene1    Gene2 PathwayA   Gene1|Gene2
#> 2: Gene11   Gene12 PathwayF Gene11|Gene12
#> 3: Gene13   Gene14 PathwayG Gene13|Gene14
#> 4: Gene15   Gene16 PathwayH Gene15|Gene16
#> 5: Gene17   Gene18 PathwayI Gene17|Gene18
#> 6: Gene19   Gene20 PathwayJ Gene19|Gene20

# Spatial edge lists carry coordinates -- they propagate to $edges
setdiff(c("x.from", "y.from", "x.to", "y.to"), colnames(result$edges))
#> character(0)
# character(0) -- all four coordinate columns are present
# }
```
