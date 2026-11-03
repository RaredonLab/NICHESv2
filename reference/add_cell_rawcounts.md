# Add Raw UMI Counts to `$cell.data`

Adds a `count` column to `niches.obj$cell.data` from a user-supplied
gene x cell raw count matrix. The `count` column records the original
sequencing depth for each (cell, gene) pair as a provenance record; it
is not used by any other NICHESv2 function.

## Usage

``` r
add_cell_rawcounts(niches.obj, count.mtx, verbose = TRUE)
```

## Arguments

- niches.obj:

  A `NICHESObject`.

- count.mtx:

  A gene x cell raw count matrix. Either a base `matrix` or a sparse
  `Matrix` object (e.g., `dgCMatrix`). Must have `rownames` (gene
  symbols) and `colnames` (cell barcodes) covering every cell in
  `niches.obj$cell.list`.

- verbose:

  Logical. If `TRUE` (default), emits messages for missing genes, new
  rows added, and an overwrite notice if a `count` column already
  exists.

## Value

Updated `NICHESObject` with a `count` column in `$cell.data`. Row count
may increase if any (cell, gene) pairs had `count > 0` and `score = 0`.

## Row expansion

`$cell.data` is sparse: it originally contains only (cell, gene) pairs
where `score > 0`. When raw counts are provided, new rows are appended
for any (cell, gene) pair where `count > 0` but `score = 0` (e.g., a
gene detected in sequencing but not expressed in the scoring layer). New
rows carry `score = 0` and `score.norm = 0`.

## Gene coverage

`count.mtx` must cover every cell in `niches.obj$cell.list` (hard error
if any are missing). Genes in `niches.obj$gene.list` absent from
`rownames(count.mtx)` receive `count = NA` with a message when
`verbose = TRUE`. Extra rows in `count.mtx` not in `$gene.list` are
silently ignored.

## Validation

Non-zero stored values are checked for validity. Warnings (not errors)
are emitted if any values are negative or non-integer; the pipeline
continues regardless. These checks help detect accidental use of a
normalized layer instead of raw counts.

## See also

[`create_NICHESObject`](https://raredonlab.github.io/NICHESv2/reference/create_NICHESObject.md)
for object construction;
[`extract_NICHESInputs_Seurat`](https://raredonlab.github.io/NICHESv2/reference/extract_NICHESInputs_Seurat.md)
and
[`extract_NICHESInputs_AnnData`](https://raredonlab.github.io/NICHESv2/reference/extract_NICHESInputs_AnnData.md)
for obtaining a count matrix from Seurat or AnnData objects.

## Examples

``` r
# \donttest{
data(NICHESv2_inputs)
data(NICHESv2_example)

# NICHESv2_inputs$count.mtx is the raw count matrix used to build
# NICHESv2_example. Attach it as the provenance count column.
obj <- add_cell_rawcounts(NICHESv2_example, NICHESv2_inputs$count.mtx)
head(obj$cell.data[, c("cell", "gene", "score", "score.norm", "count")])
#> Key: <cell, gene>
#>      cell   gene score score.norm count
#>    <char> <char> <num>      <num> <num>
#> 1: Cell_1 Gene10     7 0.10606061     7
#> 2: Cell_1 Gene11     6 0.09090909     6
#> 3: Cell_1 Gene12     2 0.03030303     2
#> 4: Cell_1 Gene25     5 0.07575758     5
#> 5: Cell_1 Gene33     6 0.09090909     6
#> 6: Cell_1  Gene4     8 0.12121212     8
# }
```
