# Build a Type-Stratified Sampled Edge List for Non-Spatial scRNA Data

Constructs a directed edge list for non-spatial single-cell
communication analysis by randomly sampling cell pairs within each
ordered cell type crossing. For every directed type pair (A -\> B), the
function:

1.  Computes a per-crossing cap: `min(n_A, n_B)` (further bounded by
    `max.cells` if supplied).

2.  Independently samples `cap` cell barcodes from type A (senders) and
    `cap` barcodes from type B (receivers) without replacement.

3.  Pairs them positionally – the k-th sampled sender is paired with the
    k-th sampled receiver.

All T^2 directed crossings are covered, including autocrine (A -\> A),
where both sender and receiver are sampled from the same type.

**Important:** Individual edges are biologically synthetic – the
specific cell pairing within a type crossing is arbitrary. Only
aggregate signals across all edges of a given type pair carry biological
meaning. Results differ across runs unless `seed` is fixed.

## Usage

``` r
create_Edgelist_SamplingByType(
  meta.data,
  cell.type.col,
  seed = 42L,
  max.cells = NULL
)
```

## Arguments

- meta.data:

  `data.frame` with one row per cell. Row names must be unique cell
  barcodes. Must contain the column specified by `cell.type.col`.

- cell.type.col:

  `character(1)`. Name of the column in `meta.data` that holds cell type
  labels.

- seed:

  `integer(1)`. Random seed for reproducibility. Default `42L`. Pass a
  different value to obtain an alternative sampling draw; pass `NULL` to
  skip seed-setting (not recommended for publication). The caller's
  random number state is restored on exit.

- max.cells:

  `integer(1)` or `NULL`. Optional upper cap on the number of edges per
  type-pair crossing. When `NULL` (default), the cap is `min(n_i, n_j)`
  for each crossing. When set, the cap becomes
  `min(n_i, n_j, max.cells)`, preventing very large types from producing
  an excessive number of edges even when paired with another large type.

## Value

A `data.table` keyed on `(from, to)` with columns:

- `from`:

  Sending cell barcode (character).

- `to`:

  Receiving cell barcode (character).

Carries metadata attributes: `mode` (`"nonspatial"`), `param` (named
list of `cell.type.col`, `seed`, `max.cells`), `n.cells` (total cells in
`meta.data`), `cell.names` (all cell barcodes).

## Details

This function is a low-level primitive: it builds the sampled
non-spatial edge list that
[`create_NICHESObject`](https://raredonlab.github.io/NICHESv2/reference/create_NICHESObject.md)
consumes internally when called with `mode = "nonspatial"`. Most users
should call `create_NICHESObject(mode = "nonspatial", ...)`, which
invokes this function and assembles the result into a full
`NICHESObject` in one step. `create_Edgelist_SamplingByType()` is
exposed for power users who want to inspect the sampled edges, fix the
random seed explicitly outside the wrapper, or feed the edges into a
custom scoring pipeline.

Cell types with zero cells after barcode-matrix intersection are
silently skipped. Type crossings where the effective cap is zero (e.g.,
one type has no cells after filtering) are also silently skipped.

## See also

[`create_NICHESObject`](https://raredonlab.github.io/NICHESv2/reference/create_NICHESObject.md),
[`create_Edgelist_Spatial`](https://raredonlab.github.io/NICHESv2/reference/create_Edgelist_Spatial.md)

## Examples

``` r
set.seed(1)
meta <- data.frame(
  cell.type = sample(c("ATI", "gCap", "Macro"), 300, replace = TRUE),
  row.names = paste0("cell_", 1:300)
)
el <- create_Edgelist_SamplingByType(meta, cell.type.col = "cell.type",
                                     seed = 42L, max.cells = 50L)
#> [create_Edgelist_SamplingByType] 3 cell types | seed = 42 | max.cells = 50 | 450 sampled directed edges.
```
