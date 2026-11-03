# Merge multiple NICHESObjects into one

Combines a list of NICHESObjects – or any number passed as individual
arguments – into a single NICHESObject by concatenating all per-cell and
per-edge data. This is the complement to the `sample.col` splitting that
[`create_NICHESObject`](https://raredonlab.github.io/NICHESv2/reference/create_NICHESObject.md)
performs: it can be used to merge objects that were produced
independently and then need to be analyzed together.

## Usage

``` r
merge_NICHESObjects(..., verbose = TRUE)
```

## Arguments

- ...:

  Either (a) multiple NICHESObjects passed directly as positional
  arguments, or (b) a single named or unnamed list of NICHESObjects,
  or (c) an explicit `list(a, b, c)` call.

- verbose:

  Logical. Whether to emit progress messages and warnings. Default
  `TRUE`.

## Value

A single `NICHESObject`.

## Details

**Mode compatibility.** All objects must share the same
`attr(., "mode")` (`"spatial"` or `"nonspatial"`). A mismatch is a hard
error because the slot schemas differ (spatial objects carry `x`/`y` in
`$cell.meta`; edge semantics differ).

**Barcode uniqueness.** Cell barcodes must be globally unique across all
objects. If the same barcode string appears in two or more objects the
function stops with a clear error. Prepend a sample-specific prefix
before constructing the objects (e.g. `paste0(sample.id, "_", barcode)`)
to avoid collisions – Seurat does this automatically via
`RenameCells()`.

**Gene / LRM union.** When gene or LRM lists differ across objects,
their union is taken. Genes absent in one object are implicitly zero in
`$cell.data` (which is sparse – zero rows are never stored).

**Slot-by-slot strategy:**

- `$cell.list`, `$cell.data`, `$edge.list`, `$edge.data` – simple
  concatenation.

- `$cell.meta`, `$edge.meta` – concatenated with `fill = TRUE`; extra
  columns in some objects become `NA` in others.

- `$gene.list`, `$LRM.list` – union.

- `$gene.meta`, `$LRM.meta` – concatenated, deduplicated by key (first
  occurrence kept).

- `$segmentation` – merged via `sf::rbind()` when all objects have it;
  `NULL` with a warning when only some do.

- Neighborhood slots – always `NULL` (reserved).

## See also

[`create_NICHESObject`](https://raredonlab.github.io/NICHESv2/reference/create_NICHESObject.md)
for the `sample.col` argument that automates splitting and merging
within one call.

## Examples

``` r
data(NICHESv2_inputs)
build_sample <- function(cells) {
  create_NICHESObject(
    count.mtx = NICHESv2_inputs$count.mtx[, cells],
    meta.data = NICHESv2_inputs$meta.data[cells, ],
    LRM.db    = NICHESv2_inputs$LRM.db,
    mode      = "spatial",
    k         = 3L,
    verbose   = FALSE
  )
}
barcodes <- rownames(NICHESv2_inputs$meta.data)
obj.s1 <- build_sample(barcodes[1:40])
#> [create_Edgelist_Spatial] Mode: kNN    | Param: 3 | Cells: 40 | Inter-cell: 74 | Autocrine: 40
#> [compute_CellToCell] 114 edges x 15 LR pairs | method = 'product' | cores = 1
#> [compute_CellToCell] Done. 72 non-zero edge signals.
obj.s2 <- build_sample(barcodes[41:80])
#> [create_Edgelist_Spatial] Mode: kNN    | Param: 3 | Cells: 40 | Inter-cell: 84 | Autocrine: 40
#> [compute_CellToCell] 124 edges x 15 LR pairs | method = 'product' | cores = 1
#> [compute_CellToCell] Done. 57 non-zero edge signals.

# Individual objects
merged <- merge_NICHESObjects(obj.s1, obj.s2, verbose = FALSE)

# A pre-built list
obj.list <- list(obj.s1, obj.s2)
merged <- merge_NICHESObjects(obj.list, verbose = FALSE)
```
