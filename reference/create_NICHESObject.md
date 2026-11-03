# Create a NICHESObject from Single-Cell Data

The primary entry point for NICHESv2. Supports two analysis modes:

- `"spatial"`:

  Edges defined by physical proximity using kNN or radius in x/y space.
  Requires `x` and `y` in `meta.data`.

- `"nonspatial"`:

  NICHES v1-style type-stratified sampling. Requires `cell.type.col`.

## Usage

``` r
create_NICHESObject(
  count.mtx,
  meta.data,
  LRM.db,
  mode = c("spatial", "nonspatial"),
  species = "human",
  sample.col = NULL,
  cell.type.col = NULL,
  ligand.col = "ligand",
  receptor.col = "receptor",
  k = NULL,
  rad = NULL,
  seed = 42L,
  max.cells = NULL,
  method = c("product", "minimum"),
  normalize.method = c("prop", "none"),
  gene.meta = NULL,
  n.cores = 1L,
  verbose = TRUE
)
```

## Arguments

- count.mtx:

  Gene x cell matrix (sparse `dgCMatrix` or dense). Row names = gene
  symbols; column names = cell barcodes. Both name sets must match
  `rownames(meta.data)` and the columns of `LRM.db` respectively.

- meta.data:

  `data.frame` with row names = cell barcodes. Must contain `x` and `y`
  columns when `mode = "spatial"`, and the column named by
  `cell.type.col` when `mode = "nonspatial"`.

- LRM.db:

  Either a character string naming a built-in database
  (`"connectomedb2025"`, `"fantom5"` or `"omnipath"`) or a `data.frame`
  with at least the columns named by `ligand.col` and `receptor.col`.
  See
  [`load_LRM_database`](https://raredonlab.github.io/NICHESv2/reference/load_LRM_database.md)
  for details and supported species.

- mode:

  Character. `"spatial"` (default, kNN/radius edges) or `"nonspatial"`
  (type-stratified random sampling).

- species:

  Species for built-in LRM databases. One of `"human"` (default),
  `"mouse"`, `"rat"`, or `"pig"` (`"pig"` is not available from
  OmniPath). Silently ignored when `LRM.db` is a `data.frame`.

- sample.col:

  Character or `NULL` (default). Name of a column in `meta.data`
  identifying the biological sample of origin for each cell. When
  supplied, the pipeline runs independently within each sample and the
  per-sample objects are merged. See the Sample-aware processing section
  above for the barcode-uniqueness requirement.

- cell.type.col:

  Character or `NULL`. Column in `meta.data` holding cell type labels.
  **Required for `mode = "nonspatial"`**; ignored otherwise.

- ligand.col:

  Character. Column name for ligands in `LRM.db`. Default `"ligand"`.

- receptor.col:

  Character. Column name for receptors in `LRM.db`. Default
  `"receptor"`.

- k:

  Integer or `NULL`. Number of nearest neighbors for kNN edges.
  **Spatial mode only; supply exactly one of `k` or `rad`.**

- rad:

  Numeric or `NULL`. Neighborhood radius for fixed-radius edges (same
  units as the `x`/`y` columns). **Spatial mode only; supply exactly one
  of `k` or `rad`.** Use
  [`Explore_Radius_Neighborhood`](https://raredonlab.github.io/NICHESv2/reference/Explore_Radius_Neighborhood.md)
  to choose a sensible value.

- seed:

  Integer. Random seed for reproducible type-stratified sampling.
  Default `42L`. Nonspatial mode only.

- max.cells:

  Integer or `NULL`. Per cell-type-pair edge cap. Default `NULL` (no
  cap). Nonspatial mode only.

- method:

  Per-edge LR scoring method. `"product"` (default) multiplies sender
  ligand expression by receiver receptor expression; `"minimum"` takes
  the elementwise minimum. Both handle gene complexes (underscore-joined
  symbols).

- normalize.method:

  Edge-level normalization for `score.norm` in `$edge.data`. `"prop"`
  (default) scales each edge's per-LRM scores to sum to 1 across LRMs;
  `"none"` leaves `score.norm` equal to `score`.

- gene.meta:

  Optional `data.frame` with a `"gene"` column (and any number of
  additional annotation columns) joined onto the `$gene.meta` slot.
  Default `NULL`.

- n.cores:

  Integer. Number of CPU cores for the inner scoring pipeline. Default
  `1L`. Always uses socket clusters
  ([`parallel::makeCluster`](https://rdrr.io/r/parallel/makeCluster.html))
  for cross-platform compatibility.

- verbose:

  Logical. Print progress messages. Default `TRUE`.

## Value

A `NICHESObject` with attribute `"mode"` set to either `"spatial"` or
`"nonspatial"`. When `sample.col` is provided, attribute `"sample.col"`
is also set to the column name. Slots `$cell.list`, `$cell.data`,
`$cell.meta`, `$edge.list`, `$edge.data`, `$edge.meta`, `$gene.list`,
`$gene.meta`, `$LRM.list`, and `$LRM.meta` are always populated. The
three neighborhood slots are always populated and
`attr(obj, "edge.filter.mode")` and `attr(obj, "remove.autocrine")` are
set on the object. `$segmentation` always remains `NULL` until
[`add_segmentation`](https://raredonlab.github.io/NICHESv2/reference/add_segmentation.md)
is called.

## Pipeline overview

Internally, `create_NICHESObject()` chains the following steps:

1.  Load the LR database
    ([`load_LRM_database`](https://raredonlab.github.io/NICHESv2/reference/load_LRM_database.md))
    when `LRM.db` is a character string; otherwise the supplied
    `data.frame` is used directly.

2.  Build the edge list –
    [`create_Edgelist_Spatial`](https://raredonlab.github.io/NICHESv2/reference/create_Edgelist_Spatial.md)
    for `mode = "spatial"` or
    [`create_Edgelist_SamplingByType`](https://raredonlab.github.io/NICHESv2/reference/create_Edgelist_SamplingByType.md)
    for `mode = "nonspatial"`.

3.  Score every edge x LRM combination via
    [`compute_CellToCell`](https://raredonlab.github.io/NICHESv2/reference/compute_CellToCell.md),
    sparsifying as it goes.

4.  Assemble the 15-slot `NICHESObject` (cell, edge, gene, LRM,
    segmentation, and aggregation slots).

5.  Extract neighborhoods via
    [`add_neighborhoods`](https://raredonlab.github.io/NICHESv2/reference/add_neighborhoods.md)
    using `edge.filter.mode = "both_in_neighborhood"` and
    `remove.autocrine = FALSE` (fixed by design).

Each step is also exported individually for users who need finer
control.

## Sample-aware processing

When `sample.col` is supplied, `meta.data` is split by that column and
the full pipeline runs *independently within each sample*. Per-sample
objects are then combined via
[`merge_NICHESObjects`](https://raredonlab.github.io/NICHESv2/reference/merge_NICHESObjects.md).
This prevents physically meaningless cross-sample edges (a kNN search on
coordinates pooled across slides will happily connect cells that were
never on the same tissue).

**Barcodes must be globally unique across samples.** A duplicate barcode
is a hard error during the merge step. Prepend a sample-specific prefix
before construction (e.g.\\ `paste0(sample.id, "_", barcode)`, which is
what
[`Seurat::RenameCells()`](https://satijalab.github.io/seurat-object/reference/RenameCells.html)
does automatically).

Outer iteration over samples is always sequential; the inner pipeline
uses `n.cores` for parallelism, so socket clusters are never nested.

## See also

[`NICHESv2_example`](https://raredonlab.github.io/NICHESv2/reference/NICHESv2_example.md)
for an example object produced by this function;
[`merge_NICHESObjects`](https://raredonlab.github.io/NICHESv2/reference/merge_NICHESObjects.md)
for combining objects produced independently;
[`compute_CellToCell`](https://raredonlab.github.io/NICHESv2/reference/compute_CellToCell.md)
for the underlying scoring engine;
[`create_Edgelist_Spatial`](https://raredonlab.github.io/NICHESv2/reference/create_Edgelist_Spatial.md),
[`create_Edgelist_SamplingByType`](https://raredonlab.github.io/NICHESv2/reference/create_Edgelist_SamplingByType.md)
for the edge-list constructors invoked internally;
[`load_LRM_database`](https://raredonlab.github.io/NICHESv2/reference/load_LRM_database.md)
for the LR-database loader;
[`aggregate_NICHESObject`](https://raredonlab.github.io/NICHESv2/reference/aggregate_NICHESObject.md)
for downstream summarization.

## Examples

``` r
# \donttest{
data(NICHESv2_inputs)

# ---- Basic spatial mode (kNN edges) ----
obj <- create_NICHESObject(
  count.mtx    = NICHESv2_inputs$count.mtx,
  meta.data    = NICHESv2_inputs$meta.data,
  LRM.db       = NICHESv2_inputs$LRM.db,
  mode         = "spatial",
  ligand.col   = "ligand",
  receptor.col = "receptor",
  k            = 5L,
  n.cores      = 1L,
  verbose      = FALSE
)
#> [create_Edgelist_Spatial] Mode: kNN    | Param: 5 | Cells: 80 | Inter-cell: 274 | Autocrine: 80
#> [compute_CellToCell] 354 edges x 15 LR pairs | method = 'product' | cores = 1
#> [compute_CellToCell] Done. 190 non-zero edge signals.
print(obj)
#> -- NICHESObject --------------------------------------------
#>   Cells        : 80  (80 non-empty)
#>   Genes        : 50  (50 valid)
#>   Edges        : 354  (141 non-empty)
#>   LRM pairs    : 15  (15 valid)
#>   Segmentation : none
#>   Neighborhoods: 80 cells  |  1364 edges  (mode: both_in_neighborhood, autocrine: included)
#>   Aggregations : none
#> 
#>   Slot summary:
#>     $cell.list     -- character vector of cell barcodes
#>     $cell.data     -- (cell, gene, score, score.norm) [sparse]
#>     $cell.meta     -- per-cell metadata
#>     $edge.list     -- (edge, SendingCell, ReceivingCell)
#>     $edge.data     -- (edge, LRM, score, score.norm) [sparse]
#>     $edge.meta     -- per-edge metadata with Sending/Receiving/Joint cols
#>     $neighborhood.cell.list -- named list: cell -> neighbor barcodes
#>     $neighborhood.edge.list -- named list: cell -> edge IDs
#>     $neighborhood.meta      -- (cell, n.neighbors, n.edges) [keyed]
#>     $gene.list     -- character vector of gene names
#>     $gene.meta     -- per-gene metadata
#>     $LRM.list      -- character vector of 'ligand|receptor' pairs
#>     $LRM.meta      -- per-LRM metadata
#>     $segmentation  -- sf object of cell polygons, or NULL
#>     $aggregations  -- named list of aggregation data.tables
#> ------------------------------------------------------------
attr(obj, "mode")           # "spatial"
#> [1] "spatial"

# ---- Nonspatial mode (NICHES v1-style type-stratified sampling) ----
obj.sampled <- create_NICHESObject(
  count.mtx     = NICHESv2_inputs$count.mtx,
  meta.data     = NICHESv2_inputs$meta.data,
  LRM.db        = NICHESv2_inputs$LRM.db,
  mode          = "nonspatial",
  cell.type.col = "celltype",
  ligand.col    = "ligand",
  receptor.col  = "receptor",
  n.cores       = 1L,
  verbose       = FALSE
)
#> [create_Edgelist_SamplingByType] 3 cell types | seed = 42 | max.cells = none | 230 sampled directed edges.
#> [compute_CellToCell] 230 edges x 15 LR pairs | method = 'product' | cores = 1
#> [compute_CellToCell] Done. 123 non-zero edge signals.
attr(obj.sampled, "mode")   # "nonspatial"
#> [1] "nonspatial"

# ---- Sample-aware processing (sample.col) ----
# Synthesise a 'sample' column to demonstrate the split-and-merge path.
# Real workflows: ensure barcodes are globally unique across samples.
md <- NICHESv2_inputs$meta.data
md$sample <- rep(c("S1", "S2"), each = nrow(md) / 2L)

obj.split <- create_NICHESObject(
  count.mtx    = NICHESv2_inputs$count.mtx,
  meta.data    = md,
  LRM.db       = NICHESv2_inputs$LRM.db,
  mode         = "spatial",
  sample.col   = "sample",
  ligand.col   = "ligand",
  receptor.col = "receptor",
  k            = 5L,
  n.cores      = 1L,
  verbose      = FALSE
)
#> [create_Edgelist_Spatial] Mode: kNN    | Param: 5 | Cells: 40 | Inter-cell: 142 | Autocrine: 40
#> [compute_CellToCell] 182 edges x 15 LR pairs | method = 'product' | cores = 1
#> [compute_CellToCell] Done. 113 non-zero edge signals.
#> [create_Edgelist_Spatial] Mode: kNN    | Param: 5 | Cells: 40 | Inter-cell: 146 | Autocrine: 40
#> [compute_CellToCell] 186 edges x 15 LR pairs | method = 'product' | cores = 1
#> [compute_CellToCell] Done. 80 non-zero edge signals.
attr(obj.split, "sample.col")   # "sample"
#> [1] "sample"
# }
```
