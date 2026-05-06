# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Package overview

NICHESv2 is an R package for spatial cell-to-cell ligand-receptor communication analysis. It computes LR signalling scores between neighbouring cells in single-cell spatial transcriptomics data and organises results in a memory-efficient S3 object (`NICHESObject`). It is in active development (v0.0.9000); a formal 0.1.0 release will accompany the methods paper.

## Development commands

All standard package development goes through `devtools`:

```r
devtools::load_all()      # load package in-session during development
devtools::document()      # regenerate NAMESPACE and man/ from roxygen2 comments
devtools::test()          # run full test suite
devtools::check()         # R CMD CHECK (includes tests, examples, vignettes)
devtools::build()         # build source tarball
```

Run a single test file:
```r
testthat::test_file("tests/testthat/test-create_NICHESObject.R")
```

Tests auto-generate fixtures on first run via `tests/testthat.R`.

## Architecture

### NICHESObject — 15-slot S3 class

The central data structure is defined in `R/NICHESObject_class.R`. It is a plain S3 list (not R6/Reference), so **it is copied on every assignment** — always reassign the result of any mutating function (`obj <- add_meta_column(obj, ...)`).

The 15 slots split into three groups:

- **Cell level**: `$cell.list`, `$cell.data`, `$cell.meta`
- **Edge level**: `$edge.list`, `$edge.data`, `$edge.meta`
- **Neighborhood level**: `$neighborhood.cell.list`, `$neighborhood.edge.list`, `$neighborhood.meta`
- **LRM / gene level**: `$gene.list`, `$gene.meta`, `$LRM.list`, `$LRM.meta`
- **Extras**: `$segmentation` (sf polygons), `$aggregations` (named list, empty at construction)

Two object attributes (not slots) track mode: `attr(obj, "mode")` (`"spatial"` or `"nonspatial"`) and `attr(obj, "sample.col")`.

### Sparse long-format storage

`$cell.data` and `$edge.data` are **keyed `data.table` objects in long format** — one row per (cell, gene) or (edge, LRM). This scales with interaction density rather than total dimensions and is fundamental to the package's memory model. Do not convert these to wide matrices without an explicit reason.

Key string separators (used throughout):
- Edge IDs: `"SendingCell|ReceivingCell"`
- LRM IDs: `"ligand|receptor"`
- Gene complexes: `"GeneA_GeneB"` (underscore; handled automatically in `compute_LRM_Score()`)

### Construction pipeline

`create_NICHESObject()` (`R/create_NICHESObject.R`) is the public entry point. It dispatches to:

- `create_NICHESObject_single()` — one sample or no `sample.col`
- `create_NICHESObject_by_sample()` — splits by `sample.col`, processes each partition independently, merges via `merge_NICHESObjects()`

**Barcodes must be globally unique across all samples** when using `sample.col`.

The single-partition pipeline (in order): validate inputs → build edge list → optional smoothing → score LR pairs via `compute_CellToCell()` → assemble slots → extract neighborhoods.

Edge list construction:
- **Spatial kNN** (`k=`): `RANN::nn2`; only mutual edges are kept
- **Spatial radius** (`rad=`): `dbscan::frNN`; inherently symmetric
- Both modes add autocrine self-loop edges for every cell

LR scoring (`compute_CellToCell()` → `compute_LRM_Score()`): supports `"product"` (L×R) and `"minimum"` (min(L,R)) methods. Gene complexes take the product or minimum across all component genes.

### Aggregation framework

`$aggregations` starts as an empty list. Aggregations are populated on demand:

- `add_neighborhood_edge_agg()` — adds `$aggregations$neighborhood.edge.agg`: 27-column table with directional per-(cell, LRM) scores (in, out, cross, self.autocrine, other.autocrine; sum/mean/cmean variants and in/out ratio)
- `add_neighborhood_composition()` — adds `$aggregations$neighborhood.composition`: per-cell neighbor cell-type proportions
- `aggregate_NICHESObject()` calls both as a convenience wrapper

### Neighborhood filtering

`add_neighborhoods()` (called internally at construction) supports three `edge.filter.mode` values: `"both_in_neighborhood"` (default at construction), `"any_in_neighborhood"`, and `"central_cell_only"`. The mode set at construction is stored in the object.

### Metadata operations — important distinction

- `add_cell_meta()` **fully replaces** `$cell.meta` and rebuilds `$edge.meta` from scratch. Use only when replacing the entire metadata table.
- `add_meta_column()` appends or overwrites a single column without touching anything else.

### Parallelism

Always use `parallel::makeCluster()` (socket clusters). **Never use `mclapply`** — it is not cross-platform. The `n.cores` parameter in user-facing functions wires directly to socket clusters.

### LR databases

Three built-in options via `load_LRM_database()`:
- `"connectomedb2025"` (default) — bundled `.rda` files, human/mouse/rat/pig, no internet required
- `"fantom5"` — bundled `.rda` files, same four species
- `"omnipath"` — requires `OmnipathR` (Suggests) and internet access

Custom data frames are also accepted; must have `ligand` and `receptor` columns (configurable via `ligand.col`/`receptor.col`).

### Data ingest

- `extract_NICHESInputs_Seurat()` — supports Seurat V4 and V5; detects assay class automatically; **spatial coordinates must be extracted manually** and added to the returned metadata before calling `create_NICHESObject()`
- `extract_NICHESInputs_AnnData()` — reads `.h5ad` via `rhdf5` only; no Python required

### Known issues

- `smooth = TRUE` is untested end-to-end in `create_NICHESObject()` (nonspatial mode only)
- `[.NICHESObject` subsetting does not propagate `mode` or `sample.col` attributes

## Repository layout notes

- `R/` — all package source; every exported function has its own file
- `data/` — bundled `.rda` files (count matrix, example object, LR databases for 4 species)
- `data-raw/` — scripts that build the bundled `.rda` files; excluded from package build
- `tests/testthat/` — one test file per source file; `helper-testutils.R` provides shared fixtures
- `vignettes/` — three vignettes: `getting-started.Rmd`, `spatial.Rmd`, `scrna-sampled.Rmd`
- `demo-scripts/` — standalone R scripts demonstrating real-dataset usage (not part of the package build)
- `working-scripts-MSBR/` and `working-scripts-NW/` — exploratory and prototyping scripts; not part of the package
- `package management/` — design documents
