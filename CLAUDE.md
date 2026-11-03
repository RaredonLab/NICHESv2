# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Package overview

NICHESv2 is an R package for spatial cell-to-cell ligand-receptor communication analysis. It computes LR signaling scores between pairs of cells in spatial and non-spatial single-cell transcriptomics data and organizes results in a memory-efficient S3 object (`NICHESObject`). Version 1.0.0 is being prepared for CRAN, alongside the methods paper. `working-scripts-NW/v1_facts_sheet.md` holds the verified facts that all documents copy.

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

The 15 slots split into five groups:

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

The single-partition pipeline (in order): validate inputs → build edge list → score LR pairs via `compute_CellToCell()` → assemble slots → extract neighborhoods.

Edge list construction:
- **Spatial kNN** (`k=`): `RANN::nn2`; only mutual edges are kept
- **Spatial radius** (`rad=`): `dbscan::frNN`; inherently symmetric
- Both modes add autocrine self-loop edges for every cell

LR scoring (`compute_CellToCell()` → `compute_LRM_Score()`): supports `"product"` (L×R) and `"minimum"` (min(L,R)) methods. Gene complexes take the product or minimum across all component genes.

### Aggregation framework

`$aggregations` starts as an empty list. Aggregations are populated on demand:

- `add_neighborhood_edge_agg()` — adds `$aggregations$neighborhood.edge.agg`: 27-column table with directional per-(cell, LRM) scores (in, out, cross, self.autocrine, other.autocrine; sum/mean/cmean variants and in/out ratio). Uses sparse matrix multiplication internally: five binary A matrices (cells × edges) are built once and multiplied by the score matrix B (edges × LRMs) in batches of `chunk.size` LRMs (default 50). This bounds peak memory while reducing iterations from n.lrms to ceil(n.lrms / chunk.size).
- `add_neighborhood_composition()` — adds `$aggregations$neighborhood.composition`: per-cell neighbor cell-type proportions
- `aggregate_NICHESObject()` calls both as a convenience wrapper; accepts `chunk.size` and passes it through to `add_neighborhood_edge_agg()`

### Neighborhood filtering

`add_neighborhoods()` (called internally at construction) supports three `edge.filter.mode` values: `"both_in_neighborhood"` (default at construction), `"any_in_neighborhood"`, and `"central_cell_only"`. The mode set at construction is stored in the object.

### Subsetting

`[.NICHESObject` filters neighborhood slots and aggregation tables by cell, leaving each retained cell's entries unchanged, so they still describe its whole neighborhood in the full tissue. Edge slots keep only edges between retained cells, so a subset's neighborhood lists can name edges that are absent from `$edge.list`. `add_neighborhood_edge_agg()` and `add_neighborhood_composition()` therefore stop on such a subset through `.stop_if_neighborhood_edges_missing()`, because recomputing would divide partial sums by full edge counts. The workflow is to aggregate first, then subset. `add_neighborhoods(obj, overwrite = TRUE)` treats a subset as a separate tissue and clears the carried-over aggregations.

### Metadata operations — important distinction

- `add_cell_meta()` **fully replaces** `$cell.meta` and rebuilds `$edge.meta` from scratch. Use only when replacing the entire metadata table.
- `add_meta_column()` appends or overwrites a single column without touching anything else.

### Parallelism

Always use `parallel::makeCluster()` (socket clusters). **Never use `mclapply`** — it is not cross-platform. The `n.cores` parameter in user-facing functions wires directly to socket clusters. Workers call `data.table::setDTthreads(1)`. The package must not change `data.table`'s thread setting on load, because each worker loads the namespace and would lose that limit.

**Exception**: `add_neighborhood_edge_agg()` uses sparse matrix multiplication (Matrix package) instead of socket clusters. Its `n.cores` parameter is retained for API compatibility only and is not used. New functions that adopt a similar matrix-based approach should follow the same pattern: keep `n.cores` in the signature, validate it, but skip cluster creation.

### LR databases

Three built-in options via `load_LRM_database()`:
- `"connectomedb2025"` (default) — bundled `.rda` files, human/mouse/rat/pig, no internet required
- `"fantom5"` — bundled `.rda` files, same four species
- `"omnipath"` — requires internet access; uses `OmnipathR` (Suggests) when available, with a base-R REST fallback

Custom data frames are also accepted; must have `ligand` and `receptor` columns (configurable via `ligand.col`/`receptor.col`).

**OmniPath resilient loading** (ported from NICHESv1 PR #72 / issue #70): newer OmnipathR (≥ ~3.17, incl. 4.1.0) resolves the `organism` argument via `ncbi_taxid()`, which downloads/coalesces species tables from Ensembl, OMA, and UniProt — when any of those ancillary services is unreachable the call aborts with `"Can't combine ..1 and ..3"` even though the OmniPath interaction server is fine. `.load_omnipath()` therefore delegates to `.niches_fetch_omnipath_ligrec()`, which (1) tries the standard OmnipathR path with its logger muted, (2) falls back to a direct OmniPath REST query (`omnipathdb.org/interactions?datasets=ligrecextra&organisms=<taxid>`) using the already-known NCBI taxon id, bypassing the organism-name lookup, and (3) raises an informative network error only if both fail. Fallback is base-R only (`utils::download.file`/`read.delim`), no new dependencies. All four OmniPath helpers live in `R/load_LRM_database.R`.

### Data ingest

- `extract_NICHESInputs_Seurat()` — supports Seurat V4 and V5; detects assay class automatically; **spatial coordinates must be extracted manually** and added to the returned metadata before calling `create_NICHESObject()`
- `extract_NICHESInputs_AnnData()` — reads `.h5ad` via `rhdf5` only; no Python required

## Repository layout notes

- `R/` — all package source; every exported function has its own file
- `data/` — bundled `.rda` files (count matrix, example object, LR databases for 4 species)
- `inst/extdata/` — `NICHESv2_example.h5ad`, a small synthetic file used by the AnnData example
- `data-raw/` — scripts that build the bundled data files; excluded from package build
- `tests/testthat/` — one test file per source file; `helper-testutils.R` provides shared fixtures
- `vignettes/` — one vignette, `NICHESv2.Rmd`, which runs on bundled data; `vignettes/articles/` holds the website pages and is excluded from the package build
- `.github/workflows/` — `R-CMD-check.yaml` runs the checks; `pkgdown.yaml` builds the website and publishes it from `main`, removing this file first because pkgdown turns top-level Markdown files into pages
- `demo-scripts/` — standalone R scripts demonstrating real-dataset usage (not part of the package build)
- `working-scripts-MSBR/` and `working-scripts-NW/` — exploratory and prototyping scripts; not part of the package
- `package management/` — design documents
