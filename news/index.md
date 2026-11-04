# Changelog

## NICHESv2 1.0.0

First public release.

### Main features

- [`create_NICHESObject()`](https://raredonlab.github.io/NICHESv2/reference/create_NICHESObject.md)
  scores ligand-receptor (LR) signaling between pairs of cells, in
  `"spatial"` mode (k-nearest neighbors or fixed radius) or
  `"nonspatial"` mode (type-stratified random sampling).
- The 15-slot `NICHESObject` stores non-zero scores in sparse
  long-format `data.table`s, with
  [`print()`](https://rdrr.io/r/base/print.html),
  [`summary()`](https://rdrr.io/r/base/summary.html),
  [`dim()`](https://rdrr.io/r/base/dim.html),
  [`length()`](https://rdrr.io/r/base/length.html) and `[` methods.
- [`aggregate_NICHESObject()`](https://raredonlab.github.io/NICHESv2/reference/aggregate_NICHESObject.md)
  summarizes each cell’s neighborhood by LR pair and, in spatial mode,
  by neighboring cell type.
- [`load_LRM_database()`](https://raredonlab.github.io/NICHESv2/reference/load_LRM_database.md)
  provides ConnectomeDB2025 (default) and FANTOM5, both bundled for
  human, mouse, rat and pig, and OmniPath for human, mouse and rat.
- [`extract_NICHESInputs_Seurat()`](https://raredonlab.github.io/NICHESv2/reference/extract_NICHESInputs_Seurat.md)
  and
  [`extract_NICHESInputs_AnnData()`](https://raredonlab.github.io/NICHESv2/reference/extract_NICHESInputs_AnnData.md)
  read Seurat V4 and V5 objects and `.h5ad` files. The AnnData reader
  doesn’t need Python.
- [`export_to_TissuePlex()`](https://raredonlab.github.io/NICHESv2/reference/export_to_TissuePlex.md)
  writes edges to a Parquet file for interactive viewing.
- [`download_NICHESv2_data()`](https://raredonlab.github.io/NICHESv2/reference/download_NICHESv2_data.md)
  downloads the two public datasets used in the website tutorials: a 10x
  Visium mouse brain section and interferon-beta stimulated PBMCs.

### Changes for users of the 0.0.9000 development version

- The license is now the GNU AGPL v3 instead of MIT. A commercial
  license is also available from the Raredon Lab.
- `[` now keeps each retained cell’s neighborhood entries and
  aggregation rows unchanged, so they still describe the full tissue.
  Previously, neighborhood slots kept the full object’s entries,
  aggregations were dropped, and aggregating a subset could give wrong
  means.
- Recomputing aggregations on a subset now stops with an error.
  Aggregate the full object before subsetting, or run
  `add_neighborhoods(overwrite = TRUE)` first to analyze the subset as a
  separate tissue.
- `[` now keeps the `edge.filter.mode` and `remove.autocrine`
  attributes.
- `add_neighborhoods(overwrite = TRUE)` now clears the stale
  `$aggregations$neighborhood.edge.agg` table.
- Loading the package no longer changes `data.table`’s thread setting,
  and socket workers now stay at one `data.table` thread.
- Functions with a `seed` argument now restore the caller’s random
  number state.
- Mode `"sampled"` is now `"nonspatial"`.
- `Aggregate_CellToCell_Edges()` has been removed. Use
  [`aggregate_NICHESObject()`](https://raredonlab.github.io/NICHESv2/reference/aggregate_NICHESObject.md).
- `compute_PseudoBulk_Connectivity()` and `create_Edgelist_PseudoBulk()`
  have been removed. To summarize by cell type, sum `$edge.data` using
  the cell types in `$edge.meta`.
- [`add_cell_meta()`](https://raredonlab.github.io/NICHESv2/reference/add_cell_meta.md)
  now replaces `$cell.meta` in full. Use
  [`add_meta_column()`](https://raredonlab.github.io/NICHESv2/reference/add_meta_column.md)
  to change one column.
- The default database in
  [`load_LRM_database()`](https://raredonlab.github.io/NICHESv2/reference/load_LRM_database.md)
  is now ConnectomeDB2025.
- In
  [`add_neighborhood_edge_agg()`](https://raredonlab.github.io/NICHESv2/reference/add_neighborhood_edge_agg.md),
  the `autocrine` category is now split into `self.autocrine` and
  `other.autocrine`.
- `load_LRM_database("omnipath")` now falls back to the OmniPath web
  service when OmnipathR can’t reach its ancillary services.
