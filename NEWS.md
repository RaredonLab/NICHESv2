# NICHESv2 1.0.0

First public release.

## Main features

- `create_NICHESObject()` scores ligand-receptor (LR) signaling between pairs of
  cells, in `"spatial"` mode (k-nearest neighbors or fixed radius) or
  `"nonspatial"` mode (type-stratified random sampling).
- The 15-slot `NICHESObject` stores non-zero scores in sparse long-format
  `data.table`s, with `print()`, `summary()`, `dim()`, `length()` and `[`
  methods.
- `aggregate_NICHESObject()` summarizes each cell's neighborhood by LR pair and,
  in spatial mode, by neighboring cell type.
- `compute_PseudoBulk_Connectivity()` scores each pair of cell types.
- `load_LRM_database()` provides ConnectomeDB2025 (default) and FANTOM5, both
  bundled for human, mouse, rat and pig, and OmniPath for human, mouse and rat.
- `extract_NICHESInputs_Seurat()` and `extract_NICHESInputs_AnnData()` read
  Seurat V4 and V5 objects and `.h5ad` files. The AnnData reader doesn't need
  Python.
- `export_to_TissuePlex()` writes edges to a Parquet file for interactive
  viewing.

## Changes for users of the 0.0.9000 development version

- `[` now keeps each retained cell's neighborhood entries and aggregation rows
  unchanged, so they still describe the full tissue. Previously, neighborhood
  slots kept the full object's entries, aggregations were dropped, and aggregating
  a subset could give wrong means.
- Recomputing aggregations on a subset now stops with an error. Aggregate the
  full object before subsetting, or run `add_neighborhoods(overwrite = TRUE)`
  first to analyze the subset as a separate tissue.
- `[` now keeps the `edge.filter.mode` and `remove.autocrine` attributes.
- `add_neighborhoods(overwrite = TRUE)` now clears the stale
  `$aggregations$neighborhood.edge.agg` table.
- Loading the package no longer changes `data.table`'s thread setting, and socket
  workers now stay at one `data.table` thread.
- Functions with a `seed` argument now restore the caller's random number state.
- Mode `"sampled"` is now `"nonspatial"`.
- `Aggregate_CellToCell_Edges()` has been removed. Use
  `aggregate_NICHESObject()`.
- `add_cell_meta()` now replaces `$cell.meta` in full. Use `add_meta_column()` to
  change one column.
- The default database in `load_LRM_database()` is now ConnectomeDB2025.
- In `add_neighborhood_edge_agg()`, the `autocrine` category is now split into
  `self.autocrine` and `other.autocrine`.
- `load_LRM_database("omnipath")` now falls back to the OmniPath web service when
  OmnipathR can't reach its ancillary services.
