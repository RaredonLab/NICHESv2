# NICHESv2 0.0.9000

Development version distributed to lab members and collaborators for
pre-release testing. The API is stable but not yet frozen; a formal 0.1.0
release will accompany the methods paper.

---

## New features

### NICHESObject -- 15-slot S3 class

- Long-format sparse `data.table` storage for cell-level (`$cell.data`) and
  edge-level (`$edge.data`) scores; memory footprint scales with dataset
  density rather than full matrix dimensions.
- `score.norm` is proportional normalisation within cell (cell slot) or
  within edge (edge slot).
- `$aggregations` slot holds named aggregation tables computed on demand;
  starts as an empty list and is populated by `aggregate_NICHESObject()` or
  the individual `add_neighborhood_edge_agg()` and
  `add_neighborhood_composition()` functions.
- S3 methods: `print()`, `summary()`, `dim()`, `length()`, `[` (edge
  subsetting).
- Object attributes track `mode` (`"spatial"` or `"nonspatial"`) and
  `sample.col` (name of the sample column, if any).
- Neighborhood slots (`$neighborhood.cell.list`, `$neighborhood.edge.list`,
  `$neighborhood.meta`) are always populated by `create_NICHESObject()` with
  fixed parameters (`edge.filter.mode = "both_in_neighborhood"`,
  `remove.autocrine = FALSE`). Call `add_neighborhoods(overwrite = TRUE)` to
  experiment with other combinations after construction.

### LR database loaders -- `load_LRM_database()`

- ConnectomeDB2025 is now the default database, bundled as pre-filtered
  `.rda` objects for human, mouse, rat, and pig; no internet connection
  required. Call `load_LRM_database()` with no arguments to load human
  ConnectomeDB2025 pairs.
- FANTOM5 bundled as pre-filtered `.rda` objects for human, mouse, rat, and
  pig; no internet connection required.
- OmniPath loaded via `import_ligrecextra_interactions()` for human, mouse,
  and rat (pig not supported by OmniPath).
- User-supplied data frames accepted directly with `ligand.col` and
  `receptor.col` arguments for non-standard column names.

### Edge list constructors

- `create_Edgelist_Spatial()` -- k-NN or radius-based spatial edges including
  autocrine self-loops; parallel via socket clusters.
- `create_Edgelist_PseudoBulk()` -- deterministic all-type-pairs edge list;
  no randomness.
- `create_Edgelist_SamplingByType()` -- type-stratified random sampling with
  reproducible `seed` and `max.cells` cap.

### Scoring

- `compute_CellToCell()` -- universal per-edge LR scorer; `"product"` and
  `"minimum"` methods; multi-gene complexes handled via underscore delimiting
  (`gene1_gene2`). Returns a plain named list (`$edges`, `$lr.meta`) rather
  than a NICHESObject, to allow flexible downstream assembly.
- `compute_PseudoBulk_Connectivity()` -- type-pair aggregated scores returned
  as a plain `data.table`; does not produce a NICHESObject.

### Object construction -- `create_NICHESObject()`

- Single dispatcher routes to per-sample or single-sample internal branches
  depending on `sample.col`; the outer sample loop is always sequential while
  inner scoring uses the full `n.cores` allocation.
- Hard barcode-uniqueness error across samples; tip: prefix barcodes with
  `paste0(sample.id, "_", barcode)` before calling.
- Neighborhood slots are always populated at construction; no separate
  `add_neighborhoods()` call required for the default configuration.
- Two bundled datasets for immediate exploration and reproducible examples:
  `data(NICHESv2_inputs)` (raw construction inputs) and
  `data(NICHESv2_example)` (pre-built NICHESObject with neighborhood slots
  populated).

### Object extensions

- `add_neighborhoods()` / `extract_neighborhoods()` -- re-compute neighborhood
  slots with three `edge.filter.mode` options: `"both_in_neighborhood"`,
  `"any_in_neighborhood"`, `"central_cell_only"`. Use `overwrite = TRUE` to
  replace the slots populated at construction.
- `add_segmentation()` -- attach an `sf` polygon object to `$segmentation`;
  hard error on any missing barcode.
- `add_cell_meta()` -- fully replaces `$cell.meta` with the supplied table
  and rebuilds `$edge.meta` by default. The replacement is complete: columns
  present in the old `$cell.meta` but absent from the new table are dropped.
  The user is responsible for supplying all columns they wish to retain,
  including spatial coordinates. Hard error if any cell in `$cell.list` is
  missing from the new table.
- `add_meta_column()` -- appends or overwrites a single column in `$cell.meta`
  from a named vector (names = cell barcodes) without disturbing other columns;
  emits a message if the column already exists; rebuilds `$edge.meta` by
  default.
- `add_cell_rawcounts()` -- adds a `count` column to `$cell.data` from a
  user-supplied raw count matrix, for provenance. Non-negative integer checks
  are performed with warnings; pipeline continues regardless.
- `merge_NICHESObjects()` -- combine two or more NICHESObjects; hard errors
  on mode mismatch or duplicate barcodes; union strategy for gene and LRM
  lists.

### Data ingest

- `extract_NICHESInputs_Seurat()` -- compatible with Seurat V4 and V5;
  version detected by class inspection (`inherits(assay, "Assay5")`), not
  package version string.
- `extract_NICHESInputs_AnnData()` -- reads `.h5ad` via `rhdf5` only; no
  Python, conda, basilisk, or zellkonverter required; requires R >= 4.2.

### Aggregation

- `aggregate_NICHESObject()` -- convenience wrapper that calls
  `add_neighborhood_edge_agg()` and optionally `add_neighborhood_composition()`
  in one step; always fully recomputes and replaces both slots (Seurat
  convention).
- `add_neighborhood_edge_agg()` -- populates
  `$aggregations$neighborhood.edge.agg`: a sparse keyed `data.table` of
  per-(cell, LRM) directional scores classified into five categories (in, out,
  cross, self.autocrine, other.autocrine) with sum, mean, and
  conditional-mean statistics. Output table has 27 columns.
- `add_neighborhood_composition()` -- populates
  `$aggregations$neighborhood.composition`: a zero-filled keyed `data.table`
  of per-cell neighbor cell-type counts and proportions; spatial mode only.

---

## Breaking changes

- **`"sampled"` mode renamed to `"nonspatial"`** across all functions,
  attributes, and documentation. Any code that checks
  `attr(obj, "mode") == "sampled"` must be updated to `"nonspatial"`.
- **`Aggregate_CellToCell_Edges()` removed.** This legacy function has been
  deleted. Use `aggregate_NICHESObject()`, `add_neighborhood_edge_agg()`, and
  `add_neighborhood_composition()` for all aggregation workflows.
- **`add_cell_meta()` now fully replaces `$cell.meta`** rather than joining
  new columns onto the existing table. Code that relied on column preservation
  must be updated to supply a complete metadata table.
- **`load_LRM_database()` default database changed** from FANTOM5 to
  ConnectomeDB2025. Existing code that calls `load_LRM_database()` with no
  arguments will now load ConnectomeDB2025 human pairs instead of FANTOM5
  human pairs. Pass `db = "fantom5"` explicitly to restore the previous
  behaviour.
- **`add_neighborhood_edge_agg()`: `autocrine` split into `self.autocrine`
  and `other.autocrine`.** The single `autocrine` category has been replaced
  by two sub-categories. `self.autocrine` captures only the focal cell's own
  self-loop (count.self.autocrine is always 0 or 1 per (cell, LRM));
  `other.autocrine` captures self-loop edges of non-focal neighbor cells in
  the neighborhood. The output table now has 27 columns instead of 23.
  The following columns have been removed and must be updated in any
  downstream code: `count.autocrine`, `score.autocrine.sum`,
  `score.autocrine.mean`, `score.autocrine.cmean`.

---

## Known limitations

- **`[.NICHESObject` does not propagate attributes.** Subsetting a
  NICHESObject with `[` does not copy `mode` or `sample.col` attributes to
  the result. Fix is known and scheduled for 0.1.1.
- **`Aggregate_CellToCell_Edges()` returning dense matrices** is no longer
  a limitation -- the function has been removed and replaced with a sparse
  design in `add_neighborhood_edge_agg()`.

---

## Acknowledgements

Development supported by the Raredon Lab, Yale University School of Medicine.
NICHESv2 is the successor to the original NICHES package (Raredon et al.).
Methods paper in preparation.
