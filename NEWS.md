# NICHESv2 0.0.9000

Development version distributed to lab members and collaborators for
pre-release testing. The API is stable but not yet frozen; a formal 0.1.0
release will accompany the methods paper.

---

## New features

### NICHESObject — 14-slot S3 class

- Long-format sparse `data.table` storage for cell-level (`$cell.data`) and
  edge-level (`$edge.data`) scores; memory footprint scales with dataset
  density rather than full matrix dimensions.
- Raw counts always stored in `$cell.data` even when expression smoothing is
  applied; `score.norm` is proportional normalisation within cell (cell slot)
  or within edge (edge slot).
- S3 methods: `print()`, `summary()`, `dim()`, `length()`, `[` (edge
  subsetting).
- Object attributes track `mode` (`"spatial"` or `"sampled"`) and
  `sample.col` (name of the sample column, if any).
- Neighborhood slots (`$neighborhood.cell.list`, `$neighborhood.edge.list`,
  `$neighborhood.meta`) are `NULL` until `add_neighborhoods()` is called;
  they store IDs only and are joined to score slots for analysis.

### LR database loaders — `load_LRM_database()`

- FANTOM5 bundled as pre-filtered `.rda` objects for human, mouse, rat, and
  pig; no internet connection required.
- OmniPath loaded via `import_ligrecextra_interactions()` for human, mouse,
  and rat.
- User-supplied data frames accepted directly with `ligand.col` and
  `receptor.col` arguments for non-standard column names.

### Edge list constructors

- `create_Edgelist_Spatial()` — k-NN or radius-based spatial edges including
  autocrine self-loops; parallel via socket clusters.
- `create_Edgelist_PseudoBulk()` — deterministic all-type-pairs edge list;
  no randomness.
- `create_Edgelist_SamplingByType()` — type-stratified random sampling with
  reproducible `seed` and `max.cells` cap.

### Scoring

- `compute_CellToCell()` — universal per-edge LR scorer; `"product"` and
  `"minimum"` methods; multi-gene complexes handled via underscore delimiting
  (`gene1_gene2`). Returns a plain named list (`$edges`, `$lr.meta`) rather
  than a NICHESObject, to allow flexible downstream assembly.
- `compute_PseudoBulk_Connectivity()` — type-pair aggregated scores returned
  as a plain `data.table`; does not produce a NICHESObject.

### Object construction — `create_NICHESObject()`

- Single dispatcher routes to per-sample or single-sample internal branches
  depending on `sample.col`; the outer sample loop is always sequential while
  inner scoring uses the full `n.cores` allocation.
- Hard barcode-uniqueness error across samples; tip: prefix barcodes with
  `paste0(sample.id, "_", barcode)` before calling.
- Optional PCA-based kNN expression smoothing (`smooth = TRUE`; see Known
  Limitations below).
- Two bundled datasets for immediate exploration and reproducible examples:
  `data(NICHESv2_inputs)` (raw construction inputs) and
  `data(NICHESv2_example)` (pre-built NICHESObject).

### Object extensions

- `add_neighborhoods()` / `extract_neighborhoods()` — populate neighborhood
  slots with three `edge.filter.mode` options: `"both_in_neighborhood"`,
  `"any_in_neighborhood"`, `"central_cell_only"`.
- `add_segmentation()` — attach an `sf` polygon object to `$segmentation`;
  hard error on any missing barcode.
- `add_cell_meta()` — join new metadata columns onto `$cell.meta` and
  optionally rebuild `$edge.meta` via `rebuild.edge.meta = TRUE`.
- `merge_NICHESObjects()` — combine two or more NICHESObjects; hard errors
  on mode mismatch or duplicate barcodes; union strategy for gene and LRM
  lists.

### Data ingest

- `extract_NICHESInputs_Seurat()` — compatible with Seurat V4 and V5;
  version detected by class inspection (`inherits(assay, "Assay5")`), not
  package version string.
- `extract_NICHESInputs_AnnData()` — reads `.h5ad` via `rhdf5` only; no
  Python, conda, basilisk, or zellkonverter required; requires R ≥ 4.2.

### Visualisation

- `CellToCell_InSitu()` — batch PNG export for a list of LR pairs via socket
  cluster parallelism; background is pre-computed once before workers start.
- `CellToCell_InSitu_Single_LRM()` — single-pair in-situ plot; accepts a
  `precomputed.bg` argument to skip redundant background computation.

### Aggregation

- `Aggregate_CellToCell_Edges()` — returns full cell × LRM grids (zero-filled)
  from both sending and receiving perspectives as dense matrices in a named
  list (`$sending`, `$receiving`).

---

## Known limitations

- **`smooth = TRUE` is untested end-to-end.** The smoothing branch in
  `create_NICHESObject()` has no integration test coverage. The underlying
  `smooth_expression()` function is unit-tested and correct, but the call site
  in the dispatcher should be verified before use on real data.
- **`[.NICHESObject` does not propagate attributes.** Subsetting a
  NICHESObject with `[` does not copy `mode` or `sample.col` attributes to
  the result. Fix is known and scheduled for 0.1.1.
- **Neighborhood integration into `create_NICHESObject()` is pending.** Call
  `add_neighborhoods()` explicitly after object construction in the meantime.
- **`Aggregate_CellToCell_Edges()` returns dense matrices.** For very large
  datasets this may be memory-intensive; a sparse redesign is under
  consideration for a future release.

---

## Acknowledgements

Development supported by the Raredon Lab, Yale University School of Medicine.
NICHESv2 is the successor to the original NICHES package (Raredon et al.).
Methods paper in preparation.
