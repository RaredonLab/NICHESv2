# NICHESv2 NEWS

## NICHESv2 0.1.0 (2026-03)

### New features

* **`NICHESObject`** — new 13-slot S3 class replacing the intermediate
  `EdgeLRResult` as the primary output structure. Slots:
  `$cell.list`, `$cell.data`, `$cell.meta`, `$edge.list`, `$edge.data`,
  `$edge.meta`, `$neighborhood.*` (×3, reserved), `$gene.list`,
  `$gene.meta`, `$LRM.list`, `$LRM.meta`.

* **`create_NICHESObject()`** — single top-level constructor. Accepts a
  gene × cell count matrix plus a `meta.data` data.frame (must contain
  `x` and `y` columns for spatial coordinates) and orchestrates the full
  pipeline internally.

* **`load_LRM_database()`** — built-in loaders for `"fantom5"` and
  `"omnipath"` databases (via `OmnipathR`), plus a passthrough for
  user-supplied data.frames.

* **`add_cell_meta()`** — add or overwrite `$cell.meta` columns
  post-construction; rebuilds `$edge.meta` automatically.

* **`is_NICHESObject()`** — public identity check.

* **S3 methods** for `NICHESObject`: `print`, `summary`, `dim`, `length`,
  `[` (cell subsetting).

### Breaking changes vs. NICHES v1

* Output class renamed `NICHESObject` (was `EdgeLRResult`).
* LRM IDs now use `|` as separator (was `.`): `"MIF|CD74"` not `"MIF.CD74"`.
* Edge IDs use `|`: `"CellA|CellB"`.
* Signal value column renamed `score` (was `edge`).
* Normalised score column added: `score.norm`.
* `Aggregate_CellToCell_Edges()` now accepts `NICHESObject` (not
  `EdgeLRResult`); aggregated columns renamed `score.sum`, `score.mean`,
  `score.max`, `score.stat`.
* `CellToCell_InSitu()` and `CellToCell_InSitu_Single_LRM()` now accept
  `NICHESObject`.

### Complex ligand/receptor support

* Ligands and receptors may be multi-gene complexes using underscore
  notation (e.g. `"COL1A1_COL1A2"`).  All component expressions are
  combined via the chosen `method` (`"product"` or `"minimum"`) across
  the full set of ligand and receptor components.

### Internal improvements

* `utils.R` centralises shared helpers: `.parse_complex()`,
  `.compute_lrm_score()`, `.genes_available()`, `.prop_norm()`,
  `.categorical_cols()`.
* `create_Edgelist_Spatial()` is now called internally by
  `create_NICHESObject()` but remains independently exported for advanced
  workflows.
* `compute_CellToCell_Spatial()` retained as a standalone export for
  users who provide custom edge lists (e.g. contact-based adjacency).
