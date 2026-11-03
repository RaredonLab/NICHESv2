# NICHESv2 -- Package Design & Status
# Updated: 2026-11-02

================================================================================
COLLABORATION & DELIVERY CONVENTIONS (read at session start)
================================================================================

File delivery modes:
  Mode 1 -- targeted edits to existing files: provide exact Old/New text
             blocks the user applies manually via str_replace or copy-paste.
             Always show the minimal change needed, not the whole file.
  Mode 2 -- new files: deliver the complete file content as a downloadable
             file via the present_files tool. User copies to the appropriate
             directory (R/, tests/testthat/, etc.).

Mode 1 formatting rules (no exceptions):
  - Every edit is presented as a separate numbered block: "Edit N of M"
  - Old and New code blocks are always SEPARATE -- never combined or diffed.
    The user copy-pastes the Old block to find the location, then replaces
    with the New block. Mixing them makes copy-paste error-prone.
  - Edits within a file are listed in TOP-TO-BOTTOM text order, not ranked
    by importance. If an edit is optional, label it "optional" with a brief
    one-line explanation immediately after the edit heading.
  - No grep verification counts are needed for new files (Mode 2). For
    Mode 1 edits, only add grep checks when the old text could plausibly
    match in multiple places -- if the old block is unique by inspection,
    skip the grep.

Before writing any code:
  - Ask for uploads of all files that will be modified; read them fully first.
  - Never guess at file content -- always read before editing.
  - State the complete plan (which files change, logic, decisions to confirm)
    and wait for explicit confirmation before writing any code or edits.
  - Decide which files are needed and ask for them all at once. Never ask
    what file do you want to upload -- determine what is needed and request
    it by name.

Verification greps:
  - Before stating any expected grep count, compute it explicitly: count
    occurrences in the new text being added PLUS occurrences already in the
    existing file. Never guess from memory -- this causes repeated errors.
  - Use python3 -c or bash to verify counts programmatically when uncertain.

CRAN compliance on every delivery:
  - No non-ASCII characters anywhere in .R files (comments included).
    Use -- not em/en dash, ... not ellipsis, straight quotes only.
    Verify with: tools::showNonASCII(readLines("R/file.R"))
  - inherits(x, "NICHESObject") always -- never .is_NICHESObject()
  - All new NSE column names declared in globalVariables() in utils.R,
    appended at the END of the existing block with a per-function comment.
  - @seealso links only reference functions that currently exist.
  - One exported function per file, file named after the function.

End-of-session .md update:
  - When the user asks to update the .md, produce a complete replacement
    file (Mode 2), not incremental edits. Carry forward ALL existing
    sections; never silently drop any section including this one.


================================================================================
CONVENTIONS (no exceptions)
================================================================================

Functions       : underscore  e.g. create_NICHESObject(), CellToCell_InSitu()
Variables/args  : dot         e.g. n.cores, bg.color, edge.list
Edge ID         : "SendingCell|ReceivingCell"
LRM ID          : "ligand|receptor"
Gene complexes  : "gene1_gene2"
Score columns   : score (raw)  |  score.norm (normalised, proportional)

PARALLELISM (CRITICAL -- no exceptions):
  socket clusters ONLY: parallel::makeCluster + parLapply + on.exit(stopCluster)
  NEVER mclapply (Unix-only, crashes on macOS with ggplot2).
  Always clusterEvalQ: data.table::setDTthreads(1)
  The package never changes data.table's thread setting on load: an earlier
  .onLoad did, which also reset every worker that loaded the namespace.
  Exception: add_neighborhood_edge_agg() uses sparse matrix multiplication
  (Matrix package) instead of socket clusters. n.cores retained for API
  compatibility only -- not used. New matrix-based functions follow same pattern.

NAMESPACE: Never hand-edit. Use @import data.table / @import ggplot2 in
  NICHESv2-package.R. Run devtools::document() TWICE after changes.

ONE EXPORTED FUNCTION PER FILE.


================================================================================
NICHESObject -- 15 SLOTS
================================================================================

$cell.list     chr vector   -- barcodes (colnames of count.mtx), all modes
$cell.data     data.table   -- SPARSE: cell|gene|score|score.norm[|count]; key(cell,gene)
                               score.norm = gene/total within cell; count column
                               added optionally via add_cell_rawcounts(); new rows
                               appended for (cell,gene) where count>0 and score=0
$cell.meta     data.table   -- cell|x|y|<meta cols>; key: cell
$edge.list     data.table   -- edge|SendingCell|ReceivingCell; key: edge
$edge.data     data.table   -- SPARSE: edge|LRM|score|score.norm; key(edge,LRM)
                               score.norm = proportional within edge
$edge.meta     data.table   -- edge|SendingCell|ReceivingCell|Sending.*|Receiving.*|Joint.*
$neighborhood.cell.list     -- named list: cell -> character vector of neighbor barcodes
$neighborhood.edge.list     -- named list: cell -> character vector of edge IDs
$neighborhood.meta          -- data.table: cell|n.neighbors|n.edges; key: cell
$gene.list     chr vector   -- rownames(count.mtx)
$gene.meta     data.table   -- keyed: gene
$LRM.list      chr vector   -- valid "ligand|receptor" strings
$LRM.meta      data.table   -- LRM|ligand|receptor|<extras>; keyed: LRM
$segmentation  sf object    -- NULL until add_segmentation(); sf in Suggests
$aggregations  list         -- named list of aggregation data.tables:
                               $neighborhood.edge.agg   data.table, keyed (cell,LRM)
                               $neighborhood.composition data.table, keyed cell
                               Starts as list() until aggregate_NICHESObject() or
                               individual add_*() functions are called.

Neighborhood slots are always populated by create_NICHESObject() using fixed
parameters: edge.filter.mode = "both_in_neighborhood", remove.autocrine = FALSE.
They store IDs only -- join to $edge.data/$edge.meta for scores/coordinates.
When add_neighborhoods(overwrite=TRUE) is called, stale $neighborhood.edge.agg
and $neighborhood.composition entries are cleared automatically with a message.

Object attributes:
  attr(obj,"mode")             = "spatial"|"nonspatial"
  attr(obj,"sample.col")       = col name string (if split by sample)
  attr(obj,"edge.filter.mode") = set by add_neighborhoods(); always
                                 "both_in_neighborhood" from create_NICHESObject()
  attr(obj,"remove.autocrine") = set by add_neighborhoods(); always FALSE
                                 from create_NICHESObject()


================================================================================
EXPORTED FUNCTIONS
================================================================================

Each lives in its own file under R/.

-- LR DATABASE --
load_LRM_database(db, species, ligand.col, receptor.col, verbose)
  db: "connectomedb2025" (default) | "fantom5" | "omnipath" | user data.frame
  ConnectomeDB2025: human/mouse/rat/pig from bundled .rda; no internet required.
    Default: load_LRM_database() with no arguments loads human ConnectomeDB2025.
  FANTOM5: human/mouse/rat/pig from bundled .rda; no internet, no OmnipathR
  OmniPath: human/mouse/rat via import_ligrecextra_interactions(organism=<code>)
    (pig not supported by OmniPath)
  -> data.frame with "ligand", "receptor"

-- EDGE LIST CONSTRUCTORS --
create_Edgelist_Spatial(coord, k, rad, n.cores)
  -> data.table(from|to|weight|x.from|y.from|x.to|y.to) + autocrine self-loops

create_Edgelist_SamplingByType(meta.data, cell.type.col, seed=42L, max.cells)
  -> data.table(from|to) with barcodes; attrs: mode="nonspatial"

-- SCORING --
compute_CellToCell(count.mtx, edge.list, LRM.db, ligand.col, receptor.col,
                   method, n.cores)
  -> plain named list: $edges | $lr.meta
  $edges:   data.table -- from|to|ligand|receptor|LRM|score (+ x.from|y.from|
            x.to|y.to for spatial edgelists); keyed (from, to, LRM)
  $lr.meta: data.table -- ligand|receptor|LRM + extras from LRM.db;
            keyed (ligand, receptor)

-- OBJECT CONSTRUCTION --
create_NICHESObject(count.mtx, meta.data, LRM.db,
                    mode="spatial", species="human", sample.col=NULL,
                    cell.type.col=NULL, ligand.col, receptor.col,
                    k, rad, seed=42L, max.cells,
                    method, normalize.method, gene.meta, n.cores, verbose)
  Dispatcher only. Routes to .create_NICHESObject_single() or
  .create_NICHESObject_by_sample() depending on sample.col. After object
  construction, always calls add_neighborhoods() with fixed parameters:
    edge.filter.mode = "both_in_neighborhood", remove.autocrine = FALSE.
  Neighborhood slots are ALWAYS populated after create_NICHESObject().
  edge.filter.mode and remove.autocrine are NOT user-facing parameters here.
  Use add_neighborhoods() directly to experiment with other combinations.
  -> NICHESObject
  FILE: R/create_NICHESObject.R (dispatcher only -- single function in this file)

-- OBJECT EXTENSIONS --
add_neighborhoods(niches.obj, edge.filter.mode="both_in_neighborhood",
                  remove.autocrine=FALSE, overwrite=FALSE, n.cores, verbose)
  Calls extract_neighborhoods() on niches.obj$edge.list, populates 3 slots.
  overwrite=FALSE hard-errors if slots already populated.
  overwrite=TRUE clears stale $neighborhood.edge.agg and
  $neighborhood.composition before recomputing, then sets both
  attr(obj,"edge.filter.mode") and attr(obj,"remove.autocrine").
  FILE: R/add_neighborhoods.R

extract_neighborhoods(edge.list, edge.filter.mode, remove.autocrine,
                      n.cores, verbose)
  edge.filter.mode: "both_in_neighborhood" (default) |
                    "any_in_neighborhood" | "central_cell_only"
  -> list($neighborhood.cell.list, $neighborhood.edge.list, $neighborhood.meta)
  FILE: R/extract_neighborhoods.R

add_segmentation(niches.obj, seg.sf)
  Hard error if any cell missing polygon.

add_cell_meta(niches.obj, new.meta, rebuild.edge.meta=TRUE, verbose=TRUE)
  FULL REPLACEMENT of $cell.meta -- not a merge. $cell.meta becomes exactly
  new.meta. Columns in the old $cell.meta absent from new.meta are dropped.
  User is responsible for supplying all columns they wish to retain, including
  spatial coordinates (x, y). Hard error if any cell in $cell.list is missing
  from new.meta. Extra rows in new.meta are silently dropped with a verbose
  message. Accepts a "cell" column or rownames as cell barcodes. Rebuilds
  $edge.meta by default. Emits a staleness message if
  $aggregations$neighborhood.composition is populated.
  FILE: R/add_cell_meta.R

add_meta_column(niches.obj, col.name, col.values,
                rebuild.edge.meta=TRUE, verbose=TRUE)
  Appends or overwrites a SINGLE column in $cell.meta from a named vector
  (names = cell barcodes). All other columns in $cell.meta are preserved.
  Hard error if any cell in $cell.list is missing from names(col.values).
  Emits a message if col.name already exists (then overwrites).
  Rebuilds $edge.meta by default. Emits a staleness message if
  $aggregations$neighborhood.composition is populated.
  FILE: R/add_meta_column.R

add_cell_rawcounts(niches.obj, count.mtx, verbose=TRUE)
  Adds a "count" column to $cell.data from a user-supplied gene x cell raw
  count matrix (base matrix or dgCMatrix). Provenance only -- not used by
  any other NICHESv2 function.
  Row expansion (Option B): new rows appended for (cell,gene) pairs where
  count>0 but score=0; new rows carry score=0, score.norm=0.
  Hard error if any cell in $cell.list is missing from colnames(count.mtx).
  Genes in $gene.list absent from rownames(count.mtx) get count=NA with
  a verbose message. Non-negative integer checks on stored values emit
  warnings (independent); pipeline continues regardless.
  Emits message if "count" column already exists (then overwrites).
  FILE: R/add_cell_rawcounts.R

merge_NICHESObjects(..., verbose=TRUE)
  Accepts: positional args | pre-built list | explicit list() call.
  Hard errors: empty list | non-NICHESObject | mode mismatch | dup barcodes.
  Slot strategy: rbindlist data slots; union gene/LRM lists; sf::rbind segs.
  Self-contained -- no internal helper dependency.
  FILE: R/merge_NICHESObjects.R

-- DATA INGEST --
extract_NICHESInputs_Seurat(seurat.obj, assay="RNA", layer="counts", verbose)
  V4/V5 detection via inherits(assay.obj,"Assay5"), NOT package version.
  -> list(count.mtx, meta.data); user calls create_NICHESObject() next.
  Suggests: SeuratObject

extract_NICHESInputs_AnnData(file, layer="counts", obs.cols=NULL, verbose)
  rhdf5 only. No Python, conda, basilisk, or zellkonverter. R >= 4.2.
  -> list(count.mtx, meta.data)
  Suggests: rhdf5

-- VISUALIZATION --
CellToCell_InSitu(niches.obj, LRM.list, folder, ..., n.cores)
  Batch PNG; socket cluster; pre-computes all.edges before workers.

CellToCell_InSitu_Single_LRM(niches.obj, LRM, ..., precomputed.bg=NULL)

-- AGGREGATION (current design) --
aggregate_NICHESObject(niches.obj, cell.type.col=NULL, n.cores=1L,
                       chunk.size=50L, verbose=TRUE)
  Convenience wrapper. Calls add_neighborhood_edge_agg() then optionally
  add_neighborhood_composition(). Both steps always fully recompute and
  replace their slots -- no incremental update mode (Seurat convention).
  chunk.size is passed through to add_neighborhood_edge_agg().
  FILE: R/aggregate_NICHESObject.R

add_neighborhood_edge_agg(niches.obj, n.cores=1L, chunk.size=50L, verbose=TRUE)
  Five-category unified edge aggregation. For each (cell, LRM) pair,
  classifies every neighborhood edge into one of:
    out             -- focal cell is SendingCell (non-autocrine)
    in              -- focal cell is ReceivingCell (non-autocrine)
    cross           -- focal cell is neither endpoint (non-autocrine)
    self.autocrine  -- SendingCell == ReceivingCell == focal cell
    other.autocrine -- SendingCell == ReceivingCell != focal cell
  Implementation: sparse matrix multiplication (Matrix package). Builds
  five binary A matrices (cells x edges) once, then multiplies by the
  score matrix B (edges x LRMs) in batches of chunk.size LRMs. Bounds
  peak memory while reducing iterations from n.lrms to
  ceil(n.lrms / chunk.size). n.cores retained for API compatibility only.
  Always computes all three score methods (sum, mean, cmean) for all five
  categories plus totals and score.inout.ratio. Always replaces slot.
  Requires neighborhood slots to be populated.
  -> $aggregations$neighborhood.edge.agg keyed (cell,LRM)
  Column structure (27 columns):
    count.in | count.out | count.cross |
    count.self.autocrine | count.other.autocrine | count.total
    score.{cat}.sum   -- total score per category
    score.{cat}.mean  -- sum / n.edges in neighborhood (LRM-agnostic denom)
    score.{cat}.cmean -- sum / count.{cat}; NA when count == 0
    score.total.*     -- sum across all five categories
    score.inout.ratio -- score.in.sum / score.out.sum; NA when out == 0
  FILE: R/add_neighborhood_edge_agg.R

add_neighborhood_composition(niches.obj, cell.type.col, include.self=TRUE,
                              normalize=TRUE, verbose=TRUE)
  Source: $edge.list + $cell.meta only. No neighborhood slot dependency.
  Neighbor = any cell sharing an edge (either direction), deduplicated.
  Spatial mode only (hard error for nonspatial mode).
  Always overwrites -- no overwrite parameter.
  Column names: gsub("[^A-Za-z0-9.]",".",type) + .neighbor.count/.neighbor.prop
  Zero-filled against full dataset type universe.
  FILE: R/add_neighborhood_composition.R


================================================================================
INTERNAL HELPERS (utils.R)
================================================================================

Small utilities only:
  .parse_complex(gene.str)        "gene1_gene2" -> c("gene1","gene2")
  .genes_available(gene.str, set) all components in set?
  .prop_norm(x)                   x/sum(x); zeros if sum==0
  .categorical_cols(dt, exclude)  character or factor column names
  .new_NICHESObject(...)          15-slot constructor (aggregations = list())
  .coerce_meta(meta.data)         -> keyed data.table with "cell" col
  .is_NICHESObject(x)             listed here but does NOT exist in utils.R;
                                   always use inherits(x,"NICHESObject") directly
  %||%                            null-coalescing operator

Internal branch implementations of create_NICHESObject() (each in its own file):
  .create_NICHESObject_single()    R/create_NICHESObject_single.R
  .create_NICHESObject_by_sample() R/create_NICHESObject_by_sample.R

h5ad internal helpers (in R/extract_NICHESInputs_AnnData.R, unexported):
  .h5ad_attrs, .h5ad_encoding, .h5ad_obs_names, .h5ad_var_names,
  .h5ad_available_layers, .h5ad_read_sparse, .h5ad_read_matrix, .h5ad_read_obs

Internal helpers in dedicated files:
  .compute_LRM_Score(lig.str, rec.str, count.mtx,
                  sending.cells, receiv.cells, method="product")
  Single LR pair -> numeric vector of length n.edges
  method: "product" | "minimum"; supports gene complexes via "_"
  FILE: R/compute_LRM_Score.R

  .build_cell_data(count.mtx, normalize.method="prop")
  Sparse count matrix -> long-format cell|gene|score|score.norm data.table.
  FILE: R/build_cell_data.R

  .build_edge_meta(edge.list.dt, cell.meta)
  Joins cell metadata onto edge list. Adds Sending.*, Receiving.*, Joint.*
  columns (Joint only for character/factor cols).
  FILE: R/build_edge_meta.R

globalVariables declarations live in utils.R (NOT NICHESv2-package.R).
Organized in per-function comment blocks at the bottom of utils.R.
Current blocks: score/LRM/cell/edge core names, add_neighborhood_edge_agg
columns, add_neighborhood_composition columns, add_cell_rawcounts columns
(count, i.count, gene).


================================================================================
PACKAGE INFRASTRUCTURE
================================================================================

Imports: data.table(>=1.14), Matrix(>=1.5), methods, RANN(>=2.6), dbscan(>=1.1),
         ggplot2(>=3.4), parallel, stats, utils
Suggests: SeuratObject, rhdf5, OmnipathR, arrow, testthat(>=3.0), withr, sf,
          knitr, rmarkdown

Bundled data (data/):
  ncomms8866_orig.rda -- FANTOM5 original unfiltered table; R object name is
    ncomms8866 (file and object names differ). data(ncomms8866) fails with
    LazyData: false. Access via load_LRM_database("fantom5", ...) only.
  ncomms8866_human/mouse/rat/pig.rda -- FANTOM5 LR pairs, species-converted
    per Raredon et al. 2019. Loaded internally via utils::data().
  connectomedb2025_human/mouse/rat/pig.rda -- ConnectomeDB2025 LR pairs,
    produced by data-raw/build_connectomedb2025.R. Default database.
    Access via load_LRM_database("connectomedb2025", species=...) or
    load_LRM_database() (no args = human ConnectomeDB2025).
    Loaded internally via utils::data(). Direct data() access discouraged.

Tests (all passing):
  helper-testutils.R              -- keep (make_grid_coord utility)
  test-compute_LRM_Score.R        test-create_Edgelist_Spatial.R
  test-extract_neighborhoods.R    test-create_Edgelist_SamplingByType.R
  test-add_neighborhoods.R        test-compute_CellToCell.R
  test-merge_NICHESObjects.R
  test-Explore_Radius_Neighborhood.R  test-create_NICHESObject.R
  test-extract_NICHESInputs_Seurat.R  test-load_LRM_database.R
  test-extract_NICHESInputs_AnnData.R test-add_cell_meta.R
  test-add_segmentation.R             test-add_meta_column.R
  test-CellToCell_InSitu.R            test-NICHESObject_methods.R
  test-add_neighborhood_edge_agg.R    test-add_neighborhood_composition.R
  test-aggregate_NICHESObject.R       test-add_cell_rawcounts.R

Fixture strategy:
  - data(NICHESv2_example): add_cell_meta, add_meta_column,
    add_cell_rawcounts, add_segmentation, CellToCell_InSitu,
    NICHESObject_methods, add_neighborhood_edge_agg,
    add_neighborhood_composition, aggregate_NICHESObject
  - Inline synthesis: all other files
  - helper-testdata.R and inst/testdata/ DELETED
  - NICHESv2_example now always has neighborhood slots populated
  - No aggregation slots pre-populated in fixture

Skips: none when all Suggests packages are installed.

R CMD check status: 0 errors | 0 warnings | 1 note
  (N1: timestamp/network issue on dev machine -- unfixable in code, ignore)


================================================================================
CHECKLIST -- EVERY NEW EXPORTED FUNCTION (CRAN compliance)
================================================================================

Before delivering any new exported function, verify ALL of the following:
  [ ] Input guard: inherits(x, "NICHESObject") -- NOT .is_NICHESObject()
  [ ] All dot-prefixed helpers called exist in utils.R (grep to confirm)
  [ ] All NSE column names in := or .() declared in globalVariables() in utils.R
  [ ] New globalVariables block appended at END of utils.R globalVariables list
  [ ] Socket cluster: clusterExport uses envir = environment() for local vars
  [ ] No mclapply -- socket clusters only (parallel::makeCluster + parLapply)
  [ ] devtools::document() run twice after any structural change
  [ ] One exported function per file, file named after the function
  [ ] No non-ASCII characters in .R file (tools::showNonASCII to verify)
  [ ] Use -- not em/en dash; ... not ellipsis; straight quotes only
  [ ] @seealso links only reference functions that already exist
  [ ] Verification greps: count occurrences in new text AND existing file
      before stating expected counts -- never guess from memory


================================================================================
KEY GOTCHAS (read before touching these areas)
================================================================================

data.table composite key join:
  NEVER dt1[dt2, on="col"] when dt2 has composite key.
  Use: idx <- match(dt1$col, dt2$col); dt1[, new := dt2$other[idx]]

data.table reference semantics in add_*() functions:
  data.table uses pass-by-reference for := assignments. When a function
  receives a NICHESObject and modifies a data.table slot with :=, it will
  silently mutate the caller's original object even though R appears to
  pass by value. Fix: always call copy() on the slot before any := assignment.
    niches.obj$cell.data <- copy(niches.obj$cell.data)
    niches.obj$cell.data[, new.col := ...]
  This was hit in add_meta_column() and add_cell_rawcounts(). Apply
  defensively to any new add_*() function that modifies a data.table slot.

socket cluster + asNamespace():
  clusterExport(cl, c("compute_LRM_Score",".parse_complex",".genes_available"),
                envir = asNamespace("NICHESv2"))
  Works under both devtools::load_all() and installed builds. Exporting a
  function from the package does NOT remove the need to clusterExport it
  -- workers still need the symbol in their global env.
  For local functions defined inside the calling function, use
  envir = environment() instead of asNamespace().

c() on named lists prepends outer names with dot separator:
  c(list("1" = list(a = 1))) -> list("1.a" = 1)
  Fix: unname() the outer list before c(). Hit this in extract_neighborhoods()
  parallel merge path.

sf dependency: keep in Suggests, NOT Imports.

Barcode uniqueness across samples: hard error in both create_NICHESObject()
  and merge_NICHESObjects(); tip: paste0(sample.id,"_",barcode).

Outer sample loop is ALWAYS sequential -- inner compute_CellToCell() already
  uses n.cores; nested socket clusters must never be created.

FANTOM5 must NEVER be loaded via OmniPath's import_intercell_network():
  That API always returns human gene symbols regardless of the datasets
  filter. FANTOM5 must always be loaded from bundled .rda objects.

OmniPath function: use import_ligrecextra_interactions(organism=<code>),
  NOT import_intercell_network().

Seurat V4 vs V5 object detection:
  NEVER check SeuratObject package version. ALWAYS check:
    inherits(seurat.obj@assays[[assay]], "Assay5")
  TRUE = V5 Assay5 structure; FALSE = V4 Assay structure.

rhdf5 logical attributes:
  rhdf5::h5writeAttribute() does NOT support logical (TRUE/FALSE) values.
  Write boolean attributes as 0L/1L (integer) instead.

h5ad matrix orientation:
  AnnData stores all matrices as obs x vars (cells x genes). rhdf5 reads
  dense HDF5 arrays with reversed dims, so [n_obs, n_vars] HDF5 -> [n_vars,
  n_obs] in R, already genes x cells. Sparse matrices must be transposed
  explicitly after reconstruction.

Style B section markers (## Title ---- and # filename ----):
  These are RStudio-foldable. Don't 'fix' them back to banner blocks during
  future edits.

[.NICHESObject keeps full-tissue values by cell:
  Neighborhood slots and aggregation tables are indexed by cell, so each
  retained cell keeps its entries unchanged and they still describe its whole
  neighborhood. Edge slots keep only edges between retained cells, so listed
  neighborhood edges can be missing from $edge.list. Recomputing would divide
  partial sums by full edge counts, so add_neighborhood_edge_agg() and
  add_neighborhood_composition() stop through
  .stop_if_neighborhood_edges_missing(). Workflow: aggregate, then subset.
  add_neighborhoods(overwrite = TRUE) treats a subset as a separate tissue.
  All four attributes carry over: mode, sample.col, edge.filter.mode,
  remove.autocrine.

Non-ASCII characters in .R files:
  CRAN prohibits non-ASCII characters anywhere in .R files, including
  comments and roxygen strings. Em dashes, en dashes, curly quotes, and
  ellipsis characters are common offenders. Always use:
    -- instead of em/en dash
    ... instead of ellipsis character
    straight quotes only
  Detect with: tools::showNonASCII(readLines("R/file.R"))

inherits() vs .is_NICHESObject():
  .is_NICHESObject() is listed in the design doc but does NOT exist in
  utils.R. Always use inherits(x, "NICHESObject") directly in new functions.

globalVariables pattern for data.table NSE:
  All data.table column names used with NSE must be declared via
  utils::globalVariables() in utils.R (NOT NICHESv2-package.R) to suppress
  R CMD check "no visible binding" notes. Similarly, base R functions used
  without explicit :: (head, setNames, as, median, quantile, sd) need
  @importFrom tags in NICHESv2-package.R.

FANTOM5 original file vs object name mismatch:
  data/ncomms8866_orig.rda contains the R object ncomms8866. Because
  LazyData: false, data("ncomms8866") will fail (looks for ncomms8866.rda).
  Never call data(ncomms8866) in examples or code. Always use
  load_LRM_database("fantom5", species = ...) as the access interface.

S3 method source vs test file naming:
  S3 methods (print, summary, dim, length, [) are defined in
  R/NICHESObject_class.R -- NOT a file called NICHESObject_methods.R.
  The test file test-NICHESObject_methods.R is named after the behavior
  under test (methods), not the source file. Both names are correct and
  intentional.

CI is now live on dev/main via .github/workflows/R-CMD-check.yaml

Neighborhood slots always populated after create_NICHESObject():
  edge.filter.mode and remove.autocrine are fixed (not user-facing) in
  create_NICHESObject(). Users who want different combinations must call
  add_neighborhoods() directly with overwrite=TRUE. This is intentional --
  the advisor decision was to give users less freedom here to avoid
  confusion about neighborhood definitions in the main constructor.

Aggregation slot design (current):
  $aggregations is a named list with at most two entries:
    $neighborhood.edge.agg    -- always replaced on each call (no overwrite param)
    $neighborhood.composition -- always replaced on each call (no overwrite param)
  Both functions follow the Seurat convention: call = recompute + replace.
  add_neighborhood_edge_agg() has no methods or overwrite parameters.
    chunk.size (default 50L) controls LRM batch size for memory management.
  aggregate_NICHESObject() has no methods, overwrite, by.direction, include.self,
  or normalize parameters -- these are all fixed internally.
    chunk.size (default 50L) is passed through to add_neighborhood_edge_agg().

add_cell_meta() / add_meta_column() staleness rule:
  Both functions emit a MESSAGE (not a warning, not an error) when
  $aggregations$neighborhood.composition is non-NULL at call time. The message
  tells the user to recompute. The pipeline is never blocked.
  add_cell_meta() always emits the message (full meta replacement invalidates
  any prior composition result unconditionally).
  add_meta_column() emits a conditional message: "may be stale if <col.name>
  was the cell.type.col used to compute it" -- because only a column that was
  used as the cell type grouping variable actually invalidates the result.

autocrine edge classification in add_neighborhood_edge_agg():
  [RESOLVED] The single "autocrine" category has been split into two:
    self.autocrine  -- SendingCell == ReceivingCell == focal cell.
                       count is always 0 or 1 per (cell, LRM).
    other.autocrine -- SendingCell == ReceivingCell != focal cell.
                       Captures self-loop edges of neighbor cells.
  Classification priority: self-loop check runs first (before out/in/cross),
  then subdivides by whether the focal cell is the self-looping cell.
  Total column count is now 27 (up from 23).

Non-spatial mode and neighborhood composition:
  add_neighborhood_composition() hard-errors for mode="nonspatial". In
  nonspatial mode, edges are statistically paired and do not reflect physical
  proximity, so composition would reflect sampling design not local tissue
  structure. add_neighborhood_edge_agg() works in nonspatial mode but
  interpretation shifts: results reflect signaling among co-occurring cell
  pairs, not a spatial microenvironment. Document this distinction clearly
  in vignettes.

print.NICHESObject neighborhood display:
  The print method now shows neighborhood slot contents explicitly:
    Neighborhoods: N cells | M edges (mode: X, autocrine: included/removed)
  And the slot legend shows all three neighborhood slots individually.
  The old "reserved (NULL)" text has been removed.

ConnectomeDB2025 access pattern:
  connectomedb2025_<species>.rda files are accessed ONLY via
  load_LRM_database("connectomedb2025", species = ...). Direct data() calls
  work but are discouraged (consistent with FANTOM5 pattern). The default
  call load_LRM_database() with no arguments loads human ConnectomeDB2025.
  Column validation (.load_connectomedb2025() checks for "ligand" and
  "receptor" columns) runs before the verbose message, so a corrupt bundled
  object will hard-error with a clear re-run instruction.


================================================================================
TEST SUITE STATUS
================================================================================

Full rewrite completed. All test files pass under both devtools::test() and
R CMD check. Coverage: ~270+ test_that blocks.

File structure (all in tests/testthat/):
  helper-testutils.R                  -- keep (make_grid_coord utility)
  test-compute_LRM_Score.R            test-create_Edgelist_Spatial.R
  test-extract_neighborhoods.R        test-create_Edgelist_SamplingByType.R
  test-add_neighborhoods.R            test-compute_CellToCell.R
  test-merge_NICHESObjects.R
  test-Explore_Radius_Neighborhood.R  test-create_NICHESObject.R
  test-extract_NICHESInputs_Seurat.R  test-load_LRM_database.R
  test-extract_NICHESInputs_AnnData.R test-add_cell_meta.R
  test-add_segmentation.R             test-add_meta_column.R
  test-CellToCell_InSitu.R            test-NICHESObject_methods.R
  test-add_neighborhood_edge_agg.R    test-add_neighborhood_composition.R
  test-aggregate_NICHESObject.R       test-add_cell_rawcounts.R

Fixture strategy:
  - data(NICHESv2_example): add_cell_meta, add_meta_column,
    add_cell_rawcounts, add_segmentation, CellToCell_InSitu,
    NICHESObject_methods, add_neighborhood_edge_agg,
    add_neighborhood_composition, aggregate_NICHESObject
  - Inline synthesis: all other files
  - helper-testdata.R and inst/testdata/ DELETED
  - NICHESv2_example now always has neighborhood slots populated
  - No aggregation slots pre-populated in fixture

Skips: none when all Suggests packages are installed.

R CMD check status: 0 errors | 0 warnings | 1 note
  (N1: timestamp/network issue on dev machine -- unfixable in code, ignore)


================================================================================
PENDING DECISIONS
================================================================================

1. [RESOLVED] extract_neighborhoods integrated into create_NICHESObject()
   with fixed parameters: edge.filter.mode="both_in_neighborhood",
   remove.autocrine=FALSE. Not user-facing in create_NICHESObject().

2. [RESOLVED] Default edge.filter.mode = "both_in_neighborhood".

3. Performance optimization for large datasets: pre-extract gene expression
   vectors once rather than per LRM pair; pre-filter silent LRM pairs.
   Identified but not yet implemented.

4. [RESOLVED] Aggregate_CellToCell_Edges() deleted. Replaced by the current
   aggregation design: aggregate_NICHESObject(), add_neighborhood_edge_agg(),
   add_neighborhood_composition().

5. [RESOLVED] PCA kNN smoothing removed from the package, along with the
   smooth, smooth.k and n.pcs arguments.

6. [RESOLVED] Rename "sampled" mode to "nonspatial": completed. All R source
   files, test files, and attributes updated via sed + manual fix.

7. compare_neighborhood_modes() future function: compare the same aggregation
   statistic across objects built with different neighborhood parameters.
   Returns a standalone data.table (not stored in object). Not yet designed.

8. [RESOLVED] Aggregation test suite complete:
     test-add_neighborhood_edge_agg.R
     test-add_neighborhood_composition.R
     test-aggregate_NICHESObject.R
   All passing.

9. [RESOLVED] vignettes/NICHESv2.Rmd is now the only package vignette. It runs
   on bundled data. Public-data tutorials (stxBrain, ifnb) belong on the website
   and must be executed, never eval=FALSE.

10. score.total.cmean formula in add_neighborhood_edge_agg(): currently
    score.total.cmean = score.total.sum / count.total which is correct.
    Confirm this is the intended formula (re-derived from total, not summed
    from directional cmeans) before closing.

11. [RESOLVED] Split autocrine category in add_neighborhood_edge_agg():
    Completed. Single "autocrine" category replaced by self.autocrine and
    other.autocrine. Column count increased from 23 to 27. Implementation
    uses sparse matrix multiplication with chunked LRM batching (chunk.size
    parameter, default 50L) to bound peak memory on large datasets.
