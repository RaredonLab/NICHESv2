# NICHESv2 — Base Knowledge Document
# For use as context in new chat sessions
# Updated: 2026-04-02

================================================================================
PART 1 — NAMING CONVENTIONS (apply everywhere, no exceptions)
================================================================================

- Function names     : underscore separator    e.g. create_NICHESObject()
- Variable / argument names : dot separator    e.g. n.cores, bg.color, edge.list
- Score columns      : score (raw), score.norm (normalised)
- Edge ID format     : "SendingCell|ReceivingCell"   (pipe separator)
- LRM ID format      : "ligand|receptor"              (pipe separator)
- Gene complexes in LRM database : "gene1_gene2"      (underscore separator)

Strict debugging protocol (user-specified):
  1. Explain the bug / design decision in full detail
  2. Wait for user confirmation of understanding
  3. Only then produce files — and only the specific files that changed

Parallelisation rule (CRITICAL — applies everywhere without exception):
  Always use socket clusters: parallel::makeCluster + parLapply + stopCluster
  NEVER use mclapply or any fork-based method.
  Rationale: socket clusters work on Windows, macOS, and Linux.
  mclapply is Unix-only and crashes on macOS with ggplot2 font rendering.
  Always add data.table::setDTthreads(1) inside clusterEvalQ to prevent
  OpenMP thread contention across workers.
  Always use on.exit(parallel::stopCluster(cl), add = TRUE) for cleanup.

Workflow for iterating on files across sessions:
  - Attach files in chat for active iteration (guarantees latest version;
    knowledge base copies may be stale if edited locally since last upload).
  - Upload to knowledge base only stable reference files (e.g. this .md).
  - Keep only ONE .md file in the knowledge base — delete old versions to
    prevent a future session loading stale context.


================================================================================
PART 2 — NICHESObject: THE CENTRAL S3 CLASS
================================================================================

Class name : NICHESObject
Storage    : named R list with S3 class attribute "NICHESObject"

Constructed by : create_NICHESObject()   (main entry point)
Identity check : is_NICHESObject(x)

Object-level attribute:
  attr(obj, "mode") — records which analysis mode was used:
    "spatial" | "sampled"
  Set at the end of create_NICHESObject(); NOT a named slot.
  NOTE: "pseudobulk" is NO LONGER a mode of create_NICHESObject().
  Pseudobulk is handled by compute_PseudoBulk_Connectivity() instead.

S3 methods registered:
  print(niches.obj)      — concise slot summary with non-empty/valid counts
  summary(niches.obj)    — cell type breakdown, edge stats, signal density,
                           segmentation status
  dim(niches.obj)        — named integer vector: cells / genes / edges / LRMs
  length(niches.obj)     — number of cells
  niches.obj[cell.vec]   — subset to a character vector of cell barcodes;
                           filters all cell-indexed and edge-indexed slots,
                           keeps only edges where BOTH endpoints are in cell.vec;
                           also filters $segmentation if populated

print.NICHESObject output format (current):
  Cells        : <n>  (<n non-empty>)
  Genes        : <n>  (<n valid>)
  Edges        : <n>  (<n non-empty>)
  LRM pairs    : <n>  (<n valid>)
  Segmentation : none | attached (<n> cells)

  Definitions:
    non-empty cells  = cells that appear at least once in $cell.data (score > 0)
    valid genes      = genes that appear at least once in $cell.data
    non-empty edges  = edges that appear at least once in $edge.data (score > 0)
    valid LRM pairs  = LRM pairs that appear at least once in $edge.data
  All four computed via uniqueN() on the respective sparse table column.

--------------------------------------------------------------------------------
SLOT REFERENCE (all 14 slots)
--------------------------------------------------------------------------------

SLOT 1  $cell.list
  Type    : character vector
  Content : cell barcodes from the original count matrix (colnames)
            Always original single-cell barcodes in ALL modes.
  Notes   : authoritative cell registry; length = n cells

SLOT 2  $cell.data
  Type    : data.table
  Columns : cell | gene | score | score.norm
  Sparsity: SPARSE — only rows where score > 0 OR score.norm > 0 are stored
  Notes   : score      = raw count value from count matrix
            score.norm = proportional normalisation WITHIN each cell
                         (gene count / total cell count); sums to 1.0 per cell
            Always built from ORIGINAL (unsmoothed) count.mtx — even when
            smooth = TRUE in sampled mode.
  Key     : (cell, gene)

SLOT 3  $cell.meta
  Type    : data.table
  Columns : cell | x | y | <any other columns from input meta.data>
  Notes   : one row per cell; always original cell barcodes.
            x and y present only when mode = "spatial".
  Key     : cell

SLOT 4  $edge.list
  Type    : data.table
  Columns : edge | SendingCell | ReceivingCell
  Notes   : edge = paste(SendingCell, ReceivingCell, sep = "|")
            Both modes use cell barcodes as edge nodes.
            Includes autocrine self-loops (SendingCell == ReceivingCell).
  Key     : edge

SLOT 5  $edge.data
  Type    : data.table
  Columns : edge | LRM | score | score.norm
  Sparsity: SPARSE
  Notes   : score.norm sums to 1.0 per edge (proportional).
            When smooth = TRUE, scores reflect smoothed expression but
            $cell.data always retains original raw counts.
  Key     : (edge, LRM)

SLOT 6  $edge.meta
  Type    : data.table
  Columns : edge | SendingCell | ReceivingCell
            | Sending.<col>  | Receiving.<col>  | Joint.<col> (categorical only)
  Key     : edge

SLOT 7  $neighborhood.cell.list   — reserved, always NULL
SLOT 8  $neighborhood.edge.list   — reserved, always NULL
SLOT 9  $neighborhood.meta        — reserved, always NULL

SLOT 10 $gene.list   — character vector; all rownames of count.mtx
SLOT 11 $gene.meta   — data.table keyed on gene
SLOT 12 $LRM.list    — character vector of valid "ligand|receptor" strings
SLOT 13 $LRM.meta    — data.table: LRM | ligand | receptor | <extras>; keyed LRM
SLOT 14 $segmentation — sf object or NULL; populated by add_segmentation() only


================================================================================
PART 3 — FUNCTION REFERENCE
================================================================================

--------------------------------------------------------------------------------
load_LRM_database(db, ligand.col, receptor.col, verbose)
--------------------------------------------------------------------------------
PURPOSE : Load and standardise an LR pair database.
INPUT   : db — "omnipath", "fantom5", or a user data.frame
OUTPUT  : data.frame with "ligand", "receptor", plus any extras.

--------------------------------------------------------------------------------
create_Edgelist_Spatial(coord, k, rad, n.cores)
--------------------------------------------------------------------------------
PURPOSE : Directed spatial neighbourhood edge list from 2D coordinates.
OUTPUT  : data.table (from|to|weight|x.from|y.from|x.to|y.to)
          Autocrine self-loops appended. Attributes: mode, param, n.cells, cell.names
PARALLELISM: socket cluster — cells chunked.

--------------------------------------------------------------------------------
create_Edgelist_PseudoBulk(meta.data, cell.type.col)
--------------------------------------------------------------------------------
PURPOSE : T² type-to-type directed edge list (fully deterministic).
OUTPUT  : data.table (from|to) where values are type labels.
          Attributes: mode="pseudobulk", param, n.cells, cell.names (sorted types)

--------------------------------------------------------------------------------
create_Edgelist_SamplingByType(meta.data, cell.type.col, seed, max.cells)
--------------------------------------------------------------------------------
PURPOSE : Type-stratified random sampling (NICHES v1 strategy).
          cap per crossing = min(n_i, n_j) [bounded by max.cells].
          Independent samples from each type; positional pairing.
INPUT   : seed default 42L; max.cells default NULL
OUTPUT  : data.table (from|to) where values are cell barcodes.
          Attributes: mode="sampled", param=list(...), n.cells, cell.names

--------------------------------------------------------------------------------
compute_CellToCell(count.mtx, edge.list, LRM.db,
                   ligand.col, receptor.col, method, n.cores)
--------------------------------------------------------------------------------
PURPOSE : Universal scoring engine for spatial and sampled modes.
OUTPUT  : EdgeLRResult S3:
    $edges     — from|to|ligand|receptor|LRM|score [+coords if spatial]
    $cell.meta — cell [|x|y]
    $lr.meta   — LRM|ligand|receptor|<extras>
    $params    — method, n.cells, n.edges, n.lr, n.nonzero, created.at
PARALLELISM: socket cluster — LRM chunks; internal helpers via asNamespace().

--------------------------------------------------------------------------------
compute_PseudoBulk_Connectivity(count.mtx, meta.data, LRM.db,
                                 cell.type.col, ligand.col, receptor.col,
                                 method, normalize.method, n.cores, verbose)
--------------------------------------------------------------------------------
PURPOSE : Pseudobulk type-to-type connectivity. Returns a plain data.table,
          NOT a NICHESObject. This is the replacement for the removed
          mode = "pseudobulk" in create_NICHESObject().
OUTPUT  : data.table keyed (SendingType, ReceivingType, LRM):
          SendingType|ReceivingType|LRM|ligand|receptor|score|score.norm
          score.norm sums to 1.0 per type-pair. Only non-zero rows stored.
INTERNALS: create_Edgelist_PseudoBulk() + rowMeans matrix + compute_CellToCell()

--------------------------------------------------------------------------------
create_NICHESObject(count.mtx, meta.data, LRM.db,
                    mode, cell.type.col,
                    ligand.col, receptor.col,
                    k, rad, seed, max.cells,
                    smooth, smooth.k, n.pcs,
                    method, normalize.method, gene.meta, n.cores, verbose)
--------------------------------------------------------------------------------
PURPOSE : Main constructor. Two modes:

  "spatial" (default):
    - meta.data must have "x" and "y" (hard error otherwise)
    - calls create_Edgelist_Spatial(); scores original count.mtx

  "sampled":
    - cell.type.col required (hard error otherwise)
    - calls create_Edgelist_SamplingByType()
    - smooth = FALSE (default): scores original count.mtx directly
    - smooth = TRUE: calls .smooth_expression() before scoring;
      smooth.k = 9L (neighbours), n.pcs = 30L (PCA components)
    - $cell.data always uses original raw counts regardless of smooth

  "pseudobulk" mode REMOVED — use compute_PseudoBulk_Connectivity() instead.

OUTPUT  : NICHESObject with attr(obj, "mode") = "spatial" or "sampled".

--------------------------------------------------------------------------------
add_segmentation(niches.obj, seg.sf)
--------------------------------------------------------------------------------
Every cell in $cell.list must have a polygon — hard error if missing (barcode
mismatch). Extra barcodes in seg.sf not in $cell.list → warning only.

--------------------------------------------------------------------------------
add_cell_meta(niches.obj, new.meta, rebuild.edge.meta = TRUE)
--------------------------------------------------------------------------------
Add/overwrite $cell.meta columns; optionally rebuild $edge.meta.

--------------------------------------------------------------------------------
Aggregate_CellToCell_Edges(niches.obj, stat)
--------------------------------------------------------------------------------
Collapse $edge.data to per-cell sending/receiving stats (sum/mean/max).
Full cell × LRM grid with zero-fill. Returns list: $sending, $receiving.

--------------------------------------------------------------------------------
CellToCell_InSitu_Single_LRM(niches.obj, LRM, value.col, scale.upper,
                              bg.color, bg.linewidth, bg.alpha,
                              sig.linewidth, sig.alpha, viridis.option,
                              title, precomputed.bg = NULL)
--------------------------------------------------------------------------------
Two-layer ggplot segment plot for one LRM pair.
precomputed.bg: pre-built background data.table (skips Section 3 joins).
  Supplied internally by CellToCell_InSitu() for batch speed.
  Direct callers leave NULL.

--------------------------------------------------------------------------------
CellToCell_InSitu(niches.obj, LRM.list, folder, value.col, scale.upper,
                  bg.color, bg.linewidth, bg.alpha, sig.linewidth, sig.alpha,
                  viridis.option, width, height, dpi, n.cores)
--------------------------------------------------------------------------------
Batch PNG export. Socket cluster — all platforms.
3-step performance design before spawning workers:
  1. $cell.data → NULL; $edge.data filtered to plotted LRMs only
  2. all.edges (coordinate joins) pre-computed once; passed via precomputed.bg
  3. setDTthreads(1) on each worker
Returns: invisible data.frame LRM|status|path.

--------------------------------------------------------------------------------
as_matrix(x)  EdgeLRResult → dgCMatrix (edges × LRMs). Not stored; on demand.
--------------------------------------------------------------------------------


================================================================================
PART 4 — INTERNAL HELPERS (not exported, defined in utils.R)
================================================================================

.parse_complex(gene.str)       — split "gene1_gene2" → c("gene1","gene2")
.genes_available(gene.str, gene.set) — TRUE if all components in gene.set
.compute_lrm_score(...)        — core scorer for one LRM pair
.categorical_cols(dt, exclude) — names of character/factor columns
.prop_norm(x)                  — x / sum(x); all zeros if sum == 0
.sparse_emin(A, B)             — element-wise min of sparse matrices
.new_NICHESObject(...)         — low-level constructor; validates 14 slots
.build_cell_data(count.mtx, normalize.method) — sparse → long data.table
.build_edge_meta(edge.list.dt, cell.meta) — Sending/Receiving/Joint columns

.smooth_expression(count.mtx, barcodes, k, n.pcs)
  PCA-based neighbourhood smoothing for sampled mode. Five steps:
  1. Library-size normalise + log1p → sparse norm.mtx (temporary)
  2. irlba::irlba(Matrix::t(norm.mtx), nv=n.pcs, fastpath=FALSE)
     fastpath=FALSE bypasses C/CHOLMOD code that calls as_cholmod_sparse
     (removed in Matrix >= 1.6-2). Pure R matrix-vector products; sparse
     matrix never densified.
  3. RANN::nn2 kNN in PC space (query = sampled barcodes only)
  4. Average RAW counts (from original count.mtx) over k+1 neighbourhood
  5. Replace sampled columns in count.mtx; return modified matrix
  PCA embedding discarded immediately after step 3 (rm()).
  $cell.data in the returned object always uses original unsmoothed counts.


================================================================================
PART 5 — PACKAGE INFRASTRUCTURE
================================================================================

Dependencies (Imports):
  data.table (>= 1.14)  — import(data.table) in NAMESPACE (not importFrom)
  Matrix     (>= 1.5)   — sparse matrix types
  RANN       (>= 2.6)   — nn2() for kNN
  dbscan     (>= 1.1)   — frNN() for radius search
  ggplot2    (>= 3.4)   — import(ggplot2) in NAMESPACE
  irlba      (>= 2.3.0) — truncated SVD in .smooth_expression()
  scales                — scales::squish; :: notation only; DESCRIPTION only
  parallel              — socket clusters

Optional (Suggests):
  OmnipathR, sf, testthat (>= 3.0), withr

Test suite (11 files, all passing):
  test-compute_CellToCell.R
  test-compute_PseudoBulk_Connectivity.R   (NEW)
  test-create_Edgelist_Spatial.R
  test-create_Edgelist_PseudoBulk.R
  test-create_Edgelist_SamplingByType.R
  test-create_NICHESObject.R               (pseudobulk tests removed;
                                            smoothing tests added)
  test-utilities.R
  test-add_segmentation.R
  test-Aggregate_CellToCell_Edges.R
  test-CellToCell_InSitu.R

Intentional skips (3, always expected):
  1. Full-batch InSitu test when fixture has > 5 LRM pairs
  2. irlba-not-installed error path (skips when irlba is present)
  3. Empty test slot in test-compute_PseudoBulk_Connectivity.R line 70
     (replaced with a valid test in current version)

Fixtures: source("inst/testdata/generate_testdata.R")
  mode = "spatial" explicit → attr(niches.obj,"mode") == "spatial" in fixture.
  Re-run only when NICHESObject SCHEMA changes (not routine code changes).

Files DELETED from R/:
  compute_CellToCell_Spatial.R → compute_CellToCell.R

Files DELETED from tests/testthat/:
  test-compute_CellToCell_Spatial.R → test-compute_CellToCell.R

NAMESPACE management:
  devtools::document() — run TWICE after export changes.
  @import data.table and @import ggplot2 in R/NICHESv2-package.R.
  :: notation functions need only DESCRIPTION Imports, not NAMESPACE.

Known irlba / Matrix incompatibility (RESOLVED):
  Matrix >= 1.6-2 removed as_cholmod_sparse(). irlba's C fastpath calls
  this regardless of centering arguments. Fix: fastpath = FALSE in
  irlba::irlba() call inside .smooth_expression(). Requires irlba >= 2.3.0.
  On R 4.3 / macOS ARM64: install.packages("irlba") gives binary 2.3.5.1
  (sufficient). Source install fails without gfortran toolchain.

Known data-table composite-key join rule:
  Never dt1[dt2, on="col"] when dt2 has composite key. Use match():
    idx <- match(dt1$col, dt2$col)
    dt1[, new_col := dt2$other_col[idx]]

Known socket cluster + asNamespace() pattern:
  parallel::clusterExport(cl,
    varlist = c(".compute_lrm_score", ".parse_complex", ".genes_available"),
    envir = asNamespace("NICHESv2"))
  Works under devtools::load_all() and installed builds.

Known database gaps (not code bugs):
  Some LRM pairs absent from Connectome2025 human LR database
  (e.g. TGFBI-integrin pairs, THBS2-CD47).


================================================================================
PART 6 — NON-SPATIAL ANALYSIS DESIGN DECISIONS
================================================================================

Two entry points for non-spatial scRNA data:

compute_PseudoBulk_Connectivity()
  Question: "What does type A communicate to type B on average?"
  Output: plain data.table (not NICHESObject) — T² type pairs × LRMs
  Reproducible: YES; Best for: quick summary, cross-dataset comparison

create_NICHESObject(mode = "sampled")
  Question: "Distribution of possible signals across the population"
  Edge unit: cell barcode pairs, type-stratified (NICHES v1 strategy)
  Optional: smooth = TRUE — PCA kNN neighbourhood averaging of raw counts
    before scoring; k+1 cells averaged; PCA never stored
  Reproducible: only with fixed seed (default 42L)
  Best for: NICHES v1 compatibility, full population signal space

REJECTED — kNN in embedding space as edge definition:
  Transcriptomic kNN wires same-type cells preferentially, suppressing
  cross-type edges (e.g. ATI → gCap) that are biologically interesting.
  Rejected by advisor. Stratified sampling covers all type crossings equally.

REJECTED — user-supplied pre-computed embedding for smoothing:
  Complex validation surface (batch correction, normalisation, barcode coverage).
  Internal PCA is simpler, consistent, and sufficient for kNN neighbour finding.


================================================================================
PART 7 — SAMPLE-AWARE SPLITTING AND MERGING (added 2026-04-02)
================================================================================

--------------------------------------------------------------------------------
create_NICHESObject() — new argument: sample.col
--------------------------------------------------------------------------------
Argument : sample.col = NULL (character or NULL)
Purpose  : Name of the column in meta.data identifying the biological sample of
           origin for each cell.
Behaviour:
  - NULL (default): unchanged, single-pass behaviour as before.
  - Non-NULL: the function splits meta.data by sample.col, runs the full
    pipeline (.create_NICHESObject_single) on each partition sequentially,
    then merges the results via .merge_NICHESObject_list().
  - attr(obj, "sample.col") is set on the merged object.
  - Verbose message (when verbose = TRUE):
      "Cell-cell connectivity will be computed per '<sample.col>'.
       Edges exist only between cells within the same '<sample.col>'."
  - Samples with fewer than 2 cells are skipped with a warning.
  - Barcodes MUST be unique across samples. Hard error on collision.
    Tip printed: prepend sample prefix.

Parallelism: outer loop over samples is ALWAYS sequential (inner
  compute_CellToCell() already uses n.cores workers). Nested socket clusters
  must never be created.

New internal helpers added (unexported):
  .create_NICHESObject_by_sample(...)  — split/loop/merge orchestration
  .create_NICHESObject_single(...)     — encapsulates original single-pass logic
  .coerce_meta(meta.data)              — converts meta.data to keyed data.table
                                         with guaranteed "cell" column

$edge.meta sample propagation:
  sample.col flows naturally through .build_edge_meta():
  $cell.meta will have a sample column, so $edge.meta automatically gets
  Sending.<sample.col>, Receiving.<sample.col>, Joint.<sample.col>.
  No extra code needed.

--------------------------------------------------------------------------------
merge_NICHESObjects(..., verbose = TRUE)  [EXPORTED]
--------------------------------------------------------------------------------
Purpose : Merge a list of NICHESObjects into one. Complements sample.col
          splitting and can also combine independently-built objects.

Calling styles (all equivalent):
  merge_NICHESObjects(a, b, c)           # direct positional args
  merge_NICHESObjects(list(a, b, c))     # pre-built list
  merge_NICHESObjects(my.list)           # variable holding a list

Disambiguation logic in the public function:
  if length(args)==1 && is.list(args[[1]]) && !is_NICHESObject(args[[1]]):
    obj.list <- args[[1]]
  else:
    obj.list <- args   (treats each ... element as one object)

Hard errors:
  - Empty obj.list
  - Any element is not a NICHESObject
  - Mode mismatch across objects (attr "mode" must be identical)
  - Duplicate barcodes across objects (names offending barcodes, prints tip)

Slot merge strategy:
  $cell.list              c(...)
  $cell.data              rbindlist(); re-key (cell, gene)
  $cell.meta              rbindlist(fill=TRUE); re-key cell
  $edge.list              rbindlist(); re-key edge
  $edge.data              rbindlist(); re-key (edge, LRM)
  $edge.meta              rbindlist(fill=TRUE); re-key edge
  $gene.list              Reduce(union, ...)  — warn if lists differ
  $gene.meta              rbindlist(fill=TRUE) → unique(by="gene") first occ
  $LRM.list               Reduce(union, ...)  — warn if lists differ
  $LRM.meta               rbindlist(fill=TRUE) → unique(by="LRM") first occ
  $segmentation           sf::rbind() if ALL non-NULL; NULL + warning if partial
  neighborhood slots      always NULL (reserved)

Attribute propagation:
  attr(., "mode")         from objects[1] (validated identical across all)
  attr(., "sample.col")   propagated iff all objects share the same value;
                          warn + drop if inconsistent or partially present

Implementation:
  Public merge_NICHESObjects() disambiguates calling style then delegates to
  .merge_NICHESObject_list() (unexported, in utils.R). Both create_NICHESObject
  (sample.col path) and merge_NICHESObjects() share this single helper.

--------------------------------------------------------------------------------
New files
--------------------------------------------------------------------------------
  R/merge_NICHESObjects.R              — exported public function
  tests/testthat/test-merge_NICHESObjects.R — 23 tests across 8 sections

Modified files
  R/create_NICHESObject.R              — sample.col argument + branching;
                                         original logic now in
                                         .create_NICHESObject_single()
  R/utils.R                            — add .merge_NICHESObject_list(),
                                         .create_NICHESObject_by_sample(),
                                         .create_NICHESObject_single(),
                                         .coerce_meta()
  tests/testthat/test-create_NICHESObject.R — 3 new sections (sample.col
                                         validation, functional, edge cases)

%||% operator: already present in utils.R (or add once if missing).
