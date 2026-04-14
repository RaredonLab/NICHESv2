# NICHESv2 — Compact Reference
# Updated: 2026-04-09

================================================================================
CONVENTIONS (no exceptions)
================================================================================

Functions       : underscore  e.g. create_NICHESObject(), CellToCell_InSitu()
Variables/args  : dot         e.g. n.cores, bg.color, edge.list
Edge ID         : "SendingCell|ReceivingCell"
LRM ID          : "ligand|receptor"
Gene complexes  : "gene1_gene2"
Score columns   : score (raw)  |  score.norm (normalised, proportional)

PARALLELISM (CRITICAL — no exceptions):
  socket clusters ONLY: parallel::makeCluster + parLapply + on.exit(stopCluster)
  NEVER mclapply (Unix-only, crashes on macOS with ggplot2).
  Always clusterEvalQ: data.table::setDTthreads(1)

DEBUGGING PROTOCOL:
  1. Explain bug/decision fully  2. Wait for confirmation  3. Then produce files
  Before proposing changes, Claude must proactively identify ALL affected files,
  including test files — never wait to be asked. Claude must also independently
  assess whether generate_testdata.R needs re-running (re-run ONLY if schema
  changes; not required for argument additions or logic fixes).

NAMESPACE: Never hand-edit. Use @import data.table / @import ggplot2 in
  NICHESv2-package.R. Run devtools::document() TWICE after changes.

HELPER FUNCTION POLICY:
  Avoid unexported dot-prefix helpers for long or crucial package logic.
  Only use helpers for genuinely small utilities (one-liners, simple coercions).
  Always ask for approval before creating any new helper function.

FILE DELIVERY RULES (no exceptions):
  Only two delivery modes are permitted:

  MODE 1 — In-dialogue edit (≤3 changes, simple enough to apply manually):
    For each change, state ALL of:
      - File path (e.g. R/utils.R)
      - Exact block/line to find
      - Exact replacement
    Present old and new code as formatted fenced code blocks, one each.
    Tables are ONLY for file/location metadata — never put code inside tables.
    If something is meant to be copy-pasted, it must be a clean code block.

  MODE 2 — Downloadable replacement file:
    Only provide files that are either:
      (a) brand-new files that do not yet exist in R/, tests/, etc., OR
      (b) complete replacements for an existing file
    NEVER provide a file to download that is a partial snippet, helper
    extract, or auxiliary — these cause confusion about where to put them.

  ONE EXPORTED FUNCTION PER FILE. Do not group multiple exported functions
  in a single .R file.

  CHANGE SUMMARY TABLE ORDER (when multiple files are affected):
    1. Delete  (files or blocks being removed)
    2. Add     (brand-new files)
    3. Change  (edits to existing files)

  Always use devtools::load_all() + devtools::test() to verify after changes.
  Always run devtools::document() TWICE when any @export or @noRd tag changes.


================================================================================
NICHESObject — 14 SLOTS
================================================================================

$cell.list     chr vector   — barcodes (colnames of count.mtx), all modes
$cell.data     data.table   — SPARSE: cell|gene|score|score.norm; key(cell,gene)
                              score.norm = gene/total within cell; raw counts always
$cell.meta     data.table   — cell|x|y|<meta cols>; key: cell
$edge.list     data.table   — edge|SendingCell|ReceivingCell; key: edge
$edge.data     data.table   — SPARSE: edge|LRM|score|score.norm; key(edge,LRM)
                              score.norm = proportional within edge
$edge.meta     data.table   — edge|SendingCell|ReceivingCell|Sending.*|Receiving.*|Joint.*
$neighborhood.cell.list     — named list: cell → character vector of neighbor barcodes
$neighborhood.edge.list     — named list: cell → character vector of edge IDs
$neighborhood.meta          — data.table: cell|n.neighbors|n.edges; key: cell
$gene.list     chr vector   — rownames(count.mtx)
$gene.meta     data.table   — keyed: gene
$LRM.list      chr vector   — valid "ligand|receptor" strings
$LRM.meta      data.table   — LRM|ligand|receptor|<extras>; keyed: LRM
$segmentation  sf object    — NULL until add_segmentation(); sf in Suggests

Neighborhood slots are NULL until add_neighborhoods() is called.
They store IDs only (not duplicated data) — join to $edge.data/$edge.meta
for scores/coordinates when needed. This avoids memory explosion from
overlapping neighborhoods on large datasets.

Object attributes:
  attr(obj,"mode")             = "spatial"|"sampled"
  attr(obj,"sample.col")       = col name string (if split by sample)
  attr(obj,"edge.filter.mode") = set by add_neighborhoods()


================================================================================
FUNCTION REFERENCE
================================================================================

load_LRM_database(db, species="human", ligand.col, receptor.col, verbose)
  db: "omnipath" | "fantom5" | user data.frame
  species: "human"|"mouse"|"rat"|"pig" — REQUIRED for built-in shorthands
    FANTOM5 supports: human, mouse, rat, pig (loaded from bundled .rda, offline)
    OmniPath supports: human, mouse, rat ONLY (pig → hard error)
    Ignored silently when db is a data.frame
  FANTOM5 loads from bundled package data objects (ncomms8866_*) via
    utils::data() — NO OmniPath call, NO internet dependency.
  OmniPath calls import_ligrecextra_interactions(organism=<code>) — NOT
    import_intercell_network(). Organism codes: human=9606, mouse=10090, rat=10116.
  → data.frame with "ligand", "receptor"

create_Edgelist_Spatial(coord, k, rad, n.cores)
  → data.table(from|to|weight|x.from|y.from|x.to|y.to) + autocrine self-loops
  attrs: mode, param, n.cells, cell.names

create_Edgelist_PseudoBulk(meta.data, cell.type.col)
  → data.table(from|to) with type labels; deterministic

create_Edgelist_SamplingByType(meta.data, cell.type.col, seed=42L, max.cells)
  → data.table(from|to) with barcodes; attrs: mode="sampled", param, n.cells

compute_CellToCell(count.mtx, edge.list, LRM.db, ligand.col, receptor.col,
                   method, n.cores)
  LRM.db here must be a data.frame (already resolved) — string shorthands are
    resolved upstream in .create_NICHESObject_single() before this call.
  → EdgeLRResult S3: $edges|$cell.meta|$lr.meta|$params
  socket cluster over LRM chunks; helpers via asNamespace("NICHESv2")

compute_PseudoBulk_Connectivity(count.mtx, meta.data, LRM.db, cell.type.col,
                                 ligand.col, receptor.col, method,
                                 normalize.method, n.cores, verbose)
  → plain data.table (NOT NICHESObject): SendingType|ReceivingType|LRM|score|score.norm
  NOTE: "pseudobulk" mode removed from create_NICHESObject(); use this instead.

create_NICHESObject(count.mtx, meta.data, LRM.db,
                    mode="spatial", species="human", sample.col=NULL,
                    cell.type.col=NULL, ligand.col, receptor.col,
                    k, rad, seed=42L, max.cells,
                    smooth=FALSE, smooth.k=9L, n.pcs=30L,
                    method, normalize.method, gene.meta, n.cores, verbose)
  mode="spatial" : needs x,y in meta.data; k or rad required
  mode="sampled" : cell.type.col required; smooth=TRUE → PCA kNN smoothing
                   $cell.data ALWAYS uses original raw counts
  species        : passed to load_LRM_database(); ignored if LRM.db is data.frame
  sample.col     : splits by sample → sequential per-sample pipeline → merge
                   sets attr(obj,"sample.col"); hard error on barcode collision
  LRM.db string shorthands are resolved to a data.frame inside
    .create_NICHESObject_single() BEFORE calling compute_CellToCell(), so
    species never needs to propagate into the scoring layer.
  → NICHESObject
  NOTE: neighborhood extraction NOT yet integrated (pending advisor discussion).
        Currently use add_neighborhoods() as a separate step.

extract_NICHESInputs_Seurat(seurat.obj, assay="RNA", layer="counts",
                             verbose=TRUE)                       [EXPORTED]
  Extracts count matrix and cell metadata from a Seurat object (V4 or V5).
  Structure detected from assay CLASS (inherits(assay.obj, "Assay5")), NOT
    from installed SeuratObject package version — a V4 object loaded in a V5
    environment keeps its Assay class and is handled correctly.
  V4 dispatch: GetAssayData(obj, slot=layer, assay=assay)
  V5 dispatch: GetAssayData(obj, layer=layer, assay=assay) + Layers() check
  meta.data: seurat.obj@meta.data returned as-is; NO coordinate fetching.
  Prints spatial coordinate hint when verbose=TRUE — user must add x/y themselves.
  Does NOT return gene.meta — rownames(count.mtx) carry gene names already.
  Dependency: SeuratObject in Suggests (NOT Seurat; lighter, sufficient).
  → list(count.mtx = dgCMatrix, meta.data = data.frame)
  FILE: R/extract_NICHESInputs_Seurat.R

extract_NICHESInputs_AnnData(file, layer="counts", obs.cols=NULL,
                              verbose=TRUE)                      [EXPORTED]
  Reads a .h5ad file using rhdf5 ONLY — no Python, no conda, no basilisk,
    no zellkonverter. Pure C HDF5 bindings. Compatible with R >= 4.2.
  layer="counts" → reads /layers/counts group
  layer="X"      → reads /X group (adata.X)
  layer=<other>  → reads /layers/<other> group; hard error if absent
  h5ad sparse matrices: reads CSR and CSC encoding-type attrs; reconstructs
    dgCMatrix. AnnData is obs x vars (cells x genes); function transposes
    to genes x cells (NICHESv2 convention).
  obs categorical columns: reads /categories + /codes arrays; 0-based codes,
    -1 sentinel for NA.
  obs_names from obs/_index (or attr "_index" if remapped).
  var_names from var/_index (or attr "_index" if remapped).
  meta.data returned as-is; NO coordinate fetching from obsm.
  Prints spatial coordinate hint (rhdf5::h5read(file,"obsm/spatial")) when
    verbose=TRUE — user must add x/y themselves.
  obs.cols: optional character vector to subset metadata columns; hard error
    if any requested column is missing (lists available columns in error).
  Dependency: rhdf5 in Suggests. Install: BiocManager::install("rhdf5")
  → list(count.mtx = dgCMatrix, meta.data = data.frame)
  FILE: R/extract_NICHESInputs_AnnData.R

  DESIGN RATIONALE (do not revisit without good reason):
  - zellkonverter was rejected: despite claiming no-Python support via
    reader="R", it requires rhdf5 + HDF5Array + zellkonverter chain and
    triggers basilisk (conda + Python anndata) by default on macOS ARM64.
  - anndataR was rejected: requires R >= 4.5 (too new for user base 4.2-4.4).
  - rhdf5 alone is the correct answer: single Bioconductor package, pure C,
    no Python dependency at any level, R >= 4.2.

extract_NICHESInputs_Seurat / extract_NICHESInputs_AnnData DESIGN RULES:
  - Both return named list(count.mtx, meta.data) ONLY — no gene.meta.
  - Neither calls create_NICHESObject() internally — they are ingredient
    extractors only; the user explicitly calls create_NICHESObject() next.
  - Neither fetches spatial coordinates automatically; they print a hint
    message explaining where to find coordinates and how to add them.
  - Two separate functions (not one with dispatch) — keeps API transparent.

extract_neighborhoods(edge.list, edge.filter.mode, remove.autocrine,
                      n.cores, verbose)                          [EXPORTED]
  edge.list      : data.table with edge|SendingCell|ReceivingCell
  edge.filter.mode: "both_in_neighborhood" (default) |
                    "any_in_neighborhood" | "central_cell_only"
  remove.autocrine: FALSE default (adjacency only; edges still subject to filter)
  Builds bidirectional adjacency (symmetric neighbors), integer-indexed lookups,
  parallelized per-cell extraction via socket cluster.
  → list with 3 elements:
      $neighborhood.cell.list  named list of chr vectors (neighbor barcodes)
      $neighborhood.edge.list  named list of chr vectors (edge IDs)
      $neighborhood.meta       data.table(cell|n.neighbors|n.edges), keyed: cell
    attrs: "edge.filter.mode", "remove.autocrine"
  FILE: R/extract_neighborhoods.R

add_neighborhoods(niches.obj, edge.filter.mode, remove.autocrine,
                  n.cores, verbose)                              [EXPORTED]
  Calls extract_neighborhoods() on niches.obj$edge.list, writes the 3 slots.
  Overwrites if slots already populated (with message).
  Sets attr(niches.obj,"edge.filter.mode").
  FILE: R/add_neighborhoods.R (to be split from extract_neighborhoods.R)

merge_NICHESObjects(..., verbose=TRUE)  [EXPORTED]
  accepts: merge_NICHESObjects(a,b,c) | merge_NICHESObjects(list(a,b,c)) | var
  hard errors: empty list | non-NICHESObject | mode mismatch | dup barcodes
  slot strategy: rbindlist data slots; union gene/LRM lists; sf::rbind segs
  → NICHESObject

add_segmentation(niches.obj, seg.sf)
  every cell in $cell.list must have polygon — hard error if missing
  extra barcodes in seg.sf → warning only

add_cell_meta(niches.obj, new.meta, rebuild.edge.meta=TRUE)

Aggregate_CellToCell_Edges(niches.obj, stat)
  → list($sending, $receiving): full cell×LRM grid, zero-filled

CellToCell_InSitu(niches.obj, LRM.list, folder, ..., n.cores)
  batch PNG; socket cluster; pre-computes all.edges before spawning workers
  → invisible data.frame(LRM|status|path)

CellToCell_InSitu_Single_LRM(niches.obj, LRM, ..., precomputed.bg=NULL)
  precomputed.bg: skip coord joins (supplied by batch caller, NULL for direct use)


================================================================================
INTERNAL HELPERS (utils.R, unexported) — REVIEW PENDING
================================================================================

Genuinely small utilities (approved):
  .parse_complex(gene.str)          "gene1_gene2" → c("gene1","gene2")
  .genes_available(gene.str, set)   all components in set?
  .prop_norm(x)                     x/sum(x); zeros if sum==0
  .sparse_emin(A,B)                 element-wise min, sparse matrices
  .new_NICHESObject(...)            low-level 14-slot constructor
  .coerce_meta(meta.data)           → keyed data.table with "cell" col
  %||%                              null-coalescing operator

Substantial helpers (flagged for review — may warrant export or restructure):
  .compute_lrm_score(...)           core scorer, one LRM pair
  .build_cell_data(count.mtx, normalize.method)
  .build_edge_meta(edge.list.dt, cell.meta)
  .smooth_expression(count.mtx, barcodes, k, n.pcs)
    irlba::irlba(..., fastpath=FALSE)  ← Matrix>=1.6-2 fix; NEVER remove
  .merge_NICHESObject_list(obj.list, verbose)   shared merge worker
  .create_NICHESObject_by_sample(...)            split/loop/merge
  .create_NICHESObject_single(...)               original single-pass logic

h5ad internal helpers (in R/extract_NICHESInputs_AnnData.R, unexported):
  .h5ad_attrs(file, path)           reads HDF5 attributes from any path
  .h5ad_encoding(file, path)        extracts encoding-type attribute
  .h5ad_obs_names(file)             cell barcodes from obs/_index
  .h5ad_var_names(file)             gene names from var/_index
  .h5ad_available_layers(file)      "X" + names of /layers children
  .h5ad_read_sparse(file, path)     reconstructs CSR/CSC → dgCMatrix (obs x vars)
  .h5ad_read_matrix(file, path, n.obs, n.vars)  dispatches sparse/dense;
                                    always returns genes x cells (vars x obs)
  .h5ad_read_obs(file, obs.names)   reads obs group; handles categoricals


================================================================================
PACKAGE INFRASTRUCTURE
================================================================================

Imports: data.table(>=1.14), Matrix(>=1.5), RANN(>=2.6), dbscan(>=1.1),
         ggplot2(>=3.4), irlba(>=2.3.0), scales(:: only), parallel
Suggests: OmnipathR, sf, testthat(>=3.0), withr, SeuratObject, rhdf5

Bundled data (data/):
  ncomms8866_human.rda, ncomms8866_mouse.rda, ncomms8866_rat.rda,
  ncomms8866_pig.rda — FANTOM5 LR pairs, species-converted per Raredon et al.
  2019. Loaded via utils::data(). No internet or OmnipathR required for FANTOM5.

Test files (all passing):
  test-compute_CellToCell.R               test-create_Edgelist_Spatial.R
  test-compute_PseudoBulk_Connectivity.R  test-create_Edgelist_PseudoBulk.R
  test-create_NICHESObject.R              test-create_Edgelist_SamplingByType.R
  test-merge_NICHESObjects.R              test-utilities.R
  test-add_segmentation.R                 test-Aggregate_CellToCell_Edges.R
  test-CellToCell_InSitu.R
  test-extract_NICHESInputs_Seurat.R      (NEW — no fixture needed)
  test-extract_NICHESInputs_AnnData.R     (NEW — fixture written inline via rhdf5)

  TODO: test-extract_neighborhoods.R, test-add_neighborhoods.R (not yet written)

Fixtures: inst/testdata/generate_testdata.R — re-run ONLY when schema changes
  attr(niches.obj,"mode") == "spatial" must be explicit in fixture.
  skip_if_not_installed("sf") in all segmentation tests.

Intentional skips (3): full-batch InSitu (>5 LRMs), irlba-absent path,
  replaced placeholder in test-compute_PseudoBulk_Connectivity.R

Test fixture pattern for rhdf5-based tests:
  .write_test_h5ad() helper in test-extract_NICHESInputs_AnnData.R writes
  a fully spec-compliant h5ad using rhdf5 only. Key gotcha: rhdf5 does NOT
  support logical attributes — write boolean attrs as 0L/1L (integer), not
  FALSE/TRUE. This affects the "ordered" attribute on categorical groups.


================================================================================
KEY GOTCHAS (read before touching these areas)
================================================================================

irlba/Matrix:  fastpath=FALSE in .smooth_expression() — NEVER change.
               Matrix>=1.6-2 removed as_cholmod_sparse(); fastpath calls it.

data.table composite key join:
  NEVER dt1[dt2, on="col"] when dt2 has composite key.
  Use: idx <- match(dt1$col, dt2$col); dt1[, new := dt2$other[idx]]

socket cluster + asNamespace():
  clusterExport(cl, c(".compute_lrm_score",".parse_complex",".genes_available"),
                envir = asNamespace("NICHESv2"))
  Works under both devtools::load_all() and installed builds.

c() on named lists prepends outer names with dot separator:
  c(list("1" = list(a = 1))) → list("1.a" = 1)
  Fix: unname() the outer list before c(). Hit this in extract_neighborhoods()
  parallel merge path — split() names chunks "1","2",... and parLapply
  preserves those names, causing "1.barcode" prefixes.

sf dependency: keep in Suggests, NOT Imports.

Barcode uniqueness across samples: hard error in both create_NICHESObject()
  and merge_NICHESObjects(); tip: paste0(sample.id,"_",barcode).

Outer sample loop is ALWAYS sequential — inner compute_CellToCell() already
  uses n.cores; nested socket clusters must never be created.

FANTOM5 must NEVER be loaded via OmniPath's import_intercell_network():
  That API always returns human gene symbols regardless of the datasets filter,
  making non-human species silently return 0 valid LR pairs. FANTOM5 must
  always be loaded from the bundled ncomms8866_* .rda objects.

OmniPath function: use import_ligrecextra_interactions(organism=<code>), NOT
  import_intercell_network(). The latter is broader, not species-aware in the
  same call path, and was the wrong choice for a targeted LR database.

Seurat V4 vs V5 object detection:
  NEVER check SeuratObject package version. A V4 object loaded in a V5
  environment keeps its Assay class internally — the package version tells you
  nothing about the object structure. ALWAYS check:
    inherits(seurat.obj@assays[[assay]], "Assay5")
  TRUE = V5 Assay5 structure; FALSE = V4 Assay structure.

rhdf5 logical attributes:
  rhdf5::h5writeAttribute() does NOT support logical (TRUE/FALSE) values.
  Write boolean attributes as 0L/1L (integer) instead. Affects h5ad
  "ordered" attribute on categorical obs columns.

h5ad matrix orientation:
  AnnData stores all matrices as obs x vars (cells x genes). rhdf5 reads
  dense HDF5 arrays with reversed dims (C-order reversal), so a [n_obs, n_vars]
  HDF5 array comes back as [n_vars, n_obs] in R — already genes x cells. Sparse
  matrices must be transposed explicitly after reconstruction. The defensive
  orientation check in .h5ad_read_matrix() uses known n.obs/n.vars to verify.


================================================================================
PENDING DECISIONS
================================================================================

1. Integrate extract_neighborhoods into create_NICHESObject() (Option C):
   add edge.filter.mode=NULL parameter; when non-NULL, call
   add_neighborhoods() at the end. PENDING advisor discussion.

2. Default edge.filter.mode: currently "both_in_neighborhood"; confirm with
   advisor.

3. Split extract_neighborhoods.R into two files (extract_neighborhoods.R +
   add_neighborhoods.R): do after validation is confirmed.

4. Helper function audit: review the 7 substantial helpers listed above for
   possible export or restructure. Nuoya will open a separate chat.

5. Performance optimization for large datasets: pre-extract gene expression
   vectors once rather than per LRM pair; pre-filter silent LRM pairs.

6. Aggregate_CellToCell_Edges future design: sparse dgCMatrix (cells × LRMs)
   plus inverted index $active — pending advisor discussion.

7. Resolve PCA kNN smoothing validation dataset decision (own Visium data vs
   public benchmark) before writing validation scripts.

8. Potential package rename from NICHESv2.
