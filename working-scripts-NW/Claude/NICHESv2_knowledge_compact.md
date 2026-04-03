# NICHESv2 — Compact Reference
# Updated: 2026-04-03

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

NAMESPACE: Never hand-edit. Use @import data.table / @import ggplot2 in
  NICHESv2-package.R. Run devtools::document() TWICE after changes.

FILE DELIVERY RULES (no exceptions):
  Only two delivery modes are permitted:

  MODE 1 — In-dialogue edit (≤3 changes, simple enough to apply manually):
    For each change, state ALL of:
      - File path (e.g. R/utils.R)
      - Exact block/line to find (copy-pasteable old code)
      - Exact replacement (copy-pasteable new code)
    Present changes in a table with columns: File | Location | Old | New

  MODE 2 — Downloadable replacement file:
    Only provide files that are either:
      (a) brand-new files that do not yet exist in R/, tests/, etc., OR
      (b) complete replacements for an existing file
    NEVER provide a file to download that is a partial snippet, helper
    extract, or auxiliary — these cause confusion about where to put them.

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
$neighborhood.cell.list     — NULL (reserved)
$neighborhood.edge.list     — NULL (reserved)
$neighborhood.meta          — NULL (reserved)
$gene.list     chr vector   — rownames(count.mtx)
$gene.meta     data.table   — keyed: gene
$LRM.list      chr vector   — valid "ligand|receptor" strings
$LRM.meta      data.table   — LRM|ligand|receptor|<extras>; keyed: LRM
$segmentation  sf object    — NULL until add_segmentation(); sf in Suggests

Object attributes: attr(obj,"mode") = "spatial"|"sampled"
                   attr(obj,"sample.col") = col name string (if split by sample)


================================================================================
FUNCTION REFERENCE
================================================================================

load_LRM_database(db, ligand.col, receptor.col, verbose)
  db: "omnipath" | "fantom5" | user data.frame
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
  → EdgeLRResult S3: $edges|$cell.meta|$lr.meta|$params
  socket cluster over LRM chunks; helpers via asNamespace("NICHESv2")

compute_PseudoBulk_Connectivity(count.mtx, meta.data, LRM.db, cell.type.col,
                                 ligand.col, receptor.col, method,
                                 normalize.method, n.cores, verbose)
  → plain data.table (NOT NICHESObject): SendingType|ReceivingType|LRM|score|score.norm
  NOTE: "pseudobulk" mode removed from create_NICHESObject(); use this instead.

create_NICHESObject(count.mtx, meta.data, LRM.db,
                    mode="spatial", sample.col=NULL, cell.type.col=NULL,
                    ligand.col, receptor.col, k, rad, seed=42L, max.cells,
                    smooth=FALSE, smooth.k=9L, n.pcs=30L,
                    method, normalize.method, gene.meta, n.cores, verbose)
  mode="spatial" : needs x,y in meta.data; k or rad required
  mode="sampled" : cell.type.col required; smooth=TRUE → PCA kNN smoothing
                   $cell.data ALWAYS uses original raw counts
  sample.col     : splits by sample → sequential per-sample pipeline → merge
                   sets attr(obj,"sample.col"); hard error on barcode collision
  → NICHESObject

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
INTERNAL HELPERS (utils.R, unexported)
================================================================================

.parse_complex(gene.str)          "gene1_gene2" → c("gene1","gene2")
.genes_available(gene.str, set)   all components in set?
.compute_lrm_score(...)           core scorer, one LRM pair
.prop_norm(x)                     x/sum(x); zeros if sum==0
.sparse_emin(A,B)                 element-wise min, sparse matrices
.new_NICHESObject(...)            low-level 14-slot constructor
.build_cell_data(count.mtx, normalize.method)
.build_edge_meta(edge.list.dt, cell.meta)
.smooth_expression(count.mtx, barcodes, k, n.pcs)
  irlba::irlba(..., fastpath=FALSE)  ← Matrix>=1.6-2 fix; NEVER remove
.merge_NICHESObject_list(obj.list, verbose)   shared merge worker
.create_NICHESObject_by_sample(...)            split/loop/merge
.create_NICHESObject_single(...)               original single-pass logic
.coerce_meta(meta.data)                        → keyed data.table with "cell" col
%||%                                           null-coalescing operator


================================================================================
PACKAGE INFRASTRUCTURE
================================================================================

Imports: data.table(>=1.14), Matrix(>=1.5), RANN(>=2.6), dbscan(>=1.1),
         ggplot2(>=3.4), irlba(>=2.3.0), scales(:: only), parallel
Suggests: OmnipathR, sf, testthat(>=3.0), withr

Test files (all passing):
  test-compute_CellToCell.R           test-create_Edgelist_Spatial.R
  test-compute_PseudoBulk_Connectivity.R  test-create_Edgelist_PseudoBulk.R
  test-create_NICHESObject.R          test-create_Edgelist_SamplingByType.R
  test-merge_NICHESObjects.R (NEW)    test-utilities.R
  test-add_segmentation.R             test-Aggregate_CellToCell_Edges.R
  test-CellToCell_InSitu.R

Fixtures: inst/testdata/generate_testdata.R — re-run ONLY when schema changes
  attr(niches.obj,"mode") == "spatial" must be explicit in fixture.
  skip_if_not_installed("sf") in all segmentation tests.

Intentional skips (3): full-batch InSitu (>5 LRMs), irlba-absent path,
  replaced placeholder in test-compute_PseudoBulk_Connectivity.R


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

sf dependency: keep in Suggests, NOT Imports.

Barcode uniqueness across samples: hard error in both create_NICHESObject()
  and merge_NICHESObjects(); tip: paste0(sample.id,"_",barcode).

Outer sample loop is ALWAYS sequential — inner compute_CellToCell() already
  uses n.cores; nested socket clusters must never be created.
