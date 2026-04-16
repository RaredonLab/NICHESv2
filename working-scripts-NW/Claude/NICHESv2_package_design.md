# NICHESv2 — Package Design & Status
# Updated: 2026-04-10

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

NAMESPACE: Never hand-edit. Use @import data.table / @import ggplot2 in
  NICHESv2-package.R. Run devtools::document() TWICE after changes.

ONE EXPORTED FUNCTION PER FILE.


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
They store IDs only — join to $edge.data/$edge.meta for scores/coordinates.

Object attributes:
  attr(obj,"mode")             = "spatial"|"sampled"
  attr(obj,"sample.col")       = col name string (if split by sample)
  attr(obj,"edge.filter.mode") = set by add_neighborhoods()


================================================================================
EXPORTED FUNCTIONS
================================================================================

Each lives in its own file under R/.

— LR DATABASE —
load_LRM_database(db, species, ligand.col, receptor.col, verbose)
  db: "omnipath" | "fantom5" | user data.frame
  FANTOM5: human/mouse/rat/pig from bundled .rda; no internet, no OmnipathR
  OmniPath: human/mouse/rat via import_ligrecextra_interactions(organism=<code>)
  → data.frame with "ligand", "receptor"

— EDGE LIST CONSTRUCTORS —
create_Edgelist_Spatial(coord, k, rad, n.cores)
  → data.table(from|to|weight|x.from|y.from|x.to|y.to) + autocrine self-loops

create_Edgelist_PseudoBulk(meta.data, cell.type.col)
  → data.table(from|to) with type labels; deterministic

create_Edgelist_SamplingByType(meta.data, cell.type.col, seed=42L, max.cells)
  → data.table(from|to) with barcodes; attrs: mode="sampled"

— SCORING —
compute_CellToCell(count.mtx, edge.list, LRM.db, ligand.col, receptor.col,
                   method, n.cores)
  → plain named list: $edges | $lr.meta
  $edges:   data.table — from|to|ligand|receptor|LRM|score (+ x.from|y.from|
            x.to|y.to for spatial edgelists); keyed (from, to, LRM)
  $lr.meta: data.table — ligand|receptor|LRM + extras from LRM.db;
            keyed (ligand, receptor)

compute_PseudoBulk_Connectivity(count.mtx, meta.data, LRM.db, cell.type.col,
                                ligand.col, receptor.col, method,
                                normalize.method, n.cores, verbose)
  → plain data.table (NOT NICHESObject):
    SendingType|ReceivingType|LRM|score|score.norm

— OBJECT CONSTRUCTION —
create_NICHESObject(count.mtx, meta.data, LRM.db,
                    mode="spatial", species="human", sample.col=NULL,
                    cell.type.col=NULL, ligand.col, receptor.col,
                    k, rad, seed=42L, max.cells,
                    smooth=FALSE, smooth.k=9L, n.pcs=30L,
                    method, normalize.method, gene.meta, n.cores, verbose)
  Dispatcher only. Routes to .create_NICHESObject_single() or
  .create_NICHESObject_by_sample() depending on sample.col.
  → NICHESObject
  FILE: R/create_NICHESObject.R (dispatcher only — single function in this file)

— OBJECT EXTENSIONS —
add_neighborhoods(niches.obj, edge.filter.mode, remove.autocrine,
                  n.cores, verbose)
  Calls extract_neighborhoods() on niches.obj$edge.list, populates 3 slots.
  FILE: R/add_neighborhoods.R

extract_neighborhoods(edge.list, edge.filter.mode, remove.autocrine,
                      n.cores, verbose)
  edge.filter.mode: "both_in_neighborhood" (default) |
                    "any_in_neighborhood" | "central_cell_only"
  → list($neighborhood.cell.list, $neighborhood.edge.list, $neighborhood.meta)
  FILE: R/extract_neighborhoods.R

add_segmentation(niches.obj, seg.sf)
  Hard error if any cell missing polygon.

add_cell_meta(niches.obj, new.meta, rebuild.edge.meta=TRUE)
  Joins new metadata onto $cell.meta; calls build_edge_meta() to refresh
  $edge.meta unless suppressed.
  FILE: R/add_cell_meta.R

merge_NICHESObjects(..., verbose=TRUE)
  Accepts: positional args | pre-built list | explicit list() call.
  Hard errors: empty list | non-NICHESObject | mode mismatch | dup barcodes.
  Slot strategy: rbindlist data slots; union gene/LRM lists; sf::rbind segs.
  Self-contained — no internal helper dependency.
  FILE: R/merge_NICHESObjects.R

— DATA INGEST —
extract_NICHESInputs_Seurat(seurat.obj, assay="RNA", layer="counts", verbose)
  V4/V5 detection via inherits(assay.obj,"Assay5"), NOT package version.
  → list(count.mtx, meta.data); user calls create_NICHESObject() next.
  Suggests: SeuratObject

extract_NICHESInputs_AnnData(file, layer="counts", obs.cols=NULL, verbose)
  rhdf5 only. No Python, conda, basilisk, or zellkonverter. R >= 4.2.
  → list(count.mtx, meta.data)
  Suggests: rhdf5

— VISUALIZATION —
CellToCell_InSitu(niches.obj, LRM.list, folder, ..., n.cores)
  Batch PNG; socket cluster; pre-computes all.edges before workers.

CellToCell_InSitu_Single_LRM(niches.obj, LRM, ..., precomputed.bg=NULL)

— AGGREGATION —
Aggregate_CellToCell_Edges(niches.obj, stat)
  → list($sending, $receiving): full cell×LRM grid, zero-filled.


================================================================================
INTERNAL HELPERS (utils.R)
================================================================================

Small utilities only:
  .parse_complex(gene.str)        "gene1_gene2" → c("gene1","gene2")
  .genes_available(gene.str, set) all components in set?
  .prop_norm(x)                   x/sum(x); zeros if sum==0
  .categorical_cols(dt, exclude)  character or factor column names
  .new_NICHESObject(...)          14-slot constructor
  .coerce_meta(meta.data)         → keyed data.table with "cell" col
  .is_NICHESObject(x)             inherits(x,"NICHESObject")
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
  Single LR pair → numeric vector of length n.edges
  method: "product" | "minimum"; supports gene complexes via "_"
  FILE: R/compute_LRM_Score.R
  
  .build_cell_data(count.mtx, normalize.method="prop")
  Sparse count matrix → long-format cell|gene|score|score.norm data.table.
  FILE: R/build_cell_data.R

  .build_edge_meta(edge.list.dt, cell.meta)
  Joins cell metadata onto edge list. Adds Sending.*, Receiving.*, Joint.*
  columns (Joint only for character/factor cols).
  FILE: R/build_edge_meta.R

  .smooth_expression(count.mtx, barcodes, k=9L, n.pcs=30L)
  PCA + kNN smoothing. Replaces target columns with neighborhood-averaged
  raw counts. Note: irlba::irlba(..., fastpath=FALSE) and NO center
  argument — never change.
  FILE: R/smooth_expression.R


================================================================================
PACKAGE INFRASTRUCTURE
================================================================================

Imports: data.table(>=1.14), Matrix(>=1.5), RANN(>=2.6), dbscan(>=1.1),
         ggplot2(>=3.4), irlba(>=2.3.0), scales(:: only), parallel
Suggests: OmnipathR, sf, testthat(>=3.0), withr, SeuratObject, rhdf5

Bundled data (data/):
  ncomms8866_human/mouse/rat/pig.rda — FANTOM5 LR pairs, species-converted
  per Raredon et al. 2019. Loaded via utils::data().

Tests (all passing):
  test-compute_CellToCell.R               test-create_Edgelist_Spatial.R
  test-compute_PseudoBulk_Connectivity.R  test-create_Edgelist_PseudoBulk.R
  test-create_NICHESObject.R              test-create_Edgelist_SamplingByType.R
  test-merge_NICHESObjects.R              test-utilities.R
  test-add_segmentation.R                 test-Aggregate_CellToCell_Edges.R
  test-CellToCell_InSitu.R
  test-extract_NICHESInputs_Seurat.R
  test-extract_NICHESInputs_AnnData.R
  test-compute_LRM_Score.R                (NEW, now use NICHESv2::: to call internal functions)
  test-smooth_expression.R                (NEW, now use NICHESv2::: to call internal functions)

  TODO: test-extract_neighborhoods.R, test-add_neighborhoods.R, test-add_cell_meta.R

Fixtures: inst/testdata/generate_testdata.R — re-run ONLY when schema changes.

Intentional skips (3): full-batch InSitu (>5 LRMs), irlba-absent path,
  replaced placeholder in test-compute_PseudoBulk_Connectivity.R


================================================================================
KEY GOTCHAS (read before touching these areas)
================================================================================

irlba/Matrix:  fastpath=FALSE in smooth_expression() AND no center argument
               — NEVER change. Matrix>=1.6-2 removed as_cholmod_sparse();
               fastpath calls it. Omitting center avoids the dependency
               entirely; kNN distances are unchanged by centering.

data.table composite key join:
  NEVER dt1[dt2, on="col"] when dt2 has composite key.
  Use: idx <- match(dt1$col, dt2$col); dt1[, new := dt2$other[idx]]

socket cluster + asNamespace():
  clusterExport(cl, c("compute_LRM_Score",".parse_complex",".genes_available"),
                envir = asNamespace("NICHESv2"))
  Works under both devtools::load_all() and installed builds. Exporting a
  function from the package does NOT remove the need to clusterExport it
  — workers still need the symbol in their global env.

c() on named lists prepends outer names with dot separator:
  c(list("1" = list(a = 1))) → list("1.a" = 1)
  Fix: unname() the outer list before c(). Hit this in extract_neighborhoods()
  parallel merge path.

sf dependency: keep in Suggests, NOT Imports.

Barcode uniqueness across samples: hard error in both create_NICHESObject()
  and merge_NICHESObjects(); tip: paste0(sample.id,"_",barcode).

Outer sample loop is ALWAYS sequential — inner compute_CellToCell() already
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
  dense HDF5 arrays with reversed dims, so [n_obs, n_vars] HDF5 → [n_vars,
  n_obs] in R, already genes x cells. Sparse matrices must be transposed
  explicitly after reconstruction.

Style B section markers (## Title ---- and # filename ----):
  These are RStudio-foldable. Don't 'fix' them back to banner blocks during future edits.

smooth = TRUE has zero test coverage:
  The smoothing branch in create_NICHESObject() is not exercised by any
  current test. A bug at the call site (line 208 of
  create_NICHESObject_single.R) sat undetected through Pass 1 because of
  this. Add at least one minimal test before relying on this code path.


================================================================================
PENDING DECISIONS
================================================================================

1. Integrate extract_neighborhoods into create_NICHESObject():
   add edge.filter.mode=NULL parameter; when non-NULL, call
   add_neighborhoods() at the end. PENDING advisor discussion.

2. Default edge.filter.mode: currently "both_in_neighborhood"; confirm with
   advisor.

3. Performance optimization for large datasets: pre-extract gene expression
   vectors once rather than per LRM pair; pre-filter silent LRM pairs.

4. Aggregate_CellToCell_Edges future design: sparse dgCMatrix (cells × LRMs)
   plus inverted index $active — pending advisor discussion.

5. PCA kNN smoothing validation PAUSED:
   Advisor hesitant to apply smoothing to actual project data. Function
   remains exported and tested, but no further validation work (the
   four-part validation pipeline previously sketched) until advisor
   decision is resolved. Dataset choice (own Visium vs. public benchmark)
   consequently also paused.

6. Potential package rename from NICHESv2.

7. Pass 2 (CRAN docs): add @keywords internal + @details paragraph to advisor-experimental exports: extract_neighborhoods, create_Edgelist_Spatial, create_Edgelist_PseudoBulk, create_Edgelist_SamplingByType. Each needs a sentence explaining why it's exported despite not being part of the main user workflow.
