# ==============================================================================
# FILE:    demo_computeCellToCellSpatial.R
# PURPOSE: Demo script for computeCellToCellSpatial() — shows how to run the
#          function, enrich the output with cell metadata, and perform common
#          aggregation queries on the EdgeLRResult object.
# ==============================================================================

source("createEdgelistSpatial.R")
source("computeCellToCellSpatial.R")

# ==============================================================================
# SECTION 1 ── Load data
# ==============================================================================

orig.data <- readRDS("~/Library/CloudStorage/OneDrive-YaleUniversity/Lung_Tissue_CK/Osteosarcoma_Lung_Tissue/OS.merged.rds")
sub.orig  <- subset(orig.data, subset = Sample == 'OS7')

LRM.list  <- read.csv("~/Library/CloudStorage/OneDrive-YaleUniversity/Lung_Tissue_CK/Connectome2025.LRMs.human.csv")

# ==============================================================================
# SECTION 2 ── Build spatial edge list
# ==============================================================================

demo_coord <- data.frame(x = sub.orig$x,
                         y = sub.orig$y)
rownames(demo_coord) <- colnames(sub.orig)

demo_rad <- createEdgelistSpatial(
  coord   = demo_coord,
  k       = NULL,
  rad     = 10,
  n_cores = 1L
)

# ==============================================================================
# SECTION 3 ── Compute cell-to-cell signals
# ==============================================================================

ctc <- computeCellToCellSpatial(
  count.mtx    = sub.orig@assays$SCT@data,
  edge.list    = demo_rad,
  LRM.list     = LRM.list,
  ligand.col   = 'Ligand',
  receptor.col = 'Receptor',
  method       = "product",
  n_cores      = 1L
)

# Quick overview of the result
print(ctc)

# Peek at each slot
head(ctc$signals)
head(ctc$cell_meta)
head(ctc$lr_meta)
ctc$params

# ==============================================================================
# SECTION 4 ── Enrich cell_meta with Seurat metadata
#
# cell_meta starts with only a 'cell' column. Join any per-cell information
# from your Seurat object here. All aggregation queries below then work
# automatically without touching $signals.
# ==============================================================================

seurat_meta <- as.data.table(sub.orig@meta.data, keep.rownames = "cell")

# Add whichever columns are relevant — adjust to your actual metadata columns
cols_to_add <- c("cell", "cell_type", "cell_class", "Sample", "seurat_clusters")
cols_to_add <- intersect(cols_to_add, colnames(seurat_meta))   # safe subset

ctc$cell_meta <- seurat_meta[, ..cols_to_add][ctc$cell_meta, on = "cell"]

# Confirm enrichment
head(ctc$cell_meta)

# ==============================================================================
# SECTION 5 ── Aggregation examples
# ==============================================================================

# test summarizing function summarizeCellSignals.R 
source('summarizeCellSignals.R')

result <- summarizeCellSignals(ctc)

# ── 5a. Total signal sent by each individual cell ─────────────────────────────
total_sent <- ctc$signals[, .(total_sent = sum(signal)), by = from]
setorder(total_sent, -total_sent)
head(total_sent, 10)

# ── 5b. Total signal received by each individual cell ─────────────────────────
total_received <- ctc$signals[, .(total_received = sum(signal)), by = to]
setorder(total_received, -total_received)
head(total_received, 10)

# ── 5c. Total signal sent, grouped by sending cell type ───────────────────────
#   Requires cell_type in cell_meta (Section 4 above)
if ("cell_type" %in% colnames(ctc$cell_meta)) {
  sent_by_type <- ctc$signals[
    ctc$cell_meta, on = .(from = cell), nomatch = 0L
  ][, .(total_sent = sum(signal), mean_sent = mean(signal), n_edges = .N),
    by = cell_type]
  setorder(sent_by_type, -total_sent)
  print(sent_by_type)
}

# ── 5d. Signal flow between cell type pairs ───────────────────────────────────
#   Join cell_meta twice: once for sending, once for receiving
if ("cell_type" %in% colnames(ctc$cell_meta)) {
  flow_mat <- ctc$signals[
    ctc$cell_meta, on = .(from = cell), nomatch = 0L
  ][
    ctc$cell_meta, on = .(to = cell), nomatch = 0L
  ][, .(total_signal = sum(signal)),
    by = .(sending_type = cell_type, receiving_type = i.cell_type)]
  setorder(flow_mat, -total_signal)
  print(head(flow_mat, 20))
}

# ── 5e. Top LR pairs by total signal across all edges ─────────────────────────
top_lr <- ctc$signals[
  , .(total_signal = sum(signal), n_edges = .N),
  by = .(ligand, receptor)
]
setorder(top_lr, -total_signal)
head(top_lr, 20)

# ── 5f. LR pairs enriched in a specific sending cell type ─────────────────────
if ("cell_type" %in% colnames(ctc$cell_meta)) {
  query_type <- ctc$cell_meta$cell_type[1]   # replace with your cell type of interest

  lr_from_type <- ctc$signals[
    ctc$cell_meta, on = .(from = cell), nomatch = 0L
  ][cell_type == query_type,
    .(total_signal = sum(signal), n_edges = .N),
    by = .(ligand, receptor)]
  setorder(lr_from_type, -total_signal)
  head(lr_from_type, 20)
}

# ── 5g. Aggregate by LR metadata (e.g., pathway class) ────────────────────────
#   Works if LRM.list has an extra column like 'Pathway' or 'Family'
lr_meta_cols <- colnames(ctc$lr_meta)
message("Available lr_meta columns: ", paste(lr_meta_cols, collapse = ", "))

# Example — replace 'Pathway' with whatever grouping column exists in your DB
if ("Pathway" %in% lr_meta_cols) {
  pathway_signal <- ctc$signals[
    ctc$lr_meta, on = c("ligand", "receptor"), nomatch = 0L
  ][, .(total_signal = sum(signal)), by = Pathway]
  setorder(pathway_signal, -total_signal)
  print(head(pathway_signal, 20))
}

# ==============================================================================
# SECTION 6 ── Reconstruct sparse matrix if needed for downstream tools
#
# The dgCMatrix is NOT stored in ctc — call as_matrix() on demand.
# This keeps memory footprint small when you don't need it.
# ==============================================================================

# Reconstruct only when a downstream tool requires matrix input
ctc_mat <- as_matrix(ctc)
dim(ctc_mat)        # n_edges x n_lr
ctc_mat[1:5, 1:5]
