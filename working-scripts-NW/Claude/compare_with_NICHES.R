#' Benchmark: Compare your edge-LR matrix against NICHES output
#'
#' NICHES format:  rows = "Ligand—Receptor" mechanisms, cols = "CellA—CellB" edges
#' Your format:    rows = "CellA->CellB" edges,         cols = "Ligand_Receptor" pairs

compare_with_niches <- function(niches.out,
                                your.out,
                                niches_lr_sep   = "--",   # separator used in NICHES LR names
                                niches_edge_sep = "--",   # separator used in NICHES edge names
                                your_lr_sep     = "_",    # separator used in your LR names
                                your_edge_sep   = "->",   # separator used in your edge names
                                tol             = 1e-6) { 
  library(Matrix)
  
  # ── 1. Align orientation: transpose NICHES to edges × LR ──────────────────
  niches.mat <- t(niches.out)   # now: rows = edges, cols = LR pairs
  
  message("NICHES matrix (after transpose): ", nrow(niches.mat), " edges × ", ncol(niches.mat), " LR pairs")
  message("Your matrix:                     ", nrow(your.out),   " edges × ", ncol(your.out),   " LR pairs")
  
  # ── 2. Normalize LR pair names to a canonical "GeneA|GeneB" format ─────────
  normalize_lr <- function(names_vec, sep) {
    parts <- strsplit(names_vec, sep, fixed = TRUE)
    sapply(parts, function(x) paste(trimws(x[1]), trimws(x[2]), sep = "|"))
  }
  
  normalize_edge <- function(names_vec, sep) {
    parts <- strsplit(names_vec, sep, fixed = TRUE)
    sapply(parts, function(x) paste(trimws(x[1]), trimws(x[2]), sep = "|"))
  }
  
  niches_lr_canon   <- normalize_lr(colnames(niches.mat), niches_lr_sep)
  your_lr_canon     <- normalize_lr(colnames(your.out),   your_lr_sep)
  niches_edge_canon <- normalize_edge(rownames(niches.mat), niches_edge_sep)
  your_edge_canon   <- normalize_edge(rownames(your.out),   your_edge_sep)
  
  # ── 3. Find overlapping LR pairs and edges ─────────────────────────────────
  common_lr    <- intersect(niches_lr_canon,   your_lr_canon)
  common_edges <- intersect(niches_edge_canon, your_edge_canon)
  
  message(sprintf("\nOverlap summary:"))
  message(sprintf("  LR pairs  — NICHES: %d | Yours: %d | Common: %d",
                  ncol(niches.mat), ncol(your.out), length(common_lr)))
  message(sprintf("  Edges     — NICHES: %d | Yours: %d | Common: %d",
                  nrow(niches.mat), nrow(your.out), length(common_edges)))
  
  if (length(common_lr) == 0 || length(common_edges) == 0) {
    stop("No overlapping LR pairs or edges found. Check separator arguments.")
  }
  
  # ── 4. Subset both matrices to common rows/cols, reordering to match ────────
  niches_sub <- niches.mat[match(common_edges, niches_edge_canon),
                           match(common_lr,    niches_lr_canon),
                           drop = FALSE]
  
  your_sub   <- your.out[match(common_edges, your_edge_canon),
                         match(common_lr,   your_lr_canon),
                         drop = FALSE]
  
  # ── 5. Numerical comparison ────────────────────────────────────────────────
  # Convert to dense for element-wise diff (subsetted matrix should be manageable)
  diff.mat   <- as.matrix(niches_sub) - as.matrix(your_sub)
  abs_diff   <- abs(diff.mat)
  
  max_diff   <- max(abs_diff)
  mean_diff  <- mean(abs_diff)
  n_mismatch <- sum(abs_diff > tol)
  n_total    <- length(abs_diff)
  
  match_rate <- (n_total - n_mismatch) / n_total * 100
  
  message(sprintf("\n── Numerical Comparison (tol = %g) ──", tol))
  message(sprintf("  Max absolute difference:    %.6e", max_diff))
  message(sprintf("  Mean absolute difference:   %.6e", mean_diff))
  message(sprintf("  Mismatched entries:         %d / %d (%.2f%% match)", 
                  n_mismatch, n_total, match_rate))
  
  # ── 6. Spot-check: show top mismatches for debugging ──────────────────────
  if (n_mismatch > 0) {
    message("\n── Top 10 Largest Mismatches ──")
    mismatch_idx <- which(abs_diff > tol, arr.ind = TRUE)
    mismatch_df  <- data.frame(
      edge      = common_edges[mismatch_idx[, 1]],
      lr_pair   = common_lr[mismatch_idx[, 2]],
      niches_val = as.matrix(niches_sub)[mismatch_idx],
      your_val   = as.matrix(your_sub)[mismatch_idx],
      abs_diff   = abs_diff[mismatch_idx]
    )
    mismatch_df <- mismatch_df[order(-mismatch_df$abs_diff), ]
    print(head(mismatch_df, 10))
  } else {
    message("  ✓ All overlapping entries match within tolerance!")
  }
  
  # ── 7. Coverage report: what's in NICHES but not yours, and vice versa ─────
  message("\n── Coverage Gaps ──")
  only_in_niches_lr <- setdiff(niches_lr_canon, your_lr_canon)
  only_in_your_lr   <- setdiff(your_lr_canon, niches_lr_canon)
  only_in_niches_e  <- setdiff(niches_edge_canon, your_edge_canon)
  only_in_your_e    <- setdiff(your_edge_canon, niches_edge_canon)
  
  message(sprintf("  LR pairs only in NICHES: %d", length(only_in_niches_lr)))
  message(sprintf("  LR pairs only in yours:  %d", length(only_in_your_lr)))
  message(sprintf("  Edges only in NICHES:    %d", length(only_in_niches_e)))
  message(sprintf("  Edges only in yours:     %d", length(only_in_your_e)))
  
  if (length(only_in_niches_lr) > 0) {
    message("  Sample LR only in NICHES: ", 
            paste(head(only_in_niches_lr, 5), collapse = ", "))
  }
  if (length(only_in_niches_e) > 0) {
    message("  Sample edges only in NICHES: ", 
            paste(head(only_in_niches_e, 5), collapse = ", "))
  }
  
  # ── 8. Return structured result for further inspection ────────────────────
  invisible(list(
    niches_sub   = niches_sub,
    your_sub     = your_sub,
    diff_mat     = diff.mat,
    common_lr    = common_lr,
    common_edges = common_edges,
    match_rate   = match_rate,
    max_diff     = max_diff
  ))
}

# ==============================================================================
# Source functions
# ==============================================================================

source("createEdgelistSpatial.R")
source("compute_edge_LR_matrix.R")

# ==============================================================================
# Data
# ==============================================================================

orig.data <- readRDS("~/Library/CloudStorage/OneDrive-YaleUniversity/Lung_Tissue_CK/Osteosarcoma_Lung_Tissue/OS.merged.rds")

sub.orig <- subset(orig.data, subset = Sample == 'OS7')

NICHES.CTC.rad <- readRDS("~/Library/CloudStorage/OneDrive-YaleUniversity/Lung_Tissue_CK/Osteosarcoma_Lung_Tissue/OS.CTC.rad10.rds")

load("~/Library/CloudStorage/OneDrive-YaleUniversity/Lung_Tissue_CK/Osteosarcoma_Lung_Tissue/Old/OS.CTC.Robj")

sub.rad <- subset(NICHES.CTC.rad, subset = Sample.Sending == 'OS7')

sub.kNN <- subset(OS.CTC, subset = Sample.Sending == 'OS7')

# create the edge list first
demo_coord <- data.frame(x = sub.orig$x,
                         y = sub.orig$y)

rownames(demo_coord) <- colnames(sub.orig)

demo_rad <- createEdgelistSpatial(
  coord   = demo_coord,
  k       = NULL,
  rad     = 10,           # fixed: was stored in a variable that was never defined
  n_cores = 1L            # fixed: demo_n_cores was never defined
)

# Load Ligand-Receptor database
LRM.list <- read.csv("~/Library/CloudStorage/OneDrive-YaleUniversity/Lung_Tissue_CK/Connectome2025.LRMs.human.csv")

our_result <- compute_edge_LR_matrix(
    count.mtx    = sub.orig@assays$SCT@data,
    edge.list    = demo_rad,
    LRM.list     = LRM.list,
    ligand.col   = 'Ligand',
    receptor.col = 'Receptor',
    method       = "product",
    n_cores      = 1L)

# Pass raw NICHES matrix directly — the function handles transpose and
# em-dash separator conversion internally
NICHES.mtx <- sub.rad@assays$CellToCellSpatial$counts

comp <- compare_with_niches(NICHES.mtx,
                            our_result,
                            niches_lr_sep = '—',
                            niches_edge_sep = '—')
