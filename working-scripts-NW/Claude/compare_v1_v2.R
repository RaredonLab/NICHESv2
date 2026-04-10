## compare_v1_v2.R
## Quick check: does NICHESv2 compute_CellToCell() match NICHESv1 output?
## ------------------------------------------------------------------

library(Seurat)
library(data.table)

new.PPLR.niches <- readRDS("/Users/nw426/Documents/GitHub/NICHESv2-NW/new.PPLR.niches.rds")
orig.CTC        <- readRDS("/Volumes/T7/Pneum/Xenium/ALLTMA_EDGE_RAD.rds")

DefaultAssay(orig.CTC) <- 'Edge_Raw'

# ── 1. Inspect v1 object (run this first, then comment out) ──────────
cat("=== V1 object class ===\n")
print(class(orig.CTC))

if (inherits(orig.CTC, "Seurat")) {
  cat("\nAssays:", Assays(orig.CTC), "\n")
  cat("Default assay:", DefaultAssay(orig.CTC), "\n")
  cat("Dims (LRMs x edges):", dim(orig.CTC), "\n")
  cat("First 5 row names (LRMs):\n")
  print(head(rownames(orig.CTC), 5))
  cat("First 3 col names (edges):\n")
  print(head(colnames(orig.CTC), 3))
  cat("\nSlots in default assay:\n")
  print(slotNames(orig.CTC[[DefaultAssay(orig.CTC)]]))
} else {
  cat("Not Seurat — printing str() for manual inspection:\n")
  str(orig.CTC, max.level = 2)
}

# ── 2. Extract v1 scores ─────────────────────────────────────────────
# Adjust the slot/layer if needed after inspecting above.
# Original NICHES typically stores raw scores in @data or @counts.
v1.mtx <- orig.CTC@assays$Edge_Raw$counts   # try "counts" if empty

# v1 uses "—" (em-dash) as separator for both edges and LRMs;
# v2 uses "|".  Detect separator automatically:
v1.edge.sep <- if (any(grepl("—", colnames(v1.mtx)))) "—" else
               if (any(grepl("\\|", colnames(v1.mtx)))) "|" else
               stop("Cannot detect edge separator in v1 colnames")

v1.lrm.sep  <- if (any(grepl("—", rownames(v1.mtx)))) "—" else
               if (any(grepl("\\|", rownames(v1.mtx)))) "|" else
               stop("Cannot detect LRM separator in v1 rownames")

cat("\nV1 separators — edge:", v1.edge.sep, " LRM:", v1.lrm.sep, "\n")

# Convert v1 sparse matrix → long data.table
v1.triplet <- summary(as(v1.mtx, "dgCMatrix"))   # data.frame with i, j, x (1-indexed)
v1.dt <- data.table(
  edge.v1  = colnames(v1.mtx)[v1.triplet$j],
  LRM.v1   = rownames(v1.mtx)[v1.triplet$i],
  score.v1 = v1.triplet$x
)

# Normalise separators to "|" so we can join with v2
if (v1.edge.sep != "|") {
  v1.dt[, edge := gsub(v1.edge.sep, "|", edge.v1, fixed = TRUE)]
} else {
  v1.dt[, edge := edge.v1]
}
if (v1.lrm.sep != "|") {
  v1.dt[, LRM := gsub(v1.lrm.sep, "|", LRM.v1, fixed = TRUE)]
} else {
  v1.dt[, LRM := LRM.v1]
}

# ── 3. Extract v2 scores ─────────────────────────────────────────────
v2.dt <- new.PPLR.niches$edge.data[, .(edge, LRM, score.v2 = score)]

# ── 4. Join & compare ────────────────────────────────────────────────
setkey(v1.dt, edge, LRM)
setkey(v2.dt, edge, LRM)

merged <- merge(v1.dt[, .(edge, LRM, score.v1)],
                v2.dt,
                by = c("edge", "LRM"),
                all = TRUE)

n.both    <- sum(!is.na(merged$score.v1) & !is.na(merged$score.v2))
n.v1.only <- sum(!is.na(merged$score.v1) &  is.na(merged$score.v2))
n.v2.only <- sum( is.na(merged$score.v1) & !is.na(merged$score.v2))

cat("\n=== Coverage ===\n")
cat("Pairs in both:    ", n.both, "\n")
cat("V1-only pairs:    ", n.v1.only, "\n")
cat("V2-only pairs:    ", n.v2.only, "\n")

# Numerical comparison on shared pairs
shared <- merged[!is.na(score.v1) & !is.na(score.v2)]
shared[, diff := score.v2 - score.v1]

cat("\n=== Score differences (v2 − v1) ===\n")
print(summary(shared$diff))

cat("\nExact matches (diff == 0): ",
    sum(shared$diff == 0), " / ", nrow(shared),
    " (", round(100 * mean(shared$diff == 0), 2), "%)\n", sep = "")

cat("Near matches (|diff| < 1e-10): ",
    sum(abs(shared$diff) < 1e-10), " / ", nrow(shared),
    " (", round(100 * mean(abs(shared$diff) < 1e-10), 2), "%)\n", sep = "")

# Show worst mismatches if any
if (any(abs(shared$diff) >= 1e-10)) {
  cat("\nTop 10 largest mismatches:\n")
  print(shared[order(-abs(diff))][1:min(10, .N)])
}

cat("\nDone.\n")
