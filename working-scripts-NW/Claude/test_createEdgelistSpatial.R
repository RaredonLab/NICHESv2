# ==============================================================================
# FILE:    test_createEdgelistSpatial.R
# PURPOSE: Unit tests + interactive demo for createEdgelistSpatial()
#
# RUN:     source("createEdgelistSpatial.R") first, then source this file.
#          All unit tests emit PASS/FAIL. The demo block at the bottom
#          is meant for you to modify and re-run interactively.
# ==============================================================================

source("createEdgelistSpatial.R")

cat("\n", strrep("=", 70), "\n")
cat("  createEdgelistSpatial() — Test Suite\n")
cat(strrep("=", 70), "\n\n")


# ── Helper: lightweight test runner ────────────────────────────────────────────
.test <- function(label, expr) {
  result <- tryCatch(
    { val <- expr; if (isTRUE(val)) "PASS" else paste0("FAIL (got: ", val, ")") },
    error   = function(e) paste0("ERROR: ", conditionMessage(e)),
    warning = function(w) paste0("WARN:  ", conditionMessage(w))
  )
  status <- if (startsWith(result, "PASS")) "\033[32mPASS\033[0m" else "\033[31mFAIL\033[0m"
  cat(sprintf("  [%s]  %s\n", status, label))
  invisible(result)
}


# ==============================================================================
# TEST GROUP 1 — Input validation
# ==============================================================================
cat("── Group 1: Input validation ─────────────────────────────────────────────\n")

# Basic valid setup reused in many tests
set.seed(42)
.small_coord <- data.frame(x = runif(20), y = runif(20),
                            row.names = paste0("c", seq_len(20)))

.test("Rejects non-data.frame coord", {
  tryCatch({ createEdgelistSpatial(matrix(1:4, 2, 2), k = 2); FALSE },
           error = function(e) TRUE)
})

.test("Rejects coord with < 2 columns", {
  df <- data.frame(x = 1:5); rownames(df) <- paste0("c", 1:5)
  tryCatch({ createEdgelistSpatial(df, k = 2); FALSE }, error = function(e) TRUE)
})

.test("Rejects coord without row names", {
  df <- data.frame(x = 1:5, y = 1:5)    # no rownames → default "1","2",...
  # R sets integer row names by default; we test the empty-string case
  df2 <- df; rownames(df2)[1] <- ""
  tryCatch({ createEdgelistSpatial(df2, k = 2); FALSE }, error = function(e) TRUE)
})

.test("Rejects duplicate row names", {
  df <- data.frame(x = 1:4, y = 1:4, row.names = c("a","a","b","c"))
  tryCatch({ createEdgelistSpatial(df, k = 2); FALSE }, error = function(e) TRUE)
})

.test("Rejects n_cells < 2", {
  df <- data.frame(x = 1, y = 1, row.names = "lone_cell")
  tryCatch({ createEdgelistSpatial(df, k = 1); FALSE }, error = function(e) TRUE)
})

.test("Rejects invalid n_cores", {
  tryCatch({ createEdgelistSpatial(.small_coord, k = 3, n_cores = -1); FALSE },
           error = function(e) TRUE)
})


# ==============================================================================
# TEST GROUP 2 — Parameter mode resolution
# ==============================================================================
cat("\n── Group 2: Parameter mode resolution ───────────────────────────────────\n")

.test("Both NULL → default k=4, mode=kNN", {
  el <- suppressMessages(createEdgelistSpatial(.small_coord, k = NULL, rad = NULL))
  attr(el, "mode") == "kNN" && attr(el, "param") == 4L
})

.test("Both valid → radius mode (k overwritten)", {
  el <- suppressMessages(createEdgelistSpatial(.small_coord, k = 3, rad = 0.3))
  attr(el, "mode") == "radius"
})

.test("Only k valid → kNN mode", {
  el <- suppressMessages(createEdgelistSpatial(.small_coord, k = 4, rad = NULL))
  attr(el, "mode") == "kNN" && attr(el, "param") == 4L
})

.test("Only rad valid → radius mode", {
  el <- suppressMessages(createEdgelistSpatial(.small_coord, k = NULL, rad = 0.2))
  attr(el, "mode") == "radius" && attr(el, "param") == 0.2
})

.test("Negative k treated as invalid → default k=4", {
  el <- suppressMessages(createEdgelistSpatial(.small_coord, k = -2, rad = NULL))
  attr(el, "mode") == "kNN" && attr(el, "param") == 4L
})

.test("Non-integer k (e.g. 2.7) treated as invalid → default k=4", {
  el <- suppressMessages(createEdgelistSpatial(.small_coord, k = 2.7, rad = NULL))
  attr(el, "mode") == "kNN" && attr(el, "param") == 4L
})

.test("Negative rad treated as invalid → kNN mode when k is valid", {
  el <- suppressMessages(createEdgelistSpatial(.small_coord, k = 3, rad = -1))
  attr(el, "mode") == "kNN" && attr(el, "param") == 3L
})


# ==============================================================================
# TEST GROUP 3 — kNN mode: correctness on a known small example
#
# Layout (1D on the x-axis for clarity):
#
#   A=0   B=1   E=2.6   C=5   D=6
#
# With k=1:
#   A's 1NN → B  (dist 1.0);   B's 1NN → A  (dist 1.0)  → mutual ✓  keep A↔B
#   C's 1NN → D  (dist 1.0);   D's 1NN → C  (dist 1.0)  → mutual ✓  keep C↔D
#   E's 1NN → B  (dist 1.6);   B's 1NN → A  (not E)     → NOT mutual ✗ discard E↔B
#
# Expected mutual edges: (A→B, B→A) and (C→D, D→C)  →  4 directed rows
# ==============================================================================
cat("\n── Group 3: kNN correctness (known geometry) ─────────────────────────────\n")

.known_coord <- data.frame(
  x = c(0, 1, 5, 6, 2.6),
  y = c(0, 0, 0, 0, 0),
  row.names = c("A", "B", "C", "D", "E")
)

.el_k1 <- suppressMessages(createEdgelistSpatial(.known_coord, k = 1))

.test("k=1: result is a data.table", is.data.table(.el_k1))

.test("k=1: exactly 4 mutual directed inter-cell edges + 5 autocrine = 9 total", nrow(.el_k1) == 9L)
.test("k=1: A→B edge present",   nrow(.el_k1[from == "A" & to == "B"]) == 1L)
.test("k=1: B→A edge present",   nrow(.el_k1[from == "B" & to == "A"]) == 1L)
.test("k=1: C→D edge present",   nrow(.el_k1[from == "C" & to == "D"]) == 1L)
.test("k=1: D→C edge present",   nrow(.el_k1[from == "D" & to == "C"]) == 1L)

.test("k=1: non-mutual edge E→B NOT present", nrow(.el_k1[from == "E" & to == "B"]) == 0L)
.test("k=1: non-mutual edge B→E NOT present", nrow(.el_k1[from == "B" & to == "E"]) == 0L)

.test("k=1: weight for A→B ≈ 1.0", {
  w <- .el_k1[from == "A" & to == "B", weight]
  abs(w - 1.0) < 1e-9
})

# Autocrine edge checks
.test("k=1: A→A autocrine edge present with weight=0", {
  r <- .el_k1[from == "A" & to == "A"]
  nrow(r) == 1L && r$weight == 0
})
.test("k=1: every cell has exactly one autocrine edge", {
  self_edges <- .el_k1[from == to]
  nrow(self_edges) == 5L && all(self_edges$weight == 0)
})

# Larger k: with k=2, E should now form mutual edges too
.el_k2 <- suppressMessages(createEdgelistSpatial(.known_coord, k = 2))
.test("k=2: E now participates in mutual edges", nrow(.el_k2[from == "E" | to == "E"]) > 0L)


# ==============================================================================
# TEST GROUP 4 — Radius mode: correctness on the same known geometry
# ==============================================================================
cat("\n── Group 4: Radius mode correctness (known geometry) ────────────────────\n")

# rad=1.5: A(0)↔B(1) [dist 1], C(5)↔D(6) [dist 1], E(2.6)↔B(1) [dist 1.6 > 1.5 ✗]
.el_r1 <- suppressMessages(createEdgelistSpatial(.known_coord, rad = 1.0))

.test("rad=1.0: A→B present",   nrow(.el_r1[from == "A" & to == "B"]) == 1L)
.test("rad=1.0: B→A present",   nrow(.el_r1[from == "B" & to == "A"]) == 1L)
.test("rad=1.0: C→D present",   nrow(.el_r1[from == "C" & to == "D"]) == 1L)
.test("rad=1.0: no E inter-cell edges", {
  nrow(.el_r1[from == "E" & to != "E"]) == 0L &&
  nrow(.el_r1[to == "E"   & from != "E"]) == 0L
})
.test("rad=1.0: E→E autocrine edge present", {
  r <- .el_r1[from == "E" & to == "E"]
  nrow(r) == 1L && r$weight == 0
})
.test("rad=1.0: every cell has exactly one autocrine edge", {
  self_edges <- .el_r1[from == to]
  nrow(self_edges) == 5L && all(self_edges$weight == 0)
})

# rad=2.0: now E↔B should appear (dist=1.6) and also E↔A (dist=2.6 > 2 ✗)
.el_r2 <- suppressMessages(createEdgelistSpatial(.known_coord, rad = 2.0))
.test("rad=2.0: E→B present",   nrow(.el_r2[from == "E" & to == "B"]) == 1L)
.test("rad=2.0: B→E present",   nrow(.el_r2[from == "B" & to == "E"]) == 1L)
.test("rad=2.0: E→A NOT present (dist=2.6)", nrow(.el_r2[from == "E" & to == "A"]) == 0L)

# rad very small → empty edge list (with warning)
.test("rad near-zero → empty edge list warning", {
  tryCatch({
    withCallingHandlers(
      suppressMessages(createEdgelistSpatial(.known_coord, rad = 1e-9)),
      warning = function(w) {
        if (grepl("empty", conditionMessage(w), ignore.case = TRUE))
          invokeRestart("muffleWarning")
      }
    )
    TRUE  # reached here = warning was raised (or edge list is empty)
  }, error = function(e) FALSE)
})


# ==============================================================================
# TEST GROUP 5 — Output structure & metadata
# ==============================================================================
cat("\n── Group 5: Output structure & metadata ──────────────────────────────────\n")

.el <- suppressMessages(createEdgelistSpatial(.small_coord, k = 3))

.test("Output is a data.table", is.data.table(.el))

.test("Columns are always (from, to, weight)", {
  identical(names(.el), c("from", "to", "weight"))
})

.test("Column 'from' is character",  is.character(.el$from))
.test("Column 'to'   is character",  is.character(.el$to))
.test("Column 'weight' is numeric",  is.numeric(.el$weight))

.test("All 'from' values are valid cell names", all(.el$from %in% rownames(.small_coord)))
.test("All 'to'   values are valid cell names", all(.el$to   %in% rownames(.small_coord)))

.test("Autocrine edges present for all 20 cells", {
  self <- .el[from == to]
  nrow(self) == 20L && all(self$weight == 0)
})

.test("All inter-cell weights are strictly positive", {
  all(.el[from != to, weight] > 0)
})

.test("data.table key is set on (from, to)", {
  identical(key(.el), c("from", "to"))
})

.test("Metadata attr 'mode' is 'kNN'",    attr(.el, "mode")    == "kNN")
.test("Metadata attr 'param' is 3L",      attr(.el, "param")   == 3L)
.test("Metadata attr 'n_cells' is 20L",   attr(.el, "n_cells") == 20L)
.test("Metadata attr 'cell_names' correct", length(attr(.el, "cell_names")) == 20)


# ==============================================================================
# TEST GROUP 6 — Symmetry & mutual edge property
# ==============================================================================
cat("\n── Group 6: Mutual edge symmetry ─────────────────────────────────────────\n")

.test("kNN: every inter-cell edge (A→B) has a matching reverse (B→A)", {
  el       <- suppressMessages(createEdgelistSpatial(.small_coord, k = 4))
  inter    <- el[from != to]           # exclude autocrine edges
  fwd      <- inter[, .(from, to)]
  rev      <- inter[, .(from = to, to = from)]
  nrow(fintersect(fwd, rev)) == nrow(fwd)
})

.test("Radius: inherently symmetric inter-cell edges", {
  el <- suppressMessages(createEdgelistSpatial(.small_coord, rad = 0.3))
  inter <- el[from != to]
  if (nrow(inter) == 0L) return(TRUE)
  fwd <- inter[, .(from, to)]
  rev <- inter[, .(from = to, to = from)]
  nrow(fintersect(fwd, rev)) == nrow(fwd)
})


# ==============================================================================
# TEST GROUP 7 — Scalability smoke test
# ==============================================================================
cat("\n── Group 7: Scalability smoke test ───────────────────────────────────────\n")

set.seed(99)
.n_large <- 2000
.big_coord <- data.frame(
  x = runif(.n_large),
  y = runif(.n_large),
  row.names = paste0("cell_", seq_len(.n_large))
)

.t_knn <- system.time(
  .el_big_knn <- suppressMessages(createEdgelistSpatial(.big_coord, k = 10))
)
cat(sprintf("  kNN  (n=%d, k=10):   %.2f sec | %d edges\n",
            .n_large, .t_knn["elapsed"], nrow(.el_big_knn)))

.t_rad <- system.time(
  .el_big_rad <- suppressMessages(createEdgelistSpatial(.big_coord, rad = 0.05))
)
cat(sprintf("  Rad  (n=%d, r=0.05): %.2f sec | %d edges\n",
            .n_large, .t_rad["elapsed"], nrow(.el_big_rad)))

.test("Large kNN run completes in < 30s",  .t_knn["elapsed"] < 30)
.test("Large rad run completes in < 30s",  .t_rad["elapsed"] < 30)
.test("Large kNN: n_cells attribute correct", attr(.el_big_knn, "n_cells") == .n_large)


# ==============================================================================
# ██████████████████████████████████████████████████████████████████████████████
#  INTERACTIVE DEMO BLOCK
#  ─────────────────────────────────────────────────────────────────────────────
#  Edit the variables below and source this section to try the function
#  on your own data or custom synthetic examples.
# ██████████████████████████████████████████████████████████████████████████████
# ==============================================================================

cat("\n", strrep("=", 70), "\n")
cat("  INTERACTIVE DEMO — edit variables below and re-run\n")
cat(strrep("=", 70), "\n\n")

# ── [A] Define your input ──────────────────────────────────────────────────────

# Option 1: Synthetic random coordinates (change n_demo to your dataset size)
set.seed(123)
n_demo <- 500
demo_coord <- data.frame(
  x = runif(n_demo, min = 0, max = 100),
  y = runif(n_demo, min = 0, max = 100),
  row.names = paste0("Cell_", seq_len(n_demo))
)

# Option 2: Load your own data (uncomment and edit path)
# demo_coord <- read.csv("my_coords.csv", row.names = 1)

orig.data <- readRDS("~/Library/CloudStorage/OneDrive-YaleUniversity/Lung_Tissue_CK/Osteosarcoma_Lung_Tissue/OS.merged.rds")

sub.orig <- subset(orig.data, subset = Sample == 'OS7')

NICHES.CTC.rad <- readRDS("~/Library/CloudStorage/OneDrive-YaleUniversity/Lung_Tissue_CK/Osteosarcoma_Lung_Tissue/OS.CTC.rad10.rds")

load("~/Library/CloudStorage/OneDrive-YaleUniversity/Lung_Tissue_CK/Osteosarcoma_Lung_Tissue/Old/OS.CTC.Robj")

sub.rad <- subset(NICHES.CTC.rad, subset = Sample.Sending == 'OS7')

sub.kNN <- subset(OS.CTC, subset = Sample.Sending == 'OS7')

demo_coord <- data.frame(x = sub.orig$x,
                         y = sub.orig$y)

rownames(demo_coord) <- colnames(sub.orig)

# ── [B] Set parameters ────────────────────────────────────────────────────────

demo_k        <- 7      # Number of nearest neighbors  (set to NULL to disable)
demo_rad      <- 10    # Neighborhood radius           (set to NULL to disable)
demo_n_cores  <- 1L      # Cores (use >1 on Unix/macOS for large datasets)

# ── [C] Run the function ──────────────────────────────────────────────────────

demo_rad <- createEdgelistSpatial(
  coord   = demo_coord,
  k       = NULL,
  rad     = demo_rad,
  n_cores = demo_n_cores
)

demo_kNN <- createEdgelistSpatial(
  coord   = demo_coord,
  k       = demo_k,
  rad     = NULL,
  n_cores = demo_n_cores
)

# ── [D] Inspect results ───────────────────────────────────────────────────────

cat("\n--- rad Edge list preview (first 10 rows) ---\n")
print(demo_rad[1:min(10, nrow(demo_rad))])

cat("\n--- Summary statistics ---\n")
cat(sprintf("  Mode              : %s\n",   attr(demo_rad, "mode")))
cat(sprintf("  Parameter         : %s\n",   attr(demo_rad, "param")))
cat(sprintf("  Total cells       : %d\n",   attr(demo_rad, "n_cells")))
cat(sprintf("  Total edges       : %d\n",   nrow(demo_rad)))
cat(sprintf("  Autocrine edges   : %d\n",   nrow(demo_rad[from == to])))
cat(sprintf("  Inter-cell edges  : %d\n",   nrow(demo_rad[from != to])))
cat(sprintf("  Avg neighbors/cell: %.2f\n", nrow(demo_rad[from != to]) / attr(demo_rad, "n_cells")))
cat(sprintf("  Weight range (inter-cell): [%.4f, %.4f]\n",
            min(demo_rad[from != to, weight]), max(demo_rad[from != to, weight])))



cat("\n--- kNN Edge list preview (first 10 rows) ---\n")
print(demo_kNN[1:min(10, nrow(demo_kNN))])

cat("\n--- Summary statistics ---\n")
cat(sprintf("  Mode              : %s\n",   attr(demo_kNN, "mode")))
cat(sprintf("  Parameter         : %s\n",   attr(demo_kNN, "param")))
cat(sprintf("  Total cells       : %d\n",   attr(demo_kNN, "n_cells")))
cat(sprintf("  Total edges       : %d\n",   nrow(demo_kNN)))
cat(sprintf("  Autocrine edges   : %d\n",   nrow(demo_kNN[from == to])))
cat(sprintf("  Inter-cell edges  : %d\n",   nrow(demo_kNN[from != to])))
cat(sprintf("  Avg neighbors/cell: %.2f\n", nrow(demo_kNN[from != to]) / attr(demo_kNN, "n_cells")))
cat(sprintf("  Weight range (inter-cell): [%.4f, %.4f]\n",
            min(demo_kNN[from != to, weight]), max(demo_kNN[from != to, weight])))

# ── [E] Fast lookup: all neighbors of a specific cell ─────────────────────────
# (data.table key makes this O(log n))

query_cell <- rownames(demo_coord)[1]    # change to any cell name you like
neighbors_of_query <- demo_rad[from == query_cell]

cat(sprintf("\n--- Neighbors of '%s' ---\n", query_cell))
print(neighbors_of_query)

cat("\n--- All tests and demo complete ---\n\n")
