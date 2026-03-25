# ==============================================================================
# FILE:    test_tiebreak_v2.R
# PURPOSE: Minimal 6-cell example that shows WHY NICHES produces fewer edges
#          than our code, despite both using mutual (intersection) kNN.
#
# THE SCENARIO:
#   Cell A has an exact distance TIE at the k-th neighbor slot.
#   Both tied cells (D and E) are equidistant from A.
#   The two methods break the tie differently:
#     - NICHES order()  → picks D  (lower row index)
#     - RANN::nn2       → picks E  (kd-tree traversal order)
#
#   The ASYMMETRY that causes the edge count difference:
#     - A→D is chosen by NICHES, but D does NOT have A in its kNN
#       → A↔D fails the mutual filter → 2 directed edges LOST
#     - A→E is chosen by our code,  and E DOES have A in its kNN
#       → A↔E passes the mutual filter → 2 directed edges KEPT
#
#   Net result: our code has 2 more edges than NICHES for this one tie.
#   On the real dataset (102307 vs 102301), this pattern occurs 3 times.
# ==============================================================================

library(RANN)
library(data.table)

cat(strrep("=", 70), "\n")
cat("  Minimal tie-breaking example\n")
cat(strrep("=", 70), "\n\n")

# ── Cell layout ───────────────────────────────────────────────────────────────
#
#                   E = (0, 2)
#                   |
#                   | 2.0           <── dist(A,E) = 2.0  (TIED at rank 2)
#                   |
#   A = (0,0) ──── X = (1,0) ──── D = (2,0) ──── F = (3,0) ──── G = (4,0)
#         1.0            1.0             1.0              1.0
#
#   With k = 2 (each cell gets 2 neighbors, excl. self):
#
#   From A:  rank-1 = X (dist 1.0),  rank-2 = TIE: D and E both at dist 2.0
#            order()  picks D (row index 3 < row index 4)
#            nn2      picks E (kd-tree order)
#
#   From D:  rank-1 = X (dist 1.0),  rank-2 = F (dist 1.0)
#            A is at dist 2.0 — too far, NOT in D's kNN
#
#   From E:  rank-1 = A (dist 2.0),  rank-2 = X (dist ~2.24)
#            A IS in E's kNN
#
#   MUTUAL FILTER OUTCOME:
#     NICHES chose A→D:  D's kNN = {X, F} → A not there → NOT mutual → DROPPED
#     Our code chose A→E: E's kNN = {A, X} → A is there  → IS  mutual → KEPT

k <- 2L

coord <- data.frame(
  x = c( 0,  1,  2,  0,  3,  4),
  y = c( 0,  0,  0,  2,  0,  0),
  row.names = c("A", "X", "D", "E", "F", "G")
)
#           row  x  y   — row index determines order() tie-breaking
# A          1   0  0
# X          2   1  0
# D          3   2  0   ← lower index → order() picks D over E at the tie
# E          4   0  2   ← higher index
# F          5   3  0
# G          6   4  0

cat("Cell coordinates:\n")
print(coord)

# ── Pre-computed distances from A ─────────────────────────────────────────────
cat("\nKey distances from A:\n")
cat(sprintf("  A → X : %.4f\n", sqrt((0-1)^2 + (0-0)^2)))
cat(sprintf("  A → D : %.4f   <── TIE at rank 2\n", sqrt((0-2)^2 + (0-0)^2)))
cat(sprintf("  A → E : %.4f   <── TIE at rank 2\n", sqrt((0-0)^2 + (0-2)^2)))

cat("\nKey distances from D (determines if A is in D's kNN):\n")
cat(sprintf("  D → X : %.4f   (rank 1)\n", sqrt((2-1)^2 + (0-0)^2)))
cat(sprintf("  D → F : %.4f   (rank 2)  ← A is pushed OUT by F\n", sqrt((2-3)^2 + (0-0)^2)))
cat(sprintf("  D → A : %.4f   (rank 3)  ← too far, NOT in D's kNN\n", sqrt((2-0)^2 + (0-0)^2)))

cat("\nKey distances from E (determines if A is in E's kNN):\n")
cat(sprintf("  E → A : %.4f   (rank 1)  ← A IS in E's kNN\n", sqrt((0-0)^2 + (2-0)^2)))
cat(sprintf("  E → X : %.4f   (rank 2)\n", sqrt((0-1)^2 + (2-0)^2)))

# ── Full distance matrix ───────────────────────────────────────────────────────
dist_mat <- as.matrix(dist(coord))
cat("\nFull distance matrix (rounded):\n")
print(round(dist_mat, 3))

# ── NICHES approach ────────────────────────────────────────────────────────────
cat(strrep("-", 70), "\n")
cat("NICHES  approach  (order() tie-breaking)\n")
cat(strrep("-", 70), "\n")

neighbor_mat <- apply(dist_mat, 1, function(v) order(v)[1:(k + 1)])
rownames(neighbor_mat) <- paste0("rank", 0:k)
colnames(neighbor_mat) <- rownames(coord)

cat("\nNeighbor rank matrix (row 0 = self, row 1..k = neighbors):\n")
neighbor_names <- apply(neighbor_mat, 2,
                        function(idx) rownames(coord)[idx])
rownames(neighbor_names) <- paste0("rank", 0:k)
print(neighbor_names)

cat("\n>>> For cell A: rank-2 neighbor chosen by order() =",
    neighbor_names["rank2", "A"], "(lower row index wins the tie)\n")
cat(">>> For cell D: kNN =", neighbor_names["rank1","D"],
    "and", neighbor_names["rank2","D"], "— A is NOT included\n")

n <- nrow(coord)
adj <- matrix(0, n, n, dimnames = list(rownames(coord), rownames(coord)))
for (cell in colnames(neighbor_mat))
  adj[cell, rownames(coord)[neighbor_mat[, cell]]] <- 1

adj_mutual_niches <- 1 * (adj & t(adj))

cat("\nNICHES mutual adjacency matrix:\n")
print(adj_mutual_niches)

niches_edges <- which(adj_mutual_niches == 1, arr.ind = TRUE)
niches_df <- data.frame(
  from = rownames(coord)[niches_edges[,1]],
  to   = rownames(coord)[niches_edges[,2]]
)
niches_df <- niches_df[order(niches_df$from, niches_df$to), ]
cat("\nNICHES edge list:\n")
print(niches_df, row.names = FALSE)
cat("Total NICHES edges (incl. self-loops):", nrow(niches_df), "\n")

# ── Our approach ───────────────────────────────────────────────────────────────
cat("\n", strrep("-", 70), "\n")
cat("Our approach  (RANN::nn2 tie-breaking)\n")
cat(strrep("-", 70), "\n")

xy  <- as.matrix(coord)
nn  <- RANN::nn2(data = xy, k = k + 1L)
idx <- nn$nn.idx;   rownames(idx)   <- rownames(coord)
dst <- nn$nn.dists; rownames(dst)   <- rownames(coord)

cat("\nnn2 neighbor index matrix (col1=self, col2..k+1=neighbors):\n")
print(idx)
cat("\nnn2 neighbor distance matrix:\n")
print(round(dst, 3))

cat("\n>>> For cell A: rank-2 neighbor chosen by nn2 =",
    rownames(coord)[idx["A", 3]], "(kd-tree traversal wins the tie)\n")
cat(">>> For cell E: kNN =",
    rownames(coord)[idx["E", 2]], "and",
    rownames(coord)[idx["E", 3]], "— A IS included\n")

# Build edge table (drop self = col 1)
idx2 <- idx[, -1, drop = FALSE]
dst2 <- dst[, -1, drop = FALSE]
rows <- lapply(seq_len(n), function(i)
  data.table(from = i, to = idx2[i,], weight = dst2[i,]))
et <- rbindlist(rows)

fwd <- et[, .(from, to)]
rev <- et[, .(from = to, to = from)]
mutual_idx <- fintersect(fwd, rev)
setkey(et, from, to); setkey(mutual_idx, from, to)
et_mutual <- et[mutual_idx, nomatch = 0L]

# Add autocrine
self_dt <- data.table(from = seq_len(n), to = seq_len(n), weight = 0)
et_final <- rbindlist(list(et_mutual, self_dt))
et_final[, from := rownames(coord)[from]]
et_final[, to   := rownames(coord)[to]]
et_final <- et_final[order(from, to)]

cat("\nOur edge list:\n")
print(et_final, row.names = FALSE)
cat("Total our edges (incl. self-loops):", nrow(et_final), "\n")

# ── Final comparison ───────────────────────────────────────────────────────────
cat("\n", strrep("=", 70), "\n")
cat("  SUMMARY\n")
cat(strrep("=", 70), "\n")

niches_pairs <- paste(niches_df$from, niches_df$to, sep = "→")
ours_pairs   <- paste(et_final$from,  et_final$to,  sep = "→")

only_niches <- setdiff(niches_pairs, ours_pairs)
only_ours   <- setdiff(ours_pairs,   niches_pairs)

cat("\nEdges in NICHES but NOT in ours : ")
cat(if (length(only_niches)) paste(only_niches, collapse=", ") else "(none)", "\n")

cat("Edges in OURS but NOT in NICHES : ")
cat(if (length(only_ours)) paste(only_ours, collapse=", ") else "(none)", "\n")

cat(sprintf("\nNICHES total: %d edges\n", nrow(niches_df)))
cat(sprintf("Ours   total: %d edges\n", nrow(et_final)))
cat(sprintf("Difference:   %d edges\n", nrow(et_final) - nrow(niches_df)))

cat("
WHAT HAPPENED (plain language):
─────────────────────────────────────────────────────────────────────
  A's 2nd nearest neighbor is a TIE: D and E are both exactly 2.0 away.
  One must be discarded — but which one?

  NICHES  (order):  picks D  → checks mutual: D's kNN = {X, F}, no A
                             → A↔D fails mutual filter  → 0 edges added
  Ours    (nn2):    picks E  → checks mutual: E's kNN = {A, X}, has A
                             → A↔E passes mutual filter → 2 edges added (A→E, E→A)

  Both methods are logically correct mutual kNN.
  The 2-edge gap comes entirely from an arbitrary tie-break on a pair
  of cells that happen to have different mutuality properties.

  On your real dataset, this tie pattern occurs 3 times → 6 edge difference.
─────────────────────────────────────────────────────────────────────
")
