# ==============================================================================
# FILE:    generate_pseudo_data.R
# PURPOSE: Generate pseudo spatial transcriptomics datasets for benchmarking
#          computeCellToCellSpatial() against NICHES. Produces Seurat objects
#          with realistic gene expression and clustered spatial coordinates.
#
# OUTPUT:  One .rds file per (n_genes, n_cells) combination, saved to
#          <out_dir>/pseudo_gXXXX_cYYYY.rds
#
# SCALING DESIGN:
#   Two independent sweeps to isolate the effect of each dimension:
#     (A) Gene sweep  — n_genes varies, n_cells fixed at FIXED_CELLS
#     (B) Cell sweep  — n_cells varies, n_genes fixed at FIXED_GENES
#   This gives a 2×5 = 10-dataset grid (no redundant center point counted twice).
#
#   n_genes levels : 500, 2000, 5000, 10000, 20000
#   n_cells levels : 100,  500, 1000,  3000,  5000
#   Fixed n_cells  : 1000   (used in gene sweep)
#   Fixed n_genes  : 5000   (used in cell sweep)
#
# GENE NAMES:
#   Downloaded live from HGNC (genenames.org) at runtime — always current,
#   no local gene file required.
#
# SPATIAL LAYOUT:
#   Gaussian Mixture Model:
#     - K cluster centers placed randomly in [0, 1000] x [0, 1000]
#     - Cells assigned to clusters with Dirichlet-distributed weights
#       (unequal cluster sizes, as in real tissue sections)
#     - Each cluster: bivariate Normal with its own spread (sigma)
#     - ~5% of cells placed uniformly as background scatter
#
# COUNT DISTRIBUTION:
#   Negative Binomial (Gamma-Poisson) per gene:
#     - Gene mean drawn from LogNormal to mimic housekeeping vs rare genes
#     - Dispersion drawn from Gamma to mimic overdispersion
#     - Result: realistic sparsity (~85-95%), integer counts
#
# DEPS: Matrix, Seurat (>= 5.0)
# ==============================================================================

suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
})

set.seed(42)

# ==============================================================================
# SECTION 0 ── User-adjustable parameters
# ==============================================================================

out_dir    <- "/Volumes/T5 EVO/NICHESv2/pseudo_data"        # output directory (created if absent)
dir.create(out_dir, showWarnings = FALSE)

# Sweep levels
GENE_LEVELS  <- c(500, 2000, 5000, 10000, 20000)
CELL_LEVELS  <- c(100,  500, 1000,  3000,  5000)

# Fixed dimension for each sweep
FIXED_CELLS  <- 1000   # held constant during gene sweep
FIXED_GENES  <- 5000   # held constant during cell sweep

# Spatial parameters
N_CLUSTERS      <- 8     # number of spatial clusters (Gaussian blobs)
BACKGROUND_FRAC <- 0.05  # fraction of cells placed uniformly as background
SPACE_SIZE      <- 1000  # coordinate space: [0, SPACE_SIZE] x [0, SPACE_SIZE]

# ==============================================================================
# SECTION 1 ── Fetch HGNC human gene symbols
#
# Downloads the approved symbol list from HGNC (genenames.org).
# Filters to clean, alphanumeric symbols (removes antisense loci,
# pseudogenes with numeric suffixes, etc.).
# A local cache is kept so repeated runs don't re-download.
# ==============================================================================

.fetch_hgnc_genes <- function(cache_file = file.path(tempdir(), "hgnc_genes.rds")) {

  if (file.exists(cache_file)) {
    message("[HGNC] Loading cached gene list.")
    return(readRDS(cache_file))
  }

  message("[HGNC] Downloading gene symbols from genenames.org ...")
  url  <- paste0("https://www.genenames.org/cgi-bin/download/custom?",
                 "col=gd_app_sym&status=Approved&hgnc_dbtag=on",
                 "&order_by=gd_app_sym_sort&format=text&submit=submit")
  tmp  <- tempfile(fileext = ".txt")
  tryCatch(
    download.file(url, tmp, quiet = TRUE),
    error = function(e) stop("[HGNC] Download failed. Check internet connection.\n", e)
  )

  raw   <- readLines(tmp, warn = FALSE)
  genes <- trimws(raw[-1])   # drop header "Approved symbol"
  genes <- genes[nchar(genes) > 0]

  # Keep clean protein-coding-like symbols:
  #   - no -AS (antisense), -DT (divergent transcript), -IT (intronic)
  #   - no trailing Pnumber (pseudogene suffix e.g. ABCB10P1)
  #   - purely uppercase alphanumeric, length 2-12
  keep  <- !grepl("-AS|-DT|-IT", genes) &
            !grepl("P[0-9]+$",   genes) &
            nchar(genes) >= 2 & nchar(genes) <= 12 &
            grepl("^[A-Z][A-Z0-9]*$", genes)

  genes <- genes[keep]
  message(sprintf("[HGNC] %d clean gene symbols retrieved.", length(genes)))

  saveRDS(genes, cache_file)
  genes
}

ALL_HGNC_GENES <- .fetch_hgnc_genes()

# ==============================================================================
# SECTION 2 ── Spatial coordinate generation (Gaussian Mixture Model)
#
# Strategy:
#   1. Place K cluster centers by sampling from a grid with jitter so
#      clusters are spread across the tissue space, not all bunched.
#   2. Assign each cell to a cluster using Dirichlet weights (alpha=0.5)
#      which creates realistically unequal cluster sizes.
#   3. Draw (x, y) from bivariate Normal; each cluster has its own sigma
#      drawn from Uniform[20, 120] mimicking dense vs diffuse regions.
#   4. Replace ~5% of cells with Uniform draws as background scatter.
#   5. Clip all coordinates to [0, SPACE_SIZE].
# ==============================================================================

.generate_spatial_coords <- function(n_cells,
                                     n_clusters    = N_CLUSTERS,
                                     bg_frac       = BACKGROUND_FRAC,
                                     space_size    = SPACE_SIZE,
                                     seed          = NULL) {
  if (!is.null(seed)) set.seed(seed)

  # ── Cluster centers: grid-jittered so they span the space ─────────────────
  grid_k   <- ceiling(sqrt(n_clusters))
  grid_pts <- seq(space_size * 0.1, space_size * 0.9, length.out = grid_k)
  grid_xy  <- expand.grid(x = grid_pts, y = grid_pts)
  chosen   <- sample(nrow(grid_xy), n_clusters)
  centers  <- grid_xy[chosen, ]
  centers$x <- centers$x + runif(n_clusters, -space_size * 0.08, space_size * 0.08)
  centers$y <- centers$y + runif(n_clusters, -space_size * 0.08, space_size * 0.08)

  # ── Per-cluster spread ─────────────────────────────────────────────────────
  sigmas <- runif(n_clusters, 20, 120)

  # ── Cluster assignment: Dirichlet(alpha=0.5) weights ──────────────────────
  # Gamma(0.5, 1) samples, normalized → Dirichlet(0.5, ..., 0.5)
  # alpha < 1 creates sparsity: some clusters much larger than others
  raw_weights <- rgamma(n_clusters, shape = 0.5, rate = 1)
  weights     <- raw_weights / sum(raw_weights)
  cluster_ids <- sample(n_clusters, n_cells, replace = TRUE, prob = weights)

  # ── Draw coordinates ───────────────────────────────────────────────────────
  x <- numeric(n_cells)
  y <- numeric(n_cells)
  for (k in seq_len(n_clusters)) {
    idx    <- which(cluster_ids == k)
    if (length(idx) == 0) next
    x[idx] <- rnorm(length(idx), mean = centers$x[k], sd = sigmas[k])
    y[idx] <- rnorm(length(idx), mean = centers$y[k], sd = sigmas[k])
  }

  # ── Replace background fraction with uniform scatter ──────────────────────
  n_bg  <- max(1L, round(n_cells * bg_frac))
  bg_idx <- sample(n_cells, n_bg)
  x[bg_idx] <- runif(n_bg, 0, space_size)
  y[bg_idx] <- runif(n_bg, 0, space_size)

  # ── Clip to valid space ───────────────────────────────────────────────────
  x <- pmax(0, pmin(space_size, x))
  y <- pmax(0, pmin(space_size, y))

  data.frame(x = x, y = y, cluster = cluster_ids)
}

# ==============================================================================
# SECTION 3 ── Count matrix generation (Negative Binomial)
#
# Per-gene parameters:
#   mu    ~ LogNormal(log(2), 1.5)   → range ~0.05 to ~200 counts/cell
#           (mimics the mixture of rare cell-type genes and housekeeping genes)
#   theta ~ Gamma(2, 0.5)            → range ~0.1 to ~20
#           (overdispersion: smaller theta = more zero-inflation)
#
# NB parameterization: var = mu + mu^2/theta
# Implemented via Poisson(lambda) where lambda ~ Gamma(theta, theta/mu)
# This is exact and avoids dense intermediate matrices by processing genes
# in chunks and writing directly to sparse triplet format.
# ==============================================================================

.generate_count_matrix <- function(n_genes, n_cells, gene_names, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  stopifnot(length(gene_names) >= n_genes)
  genes <- sample(gene_names, n_genes)   # random subset of HGNC symbols
  cells <- paste0("Cell_", seq_len(n_cells))

  # Per-gene mean and dispersion
  mu_vec    <- rlnorm(n_genes, meanlog = log(2), sdlog = 1.5)
  theta_vec <- rgamma(n_genes, shape = 2, rate = 0.5)

  # Build sparse matrix in chunks to limit peak memory
  chunk_size <- 500L
  n_chunks   <- ceiling(n_genes / chunk_size)

  all_i <- integer(0)
  all_j <- integer(0)
  all_x <- numeric(0)

  for (ch in seq_len(n_chunks)) {
    idx   <- ((ch - 1) * chunk_size + 1):min(ch * chunk_size, n_genes)
    mu_ch <- mu_vec[idx]
    th_ch <- theta_vec[idx]

    # Draw lambda from Gamma, then counts from Poisson(lambda)
    # Gamma(shape=theta, rate=theta/mu) → E[lambda]=mu, Var[lambda]=mu^2/theta
    lambda_mat <- matrix(
      rgamma(length(idx) * n_cells,
             shape = rep(th_ch, each = n_cells),
             rate  = rep(th_ch / mu_ch, each = n_cells)),
      nrow = length(idx), ncol = n_cells
    )
    count_mat <- matrix(
      rpois(length(lambda_mat), lambda = lambda_mat),
      nrow = length(idx), ncol = n_cells
    )

    # Extract non-zero entries as triplets
    nz        <- which(count_mat > 0, arr.ind = TRUE)
    if (nrow(nz) > 0) {
      all_i <- c(all_i, idx[nz[, 1]])      # global gene index
      all_j <- c(all_j, nz[, 2])           # cell index
      all_x <- c(all_x, count_mat[nz])
    }
    rm(lambda_mat, count_mat, nz)
  }

  mtx <- Matrix::sparseMatrix(
    i        = all_i,
    j        = all_j,
    x        = all_x,
    dims     = c(n_genes, n_cells),
    dimnames = list(genes, cells)
  )

  sparsity <- 1 - nnzero(mtx) / prod(dim(mtx))
  message(sprintf("    Count matrix: %d genes x %d cells | Sparsity: %.1f%%",
                  n_genes, n_cells, sparsity * 100))
  mtx
}

# ==============================================================================
# SECTION 4 ── Assemble Seurat object
#
# Stores:
#   - RNA assay (counts layer; data layer set equal to counts)
#   - x, y, cluster in meta.data
#   - dataset provenance in misc slot
# ==============================================================================

.make_seurat <- function(count_mtx, coords, tag) {

  cells <- colnames(count_mtx)
  stopifnot(nrow(coords) == length(cells))

  meta <- data.frame(
    x       = coords$x,
    y       = coords$y,
    cluster = paste0("C", coords$cluster),
    row.names = cells
  )

  obj <- CreateSeuratObject(
    counts = count_mtx,
    assay  = "RNA",
    meta.data = meta,
    project   = tag
  )

  # Set the data layer equal to counts so that NICHES (which reads @data)
  # and our functions operate on identical values during benchmarking.
  LayerData(obj, assay = "RNA", layer = "data") <- count_mtx

  obj@misc$pseudo_params <- list(
    tag        = tag,
    n_genes    = nrow(count_mtx),
    n_cells    = ncol(count_mtx),
    created_at = Sys.time()
  )
  obj
}

# ==============================================================================
# SECTION 5 ── Build scaling grid and generate all datasets
#
# Grid:
#   Sweep A (gene sweep):  n_genes in GENE_LEVELS, n_cells = FIXED_CELLS
#   Sweep B (cell sweep):  n_cells in CELL_LEVELS, n_genes = FIXED_GENES
#
# The (FIXED_GENES, FIXED_CELLS) = (5000, 1000) combination appears in
# both sweeps. It is generated only once and reused.
# ==============================================================================

# Build grid, marking which sweep each row belongs to (for later analysis)
grid_A <- data.frame(n_genes = GENE_LEVELS, n_cells = FIXED_CELLS, sweep = "gene")
grid_B <- data.frame(n_genes = FIXED_GENES, n_cells = CELL_LEVELS, sweep = "cell")

# Deduplicate the shared center point (5000 genes, 1000 cells)
full_grid <- rbind(grid_A, grid_B)
full_grid <- full_grid[!duplicated(full_grid[, c("n_genes", "n_cells")]), ]
rownames(full_grid) <- NULL

message(sprintf("\n[generate_pseudo_data] Generating %d datasets:\n", nrow(full_grid)))
print(full_grid)
cat("\n")

for (i in seq_len(nrow(full_grid))) {

  ng   <- full_grid$n_genes[i]
  nc   <- full_grid$n_cells[i]
  tag  <- sprintf("pseudo_g%05d_c%05d", ng, nc)
  path <- file.path(out_dir, paste0(tag, ".rds"))

  if (file.exists(path)) {
    message(sprintf("[%d/%d] Skipping %s — already exists.", i, nrow(full_grid), tag))
    next
  }

  message(sprintf("[%d/%d] Generating %s ...", i, nrow(full_grid), tag))
  t0 <- proc.time()

  # Reproducible per-dataset seed derived from dimensions
  ds_seed <- ng * 100000L + nc

  coords    <- .generate_spatial_coords(nc, seed = ds_seed)
  count_mtx <- .generate_count_matrix(ng, nc, ALL_HGNC_GENES, seed = ds_seed + 1L)
  obj       <- .make_seurat(count_mtx, coords, tag)

  saveRDS(obj, path)
  elapsed <- (proc.time() - t0)[["elapsed"]]
  message(sprintf("    Saved: %s  (%.1f sec)\n", path, elapsed))
}

message("[generate_pseudo_data] All datasets complete.")
message(sprintf("Files written to: %s/", out_dir))
message("\nDataset summary:")
print(full_grid)
