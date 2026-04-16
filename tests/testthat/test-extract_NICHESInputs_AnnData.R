# ==============================================================================
# FILE:    tests/testthat/test-extract_NICHESInputs_AnnData.R
# PURPOSE: Unit tests for extract_NICHESInputs_AnnData().
#          All h5ad files are written inline via .write_test_h5ad() using rhdf5
#          only — no Python, no conda, no zellkonverter.
#          All blocks skip automatically when rhdf5 is not installed.
# ==============================================================================

# ── File-level helpers ─────────────────────────────────────────────────────────

# Writes a minimal spec-compliant .h5ad file using rhdf5 only.
# counts : dgCMatrix, genes x cells (n_vars x n_obs)
# meta   : data.frame, rownames = barcodes, nrow = ncol(counts)
.write_test_h5ad <- function(file, counts, meta, with.counts.layer = TRUE) {
  n.obs      <- ncol(counts)
  n.vars     <- nrow(counts)
  barcodes   <- colnames(counts)
  gene.names <- rownames(counts)
  
  # AnnData stores obs x vars; transpose to cells x genes for X
  X <- Matrix::t(counts)
  
  write_group_attrs <- function(path, attrs) {
    fid <- rhdf5::H5Fopen(file)
    on.exit(rhdf5::H5Fclose(fid), add = TRUE)
    gid <- rhdf5::H5Gopen(fid, path)
    on.exit(rhdf5::H5Gclose(gid), add = TRUE)
    for (nm in names(attrs))
      rhdf5::h5writeAttribute(attrs[[nm]], gid, nm)
  }
  
  rhdf5::h5createFile(file)
  
  # X (CSC: n_obs x n_vars)
  rhdf5::h5createGroup(file, "X")
  rhdf5::h5write(X@x, file, "X/data")
  rhdf5::h5write(X@i, file, "X/indices")
  rhdf5::h5write(X@p, file, "X/indptr")
  write_group_attrs("X", list(
    "encoding-type"    = "csc_matrix",
    "encoding-version" = "0.1.0",
    "shape"            = c(n.obs, n.vars)
  ))
  
  # obs
  rhdf5::h5createGroup(file, "obs")
  rhdf5::h5write(barcodes, file, "obs/_index")
  for (cn in colnames(meta)) {
    col  <- meta[[cn]]
    path <- paste0("obs/", cn)
    if (is.character(col) || is.factor(col)) {
      col   <- as.character(col)
      cats  <- sort(unique(col[!is.na(col)]))
      codes <- match(col, cats) - 1L
      codes[is.na(col)] <- -1L
      rhdf5::h5createGroup(file, path)
      rhdf5::h5write(cats,              file, paste0(path, "/categories"))
      rhdf5::h5write(as.integer(codes), file, paste0(path, "/codes"))
      write_group_attrs(path, list(
        "encoding-type"    = "categorical",
        "encoding-version" = "0.2.0",
        "ordered"          = 0L
      ))
    } else {
      rhdf5::h5write(as.vector(col), file, path)
    }
  }
  write_group_attrs("obs", list(
    "_index"           = "_index",
    "encoding-type"    = "dataframe",
    "encoding-version" = "0.2.0",
    "column-order"     = colnames(meta)
  ))
  
  # var
  rhdf5::h5createGroup(file, "var")
  rhdf5::h5write(gene.names, file, "var/_index")
  write_group_attrs("var", list(
    "_index"           = "_index",
    "encoding-type"    = "dataframe",
    "encoding-version" = "0.2.0"
  ))
  
  # layers
  rhdf5::h5createGroup(file, "layers")
  if (with.counts.layer) {
    rhdf5::h5createGroup(file, "layers/counts")
    rhdf5::h5write(X@x, file, "layers/counts/data")
    rhdf5::h5write(X@i, file, "layers/counts/indices")
    rhdf5::h5write(X@p, file, "layers/counts/indptr")
    write_group_attrs("layers/counts", list(
      "encoding-type"    = "csc_matrix",
      "encoding-version" = "0.1.0",
      "shape"            = c(n.obs, n.vars)
    ))
  }
  
  # uns (required by spec)
  rhdf5::h5createGroup(file, "uns")
  
  invisible(file)
}

.make_test_counts <- function(n.genes = 20L, n.cells = 10L) {
  set.seed(42L)
  m <- abs(Matrix::rsparsematrix(n.genes, n.cells, density = 0.6))
  rownames(m) <- paste0("gene", seq_len(n.genes))
  colnames(m) <- paste0("cell", seq_len(n.cells))
  m
}

.make_test_meta <- function(n.cells = 10L, extra.cols = NULL) {
  df <- data.frame(
    cell.type = rep(c("A", "B"), length.out = n.cells),
    row.names = paste0("cell", seq_len(n.cells)),
    stringsAsFactors = FALSE
  )
  if (!is.null(extra.cols)) df <- cbind(df, extra.cols)
  df
}

# ── Output structure ───────────────────────────────────────────────────────────

test_that("extract_NICHESInputs_AnnData: returns list with count.mtx and meta.data", {
  skip_if_not_installed("rhdf5")
  skip_if_not_installed("withr")
  
  counts   <- .make_test_counts()
  meta     <- .make_test_meta()
  tmp.file <- withr::local_tempfile(fileext = ".h5ad")
  .write_test_h5ad(tmp.file, counts, meta)
  
  result <- extract_NICHESInputs_AnnData(tmp.file, layer = "counts",
                                         verbose = FALSE)
  expect_type(result, "list")
  expect_named(result, c("count.mtx", "meta.data"))
})

test_that("extract_NICHESInputs_AnnData: count.mtx is genes x cells (dgCMatrix)", {
  skip_if_not_installed("rhdf5")
  skip_if_not_installed("withr")
  
  counts   <- .make_test_counts(n.genes = 20L, n.cells = 10L)
  meta     <- .make_test_meta(n.cells = 10L)
  tmp.file <- withr::local_tempfile(fileext = ".h5ad")
  .write_test_h5ad(tmp.file, counts, meta)
  
  result <- extract_NICHESInputs_AnnData(tmp.file, layer = "counts",
                                         verbose = FALSE)
  
  expect_true(inherits(result$count.mtx, "dgCMatrix"))
  expect_equal(nrow(result$count.mtx), 20L)
  expect_equal(ncol(result$count.mtx), 10L)
  expect_equal(rownames(result$count.mtx), paste0("gene", seq_len(20L)))
  expect_equal(colnames(result$count.mtx), paste0("cell", seq_len(10L)))
})

test_that("extract_NICHESInputs_AnnData: layer='X' reads adata.X", {
  skip_if_not_installed("rhdf5")
  skip_if_not_installed("withr")
  
  counts   <- .make_test_counts()
  meta     <- .make_test_meta()
  tmp.file <- withr::local_tempfile(fileext = ".h5ad")
  .write_test_h5ad(tmp.file, counts, meta, with.counts.layer = FALSE)
  
  result <- extract_NICHESInputs_AnnData(tmp.file, layer = "X",
                                         verbose = FALSE)
  expect_equal(dim(result$count.mtx), c(20L, 10L))
})

test_that("extract_NICHESInputs_AnnData: matrix values survive roundtrip", {
  skip_if_not_installed("rhdf5")
  skip_if_not_installed("withr")
  
  counts   <- .make_test_counts()
  meta     <- .make_test_meta()
  tmp.file <- withr::local_tempfile(fileext = ".h5ad")
  .write_test_h5ad(tmp.file, counts, meta)
  
  result <- extract_NICHESInputs_AnnData(tmp.file, layer = "counts",
                                         verbose = FALSE)
  expect_equal(as.matrix(result$count.mtx), as.matrix(counts), tolerance = 1e-6)
})

# ── Metadata correctness ───────────────────────────────────────────────────────

test_that("extract_NICHESInputs_AnnData: meta.data rownames are cell barcodes", {
  skip_if_not_installed("rhdf5")
  skip_if_not_installed("withr")
  
  counts   <- .make_test_counts()
  meta     <- .make_test_meta()
  tmp.file <- withr::local_tempfile(fileext = ".h5ad")
  .write_test_h5ad(tmp.file, counts, meta)
  
  result <- extract_NICHESInputs_AnnData(tmp.file, layer = "counts",
                                         verbose = FALSE)
  expect_true(is.data.frame(result$meta.data))
  expect_equal(rownames(result$meta.data), paste0("cell", seq_len(10L)))
  expect_true("cell.type" %in% colnames(result$meta.data))
})

test_that("extract_NICHESInputs_AnnData: categorical column roundtrip is correct", {
  skip_if_not_installed("rhdf5")
  skip_if_not_installed("withr")
  
  counts   <- .make_test_counts()
  meta     <- .make_test_meta()
  tmp.file <- withr::local_tempfile(fileext = ".h5ad")
  .write_test_h5ad(tmp.file, counts, meta)
  
  result <- extract_NICHESInputs_AnnData(tmp.file, layer = "counts",
                                         verbose = FALSE)
  expect_equal(result$meta.data$cell.type,
               rep(c("A", "B"), length.out = 10L))
})

test_that("extract_NICHESInputs_AnnData: numeric obs column roundtrip is correct", {
  skip_if_not_installed("rhdf5")
  skip_if_not_installed("withr")
  
  counts   <- .make_test_counts()
  extra    <- data.frame(score = seq_len(10L))
  meta     <- .make_test_meta(extra.cols = extra)
  tmp.file <- withr::local_tempfile(fileext = ".h5ad")
  .write_test_h5ad(tmp.file, counts, meta)
  
  result <- extract_NICHESInputs_AnnData(tmp.file, layer = "counts",
                                         verbose = FALSE)
  expect_true("score" %in% colnames(result$meta.data))
  expect_equal(as.integer(result$meta.data$score), seq_len(10L))
})

test_that("extract_NICHESInputs_AnnData: obs.cols subsets metadata correctly", {
  skip_if_not_installed("rhdf5")
  skip_if_not_installed("withr")
  
  counts   <- .make_test_counts()
  extra    <- data.frame(batch = rep(1L:2L, 5L))
  meta     <- .make_test_meta(extra.cols = extra)
  tmp.file <- withr::local_tempfile(fileext = ".h5ad")
  .write_test_h5ad(tmp.file, counts, meta)
  
  result <- extract_NICHESInputs_AnnData(tmp.file, layer = "counts",
                                         obs.cols = "cell.type",
                                         verbose = FALSE)
  expect_equal(colnames(result$meta.data), "cell.type")
  expect_false("batch" %in% colnames(result$meta.data))
})

# ── Input validation ───────────────────────────────────────────────────────────

test_that("extract_NICHESInputs_AnnData: errors on missing file", {
  skip_if_not_installed("rhdf5")
  
  expect_error(
    extract_NICHESInputs_AnnData("/nonexistent/path/data.h5ad", verbose = FALSE),
    "File not found"
  )
})

test_that("extract_NICHESInputs_AnnData: errors when obs.cols are absent", {
  skip_if_not_installed("rhdf5")
  skip_if_not_installed("withr")
  
  counts   <- .make_test_counts()
  meta     <- .make_test_meta()
  tmp.file <- withr::local_tempfile(fileext = ".h5ad")
  .write_test_h5ad(tmp.file, counts, meta)
  
  expect_error(
    extract_NICHESInputs_AnnData(tmp.file, layer = "counts",
                                 obs.cols = "nonexistent.col",
                                 verbose = FALSE),
    "not found in adata\\.obs"
  )
})

test_that("extract_NICHESInputs_AnnData: errors on missing layer with available names listed", {
  skip_if_not_installed("rhdf5")
  skip_if_not_installed("withr")
  
  counts   <- .make_test_counts()
  meta     <- .make_test_meta()
  tmp.file <- withr::local_tempfile(fileext = ".h5ad")
  .write_test_h5ad(tmp.file, counts, meta, with.counts.layer = TRUE)
  
  err.msg <- tryCatch(
    extract_NICHESInputs_AnnData(tmp.file, layer = "NOTLAYER", verbose = FALSE),
    error = function(e) conditionMessage(e)
  )
  expect_match(err.msg, "Available layers")
  expect_match(err.msg, "counts")
})

test_that("extract_NICHESInputs_AnnData: warns on non-.h5ad file extension", {
  skip_if_not_installed("rhdf5")
  skip_if_not_installed("withr")
  
  counts   <- .make_test_counts()
  meta     <- .make_test_meta()
  tmp.file <- withr::local_tempfile(fileext = ".txt")
  .write_test_h5ad(tmp.file, counts, meta)
  
  expect_warning(
    extract_NICHESInputs_AnnData(tmp.file, layer = "counts", verbose = FALSE),
    "h5ad"
  )
})

# ── Verbosity ──────────────────────────────────────────────────────────────────

test_that("extract_NICHESInputs_AnnData: verbose=TRUE emits spatial hint message", {
  skip_if_not_installed("rhdf5")
  skip_if_not_installed("withr")
  
  counts   <- .make_test_counts()
  meta     <- .make_test_meta()
  tmp.file <- withr::local_tempfile(fileext = ".h5ad")
  .write_test_h5ad(tmp.file, counts, meta)
  
  expect_message(
    extract_NICHESInputs_AnnData(tmp.file, layer = "counts", verbose = TRUE),
    "spatial"
  )
})

test_that("extract_NICHESInputs_AnnData: verbose=FALSE suppresses all messages", {
  skip_if_not_installed("rhdf5")
  skip_if_not_installed("withr")
  
  counts   <- .make_test_counts()
  meta     <- .make_test_meta()
  tmp.file <- withr::local_tempfile(fileext = ".h5ad")
  .write_test_h5ad(tmp.file, counts, meta)
  
  expect_no_message(
    extract_NICHESInputs_AnnData(tmp.file, layer = "counts", verbose = FALSE)
  )
})