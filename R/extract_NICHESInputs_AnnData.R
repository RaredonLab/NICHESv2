# extract_NICHESInputs_AnnData.R ------------------------------------------
# Internal helpers for parsing .h5ad (AnnData HDF5) files via rhdf5.
# No Python, no conda, no basilisk -- pure C HDF5 bindings.
#
# h5ad layout (AnnData spec):
#   /X                      -- main matrix (dense or sparse CSR/CSC)
#   /layers/<name>          -- named layers (same format as X)
#   /obs/_index             -- cell barcodes (obs_names)
#   /obs/<col>              -- cell metadata columns
#   /var/_index             -- gene names (var_names)
#   /uns                    -- unstructured metadata (ignored)
#
# All matrices in h5ad are stored obs x vars = cells x genes.
# All functions below return data in genes x cells convention (NICHESv2 standard).
#
# Sparse matrix HDF5 layout (CSR or CSC):
#   <group>/data    -- non-zero values
#   <group>/indices -- row (CSC) or col (CSR) indices, 0-based
#   <group>/indptr  -- col (CSC) or row (CSR) pointers, 0-based
#   attrs: encoding-type = "csr_matrix" | "csc_matrix"
#          shape          = [n_obs, n_vars]
#
# Dense matrix HDF5 layout:
#   <dataset>       -- 2D array stored C-order (row-major) as [n_obs, n_vars]
#   attrs: encoding-type = "array"  (may be absent in old files)
#   rhdf5 reverses dim order on read -> R matrix is [n_vars, n_obs] = genes x cells

#' @noRd
.h5ad_attrs <- function(file, path) {
  tryCatch(
    rhdf5::h5readAttributes(file, path),
    error = function(e) list()
  )
}

#' @noRd
.h5ad_encoding <- function(file, path) {
  attrs <- .h5ad_attrs(file, path)
  attrs[["encoding-type"]] %||% ""
}

# Read cell barcodes from obs/_index (or whichever column the _index attr points to).
#' @noRd
.h5ad_obs_names <- function(file) {
  obs.attrs <- .h5ad_attrs(file, "obs")
  idx.col   <- obs.attrs[["_index"]] %||% "_index"
  as.character(rhdf5::h5read(file, paste0("obs/", idx.col)))
}

# Read gene names from var/_index.
#' @noRd
.h5ad_var_names <- function(file) {
  var.attrs <- .h5ad_attrs(file, "var")
  idx.col   <- var.attrs[["_index"]] %||% "_index"
  as.character(rhdf5::h5read(file, paste0("var/", idx.col)))
}

# List all available layer names: "X" plus anything stored in /layers.
#' @noRd
.h5ad_available_layers <- function(file) {
  h5.info     <- rhdf5::h5ls(file)
  layer.names <- h5.info$name[h5.info$group == "/layers"]
  c("X", layer.names)
}

# Read a CSR or CSC sparse matrix group.
# Returns dgCMatrix in obs x vars orientation (cells x genes).
#' @noRd
.h5ad_read_sparse <- function(file, path) {
  attrs  <- .h5ad_attrs(file, path)
  enc    <- attrs[["encoding-type"]]
  shape  <- as.integer(attrs[["shape"]])  # [n_obs, n_vars]
  n.obs  <- shape[1L]
  n.vars <- shape[2L]

  vals    <- as.numeric(rhdf5::h5read(file, paste0(path, "/data")))
  indices <- as.integer(rhdf5::h5read(file, paste0(path, "/indices"))) + 1L
  indptr  <- as.integer(rhdf5::h5read(file, paste0(path, "/indptr")))

  if (enc == "csr_matrix") {
    # indptr length = n_obs + 1; indices are column positions
    Matrix::sparseMatrix(
      i    = rep.int(seq_len(n.obs), diff(indptr)),
      j    = indices,
      x    = vals,
      dims = c(n.obs, n.vars)
    )
  } else {
    # csc_matrix: indptr length = n_vars + 1; indices are row positions
    Matrix::sparseMatrix(
      i    = indices,
      j    = rep.int(seq_len(n.vars), diff(indptr)),
      x    = vals,
      dims = c(n.obs, n.vars)
    )
  }
}

# Read any matrix (sparse or dense) from an h5ad path.
# Always returns dgCMatrix in vars x obs = genes x cells orientation.
#' @noRd
.h5ad_read_matrix <- function(file, path, n.obs, n.vars) {
  enc <- .h5ad_encoding(file, path)

  if (enc %in% c("csr_matrix", "csc_matrix")) {
    mat <- .h5ad_read_sparse(file, path)   # obs x vars (cells x genes)
    return(methods::as(Matrix::t(mat), "dgCMatrix"))  # -> vars x obs
  }

  # Dense array. rhdf5 reverses C-order HDF5 dims on read:
  # a [n_obs, n_vars] HDF5 array becomes [n_vars, n_obs] in R = genes x cells.
  raw <- rhdf5::h5read(file, path)
  mat <- if (is.matrix(raw)) raw else as.matrix(raw)

  # Defensive orientation check: if dims are obs x vars, transpose.
  if (nrow(mat) == n.obs && ncol(mat) == n.vars) mat <- t(mat)

  methods::as(mat, "dgCMatrix")
}

# Read obs group as a data.frame. Handles categorical and plain columns.
#' @noRd
.h5ad_read_obs <- function(file, obs.names) {
  h5.info   <- rhdf5::h5ls(file)
  obs.attrs <- .h5ad_attrs(file, "obs")
  idx.col   <- obs.attrs[["_index"]] %||% "_index"

  # Direct children of /obs; drop the index column and legacy __categories group
  all.cols  <- h5.info$name[h5.info$group == "/obs"]
  col.names <- all.cols[!all.cols %in% c(idx.col, "__categories")]

  if (length(col.names) == 0L) {
    return(data.frame(row.names = obs.names))
  }

  result <- vector("list", length(col.names))
  names(result) <- col.names

  for (cn in col.names) {
    path <- paste0("obs/", cn)
    enc  <- .h5ad_encoding(file, path)

    if (enc == "categorical") {
      # Newer AnnData encodes categoricals as categories array + codes array.
      # codes are 0-based integers; -1 means NA.
      cats    <- as.character(rhdf5::h5read(file, paste0(path, "/categories")))
      codes.0 <- as.integer(rhdf5::h5read(file, paste0(path, "/codes")))
      na.mask <- (codes.0 == -1L)
      codes.0[na.mask] <- 0L          # prevent 0-index into cats
      vals    <- cats[codes.0 + 1L]   # 0-based -> 1-based
      vals[na.mask] <- NA_character_
      result[[cn]] <- vals
    } else {
      result[[cn]] <- tryCatch(
        as.vector(rhdf5::h5read(file, path)),
        error = function(e) {
          warning("Could not read obs column '", cn, "': ",
                  conditionMessage(e), call. = FALSE)
          rep(NA_character_, length(obs.names))
        }
      )
    }
  }

  df <- as.data.frame(result, stringsAsFactors = FALSE)
  rownames(df) <- obs.names
  df
}

## Exported function ----

#' Extract inputs for create_NICHESObject from an AnnData (.h5ad) file
#'
#' Reads a \code{.h5ad} file using \pkg{rhdf5} (no Python, no conda, no
#' basilisk -- pure C HDF5 bindings, compatible with R \eqn{\geq} 4.2) and
#' returns the count matrix and cell metadata as a named list ready to pass
#' directly to \code{\link{create_NICHESObject}}.
#'
#' The AnnData layer structure maps to the \code{layer} argument as follows:
#' \itemize{
#'   \item \code{adata.X}                    \eqn{\rightarrow} \code{layer = "X"}
#'   \item \code{adata.layers[["counts"]]}   \eqn{\rightarrow} \code{layer = "counts"}
#'   \item any other named layer             \eqn{\rightarrow} \code{layer = "<name>"}
#' }
#'
#' @param file Character. Full path to the \code{.h5ad} file, including the
#'   file name and extension (e.g. \code{"/data/samples/sample1.h5ad"}).
#' @param layer Character. Which data to extract. Default \code{"counts"}
#'   (raw UMI counts stored in \code{adata.layers[["counts"]]}, recommended
#'   for NICHESv2). Use \code{"X"} to extract \code{adata.X} directly. A
#'   hard error with available layer names is thrown if the layer is absent.
#' @param obs.cols Character vector or \code{NULL}. Subset of
#'   \code{adata.obs} columns to keep. \code{NULL} (default) keeps all.
#' @param verbose Logical. If \code{TRUE} (default), prints a summary and a
#'   spatial coordinate reminder.
#'
#' @return A named list:
#'   \describe{
#'     \item{\code{count.mtx}}{Sparse \code{dgCMatrix}, genes \eqn{\times}
#'       cells. Row names are gene names; column names are cell barcodes.}
#'     \item{\code{meta.data}}{Data frame from \code{adata.obs} (optionally
#'       column-subsetted). Row names are cell barcodes.}
#'   }
#'
#' @details
#' Install \pkg{rhdf5} with:
#' \preformatted{
#'   BiocManager::install("rhdf5")
#' }
#'
#' The function does \strong{not} fetch spatial coordinates from
#' \code{adata.obsm}. For spatial mode, add \code{x} and \code{y} columns
#' to \code{meta.data} yourself:
#' \preformatted{
#'   library(rhdf5)
#'   coords <- t(rhdf5::h5read(file, "obsm/spatial"))
#'   inputs$meta.data$x <- coords[, 1]
#'   inputs$meta.data$y <- coords[, 2]
#' }
#'
#' @examples
#' \dontrun{
#'   inputs <- extract_NICHESInputs_AnnData(
#'     file  = "/path/to/sample.h5ad",
#'     layer = "counts"
#'   )
#'   niches.obj <- create_NICHESObject(
#'     count.mtx = inputs$count.mtx,
#'     meta.data  = inputs$meta.data,
#'     LRM.db     = lrm.db,
#'     mode       = "spatial",
#'     k          = 6
#'   )
#' }
#'
#' @seealso \code{\link{create_NICHESObject}},
#'   \code{\link{extract_NICHESInputs_Seurat}}
#'
#' @export
extract_NICHESInputs_AnnData <- function(file,
                                          layer    = "counts",
                                          obs.cols = NULL,
                                          verbose  = TRUE) {

  # -- dependency check ---------------------------------------------------------
  if (!requireNamespace("rhdf5", quietly = TRUE)) {
    stop(
      "Package 'rhdf5' is required but is not installed.\n",
      "  Install with:\n",
      "    BiocManager::install(\"rhdf5\")"
    )
  }

  # -- validate file -------------------------------------------------------------
  if (!file.exists(file)) {
    stop("File not found: ", file)
  }
  if (!grepl("\\.h5ad$", file, ignore.case = TRUE)) {
    warning(
      "File '", basename(file), "' does not have a .h5ad extension. ",
      "Proceeding, but verify this is a valid AnnData file."
    )
  }

  # -- read dimension names (used for orientation check and naming) -----------
  obs.names <- .h5ad_obs_names(file)   # cell barcodes
  var.names <- .h5ad_var_names(file)   # gene names
  n.obs     <- length(obs.names)
  n.vars    <- length(var.names)

  # -- validate requested layer -------------------------------------------------
  available.layers <- .h5ad_available_layers(file)
  if (!layer %in% available.layers) {
    stop(
      "Layer '", layer, "' not found in '", basename(file), "'.\n",
      "  Available layers: ", paste(available.layers, collapse = ", "),
      "\n  To use adata.X directly, set layer = \"X\"."
    )
  }

  # -- read matrix -> genes x cells ----------------------------------------------
  mat.path  <- if (layer == "X") "X" else paste0("layers/", layer)
  count.mtx <- .h5ad_read_matrix(file, mat.path, n.obs, n.vars)
  rownames(count.mtx) <- var.names
  colnames(count.mtx) <- obs.names

  # -- read cell metadata --------------------------------------------------------
  meta.data <- .h5ad_read_obs(file, obs.names)

  if (!is.null(obs.cols)) {
    missing.cols <- setdiff(obs.cols, colnames(meta.data))
    if (length(missing.cols) > 0L) {
      stop(
        "Requested obs.cols not found in adata.obs: ",
        paste(missing.cols, collapse = ", "), "\n",
        "  Available columns: ", paste(colnames(meta.data), collapse = ", ")
      )
    }
    meta.data <- meta.data[, obs.cols, drop = FALSE]
  }

  # -- messages -----------------------------------------------------------------
  if (verbose) {
    n.cols <- ncol(meta.data)
    shown  <- head(colnames(meta.data), 5L)
    extra  <- if (n.cols > 5L) paste0(", ... (", n.cols - 5L, " more)") else ""
    message(
      "[NICHESv2] Extracted ", n.vars, " genes x ", n.obs, " cells",
      " from '", basename(file), "', layer '", layer, "'."
    )
    message(
      "[NICHESv2] Cell metadata: ", n.cols, " column(s): ",
      paste(shown, collapse = ", "), extra
    )
    message(
      "[NICHESv2] NOTE: For spatial mode, create_NICHESObject() requires\n",
      "  columns named 'x' and 'y' in meta.data. If your spatial coordinates\n",
      "  are in adata.obsm[[\"spatial\"]], extract them with:\n",
      "    library(rhdf5)\n",
      "    coords <- t(rhdf5::h5read(file, \"obsm/spatial\"))\n",
      "  then add them to inputs$meta.data before calling create_NICHESObject()."
    )
  }

  list(count.mtx = count.mtx,
       meta.data  = meta.data)
}
