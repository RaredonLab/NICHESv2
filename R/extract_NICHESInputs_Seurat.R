#' Extract inputs for create_NICHESObject from a Seurat object
#'
#' Extracts the count matrix and cell metadata from a Seurat object (V4 or V5
#' internal structure) and returns them as a named list ready to pass directly
#' to \code{\link{create_NICHESObject}}. The object structure is detected from
#' the assay class, not the installed package version, so V4 objects loaded in
#' a V5 environment are handled correctly.
#'
#' @param seurat.obj A Seurat object. Both V4 (\code{Assay}) and V5
#'   (\code{Assay5}) internal structures are supported.
#' @param assay Character. The assay to extract from. Default \code{"RNA"}.
#'   Common alternatives: \code{"SCT"}, \code{"Spatial"}, \code{"Xenium"}.
#' @param layer Character. The data layer (V5) or slot (V4) to extract.
#'   Default \code{"counts"} (raw UMI counts), a starting point rather than a
#'   recommendation: this choice affects which edges score highest.
#'   Other options: \code{"data"} (log-normalized), \code{"scale.data"}
#'   (scaled; V4 only).
#' @param verbose Logical. If \code{TRUE} (default), prints a summary of what
#'   was extracted and a reminder about spatial coordinate requirements.
#'
#' @return A named list with two elements:
#'   \describe{
#'     \item{\code{count.mtx}}{Sparse matrix (genes \eqn{\times} cells),
#'       class \code{dgCMatrix}. Row names are gene names; column names are
#'       cell barcodes.}
#'     \item{\code{meta.data}}{Data frame of cell metadata extracted directly
#'       from \code{seurat.obj@meta.data}, unmodified. Row names are cell
#'       barcodes.}
#'   }
#'
#' @details
#' The function does \strong{not} fetch spatial coordinates from
#' \code{seurat.obj@images}. For spatial mode in \code{create_NICHESObject()},
#' columns named \code{x} and \code{y} must be present in \code{meta.data}.
#' If your coordinates are in \code{seurat.obj@images}, add them to the
#' returned \code{meta.data} yourself before calling
#' \code{create_NICHESObject()}:
#' \preformatted{
#'   coords <- SeuratObject::GetTissueCoordinates(seurat.obj)
#'   inputs$meta.data$x <- coords$imagerow[match(rownames(inputs$meta.data),
#'                                                rownames(coords))]
#'   inputs$meta.data$y <- coords$imagecol[match(rownames(inputs$meta.data),
#'                                                rownames(coords))]
#' }
#'
#' @examples
#' if (requireNamespace("SeuratObject", quietly = TRUE)) {
#'   data(NICHESv2_inputs)
#'   seurat.obj <- SeuratObject::CreateSeuratObject(
#'     counts    = NICHESv2_inputs$count.mtx,
#'     meta.data = NICHESv2_inputs$meta.data
#'   )
#'   inputs <- extract_NICHESInputs_Seurat(seurat.obj, assay = "RNA",
#'                                         layer = "counts", verbose = FALSE)
#'   niches.obj <- create_NICHESObject(
#'     count.mtx = inputs$count.mtx,
#'     meta.data = inputs$meta.data,
#'     LRM.db    = NICHESv2_inputs$LRM.db,
#'     mode      = "spatial",
#'     k         = 5L,
#'     verbose   = FALSE
#'   )
#'   niches.obj
#' }
#'
#' @seealso \code{\link{create_NICHESObject}},
#'   \code{\link{extract_NICHESInputs_AnnData}}
#'
#' @export
extract_NICHESInputs_Seurat <- function(seurat.obj,
                                         assay   = "RNA",
                                         layer   = "counts",
                                         verbose = TRUE) {

  # -- dependency check --------------------------------------------------------
  if (!requireNamespace("SeuratObject", quietly = TRUE)) {
    stop(
      "Package 'SeuratObject' is required but is not installed.\n",
      "  Install with: install.packages('SeuratObject')"
    )
  }

  # -- validate assay -----------------------------------------------------------
  available.assays <- names(seurat.obj@assays)
  if (!assay %in% available.assays) {
    stop(
      "Assay '", assay, "' not found in seurat.obj.\n",
      "  Available assays: ", paste(available.assays, collapse = ", ")
    )
  }

  assay.obj <- seurat.obj@assays[[assay]]

  # -- detect object structure from assay class, not package version ------------
  # A V4 Seurat object loaded under SeuratObject >= 5 keeps its Assay class.
  # A V5 object uses Assay5. Always check the object, not the package.
  is.v5 <- inherits(assay.obj, "Assay5")

  # -- extract count matrix -----------------------------------------------------
  if (is.v5) {
    available.layers <- SeuratObject::Layers(assay.obj)
    if (!layer %in% available.layers) {
      stop(
        "Layer '", layer, "' not found in assay '", assay, "' (V5 Assay5).\n",
        "  Available layers: ", paste(available.layers, collapse = ", ")
      )
    }
    count.mtx <- SeuratObject::GetAssayData(seurat.obj, layer = layer,
                                             assay = assay)
  } else {
    # V4 Assay: fixed slot names. Read the slot directly, because the 'slot'
    # argument of GetAssayData() is deprecated and SeuratObject < 5 has no
    # 'layer' argument. Slots that were never populated hold an empty matrix
    # (e.g. scale.data before ScaleData()).
    available.slots <- c("counts", "data", "scale.data")
    if (!layer %in% available.slots) {
      stop(
        "Slot '", layer, "' not found in assay '", assay, "' (V4 Assay).\n",
        "  Available slots: ", paste(available.slots, collapse = ", ")
      )
    }
    count.mtx <- methods::slot(assay.obj, layer)
    if (prod(dim(count.mtx)) == 0L) {
      stop(
        "Slot '", layer, "' in assay '", assay, "' (V4 Assay) is empty.\n",
        "  'counts' is populated at object creation.\n",
        "  'data' requires NormalizeData(); 'scale.data' requires ScaleData()."
      )
    }
  }

  # -- coerce to dgCMatrix (defensive; GetAssayData usually returns dgCMatrix) -
  if (!inherits(count.mtx, "dgCMatrix")) {
    count.mtx <- methods::as(count.mtx, "dgCMatrix")
  }

  # -- extract meta.data exactly as stored; no column additions -----------------
  meta.data <- seurat.obj@meta.data

  # -- messages -----------------------------------------------------------------
  if (verbose) {
    n.genes <- nrow(count.mtx)
    n.cells <- ncol(count.mtx)
    n.cols  <- ncol(meta.data)
    shown   <- head(colnames(meta.data), 5L)
    extra   <- if (n.cols > 5L) paste0(", ... (", n.cols - 5L, " more)") else ""

    message(
      "[NICHESv2] Extracted ", n.genes, " genes x ", n.cells, " cells",
      " from assay '", assay, "', ",
      if (is.v5) "layer '" else "slot '", layer, "' ",
      "(", if (is.v5) "V5 Assay5" else "V4 Assay", ")."
    )
    message(
      "[NICHESv2] Scoring representation is set by this ",
      if (is.v5) "layer" else "slot", ": '", layer, "' is passed to\n",
      "  create_NICHESObject() as-is. Raw counts (e.g. 'counts') give ",
      "absolute-magnitude\n",
      "  scores; normalized expression (e.g. 'data') gives relative scores. ",
      "The choice\n",
      "  affects which edges score highest; re-extract with a different '",
      if (is.v5) "layer" else "slot", "' to change it."
    )
    message(
      "[NICHESv2] Cell metadata: ", n.cols, " column(s): ",
      paste(shown, collapse = ", "), extra
    )
    message(
      "[NICHESv2] NOTE: For spatial mode, create_NICHESObject() requires\n",
      "  columns named 'x' and 'y' in meta.data. If your spatial coordinates\n",
      "  are stored in seurat.obj@images, extract them with:\n",
      "    coords <- SeuratObject::GetTissueCoordinates(seurat.obj)\n",
      "  then add them to inputs$meta.data before calling create_NICHESObject()."
    )
  }

  return(list(count.mtx = count.mtx,
              meta.data  = meta.data))
}
