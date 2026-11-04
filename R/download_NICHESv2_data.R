# download_NICHESv2_data.R ------------------------------------------------
# Download one of the public datasets used by the website tutorials from a
# release of the NICHESv2 GitHub repository, check it, cache it, and return
# inputs ready for create_NICHESObject(). The files are built by
# data-raw/build_tutorial_data.R.
#
# FUNCTIONS:
#   download_NICHESv2_data()  -- main function
#   .read_inputs()            -- read inputs.rds with Matrix loaded
#   .niches_data_url()        -- release folder the files are attached to
#   .niches_datasets()        -- file name, MD5 sum and size of each dataset
#
# DEPS: tools, utils


# ==============================================================================
#' Download a Public Tutorial Dataset
#'
#' @description
#' Downloads one of the two public datasets used in the website tutorials,
#' checks its MD5 sum, and returns a count matrix and metadata ready for
#' \code{create_NICHESObject()}. The files are attached to a release of the
#' NICHESv2 GitHub repository, so no other package is needed. Each dataset is
#' downloaded once and then read from \code{cache.dir}.
#'
#' @section Datasets:
#' \describe{
#'   \item{\code{"stxBrain"}}{10x Genomics Visium mouse brain, sagittal-anterior
#'     section 1: 2,696 spots and 31,053 genes, with spot coordinates in
#'     \code{x} and \code{y}. About 14 MB. CC BY 4.0, credit 10x Genomics.}
#'   \item{\code{"ifnb"}}{Human PBMCs, untreated or stimulated with
#'     interferon-beta (Kang et al., 2018): 13,999 cells and 14,053 genes, with
#'     the condition in \code{stim} and cell-type labels from the Satija Lab in
#'     \code{seurat_annotations}. About 11 MB. CC BY 4.0.}
#' }
#' Each dataset's folder in \code{cache.dir} also holds a \code{PROVENANCE.md}
#' file stating its source, license and the changes made.
#'
#' @param dataset   \code{"stxBrain"} or \code{"ifnb"}.
#' @param cache.dir Folder to keep downloaded files in. Default: the NICHESv2
#'   user cache folder, \code{tools::R_user_dir("NICHESv2", "cache")}. Delete
#'   the folder to free the space.
#' @param force     Download again even if the dataset is already cached.
#'   Default \code{FALSE}.
#' @param verbose   Print progress messages. Default \code{TRUE}.
#'
#' @return A list with \code{count.mtx} (a sparse genes x cells matrix of raw
#'   counts) and \code{meta.data} (a data frame with one row per cell), the
#'   same shape that \code{extract_NICHESInputs_Seurat()} returns. If the
#'   download fails, or the file doesn't match its MD5 sum, a message explains
#'   why and \code{NULL} is returned invisibly.
#'
#' @references
#' Kang HM, Subramaniam M, Targ S, et al. (2018). Multiplexed droplet
#' single-cell RNA-sequencing using natural genetic variation.
#' \emph{Nature Biotechnology}, 36(1), 89--94. \doi{10.1038/nbt.4042}
#'
#' @examples
#' \donttest{
#' # Downloads about 11 MB, kept in a temporary folder here.
#' inputs <- download_NICHESv2_data("ifnb", cache.dir = tempdir())
#' if (!is.null(inputs)) dim(inputs$count.mtx)
#' }
#'
#' @export
download_NICHESv2_data <- function(dataset,
                                   cache.dir = tools::R_user_dir("NICHESv2", "cache"),
                                   force     = FALSE,
                                   verbose   = TRUE) {

  datasets <- .niches_datasets()
  if (missing(dataset) || !is.character(dataset) || length(dataset) != 1L ||
      !dataset %in% names(datasets))
    stop(
      "[download_NICHESv2_data] 'dataset' must be one of: ",
      paste0("'", names(datasets), "'", collapse = ", "), ".",
      call. = FALSE
    )

  info     <- datasets[[dataset]]
  rds.file <- file.path(cache.dir, dataset, "inputs.rds")

  if (!force && file.exists(rds.file)) {
    if (verbose)
      message("[download_NICHESv2_data] Reading cached '", dataset, "' from ",
              dirname(rds.file))
    return(.read_inputs(rds.file))
  }

  # -- Download and check the zip ---------------------------------------------
  dir.create(cache.dir, recursive = TRUE, showWarnings = FALSE)
  zip.file <- file.path(cache.dir, info$file)
  url      <- paste(.niches_data_url(), info$file, sep = "/")

  if (verbose)
    message("[download_NICHESv2_data] Downloading '", dataset, "' (",
            info$size, ") ...")
  old.options <- options(timeout = max(600, getOption("timeout")))
  on.exit(options(old.options), add = TRUE)

  downloaded <- tryCatch(
    utils::download.file(url, zip.file, mode = "wb", quiet = TRUE) == 0L,
    error   = function(e) FALSE,
    warning = function(w) FALSE
  )
  if (!downloaded || !file.exists(zip.file)) {
    unlink(zip.file)
    message(
      "[download_NICHESv2_data] Could not download '", dataset, "' from\n  ",
      url, "\n  Check your internet connection and try again."
    )
    return(invisible(NULL))
  }

  if (!identical(unname(tools::md5sum(zip.file)), info$md5)) {
    unlink(zip.file)
    message(
      "[download_NICHESv2_data] The downloaded file for '", dataset,
      "' doesn't match its expected MD5 sum, so it was deleted.\n  ",
      "Try again with force = TRUE. If this keeps happening, please open an ",
      "issue at https://github.com/RaredonLab/NICHESv2/issues"
    )
    return(invisible(NULL))
  }

  # -- Unpack into cache.dir/<dataset>/ ---------------------------------------
  unlink(file.path(cache.dir, dataset), recursive = TRUE)
  utils::unzip(zip.file, exdir = cache.dir)
  unlink(zip.file)

  .read_inputs(rds.file)
}


## INTERNAL HELPERS ----

#' Read an inputs.rds file. Matrix is loaded first: otherwise readRDS() loads it
#' while unserializing, and the first dim() call on the matrix returns NULL.
#' @noRd
.read_inputs <- function(rds.file) {
  loadNamespace("Matrix")
  readRDS(rds.file)
}

#' Release folder that holds the dataset files. The NICHESv2.data.url option
#' points elsewhere for tests and for building the site before an upload.
#' @noRd
.niches_data_url <- function() {
  getOption(
    "NICHESv2.data.url",
    "https://github.com/RaredonLab/NICHESv2/releases/download/data-v1"
  )
}

#' File name, MD5 sum and download size of each dataset. The sums come from
#' data-raw/build_tutorial_data.R and must match the uploaded files.
#' @noRd
.niches_datasets <- function() {
  list(
    stxBrain = list(file = "stxBrain.zip",
                    md5  = "f2085cbb105198f29e94456eeebdf261",
                    size = "14 MB"),
    ifnb     = list(file = "ifnb.zip",
                    md5  = "bd4d47cac9dc110aa253f4c6e95fcc43",
                    size = "11 MB")
  )
}
