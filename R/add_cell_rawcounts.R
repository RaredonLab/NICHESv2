# add_cell_rawcounts.R --------------------------------------------------------
# Add a raw UMI count column to $cell.data of a NICHESObject from a
# user-supplied gene x cell count matrix.
#
# Row expansion (Option B): new rows are added to $cell.data for (cell, gene)
# pairs where count > 0 but score = 0. These rows carry score = 0,
# score.norm = 0, and the observed count value.
#
# FUNCTIONS:
#   add_cell_rawcounts()  -- exported
#
# DEPS: data.table (>= 1.14), Matrix (>= 1.5)


# ==============================================================================
#' Add Raw UMI Counts to \code{$cell.data}
#'
#' @description
#' Adds a \code{count} column to \code{niches.obj$cell.data} from a
#' user-supplied gene x cell raw count matrix. The \code{count} column records
#' the original sequencing depth for each (cell, gene) pair as a provenance
#' record; it is not used by any other NICHESv2 function.
#'
#' @section Row expansion:
#' \code{$cell.data} is sparse: it originally contains only (cell, gene) pairs
#' where \code{score > 0}. When raw counts are provided, new rows are appended
#' for any (cell, gene) pair where \code{count > 0} but \code{score = 0}
#' (e.g., a gene detected in sequencing but not expressed in the scoring
#' layer). New rows carry \code{score = 0} and \code{score.norm = 0}.
#'
#' @section Gene coverage:
#' \code{count.mtx} must cover every cell in \code{niches.obj$cell.list}
#' (hard error if any are missing). Genes in \code{niches.obj$gene.list}
#' absent from \code{rownames(count.mtx)} receive \code{count = NA} with a
#' message when \code{verbose = TRUE}. Extra rows in \code{count.mtx} not in
#' \code{$gene.list} are silently ignored.
#'
#' @section Validation:
#' Non-zero stored values are checked for validity. Warnings (not errors) are
#' emitted if any values are negative or non-integer; the pipeline continues
#' regardless. These checks help detect accidental use of a normalized layer
#' instead of raw counts.
#'
#' @param niches.obj A \code{NICHESObject}.
#' @param count.mtx  A gene x cell raw count matrix. Either a base
#'   \code{matrix} or a sparse \code{Matrix} object (e.g., \code{dgCMatrix}).
#'   Must have \code{rownames} (gene symbols) and \code{colnames} (cell
#'   barcodes) covering every cell in \code{niches.obj$cell.list}.
#' @param verbose    Logical. If \code{TRUE} (default), emits messages for
#'   missing genes, new rows added, and an overwrite notice if a \code{count}
#'   column already exists.
#'
#' @return Updated \code{NICHESObject} with a \code{count} column in
#'   \code{$cell.data}. Row count may increase if any (cell, gene) pairs had
#'   \code{count > 0} and \code{score = 0}.
#'
#' @seealso \code{\link{create_NICHESObject}} for object construction;
#'   \code{\link{extract_NICHESInputs_Seurat}} and
#'   \code{\link{extract_NICHESInputs_AnnData}} for obtaining a count matrix
#'   from Seurat or AnnData objects.
#'
#' @examples
#' \donttest{
#' data(NICHESv2_inputs)
#' data(NICHESv2_example)
#'
#' # NICHESv2_inputs$count.mtx is the raw count matrix used to build
#' # NICHESv2_example. Attach it as the provenance count column.
#' obj <- add_cell_rawcounts(NICHESv2_example, NICHESv2_inputs$count.mtx)
#' head(obj$cell.data[, c("cell", "gene", "score", "score.norm", "count")])
#' }
#'
#' @export
#' @import data.table
add_cell_rawcounts <- function(niches.obj, count.mtx, verbose = TRUE) {

  ## Section 1: input guards --------------------------------------------------

  if (!inherits(niches.obj, "NICHESObject"))
    stop("[add_cell_rawcounts] 'niches.obj' must be a NICHESObject.")

  if (!is.matrix(count.mtx) && !inherits(count.mtx, "Matrix"))
    stop(
      "[add_cell_rawcounts] 'count.mtx' must be a base matrix or a sparse ",
      "Matrix object (e.g., dgCMatrix)."
    )

  if (is.null(rownames(count.mtx)) || is.null(colnames(count.mtx)))
    stop(
      "[add_cell_rawcounts] 'count.mtx' must have rownames (genes) and ",
      "colnames (cells)."
    )

  ## Section 2: cell coverage -------------------------------------------------

  missing.cells <- setdiff(niches.obj$cell.list, colnames(count.mtx))
  if (length(missing.cells) > 0L)
    stop(sprintf(
      "[add_cell_rawcounts] 'count.mtx' is missing %d cell(s) from $cell.list.\n  First missing (up to 3): %s",
      length(missing.cells),
      paste(head(missing.cells, 3L), collapse = ", ")
    ))

  ## Section 3: gene coverage -------------------------------------------------

  genes.available <- intersect(niches.obj$gene.list, rownames(count.mtx))
  genes.missing   <- setdiff(niches.obj$gene.list,  rownames(count.mtx))

  if (length(genes.missing) > 0L && verbose)
    message(sprintf(
      "[add_cell_rawcounts] %d gene(s) in $gene.list are absent from 'count.mtx' -- their 'count' will be NA.\n  First absent (up to 3): %s",
      length(genes.missing),
      paste(head(genes.missing, 3L), collapse = ", ")
    ))

  ## Section 4: subset to relevant cells and genes ----------------------------
  # Align to $cell.list column order; extra cells in count.mtx are dropped.

  sub.mtx <- count.mtx[genes.available, niches.obj$cell.list, drop = FALSE]

  ## Section 5: non-negative integer validation on stored non-zero values -----

  if (inherits(sub.mtx, "Matrix")) {
    stored.vals <- sub.mtx@x
  } else {
    stored.vals <- sub.mtx[sub.mtx != 0]
  }

  if (length(stored.vals) > 0L) {
    if (any(stored.vals < 0))
      warning(
        "[add_cell_rawcounts] 'count.mtx' contains negative values. ",
        "Raw UMI counts should be non-negative. ",
        "Please verify the correct matrix was supplied."
      )
    if (any(stored.vals != floor(stored.vals)))
      warning(
        "[add_cell_rawcounts] 'count.mtx' contains non-integer values. ",
        "Raw UMI counts should be whole numbers. ",
        "Please verify the correct matrix was supplied ",
        "(e.g., not a normalized layer)."
      )
  }

  ## Section 6: existing count column notice ----------------------------------

  if ("count" %in% colnames(niches.obj$cell.data) && verbose)
    message(
      "[add_cell_rawcounts] 'count' column already exists in $cell.data ",
      "and will be overwritten."
    )

  ## Section 7: convert to long-format data.table of non-zero entries ---------
  #
  # Extract directly from dgCMatrix internals for reliability across Matrix
  # versions:
  #   @p  -- column pointer vector (length = ncol + 1); diff(@p) = nnz per col
  #   @i  -- 0-based row indices of stored non-zero entries
  #   @x  -- stored non-zero values
  #
  # This avoids coercion-related deprecation warnings in Matrix >= 1.5.

  if (!inherits(sub.mtx, "Matrix"))
    sub.mtx <- Matrix::Matrix(sub.mtx, sparse = TRUE)
  sub.mtx <- as(sub.mtx, "dgCMatrix")

  n.col   <- ncol(sub.mtx)
  col.idx <- rep(seq_len(n.col), diff(sub.mtx@p))   # 1-based column index per nnz entry
  row.idx <- sub.mtx@i + 1L                           # 0-based to 1-based row index

  count.long <- data.table(
    cell  = colnames(sub.mtx)[col.idx],
    gene  = rownames(sub.mtx)[row.idx],
    count = as.numeric(sub.mtx@x)
  )
  setkeyv(count.long, c("cell", "gene"))

  ## Section 8: add count to existing $cell.data rows ------------------------
  #
  # copy() is required to break data.table reference semantics before :=
  # assignment; without it, modifying the copy would silently mutate the
  # original object in the calling environment.

  niches.obj$cell.data         <- copy(niches.obj$cell.data)
  niches.obj$cell.data[, count := NA_real_]

  # Non-zero count entries: join onto matching (cell, gene) rows
  niches.obj$cell.data[count.long, count := i.count, on = c("cell", "gene")]

  # Genes present in count.mtx but with stored count = 0 (not in count.long):
  # set to 0 explicitly to distinguish from genes absent from count.mtx (NA).
  niches.obj$cell.data[
    is.na(count) & gene %in% genes.available,
    count := 0
  ]
  # Rows where gene is in genes.missing retain count = NA.

  ## Section 9: Option B row expansion ----------------------------------------
  #
  # Add new rows for (cell, gene) pairs where count > 0 but score = 0.
  # These pairs are in count.long but not yet in $cell.data.

  existing.key <- niches.obj$cell.data[, .(cell, gene)]
  new.rows     <- count.long[!existing.key, on = c("cell", "gene")]

  if (nrow(new.rows) > 0L) {
    new.rows[, score      := 0]
    new.rows[, score.norm := 0]
    # Align column order to match $cell.data before binding
    new.rows <- new.rows[, colnames(niches.obj$cell.data), with = FALSE]
    niches.obj$cell.data <- rbindlist(
      list(niches.obj$cell.data, new.rows),
      use.names = TRUE
    )
    if (verbose)
      message(sprintf(
        "[add_cell_rawcounts] Added %d new (cell, gene) row(s) to $cell.data where count > 0 and score = 0.",
        nrow(new.rows)
      ))
  }

  setkeyv(niches.obj$cell.data, c("cell", "gene"))
  niches.obj
}
