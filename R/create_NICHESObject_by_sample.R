# create_NICHESObject_by_sample.R -----------------------------------------
# Internal branch implementation of create_NICHESObject() for the
# sample.col split/loop/merge path.
#
# This function is NOT exported. It is one of two branch implementations of
# the public create_NICHESObject() dispatcher; the other is
# .create_NICHESObject_single(), which lives in create_NICHESObject_single.R.


## .create_NICHESObject_by_sample() ----
# Unexported. Splits meta.data by sample.col, calls .create_NICHESObject_single
# on each partition sequentially (inner function already parallelises via
# n.cores socket cluster), then merges via merge_NICHESObjects().

#' @noRd
.create_NICHESObject_by_sample <- function(count.mtx,
                                           meta.data,
                                           LRM.db,
                                           mode,
                                           species,
                                           sample.col,
                                           cell.type.col,
                                           ligand.col,
                                           receptor.col,
                                           k,
                                           rad,
                                           seed,
                                           max.cells,
                                           smooth,
                                           smooth.k,
                                           n.pcs,
                                           method,
                                           normalize.method,
                                           gene.meta,
                                           n.cores,
                                           verbose) {

  sample.values <- unique(meta.data[[sample.col]])
  n.samples     <- length(sample.values)

  if (verbose) {
    message(sprintf(
      "[create_NICHESObject] Found %d unique value(s) in '%s': %s",
      n.samples, sample.col,
      paste(sample.values, collapse = ", ")
    ))
    message(sprintf(
      "[create_NICHESObject] Cell-cell connectivity will be computed per '%s'. Edges exist only between cells within the same '%s'.",
      sample.col, sample.col
    ))
  }

  # Sequential outer loop -- .create_NICHESObject_single() already uses n.cores
  obj.list <- vector("list", n.samples)
  names(obj.list) <- as.character(sample.values)

  for (i in seq_along(sample.values)) {
    s <- sample.values[[i]]

    if (verbose)
      message(sprintf("[create_NICHESObject] [%d/%d] Processing %s = \"%s\" ...",
                      i, n.samples, sample.col, s))

    # Subset using rownames -- meta.data stays a plain data.frame throughout
    sample.barcodes  <- rownames(meta.data)[meta.data[[sample.col]] == s]
    sample.meta.data <- meta.data[sample.barcodes, , drop = FALSE]

    if (length(sample.barcodes) < 2L) {
      warning(sprintf(
        "[create_NICHESObject] Sample \"%s\" has fewer than 2 cells -- skipping.",
        s
      ), call. = FALSE)
      obj.list[[i]] <- NULL
      next
    }

    # Subset count matrix to this sample's cells only
    sample.count.mtx <- count.mtx[, sample.barcodes, drop = FALSE]

    obj.list[[i]] <- .create_NICHESObject_single(
      count.mtx        = sample.count.mtx,
      meta.data        = sample.meta.data,
      LRM.db           = LRM.db,
      mode             = mode,
      species          = species,
      cell.type.col    = cell.type.col,
      ligand.col       = ligand.col,
      receptor.col     = receptor.col,
      k                = k,
      rad              = rad,
      seed             = seed,
      max.cells        = max.cells,
      smooth           = smooth,
      smooth.k         = smooth.k,
      n.pcs            = n.pcs,
      method           = method,
      normalize.method = normalize.method,
      gene.meta        = gene.meta,
      n.cores          = n.cores,
      verbose          = FALSE    # suppress per-sample verbosity
    )
  }

  # Drop NULLs from skipped samples
  obj.list <- Filter(Negate(is.null), obj.list)

  if (length(obj.list) == 0L)
    stop("[create_NICHESObject] All samples were skipped (each had fewer than 2 cells).",
         call. = FALSE)

  merged.obj <- merge_NICHESObjects(obj.list, verbose = verbose)
  attr(merged.obj, "sample.col") <- sample.col

  merged.obj
}
