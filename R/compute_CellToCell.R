# compute_CellToCell.R ----------------------------------------------------
# Universal cell-to-cell ligand-receptor signal scoring engine.
# Accepts any edge list produced by create_Edgelist_Spatial() or
# create_Edgelist_SamplingByType().
# Supports single-gene and multi-gene complex ligands / receptors
# (underscore-delimited, e.g. "COL1A1_COL1A2"). Returns a plain
# 2-element list suitable for assembly into a NICHESObject via
# create_NICHESObject().
#
# FUNCTIONS:
#   .sparse_emin()        (internal) element-wise min of sparse matrices
#   compute_CellToCell()  main function
#
# OUTPUT STRUCTURE (plain named list):
#   $edges   - data.table, one row per non-zero (edge x LRM)
#   $lr.meta - data.table, one row per valid LRM pair
#
# DEPS: Matrix (>= 1.5), data.table (>= 1.14), parallel


## INTERNAL: Element-wise minimum of two non-negative sparse matrices ----

#' @noRd
.sparse_emin <- function(A, B) {

  stopifnot(identical(dim(A), dim(B)))

  At <- as(A, "TsparseMatrix")
  Bt <- as(B, "TsparseMatrix")

  dtA <- data.table(i = At@i, j = At@j, vA = At@x, key = c("i", "j"))
  dtB <- data.table(i = Bt@i, j = Bt@j, vB = Bt@x, key = c("i", "j"))

  common <- dtA[dtB, nomatch = 0L]

  if (nrow(common) == 0L)
    return(Matrix::sparseMatrix(
      i = integer(0), j = integer(0), x = numeric(0),
      dims = dim(A), dimnames = dimnames(A)
    ))

  Matrix::sparseMatrix(
    i        = common$i + 1L,
    j        = common$j + 1L,
    x        = pmin(common$vA, common$vB),
    dims     = dim(A),
    dimnames = dimnames(A)
  )
}


# ==============================================================================
#' Compute Cell-to-Cell Ligand-Receptor Signals
#'
#' @description
#' Universal scoring engine for all NICHESv2 edge modes. For each edge and each
#' valid LR pair, computes a signal score. The edge list may come from either
#' edgelist constructor -- spatial or type-sampled -- and the count matrix
#' column names must be the cell barcodes it uses.
#'
#' Supports \strong{multi-gene complexes}: a ligand or receptor represented as
#' \code{"gene1_gene2"} is split on \code{"_"} and all component expressions
#' are combined via the chosen method (product or minimum) \emph{across all
#' ligand and receptor components together}.
#'
#' Built-in LRM databases can be loaded by name (see
#' \code{\link{load_LRM_database}}):
#' \itemize{
#'   \item \code{"fantom5"}  -- Ramilowski et al. 2015. Bundled with
#'     NICHESv2 (human, mouse, rat, pig). No internet or extra packages
#'     required.
#'   \item \code{"omnipath"} -- OmniPath intercell network (Turei et al.
#'     2021). Requires the \pkg{OmnipathR} package and internet access.
#' }
#'
#' @section Complex scoring example:
#' LRM \code{"COL1A1_COL1A2 -> ITGB1"} with \code{method = "product"}:
#' \deqn{score = \text{expr}_{send}(COL1A1) \times \text{expr}_{send}(COL1A2) \times \text{expr}_{recv}(ITGB1)}
#' With \code{method = "minimum"}:
#' \deqn{score = \min(\text{expr}_{send}(COL1A1),\ \text{expr}_{send}(COL1A2),\ \text{expr}_{recv}(ITGB1))}
#'
#' @section Direct vs wrapped use:
#' Most users should call \code{\link{create_NICHESObject}}, which invokes
#' \code{compute_CellToCell()} internally and assembles the result into a
#' full \code{NICHESObject}. Call this function directly when you only need
#' the raw \code{(edge, LRM, score)} table -- for example, to score a
#' custom edge list (e.g.\ a CCC graph from another method) without
#' building the surrounding object, or to integrate the score table into
#' a non-NICHESv2 downstream pipeline.
#'
#' @section Parallelization:
#' Work is partitioned across LRM pairs into \code{n.cores} chunks and
#' processed by a socket cluster (\code{parallel::makeCluster}), which
#' works identically on Windows, macOS, and Linux. Each worker has its
#' \code{data.table} thread count throttled to 1 to prevent OpenMP
#' contention with the outer parallelism. The cluster is torn down via
#' \code{on.exit(stopCluster(cl))} so it is released even if scoring
#' errors midway.
#'
#' @param count.mtx    Gene x cell expression matrix (sparse \code{dgCMatrix}
#'   or dense). Row names = gene symbols; column names = cell barcodes matching
#'   \code{edge.list$from} / \code{$to}.
#' @param edge.list    \code{data.table} from
#'   \code{\link{create_Edgelist_Spatial}} or
#'   \code{\link{create_Edgelist_SamplingByType}}; or any \code{data.frame}
#'   with \code{from} and \code{to} columns whose values match
#'   \code{colnames(count.mtx)}. Coordinate columns
#'   (\code{x.from}, \code{y.from}, \code{x.to}, \code{y.to}), when
#'   present, are propagated to the output.
#' @param LRM.db       LR pair table (\code{data.frame} or \code{data.table})
#'   \strong{or} a character string (\code{"fantom5"} or \code{"omnipath"})
#'   loaded via \code{\link{load_LRM_database}}. When a string is passed,
#'   \code{ligand.col} and \code{receptor.col} are silently overridden to
#'   the canonical \code{"ligand"} and \code{"receptor"}.
#' @param ligand.col   Character. Column name for ligand symbols in
#'   \code{LRM.db}. Default \code{"ligand"}. Ignored when \code{LRM.db} is
#'   a built-in database name.
#' @param receptor.col Character. Column name for receptor symbols in
#'   \code{LRM.db}. Default \code{"receptor"}. Ignored when \code{LRM.db}
#'   is a built-in database name.
#' @param method       Per-edge LR scoring method. \code{"product"}
#'   (default) multiplies sender ligand expression by receiver receptor
#'   expression; \code{"minimum"} takes the elementwise minimum across all
#'   ligand and receptor components. See the Complex scoring section.
#' @param n.cores      Integer. CPU cores for parallel LRM-chunk
#'   processing. Default \code{1L}. See the Parallelization section.
#'
#' @return A plain named list with two elements:
#' \describe{
#'   \item{\code{$edges}}{\code{data.table}: \code{from}, \code{to},
#'     \code{ligand}, \code{receptor}, \code{LRM}
#'     (\code{"ligand|receptor"}), \code{score}. Additional columns
#'     \code{x.from}, \code{y.from}, \code{x.to}, \code{y.to} appear
#'     \strong{only when} the input \code{edge.list} carries coordinates
#'     (i.e.\ for spatial edge lists). One row per non-zero
#'     (edge x LRM) combination -- silent pairs are omitted for
#'     sparsity. Keyed on \code{(from, to, LRM)}.}
#'   \item{\code{$lr.meta}}{\code{data.table}: \code{ligand},
#'     \code{receptor}, \code{LRM}, plus any extra annotation columns
#'     present in \code{LRM.db} (e.g.\ \code{pathway}). One row per valid
#'     LRM pair after gene-availability filtering. Keyed on
#'     \code{(ligand, receptor)}.}
#' }
#' Both elements are \code{data.table}s so downstream joins on
#' \code{from}/\code{to}/\code{LRM} are O(log n).
#'
#' @seealso
#'   \code{\link{create_NICHESObject}} for the wrapped pipeline that
#'     assembles a full \code{NICHESObject};
#'   \code{\link{create_Edgelist_Spatial}},
#'   \code{\link{create_Edgelist_SamplingByType}} for the edge-list
#'     constructors;
#'   \code{\link{load_LRM_database}} for the LR-database loader.
#'
#' @examples
#' \donttest{
#' data(NICHESv2_inputs)
#'
#' # Build a spatial kNN edge list (coord must have row names = barcodes)
#' coord <- data.frame(
#'   x = NICHESv2_inputs$meta.data$x,
#'   y = NICHESv2_inputs$meta.data$y,
#'   row.names = rownames(NICHESv2_inputs$meta.data)
#' )
#' edge.list <- create_Edgelist_Spatial(coord = coord, k = 5L, n.cores = 1L)
#'
#' # Score every edge x LRM pair (returns a plain 2-element list)
#' result <- compute_CellToCell(
#'   count.mtx    = NICHESv2_inputs$count.mtx,
#'   edge.list    = edge.list,
#'   LRM.db       = NICHESv2_inputs$LRM.db,
#'   ligand.col   = "ligand",
#'   receptor.col = "receptor",
#'   method       = "product",
#'   n.cores      = 1L
#' )
#'
#' str(result, max.level = 1)   # $edges and $lr.meta
#' head(result$edges)           # one row per non-zero (edge, LRM)
#' head(result$lr.meta)         # one row per valid LRM pair
#'
#' # Spatial edge lists carry coordinates -- they propagate to $edges
#' setdiff(c("x.from", "y.from", "x.to", "y.to"), colnames(result$edges))
#' # character(0) -- all four coordinate columns are present
#' }
#'
#' @export
compute_CellToCell <- function(count.mtx,
                                edge.list,
                                LRM.db,
                                ligand.col   = "ligand",
                                receptor.col = "receptor",
                                method       = c("product", "minimum"),
                                n.cores      = 1L) {

  method  <- match.arg(method)
  n.cores <- as.integer(n.cores)

  # -- Resolve LRM.db (accept name string or data.frame) --------------------
  if (is.character(LRM.db) && length(LRM.db) == 1L) {
    LRM.db <- load_LRM_database(LRM.db, verbose = TRUE)
    ligand.col   <- "ligand"
    receptor.col <- "receptor"
  }

  if (!is.data.frame(LRM.db) && !is.data.table(LRM.db))
    stop("[compute_CellToCell] 'LRM.db' must be a data.frame or a ",
         "character string ('fantom5' or 'omnipath').")

  ## SECTION 1: Validate & filter LRM pairs ----
  #
  # A pair is valid only when ALL gene components of BOTH ligand and receptor
  # are present as row names in count.mtx.

  ligands   <- as.character(LRM.db[[ligand.col]])
  receptors <- as.character(LRM.db[[receptor.col]])

  if (is.null(ligands) || is.null(receptors))
    stop("[compute_CellToCell] 'ligand.col' or 'receptor.col' not found in LRM.db.")

  genes.in.mtx <- rownames(count.mtx)

  valid.lr <- mapply(
    function(lig, rec) .genes_available(lig, genes.in.mtx) && .genes_available(rec, genes.in.mtx),
    ligands, receptors
  )

  n.dropped <- sum(!valid.lr)
  if (n.dropped > 0L)
    message(sprintf(
      "[compute_CellToCell] Dropping %d LR pair(s): gene(s) absent from count matrix.",
      n.dropped
    ))

  LRM.valid  <- LRM.db[valid.lr, , drop = FALSE]
  ligands    <- ligands[valid.lr]
  receptors  <- receptors[valid.lr]
  n.lr       <- length(ligands)

  if (n.lr == 0L)
    stop("[compute_CellToCell] No valid LR pairs remain after gene filtering.")

  ## SECTION 2: Validate & filter edges ----

  cells.in.mtx   <- colnames(count.mtx)
  valid.edges    <- edge.list$from %in% cells.in.mtx & edge.list$to %in% cells.in.mtx
  n.edge.dropped <- sum(!valid.edges)

  if (n.edge.dropped > 0L)
    message(sprintf(
      "[compute_CellToCell] Dropping %d edge(s): identifiers absent from count matrix.",
      n.edge.dropped
    ))

  edge.list     <- edge.list[valid.edges, ]
  sending.cells <- edge.list$from
  receiv.cells  <- edge.list$to
  n.edges       <- length(sending.cells)

  if (n.edges == 0L)
    stop("[compute_CellToCell] No valid edges remain after identifier filtering.")

  message(sprintf(
    "[compute_CellToCell] %d edges x %d LR pairs | method = '%s' | cores = %d",
    n.edges, n.lr, method, n.cores
  ))

  ## SECTION 3: Parallel computation over LRM chunks (socket cluster) ----
  #
  # LRM pairs are split into n.cores chunks; each worker processes one chunk,
  # calling .compute_LRM_Score() per pair and returning non-zero (edge, LRM,
  # score) rows for its chunk.
  #
  # Internal scoring helpers are exported from the NICHESv2 namespace via
  # asNamespace(), which works both with devtools::load_all() and installed
  # builds. on.exit() guarantees stopCluster() even if an error occurs mid-run.
  # data.table threads are throttled to 1 per worker to prevent OpenMP
  # contention across socket-cluster workers.

  lr.idx.chunks <- if (n.cores == 1L) {
    list(seq_len(n.lr))
  } else {
    split(seq_len(n.lr),
          cut(seq_len(n.lr), breaks = n.cores, labels = FALSE))
  }

  .process_chunk <- function(chunk.idx) {

    chunk.results <- vector("list", length(chunk.idx))

    for (k in seq_along(chunk.idx)) {
      i <- chunk.idx[k]

      score.vec <- .compute_LRM_Score(
        lig.str       = ligands[i],
        rec.str       = receptors[i],
        count.mtx     = count.mtx,
        sending.cells = sending.cells,
        receiv.cells  = receiv.cells,
        method        = method
      )

      nonzero <- which(score.vec != 0)
      if (length(nonzero) == 0L) next

      chunk.results[[k]] <- data.table(
        from     = sending.cells[nonzero],
        to       = receiv.cells[nonzero],
        ligand   = ligands[i],
        receptor = receptors[i],
        score    = score.vec[nonzero]
      )
    }

    filtered <- Filter(Negate(is.null), chunk.results)
    if (length(filtered) == 0L) NULL else rbindlist(filtered)
  }

  if (n.cores > 1L) {
    cl <- parallel::makeCluster(n.cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterEvalQ(cl, {
      library(data.table)
      library(Matrix)
      data.table::setDTthreads(1)
    })
    parallel::clusterExport(cl,
      varlist = c("ligands", "receptors", "count.mtx",
                  "sending.cells", "receiv.cells", "method"),
      envir = environment())
    parallel::clusterExport(cl,
      varlist = c(".compute_LRM_Score", ".parse_complex", ".genes_available"),
      envir = asNamespace("NICHESv2"))
    chunk.results <- parallel::parLapply(cl, lr.idx.chunks, .process_chunk)
  } else {
    chunk.results <- lapply(lr.idx.chunks, .process_chunk)
  }

  ## SECTION 4: Assemble edge data.table ----

  raw.chunks <- Filter(Negate(is.null), chunk.results)
  rm(chunk.results)

  has.coords <- all(c("x.from", "y.from", "x.to", "y.to") %in% colnames(edge.list))

  # When every score is zero, every chunk worker returns NULL and raw.chunks
  # is empty. rbindlist on an empty list produces a zero-column data.table,
  # so column operations like LRM := paste(ligand, ...) would crash because
  # the columns do not exist. We build the correct zero-row schema explicitly.
  if (length(raw.chunks) == 0L) {
    coord.cols <- if (has.coords)
      list(x.from = numeric(0), y.from = numeric(0),
           x.to   = numeric(0), y.to   = numeric(0))
    else list()
    edges <- do.call(data.table, c(
      list(from = character(0), to = character(0),
           ligand = character(0), receptor = character(0),
           LRM = character(0), score = numeric(0)),
      coord.cols
    ))
  } else {
    edges <- rbindlist(raw.chunks)

    # LRM column: "ligand|receptor"
    edges[, LRM := paste(ligand, receptor, sep = "|")]

    # Attach sending/receiving coordinates if present in edge.list (spatial only)
    if (has.coords) {
      coord.lut <- unique(as.data.table(edge.list)[, .(from, to, x.from, y.from, x.to, y.to)])
      setkeyv(coord.lut, c("from", "to"))
      edges <- coord.lut[edges, on = c("from", "to")]
    }
  }

  setkeyv(edges, c("from", "to", "LRM"))
  n.nonzero <- nrow(edges)

  message(sprintf("[compute_CellToCell] Done. %d non-zero edge signals.", n.nonzero))

  ## SECTION 5: Build lr.meta ----

  # lr.meta: one row per valid LRM pair, standardise column names
  lr.meta <- as.data.table(LRM.valid)
  setnames(lr.meta,
           old = c(ligand.col, receptor.col),
           new = c("ligand", "receptor"),
           skip_absent = FALSE)
  lr.meta <- unique(lr.meta, by = c("ligand", "receptor"))
  lr.meta[, LRM := paste(ligand, receptor, sep = "|")]
  setkeyv(lr.meta, c("ligand", "receptor"))

  ## SECTION 6: Return plain 2-element list ----

  list(
    edges   = edges,
    lr.meta = lr.meta
  )
}
