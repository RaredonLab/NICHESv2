# create_Edgelist_SamplingByType.R ----------------------------------------
# Build a cell-barcode-level directed edge list for non-spatial
# single-cell communication analysis by stratified type-pair sampling.
# For each ordered cell type pair (i -> j), independently samples up
# to min(n_i, n_j) cell barcodes from each side and pairs them
# positionally, following the NICHES v1 sampling strategy.
#
# FUNCTIONS:
#   create_Edgelist_SamplingByType()  -- main function
#
# DEPS: data.table (>= 1.14)


# ==============================================================================
#' Build a Type-Stratified Sampled Edge List for Non-Spatial scRNA Data
#'
#' @description
#' Constructs a directed edge list for non-spatial single-cell communication
#' analysis by randomly sampling cell pairs within each ordered cell type
#' crossing. For every directed type pair (A -> B), the function:
#' \enumerate{
#'   \item Computes a per-crossing cap: \code{min(n_A, n_B)} (further bounded
#'     by \code{max.cells} if supplied).
#'   \item Independently samples \code{cap} cell barcodes from type A (senders)
#'     and \code{cap} barcodes from type B (receivers) without replacement.
#'   \item Pairs them positionally -- the k-th sampled sender is paired with
#'     the k-th sampled receiver.
#' }
#' All T^2 directed crossings are covered, including autocrine (A -> A), where
#' both sender and receiver are sampled from the same type.
#'
#' \strong{Important:} Individual edges are biologically synthetic -- the
#' specific cell pairing within a type crossing is arbitrary. Only aggregate
#' signals across all edges of a given type pair carry biological meaning.
#' Results differ across runs unless \code{seed} is fixed.
#'
#' @param meta.data    \code{data.frame} with one row per cell. Row names must
#'   be unique cell barcodes. Must contain the column specified by
#'   \code{cell.type.col}.
#' @param cell.type.col \code{character(1)}. Name of the column in
#'   \code{meta.data} that holds cell type labels.
#' @param seed         \code{integer(1)}. Random seed for reproducibility.
#'   Default \code{42L}. Pass a different value to obtain an alternative
#'   sampling draw; pass \code{NULL} to skip seed-setting (not recommended
#'   for publication). The caller's random number state is restored on exit.
#' @param max.cells    \code{integer(1)} or \code{NULL}. Optional upper cap
#'   on the number of edges per type-pair crossing. When \code{NULL} (default),
#'   the cap is \code{min(n_i, n_j)} for each crossing. When set, the cap
#'   becomes \code{min(n_i, n_j, max.cells)}, preventing very large types from
#'   producing an excessive number of edges even when paired with another large
#'   type.
#'
#' @return A \code{data.table} keyed on \code{(from, to)} with columns:
#' \describe{
#'   \item{\code{from}}{Sending cell barcode (character).}
#'   \item{\code{to}}{Receiving cell barcode (character).}
#' }
#' Carries metadata attributes: \code{mode} (\code{"nonspatial"}), \code{param}
#' (named list of \code{cell.type.col}, \code{seed}, \code{max.cells}),
#' \code{n.cells} (total cells in \code{meta.data}), \code{cell.names} (all
#' cell barcodes).
#'
#' @details
#' This function is a low-level primitive: it builds the sampled
#' non-spatial edge list that \code{\link{create_NICHESObject}} consumes
#' internally when called with \code{mode = "nonspatial"}. Most users should
#' call \code{create_NICHESObject(mode = "nonspatial", ...)}, which invokes
#' this function and assembles the result into a full \code{NICHESObject}
#' in one step. \code{create_Edgelist_SamplingByType()} is exposed for
#' power users who want to inspect the sampled edges, fix the random
#' seed explicitly outside the wrapper, or feed the edges into a custom
#' scoring pipeline.
#'
#' Cell types with zero cells after barcode-matrix intersection are
#' silently skipped. Type crossings where the effective cap is zero
#' (e.g., one type has no cells after filtering) are also silently
#' skipped.
#'
#' @examples
#' set.seed(1)
#' meta <- data.frame(
#'   cell.type = sample(c("ATI", "gCap", "Macro"), 300, replace = TRUE),
#'   row.names = paste0("cell_", 1:300)
#' )
#' el <- create_Edgelist_SamplingByType(meta, cell.type.col = "cell.type",
#'                                      seed = 42L, max.cells = 50L)
#'
#' @seealso \code{\link{create_NICHESObject}},
#'   \code{\link{create_Edgelist_Spatial}}
#'
#' @keywords internal
#' @export
create_Edgelist_SamplingByType <- function(meta.data,
                                            cell.type.col,
                                            seed      = 42L,
                                            max.cells = NULL) {

  ## SECTION 1: Input validation ----

  if (!is.data.frame(meta.data))
    stop("[create_Edgelist_SamplingByType] 'meta.data' must be a data.frame.")
  if (is.null(rownames(meta.data)) || any(rownames(meta.data) == ""))
    stop("[create_Edgelist_SamplingByType] 'meta.data' must have unique, non-empty row names (cell barcodes).")
  if (!is.character(cell.type.col) || length(cell.type.col) != 1L)
    stop("[create_Edgelist_SamplingByType] 'cell.type.col' must be a single character string.")
  if (!cell.type.col %in% colnames(meta.data))
    stop(sprintf(
      "[create_Edgelist_SamplingByType] Column '%s' not found in 'meta.data'.",
      cell.type.col
    ))
  if (!is.null(max.cells)) {
    if (!is.numeric(max.cells) || length(max.cells) != 1L ||
        is.na(max.cells) || max.cells < 1L || max.cells != floor(max.cells))
      stop("[create_Edgelist_SamplingByType] 'max.cells' must be a single positive integer or NULL.")
    max.cells <- as.integer(max.cells)
  }

  ## SECTION 2: Prepare per-type barcode lists ----

  if (!is.null(seed)) {
    rng.state <- .save_rng()
    on.exit(.restore_rng(rng.state), add = TRUE)
    set.seed(seed)
  }

  cell.types  <- sort(unique(as.character(meta.data[[cell.type.col]])))
  n.types     <- length(cell.types)
  all.barcodes <- rownames(meta.data)

  if (n.types < 2L)
    stop("[create_Edgelist_SamplingByType] At least 2 distinct cell types are required.")

  # Build a named list: type label -> character vector of barcodes
  type.cells <- lapply(
    setNames(cell.types, cell.types),
    function(ct) all.barcodes[meta.data[[cell.type.col]] == ct]
  )
  type.counts <- vapply(type.cells, length, integer(1L))

  ## SECTION 3: Sample edges for every directed type-pair crossing ----
  #
  # For each ordered pair (i -> j):
  #   cap = min(n_i, n_j)          -- natural balance cap
  #   cap = min(cap, max.cells)     -- optional user-supplied ceiling
  #
  # Senders and receivers are sampled independently from their respective types
  # and paired positionally. This matches the NICHES v1 strategy.
  # The inner loop is fast (T^2 iterations of sample()), so no parallelism
  # is needed here; the bottleneck is downstream scoring, not edge construction.

  pair.list <- vector("list", n.types * n.types)
  pair.idx  <- 0L

  for (i in seq_along(cell.types)) {
    for (j in seq_along(cell.types)) {
      pair.idx <- pair.idx + 1L

      cap <- min(type.counts[i], type.counts[j])
      if (!is.null(max.cells)) cap <- min(cap, max.cells)
      if (cap < 1L) next

      senders   <- sample(type.cells[[i]], size = cap, replace = FALSE)
      receivers <- sample(type.cells[[j]], size = cap, replace = FALSE)

      pair.list[[pair.idx]] <- data.table(from = senders, to = receivers)
    }
  }

  edge.dt <- rbindlist(Filter(Negate(is.null), pair.list))

  if (nrow(edge.dt) == 0L)
    stop("[create_Edgelist_SamplingByType] No edges were produced. ",
         "Check that cell type labels are non-empty and types have sufficient cells.")

  setkeyv(edge.dt, c("from", "to"))

  ## SECTION 4: Attach metadata attributes and report ----

  attr(edge.dt, "mode")    <- "nonspatial"
  attr(edge.dt, "param")   <- list(cell.type.col = cell.type.col,
                                    seed          = seed,
                                    max.cells     = max.cells)
  attr(edge.dt, "n.cells")    <- nrow(meta.data)
  attr(edge.dt, "cell.names") <- all.barcodes

  message(sprintf(
    paste0("[create_Edgelist_SamplingByType] %d cell types | seed = %s | ",
           "max.cells = %s | %d sampled directed edges."),
    n.types,
    if (is.null(seed)) "not set" else as.character(seed),
    if (is.null(max.cells)) "none" else as.character(max.cells),
    nrow(edge.dt)
  ))

  invisible(edge.dt)
}
