# create_Edgelist_PseudoBulk.R --------------------------------------------
# Build a type-to-type directed edge list for pseudobulk non-spatial
# cell-to-cell communication analysis. Enumerates all T^2 directed
# ordered pairs of cell types (including autocrine type -> same type).
# No sampling or randomness -- the output is fully deterministic.
#
# FUNCTIONS:
#   create_Edgelist_PseudoBulk()  -- main function
#
# DEPS: data.table (>= 1.14)


# ==============================================================================
#' Build a Pseudobulk Type-to-Type Edge List
#'
#' @description
#' Constructs a directed edge list covering all ordered cell type pairs for
#' pseudobulk non-spatial analysis. Each edge represents a directed
#' type-to-type communication channel (e.g., \code{"ATI"} -> \code{"gCap"}).
#' All T^2 directed pairs are included, encompassing both cross-type and
#' autocrine (same-type) edges.
#'
#' The output edge list is intended to be passed to \code{create_NICHESObject()}
#' with \code{mode = "pseudobulk"}, where the per-type mean expression matrix is
#' computed and paired with this edge list for scoring. Do not pass individual
#' cell barcodes as \code{count.mtx} column names when using this edge list --
#' the column names must be cell type labels.
#'
#' @param meta.data    \code{data.frame} with one row per cell. Row names must
#'   be unique cell barcodes. Must contain the column specified by
#'   \code{cell.type.col}.
#' @param cell.type.col \code{character(1)}. Name of the column in
#'   \code{meta.data} that holds cell type labels.
#'
#' @return A \code{data.table} keyed on \code{(from, to)} with columns:
#' \describe{
#'   \item{\code{from}}{Sending cell type label (character).}
#'   \item{\code{to}}{Receiving cell type label (character).}
#' }
#' Carries metadata attributes: \code{mode} (\code{"pseudobulk"}),
#' \code{param} (the \code{cell.type.col} value), \code{n.cells} (total
#' number of cells in \code{meta.data}), \code{cell.names} (sorted unique
#' cell type labels).
#'
#' @details
#' Cell types are sorted alphabetically before enumeration so the output order
#' is deterministic regardless of the order they appear in \code{meta.data}.
#'
#' @examples
#' meta <- data.frame(
#'   cell.type = c("ATI", "ATI", "gCap", "gCap", "Macro"),
#'   row.names = paste0("cell_", 1:5)
#' )
#' el <- create_Edgelist_PseudoBulk(meta, cell.type.col = "cell.type")
#'
#' @seealso \code{\link{create_NICHESObject}},
#'   \code{\link{create_Edgelist_SamplingByType}},
#'   \code{\link{create_Edgelist_Spatial}}
#'
#' @export
create_Edgelist_PseudoBulk <- function(meta.data, cell.type.col) {

  ## SECTION 1: Input validation ----

  if (!is.data.frame(meta.data))
    stop("[create_Edgelist_PseudoBulk] 'meta.data' must be a data.frame.")
  if (is.null(rownames(meta.data)) || any(rownames(meta.data) == ""))
    stop("[create_Edgelist_PseudoBulk] 'meta.data' must have unique, non-empty row names (cell barcodes).")
  if (!is.character(cell.type.col) || length(cell.type.col) != 1L)
    stop("[create_Edgelist_PseudoBulk] 'cell.type.col' must be a single character string.")
  if (!cell.type.col %in% colnames(meta.data))
    stop(sprintf(
      "[create_Edgelist_PseudoBulk] Column '%s' not found in 'meta.data'.",
      cell.type.col
    ))

  cell.types <- sort(unique(as.character(meta.data[[cell.type.col]])))
  n.types    <- length(cell.types)
  n.cells    <- nrow(meta.data)

  if (n.types < 2L)
    stop("[create_Edgelist_PseudoBulk] At least 2 distinct cell types are required.")

  ## SECTION 2: Build T^2 directed type-pair edge list ----
  #
  # All ordered pairs (i, j) for i in cell.types, j in cell.types -- including
  # i == j (autocrine: type communicates with itself). This covers all possible
  # type-to-type signalling channels without any sampling or randomness.

  edge.dt <- data.table(
    from = rep(cell.types, each  = n.types),
    to   = rep(cell.types, times = n.types)
  )

  setkeyv(edge.dt, c("from", "to"))

  ## SECTION 3: Attach metadata attributes and report ----

  attr(edge.dt, "mode")       <- "pseudobulk"
  attr(edge.dt, "param")      <- cell.type.col
  attr(edge.dt, "n.cells")    <- n.cells
  attr(edge.dt, "cell.names") <- cell.types   # sorted unique type labels

  message(sprintf(
    "[create_Edgelist_PseudoBulk] %d cell types | %d directed type-pair edges (including autocrine).",
    n.types, nrow(edge.dt)
  ))

  invisible(edge.dt)
}
