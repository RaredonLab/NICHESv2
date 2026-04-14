#' Extract Neighborhoods from an Edge List
#'
#' For each cell in the edge list, identifies its spatial neighbors and the
#' edges that fall within its neighborhood. Returns three structures suitable
#' for populating the neighborhood slots of a \code{NICHESObject}.
#'
#' Neighborhoods are defined topologically from the edge list: cell B is a
#' neighbor of cell A if any edge connects them (in either direction). The
#' \code{edge.filter.mode} parameter controls which edges are assigned to
#' each cell's neighborhood.
#'
#' @param edge.list A \code{data.table} with columns \code{edge},
#'   \code{SendingCell}, and \code{ReceivingCell}. Typically
#'   \code{niches.obj$edge.list} or output of \code{create_Edgelist_Spatial()}.
#' @param edge.filter.mode Character. One of:
#'   \describe{
#'     \item{\code{"both_in_neighborhood"}}{(default) Both endpoints of an edge
#'       must be in the cell's neighborhood (cell itself + its neighbors).}
#'     \item{\code{"any_in_neighborhood"}}{At least one endpoint must be in the
#'       neighborhood.}
#'     \item{\code{"central_cell_only"}}{The central cell must be an endpoint.}
#'   }
#' @param remove.autocrine Logical. If \code{TRUE}, removes self-loops from the
#'   adjacency (neighbor) list. Autocrine edges may still appear in
#'   \code{neighborhood.edge.list} if they pass the edge filter. Default
#'   \code{FALSE}.
#' @param n.cores Integer. Number of cores for parallel processing. Default 1
#'   (serial).
#' @param verbose Logical. Print progress messages. Default \code{TRUE}.
#'
#' @return A list with three elements:
#'   \describe{
#'     \item{\code{neighborhood.cell.list}}{Named list of character vectors.
#'       Each element is keyed by a cell barcode and contains the barcodes of
#'       its neighbors.}
#'     \item{\code{neighborhood.edge.list}}{Named list of character vectors.
#'       Each element is keyed by a cell barcode and contains the edge IDs
#'       (format \code{"SendingCell|ReceivingCell"}) belonging to that cell's
#'       neighborhood.}
#'     \item{\code{neighborhood.meta}}{A \code{data.table} with columns
#'       \code{cell}, \code{n.neighbors}, \code{n.edges}, keyed on
#'       \code{cell}.}
#'   }
#'   Attributes: \code{"edge.filter.mode"}, \code{"remove.autocrine"}.
#'
#' @export
extract_neighborhoods <- function(edge.list,
                                  edge.filter.mode = "both_in_neighborhood",
                                  remove.autocrine = FALSE,
                                  n.cores = 1,
                                  verbose = TRUE) {

  ## -----------------------------------------------------------------------
  ## Validate inputs
  ## -----------------------------------------------------------------------
  if (!data.table::is.data.table(edge.list)) {
    stop("edge.list must be a data.table")
  }
  required.cols <- c("edge", "SendingCell", "ReceivingCell")
  missing.cols <- setdiff(required.cols, names(edge.list))
  if (length(missing.cols) > 0) {
    stop("edge.list is missing required columns: ",
         paste(missing.cols, collapse = ", "))
  }
  if (nrow(edge.list) == 0) {
    stop("edge.list has zero rows")
  }

  valid.modes <- c("both_in_neighborhood", "any_in_neighborhood",
                    "central_cell_only")
  if (!edge.filter.mode %in% valid.modes) {
    stop("edge.filter.mode must be one of: ",
         paste(valid.modes, collapse = ", "))
  }

  if (n.cores < 1) stop("n.cores must be at least 1")
  available.cores <- parallel::detectCores()
  if (is.null(available.cores)) available.cores <- 1L
  if (n.cores > available.cores) {
    warning(sprintf("Requested %d cores but only %d available. Using %d.",
                    n.cores, available.cores, available.cores))
    n.cores <- available.cores
  }

  ## -----------------------------------------------------------------------
  ## Step 1: Integer index mapping
  ## -----------------------------------------------------------------------
  if (verbose) message("Step 1/4: Building integer index mapping...")

  all.cells <- unique(c(edge.list$SendingCell, edge.list$ReceivingCell))
  n.cells   <- length(all.cells)
  cell.to.idx <- setNames(seq_len(n.cells), all.cells)
  idx.to.cell <- all.cells   # positional: idx.to.cell[i] = barcode for index i

  send.idx <- unname(cell.to.idx[edge.list$SendingCell])
  recv.idx <- unname(cell.to.idx[edge.list$ReceivingCell])
  edge.ids <- edge.list$edge  # "SendingCell|ReceivingCell" strings
  n.edges  <- length(edge.ids)

  if (verbose) message(sprintf("  %s cells, %s edges",
                               format(n.cells, big.mark = ","),
                               format(n.edges, big.mark = ",")))

  ## -----------------------------------------------------------------------
  ## Step 2: Build bidirectional adjacency
  ## -----------------------------------------------------------------------
  if (verbose) message("Step 2/4: Building adjacency lists...")

  adj.dt <- unique(data.table::rbindlist(list(
    data.table::data.table(cell = send.idx, neighbor = recv.idx),
    data.table::data.table(cell = recv.idx, neighbor = send.idx)
  )))

  if (remove.autocrine) {
    n.autocrine <- sum(adj.dt$cell == adj.dt$neighbor)
    adj.dt <- adj.dt[cell != neighbor]
    if (verbose) message(sprintf("  Removed %d autocrine adjacency entries",
                                 n.autocrine))
  }

  # Plain R list for O(1) positional lookup: adj.lookup[[i]] = integer vector
  adj.lookup <- vector("list", n.cells)
  adj.grouped <- adj.dt[, .(nbrs = list(neighbor)), by = cell]
  for (r in seq_len(nrow(adj.grouped))) {
    adj.lookup[[ adj.grouped$cell[r] ]] <- adj.grouped$nbrs[[r]]
  }
  rm(adj.dt, adj.grouped)

  if (verbose) message("  Adjacency lists built")

  ## -----------------------------------------------------------------------
  ## Step 3: Build edge lookup by sender / receiver
  ## -----------------------------------------------------------------------
  if (verbose) message("Step 3/4: Building edge lookup indices...")

  edge.by.sender   <- vector("list", n.cells)
  edge.by.receiver <- vector("list", n.cells)

  sender.grouped <- data.table::data.table(
    cell = send.idx, row = seq_len(n.edges)
  )[, .(rows = list(row)), by = cell]

  for (r in seq_len(nrow(sender.grouped))) {
    edge.by.sender[[ sender.grouped$cell[r] ]] <- sender.grouped$rows[[r]]
  }

  receiver.grouped <- data.table::data.table(
    cell = recv.idx, row = seq_len(n.edges)
  )[, .(rows = list(row)), by = cell]

  for (r in seq_len(nrow(receiver.grouped))) {
    edge.by.receiver[[ receiver.grouped$cell[r] ]] <- receiver.grouped$rows[[r]]
  }
  rm(sender.grouped, receiver.grouped)

  if (verbose) message("  Edge lookup indices built")

  ## -----------------------------------------------------------------------
  ## Step 4: Per-cell neighborhood extraction
  ## -----------------------------------------------------------------------
  if (verbose) message(sprintf(
    "Step 4/4: Extracting neighborhoods [mode: %s]...", edge.filter.mode))

  start.time <- Sys.time()

  # Worker function: processes a chunk of integer cell indices.
  # All lookup objects are accessed from the enclosing environment (serial)
  # or from the exported cluster environment (parallel).
  process.chunk <- function(cell.indices) {
    n.chunk <- length(cell.indices)
    chunk.cell.list <- vector("list", n.chunk)
    chunk.edge.list <- vector("list", n.chunk)
    chunk.names     <- character(n.chunk)
    chunk.n.nbrs    <- integer(n.chunk)
    chunk.n.edges   <- integer(n.chunk)

    for (i in seq_len(n.chunk)) {
      ci <- cell.indices[i]
      nbr.idx <- adj.lookup[[ci]]
      if (is.null(nbr.idx)) nbr.idx <- integer(0)
      hood.idx <- c(ci, nbr.idx)

      # Collect candidate edge row indices from all hood members
      candidate.rows <- unique(c(
        unlist(edge.by.sender[hood.idx],   use.names = FALSE),
        unlist(edge.by.receiver[hood.idx], use.names = FALSE)
      ))

      if (length(candidate.rows) == 0L) {
        chunk.names[i]       <- idx.to.cell[ci]
        chunk.cell.list[[i]] <- idx.to.cell[nbr.idx]
        chunk.edge.list[[i]] <- character(0)
        chunk.n.nbrs[i]      <- length(nbr.idx)
        chunk.n.edges[i]     <- 0L
        next
      }

      # Filter by mode
      s <- send.idx[candidate.rows]
      r <- recv.idx[candidate.rows]

      if (edge.filter.mode == "both_in_neighborhood") {
        keep <- (s %in% hood.idx) & (r %in% hood.idx)
      } else if (edge.filter.mode == "any_in_neighborhood") {
        keep <- (s %in% hood.idx) | (r %in% hood.idx)
      } else {
        # central_cell_only
        keep <- (s == ci) | (r == ci)
      }

      local.rows <- candidate.rows[keep]

      chunk.names[i]       <- idx.to.cell[ci]
      chunk.cell.list[[i]] <- idx.to.cell[nbr.idx]
      chunk.edge.list[[i]] <- edge.ids[local.rows]
      chunk.n.nbrs[i]      <- length(nbr.idx)
      chunk.n.edges[i]     <- length(local.rows)
    }

    names(chunk.cell.list) <- chunk.names
    names(chunk.edge.list) <- chunk.names

    list(
      cell.list = chunk.cell.list,
      edge.list = chunk.edge.list,
      meta = data.table::data.table(
        cell        = chunk.names,
        n.neighbors = chunk.n.nbrs,
        n.edges     = chunk.n.edges
      )
    )
  }

  all.indices <- seq_len(n.cells)

  if (n.cores == 1) {
    ## --- Serial path ---
    if (verbose) message("  Running in serial mode...")
    result <- process.chunk(all.indices)

  } else {
    ## --- Parallel path (socket cluster) ---
    if (verbose) message(sprintf("  Running in parallel with %d cores...",
                                 n.cores))

    chunk.size <- ceiling(n.cells / n.cores)
    chunks <- split(all.indices,
                    ceiling(seq_along(all.indices) / chunk.size))

    cl <- parallel::makeCluster(n.cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)

    parallel::clusterEvalQ(cl, {
      library(data.table)
      data.table::setDTthreads(1)
    })

    parallel::clusterExport(cl,
      varlist = c("adj.lookup", "edge.by.sender", "edge.by.receiver",
                  "send.idx", "recv.idx", "edge.ids", "idx.to.cell",
                  "edge.filter.mode"),
      envir = environment()
    )

    chunk.results <- parallel::parLapply(cl, chunks, process.chunk)

    # Merge chunk results (unname to prevent c() from prefixing chunk names)
    result <- list(
      cell.list = do.call(c, unname(lapply(chunk.results, `[[`, "cell.list"))),
      edge.list = do.call(c, unname(lapply(chunk.results, `[[`, "edge.list"))),
      meta = data.table::rbindlist(lapply(chunk.results, `[[`, "meta"))
    )
  }

  elapsed <- as.numeric(difftime(Sys.time(), start.time, units = "secs"))
  if (verbose) {
    message(sprintf("  Complete: %.1f seconds (%.0f cells/sec)",
                    elapsed, n.cells / elapsed))
    message(sprintf("  Avg neighbors/cell: %.1f",
                    mean(result$meta$n.neighbors)))
    message(sprintf("  Avg edges/cell: %.1f",
                    mean(result$meta$n.edges)))
  }

  ## -----------------------------------------------------------------------
  ## Finalize output
  ## -----------------------------------------------------------------------
  neighborhood.cell.list <- result$cell.list
  neighborhood.edge.list <- result$edge.list
  neighborhood.meta      <- result$meta
  data.table::setkey(neighborhood.meta, cell)

  output <- list(
    neighborhood.cell.list = neighborhood.cell.list,
    neighborhood.edge.list = neighborhood.edge.list,
    neighborhood.meta      = neighborhood.meta
  )

  attr(output, "edge.filter.mode") <- edge.filter.mode
  attr(output, "remove.autocrine") <- remove.autocrine

  return(output)
}


#' Add Neighborhood Slots to a NICHESObject
#'
#' Extracts neighborhoods from the edge list of a \code{NICHESObject} and
#' populates the three neighborhood slots: \code{$neighborhood.cell.list},
#' \code{$neighborhood.edge.list}, and \code{$neighborhood.meta}.
#'
#' This function calls \code{\link{extract_neighborhoods}} internally.
#' If the neighborhood slots are already populated, they are overwritten
#' with a message.
#'
#' @param niches.obj A \code{NICHESObject} with a non-empty \code{$edge.list}.
#' @param edge.filter.mode Character. Passed to
#'   \code{\link{extract_neighborhoods}}. Default
#'   \code{"both_in_neighborhood"}.
#' @param remove.autocrine Logical. Passed to
#'   \code{\link{extract_neighborhoods}}. Default \code{FALSE}.
#' @param n.cores Integer. Number of cores for parallel processing. Default 1.
#' @param verbose Logical. Print progress messages. Default \code{TRUE}.
#'
#' @return The input \code{NICHESObject} with the three neighborhood slots
#'   populated. Attribute \code{"edge.filter.mode"} is set on the object.
#'
#' @export
add_neighborhoods <- function(niches.obj,
                              edge.filter.mode = "both_in_neighborhood",
                              remove.autocrine = FALSE,
                              n.cores = 1,
                              verbose = TRUE) {

  if (!inherits(niches.obj, "NICHESObject")) {
    stop("niches.obj must be a NICHESObject")
  }
  if (is.null(niches.obj$edge.list) || nrow(niches.obj$edge.list) == 0) {
    stop("niches.obj$edge.list is empty; cannot extract neighborhoods")
  }

  # Warn if overwriting
  if (!is.null(niches.obj$neighborhood.cell.list)) {
    if (verbose) message("Overwriting existing neighborhood slots")
  }

  # Extract neighborhoods
  result <- extract_neighborhoods(
    edge.list        = niches.obj$edge.list,
    edge.filter.mode = edge.filter.mode,
    remove.autocrine = remove.autocrine,
    n.cores          = n.cores,
    verbose          = verbose
  )

  # Write to slots
  niches.obj$neighborhood.cell.list <- result$neighborhood.cell.list
  niches.obj$neighborhood.edge.list <- result$neighborhood.edge.list
  niches.obj$neighborhood.meta      <- result$neighborhood.meta

  # Store parameter on the object
  attr(niches.obj, "edge.filter.mode") <- attr(result, "edge.filter.mode")

  return(niches.obj)
}
