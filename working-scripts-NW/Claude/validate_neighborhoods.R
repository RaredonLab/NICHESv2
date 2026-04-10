# =============================================================================
# validate_neighborhoods.R
#
# Validation script to verify that extract_neighborhoods() produces results
# identical to the advisor's extract_neighborhoods_optimized().
#
# Prerequisites:
#   1. Source the advisor's function:
#        source("path/to/extract_neighborhoods_optimized.R")
#   2. Source (or load_all) the NICHESv2 package so extract_neighborhoods()
#      is available.
#   3. Have an edge.list data.table (from a NICHESObject or from
#      create_Edgelist_Spatial()) ready.
#
# Usage:
#   result <- validate_neighborhoods(edge.list,
#                                    edge.filter.mode = "both_in_neighborhood",
#                                    n.cores = 4)
# =============================================================================

library(data.table)

# -----------------------------------------------------------------------------
# Converter 1: Advisor's nested list -> canonical form
# -----------------------------------------------------------------------------
# The advisor's output is a named list, one entry per cell. Each entry:
#   list(node_id, neighbors, n_neighbors, local_edges <data.table>,
#        n_local_edges)
#
# Canonical form per cell:
#   list(neighbors = sorted character vector,
#        edges     = sorted character vector of "Sending|Receiving" IDs)
# -----------------------------------------------------------------------------
advisor_to_canonical <- function(advisor.result, verbose = TRUE) {
  if (verbose) message("Converting advisor result to canonical form...")
  start.time <- Sys.time()

  cell.names <- names(advisor.result)
  n.cells    <- length(cell.names)
  canonical  <- vector("list", n.cells)
  names(canonical) <- cell.names

  for (i in seq_len(n.cells)) {
    x <- advisor.result[[i]]

    # Neighbors: sort for deterministic comparison
    nbrs <- sort(as.character(x$neighbors))

    # Edges: reconstruct "Sending|Receiving" IDs from local_edges columns,
    # then sort. The advisor stores original column names sending.barcode
    # and receiving.barcode in the local_edges data.table.
    if (!is.null(x$local_edges) && nrow(x$local_edges) > 0) {
      edges <- sort(paste0(x$local_edges$sending.barcode, "|",
                           x$local_edges$receiving.barcode))
    } else {
      edges <- character(0)
    }

    canonical[[i]] <- list(neighbors = nbrs, edges = edges)

    if (verbose && (i %% 50000 == 0)) {
      message(sprintf("  Converted %s / %s cells (%.0f%%)",
                      format(i, big.mark = ","),
                      format(n.cells, big.mark = ","),
                      100 * i / n.cells))
    }
  }

  elapsed <- as.numeric(difftime(Sys.time(), start.time, units = "secs"))
  if (verbose) message(sprintf("  Done: %.1f seconds", elapsed))

  return(canonical)
}


# -----------------------------------------------------------------------------
# Converter 2: Our extract_neighborhoods() output -> canonical form
# -----------------------------------------------------------------------------
# Our output has:
#   $neighborhood.cell.list  — named list of character vectors (neighbors)
#   $neighborhood.edge.list  — named list of character vectors (edge IDs)
#
# Just sort each vector for deterministic comparison.
# -----------------------------------------------------------------------------
niches_to_canonical <- function(our.result, verbose = TRUE) {
  if (verbose) message("Converting NICHESv2 result to canonical form...")
  start.time <- Sys.time()

  neighborhood.cell.list <- our.result$neighborhood.cell.list
  neighborhood.edge.list <- our.result$neighborhood.edge.list

  cell.names <- names(neighborhood.cell.list)
  n.cells    <- length(cell.names)
  canonical  <- vector("list", n.cells)
  names(canonical) <- cell.names

  for (i in seq_len(n.cells)) {
    cell <- cell.names[i]
    canonical[[i]] <- list(
      neighbors = sort(as.character(neighborhood.cell.list[[cell]])),
      edges     = sort(as.character(neighborhood.edge.list[[cell]]))
    )

    if (verbose && (i %% 50000 == 0)) {
      message(sprintf("  Converted %s / %s cells (%.0f%%)",
                      format(i, big.mark = ","),
                      format(n.cells, big.mark = ","),
                      100 * i / n.cells))
    }
  }

  elapsed <- as.numeric(difftime(Sys.time(), start.time, units = "secs"))
  if (verbose) message(sprintf("  Done: %.1f seconds", elapsed))

  return(canonical)
}


# -----------------------------------------------------------------------------
# Comparison function
# -----------------------------------------------------------------------------
compare_neighborhoods <- function(canonical.a, canonical.b,
                                  label.a = "Advisor",
                                  label.b = "NICHESv2",
                                  verbose = TRUE) {
  if (verbose) message("Comparing canonical neighborhood representations...")

  # Check cell sets
  cells.a <- sort(names(canonical.a))
  cells.b <- sort(names(canonical.b))

  if (!identical(cells.a, cells.b)) {
    only.a <- setdiff(cells.a, cells.b)
    only.b <- setdiff(cells.b, cells.a)
    message(sprintf("FAIL: Cell sets differ!"))
    message(sprintf("  Only in %s: %d cells", label.a, length(only.a)))
    message(sprintf("  Only in %s: %d cells", label.b, length(only.b)))
    if (length(only.a) > 0) {
      message(sprintf("  First 5 in %s only: %s", label.a,
                      paste(head(only.a, 5), collapse = ", ")))
    }
    if (length(only.b) > 0) {
      message(sprintf("  First 5 in %s only: %s", label.b,
                      paste(head(only.b, 5), collapse = ", ")))
    }
    return(list(identical = FALSE, reason = "cell_set_mismatch",
                only.a = only.a, only.b = only.b))
  }

  # Compare per cell
  n.cells <- length(cells.a)
  mismatch.nbrs  <- character(0)
  mismatch.edges <- character(0)

  for (i in seq_len(n.cells)) {
    cell <- cells.a[i]
    if (!identical(canonical.a[[cell]]$neighbors,
                   canonical.b[[cell]]$neighbors)) {
      mismatch.nbrs <- c(mismatch.nbrs, cell)
    }
    if (!identical(canonical.a[[cell]]$edges,
                   canonical.b[[cell]]$edges)) {
      mismatch.edges <- c(mismatch.edges, cell)
    }
  }

  # Report
  if (verbose) {
    message(sprintf("  Cells compared: %s",
                    format(n.cells, big.mark = ",")))
    message(sprintf("  Neighbor mismatches: %d", length(mismatch.nbrs)))
    message(sprintf("  Edge mismatches:     %d", length(mismatch.edges)))

    if (length(mismatch.nbrs) > 0) {
      message(sprintf("\n  First 5 neighbor mismatches (%s vs %s):",
                      label.a, label.b))
      for (cell in head(mismatch.nbrs, 5)) {
        n.a <- length(canonical.a[[cell]]$neighbors)
        n.b <- length(canonical.b[[cell]]$neighbors)
        message(sprintf("    %s: %s=%d, %s=%d",
                        cell, label.a, n.a, label.b, n.b))
        # Show symmetric difference
        only.in.a <- setdiff(canonical.a[[cell]]$neighbors,
                             canonical.b[[cell]]$neighbors)
        only.in.b <- setdiff(canonical.b[[cell]]$neighbors,
                             canonical.a[[cell]]$neighbors)
        if (length(only.in.a) > 0) {
          message(sprintf("      Only in %s: %s", label.a,
                          paste(head(only.in.a, 3), collapse = ", ")))
        }
        if (length(only.in.b) > 0) {
          message(sprintf("      Only in %s: %s", label.b,
                          paste(head(only.in.b, 3), collapse = ", ")))
        }
      }
    }

    if (length(mismatch.edges) > 0) {
      message(sprintf("\n  First 5 edge mismatches (%s vs %s):",
                      label.a, label.b))
      for (cell in head(mismatch.edges, 5)) {
        n.a <- length(canonical.a[[cell]]$edges)
        n.b <- length(canonical.b[[cell]]$edges)
        message(sprintf("    %s: %s=%d, %s=%d",
                        cell, label.a, n.a, label.b, n.b))
        only.in.a <- setdiff(canonical.a[[cell]]$edges,
                             canonical.b[[cell]]$edges)
        only.in.b <- setdiff(canonical.b[[cell]]$edges,
                             canonical.a[[cell]]$edges)
        if (length(only.in.a) > 0) {
          message(sprintf("      Only in %s: %s", label.a,
                          paste(head(only.in.a, 3), collapse = ", ")))
        }
        if (length(only.in.b) > 0) {
          message(sprintf("      Only in %s: %s", label.b,
                          paste(head(only.in.b, 3), collapse = ", ")))
        }
      }
    }

    if (length(mismatch.nbrs) == 0 && length(mismatch.edges) == 0) {
      message("\n  PASS: All neighborhoods are identical!")
    } else {
      message("\n  FAIL: Mismatches found (see details above)")
    }
  }

  return(list(
    identical               = length(mismatch.nbrs) == 0 &&
                              length(mismatch.edges) == 0,
    n.cells                 = n.cells,
    n.neighbor.mismatches   = length(mismatch.nbrs),
    n.edge.mismatches       = length(mismatch.edges),
    neighbor.mismatch.cells = mismatch.nbrs,
    edge.mismatch.cells     = mismatch.edges
  ))
}


# =============================================================================
# Main validation driver
# =============================================================================
#' Run both functions on the same edge.list and compare outputs.
#'
#' @param edge.list data.table with edge|SendingCell|ReceivingCell (NICHESv2
#'   format). The advisor's function will receive a reformatted copy with
#'   sending.barcode and receiving.barcode columns.
#' @param edge.filter.mode Character. Filter mode to test. Both functions
#'   receive the same mode. Default "both_in_neighborhood".
#' @param remove.autocrine Logical. Default FALSE.
#' @param n.cores Integer. Used for BOTH functions. Default 1.
#' @param verbose Logical. Default TRUE.
validate_neighborhoods <- function(edge.list,
                                   edge.filter.mode = "both_in_neighborhood",
                                   remove.autocrine = FALSE,
                                   n.cores = 1,
                                   verbose = TRUE) {

  message("=" , strrep("=", 59))
  message("  NEIGHBORHOOD VALIDATION")
  message("  Mode: ", edge.filter.mode,
          " | Autocrine: ", remove.autocrine,
          " | Cores: ", n.cores)
  message(strrep("=", 60))

  # --- 1. Run advisor's function -------------------------------------------
  message("\n[1/4] Running advisor's extract_neighborhoods_optimized()...")

  # Reformat edge.list to advisor's expected column names
  advisor.input <- data.frame(
    sending.barcode   = edge.list$SendingCell,
    receiving.barcode = edge.list$ReceivingCell,
    stringsAsFactors  = FALSE
  )

  # Map our edge.filter.mode to advisor's edge_filter_mode (same strings)
  t1 <- Sys.time()
  advisor.result <- extract_neighborhoods_optimized(
    edges_df          = advisor.input,
    node_ids          = NULL,
    remove_autocrine  = remove.autocrine,
    edge_filter_mode  = edge.filter.mode,
    n_cores           = n.cores
  )
  t.advisor <- as.numeric(difftime(Sys.time(), t1, units = "secs"))
  message(sprintf("  Advisor function: %.1f seconds", t.advisor))

  # --- 2. Run our function -------------------------------------------------
  message("\n[2/4] Running NICHESv2 extract_neighborhoods()...")

  t2 <- Sys.time()
  our.result <- extract_neighborhoods(
    edge.list        = edge.list,
    edge.filter.mode = edge.filter.mode,
    remove.autocrine = remove.autocrine,
    n.cores          = n.cores,
    verbose          = verbose
  )
  t.ours <- as.numeric(difftime(Sys.time(), t2, units = "secs"))
  message(sprintf("  NICHESv2 function: %.1f seconds", t.ours))

  # --- 3. Convert to canonical form ----------------------------------------
  message("\n[3/4] Converting to canonical form for comparison...")

  canonical.advisor <- advisor_to_canonical(advisor.result, verbose = verbose)
  canonical.ours    <- niches_to_canonical(our.result, verbose = verbose)

  # --- 4. Compare -----------------------------------------------------------
  message("\n[4/4] Comparing results...")

  comparison <- compare_neighborhoods(canonical.advisor, canonical.ours,
                                      label.a = "Advisor",
                                      label.b = "NICHESv2",
                                      verbose = verbose)

  # Summary
  message("\n", strrep("=", 60))
  message("  VALIDATION SUMMARY")
  message(strrep("=", 60))
  message(sprintf("  Advisor time:  %.1f sec", t.advisor))
  message(sprintf("  NICHESv2 time: %.1f sec", t.ours))
  message(sprintf("  Speedup:       %.2fx", t.advisor / t.ours))
  message(sprintf("  Result:        %s",
                  if (comparison$identical) "PASS" else "FAIL"))
  message(strrep("=", 60))

  return(invisible(list(
    comparison    = comparison,
    time.advisor  = t.advisor,
    time.ours     = t.ours,
    advisor.result = advisor.result,
    our.result     = our.result
  )))
}
