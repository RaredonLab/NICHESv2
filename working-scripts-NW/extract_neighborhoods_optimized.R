extract_neighborhoods_optimized <- function(edges_df, 
                                            node_ids = NULL,
                                            remove_autocrine = FALSE,
                                            edge_filter_mode = "both_in_neighborhood",
                                            n_cores = 1) {
  
  cat("=== OPTIMIZED NEIGHBORHOOD EXTRACTION (Integer Indexing) ===\n\n")
  
  # [Steps 0-3 remain the same until Step 4]
  
  available_cores <- parallel::detectCores()
  if (is.null(available_cores)) available_cores <- 1
  
  if (n_cores < 1) {
    stop("n_cores must be at least 1")
  }
  
  if (n_cores > available_cores) {
    warning(sprintf("Requested %d cores but only %d available. Using %d cores.",
                    n_cores, available_cores, available_cores))
    n_cores <- available_cores
  }
  
  has_pbapply <- requireNamespace("pbapply", quietly = TRUE)
  
  if (n_cores > 1) {
    cat(sprintf("Parallel processing enabled: using %d cores\n", n_cores))
    cat(sprintf("  (System has %d cores available)\n", available_cores))
    if (has_pbapply) {
      cat("  Progress bar: ENABLED (via pbapply package)\n")
    } else {
      cat("  Progress bar: DISABLED (install 'pbapply' package for progress bars)\n")
      cat("    Install with: install.packages('pbapply')\n")
    }
  } else {
    cat("Running in serial mode (single core)\n")
    cat("  Progress bar: ENABLED\n")
  }
  
  cat("\nStep 1: Validating inputs and converting to data.table...\n")
  
  edges_dt <- as.data.table(edges_df)
  
  if (!all(c("sending.barcode", "receiving.barcode") %in% names(edges_dt))) {
    stop("edges_df must contain 'sending.barcode' and 'receiving.barcode' columns")
  }
  
  valid_modes <- c("both_in_neighborhood", "any_in_neighborhood", "central_cell_only")
  if (!edge_filter_mode %in% valid_modes) {
    stop(sprintf("edge_filter_mode must be one of: %s", 
                 paste(valid_modes, collapse = ", ")))
  }
  
  cat(sprintf("  Edge filter mode: %s\n", edge_filter_mode))
  cat(sprintf("  Remove autocrine edges: %s\n", remove_autocrine))
  
  if (is.null(node_ids)) {
    node_ids <- unique(c(edges_dt$sending.barcode, edges_dt$receiving.barcode))
  }
  
  cat(sprintf("  Total nodes: %s\n", format(length(node_ids), big.mark = ",")))
  cat(sprintf("  Total edges: %s\n", format(nrow(edges_dt), big.mark = ",")))
  
  cat("\nStep 1.5: Creating integer index mapping...\n")
  cat("  (Integer operations are faster than string operations)\n")
  
  node_to_idx <- setNames(seq_along(node_ids), as.character(node_ids))
  idx_to_node <- setNames(as.character(node_ids), seq_along(node_ids))
  
  edges_dt[, send_idx := node_to_idx[as.character(sending.barcode)]]
  edges_dt[, recv_idx := node_to_idx[as.character(receiving.barcode)]]
  edges_dt[, edge_idx := .I]
  
  cat("  Integer indexing complete\n")
  
  cat("\nStep 2: Building adjacency lists with integer indices...\n")
  
  edges_fwd <- edges_dt[, .(node_idx = send_idx, neighbor_idx = recv_idx)]
  edges_bwd <- edges_dt[, .(node_idx = recv_idx, neighbor_idx = send_idx)]
  all_conn <- unique(rbind(edges_fwd, edges_bwd))
  
  if (remove_autocrine) {
    n_autocrine <- sum(all_conn$node_idx == all_conn$neighbor_idx)
    all_conn <- all_conn[node_idx != neighbor_idx]
    cat(sprintf("  Removed %s autocrine edges\n", format(n_autocrine, big.mark = ",")))
  } else {
    n_autocrine <- sum(all_conn$node_idx == all_conn$neighbor_idx)
    cat(sprintf("  Keeping autocrine edges (found %s)\n", format(n_autocrine, big.mark = ",")))
  }
  
  adj_list <- all_conn[, .(neighbors = list(unique(neighbor_idx))), by = node_idx]
  
  cat(sprintf("  Built adjacency lists for %s nodes\n", 
              format(nrow(adj_list), big.mark = ",")))
  
  cat("\nStep 3: Building edge lookup indices with integer indices...\n")
  
  edge_by_sender <- edges_dt[, .(edge_indices = list(edge_idx)), by = send_idx]
  setnames(edge_by_sender, "send_idx", "node_idx")
  
  edge_by_receiver <- edges_dt[, .(edge_indices = list(edge_idx)), by = recv_idx]
  setnames(edge_by_receiver, "recv_idx", "node_idx")
  
  cat("  Edge indices built\n")
  
  # --------------------------------------------------------------------------
  # Step 4: Define the worker function - FIXED VERSION
  # --------------------------------------------------------------------------
  
  process_node_chunk <- function(node_idx_chunk, 
                                 adj_list, 
                                 edge_by_sender, 
                                 edge_by_receiver,
                                 edges_dt,
                                 edge_filter_mode,
                                 idx_to_node) {
    
    # FIXED: Use different variable name to avoid shadowing
    process_single_node <- function(current_node_idx) {
      
      # FIXED: Use explicit filtering with the correct variable
      neighbor_result <- adj_list[node_idx == current_node_idx]$neighbors
      neighbor_indices <- if (length(neighbor_result) > 0) neighbor_result[[1]] else integer(0)
      
      # Define the neighborhood
      hood_indices <- c(current_node_idx, neighbor_indices)
      
      # Collect edges - FIXED: Use current_node_idx consistently
      sender_edges <- edge_by_sender[node_idx %in% hood_indices]$edge_indices
      receiver_edges <- edge_by_receiver[node_idx %in% hood_indices]$edge_indices
      
      relevant_edge_indices <- unique(c(
        unlist(sender_edges),
        unlist(receiver_edges)
      ))
      
      # Apply edge filtering - FIXED: Use current_node_idx
      if (edge_filter_mode == "both_in_neighborhood") {
        local_edges <- edges_dt[relevant_edge_indices][
          (send_idx %in% hood_indices) & (recv_idx %in% hood_indices)
        ]
        
      } else if (edge_filter_mode == "any_in_neighborhood") {
        local_edges <- edges_dt[relevant_edge_indices][
          (send_idx %in% hood_indices) | (recv_idx %in% hood_indices)
        ]
        
      } else if (edge_filter_mode == "central_cell_only") {
        local_edges <- edges_dt[relevant_edge_indices][
          (send_idx == current_node_idx) | (recv_idx == current_node_idx)
        ]
      }
      
      # Convert back to original barcodes
      node_barcode <- idx_to_node[as.character(current_node_idx)]
      neighbor_barcodes <- idx_to_node[as.character(neighbor_indices)]
      
      return(list(
        node_id = node_barcode,
        neighbors = neighbor_barcodes,
        n_neighbors = length(neighbor_barcodes),
        local_edges = local_edges,
        n_local_edges = nrow(local_edges)
      ))
    }
    
    chunk_results <- lapply(node_idx_chunk, process_single_node)
    result_names <- sapply(chunk_results, function(x) as.character(x$node_id))
    names(chunk_results) <- result_names
    
    return(chunk_results)
  }
  
  # --------------------------------------------------------------------------
  # Step 5: Extract neighborhoods - FIXED VERSION
  # --------------------------------------------------------------------------
  cat("\nStep 4: Extracting neighborhoods for all nodes...\n")
  cat(sprintf("  Processing %s nodes\n", format(length(node_ids), big.mark = ",")))
  
  start_time <- Sys.time()
  
  node_indices <- node_to_idx[as.character(node_ids)]
  
  if (n_cores == 1) {
    cat("  Running in serial mode with progress updates...\n\n")
    
    progress_every <- max(1, floor(length(node_indices) / 5000))
    neighborhoods <- vector("list", length(node_indices))
    
    for (i in seq_along(node_indices)) {
      if (i %% progress_every == 0) {
        cat(sprintf("    Progress: %s/%s (%.1f%%)\n", 
                    format(i, big.mark = ","),
                    format(length(node_indices), big.mark = ","),
                    100 * i / length(node_indices)))
      }
      
      # FIXED: Use explicit variable name
      current_node_idx <- node_indices[i]
      
      # FIXED: Use explicit filtering
      neighbor_result <- adj_list[node_idx == current_node_idx]$neighbors
      neighbor_indices <- if (length(neighbor_result) > 0) neighbor_result[[1]] else integer(0)
      
      hood_indices <- c(current_node_idx, neighbor_indices)
      
      # FIXED: Use current_node_idx consistently
      sender_edges <- edge_by_sender[node_idx %in% hood_indices]$edge_indices
      receiver_edges <- edge_by_receiver[node_idx %in% hood_indices]$edge_indices
      relevant_edge_indices <- unique(c(unlist(sender_edges), unlist(receiver_edges)))
      
      # FIXED: Use current_node_idx in filtering
      if (edge_filter_mode == "both_in_neighborhood") {
        local_edges <- edges_dt[relevant_edge_indices][
          (send_idx %in% hood_indices) & (recv_idx %in% hood_indices)
        ]
      } else if (edge_filter_mode == "any_in_neighborhood") {
        local_edges <- edges_dt[relevant_edge_indices][
          (send_idx %in% hood_indices) | (recv_idx %in% hood_indices)
        ]
      } else if (edge_filter_mode == "central_cell_only") {
        local_edges <- edges_dt[relevant_edge_indices][
          (send_idx == current_node_idx) | (recv_idx == current_node_idx)
        ]
      }
      
      node_barcode <- idx_to_node[as.character(current_node_idx)]
      neighbor_barcodes <- idx_to_node[as.character(neighbor_indices)]
      
      neighborhoods[[i]] <- list(
        node_id = node_barcode,
        neighbors = neighbor_barcodes,
        n_neighbors = length(neighbor_barcodes),
        local_edges = local_edges,
        n_local_edges = nrow(local_edges)
      )
    }
    
    names(neighborhoods) <- sapply(neighborhoods, function(x) as.character(x$node_id))
    
  } else {
    # Parallel execution
    cat(sprintf("  Running in parallel mode with %d cores...\n", n_cores))
    
    chunk_size <- ceiling(length(node_indices) / n_cores)
    node_idx_chunks <- split(node_indices, ceiling(seq_along(node_indices) / chunk_size))
    
    cat(sprintf("  Split into %d chunks of ~%s nodes each\n\n", 
                length(node_idx_chunks), 
                format(chunk_size, big.mark = ",")))
    
    cl <- parallel::makeCluster(n_cores)
    
    parallel::clusterExport(cl, 
                            varlist = c("adj_list", "edge_by_sender", "edge_by_receiver",
                                        "edges_dt", "edge_filter_mode", "idx_to_node"),
                            envir = environment())
    
    parallel::clusterEvalQ(cl, library(data.table))
    
    if (has_pbapply) {
      cat("  Processing with progress bar...\n")
      chunk_results <- pbapply::pblapply(
        X = node_idx_chunks,
        FUN = process_node_chunk,
        adj_list = adj_list,
        edge_by_sender = edge_by_sender,
        edge_by_receiver = edge_by_receiver,
        edges_dt = edges_dt,
        edge_filter_mode = edge_filter_mode,
        idx_to_node = idx_to_node,
        cl = cl
      )
    } else {
      cat("  Processing chunks (no progress bar available)...\n")
      chunk_results <- parallel::parLapply(
        cl = cl,
        X = node_idx_chunks,
        fun = process_node_chunk,
        adj_list = adj_list,
        edge_by_sender = edge_by_sender,
        edge_by_receiver = edge_by_receiver,
        edges_dt = edges_dt,
        edge_filter_mode = edge_filter_mode,
        idx_to_node = idx_to_node
      )
    }
    
    parallel::stopCluster(cl)
    
    cat("\n  Combining results from all chunks (using stack-safe method)...\n")
    
    neighborhoods <- list()
    for (i in seq_along(chunk_results)) {
      if (i %% 5 == 0) {
        cat(sprintf("    Merging chunk %d/%d\n", i, length(chunk_results)))
      }
      neighborhoods <- c(neighborhoods, chunk_results[[i]])
    }
    
    cat("  Merge complete\n")
  }
  
  end_time <- Sys.time()
  elapsed <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  cat(sprintf("\nStep 4: Complete! (%.2f seconds)\n", elapsed))
  cat(sprintf("  Processing rate: %.0f nodes/second\n", 
              length(node_ids) / elapsed))
  
  cat("\n=== EXTRACTION SUMMARY ===\n")
  
  total_neighbors <- sum(sapply(neighborhoods, function(x) x$n_neighbors))
  total_local_edges <- sum(sapply(neighborhoods, function(x) x$n_local_edges))
  
  cat(sprintf("Total nodes processed: %s\n", 
              format(length(neighborhoods), big.mark = ",")))
  cat(sprintf("Total neighbor relationships: %s\n", 
              format(total_neighbors, big.mark = ",")))
  cat(sprintf("Average neighbors per node: %.2f\n", 
              total_neighbors / length(neighborhoods)))
  cat(sprintf("Total local edges (summed across all neighborhoods): %s\n", 
              format(total_local_edges, big.mark = ",")))
  cat(sprintf("Average local edges per node: %.2f\n", 
              total_local_edges / length(neighborhoods)))
  cat(sprintf("\nTotal processing time: %.2f seconds (%.2f minutes)\n", 
              elapsed, elapsed / 60))
  
  if (n_cores > 1) {
    speedup_estimate <- min(n_cores * 0.7, n_cores - 0.5)
    cat(sprintf("Estimated speedup from parallelization: ~%.1fx\n", speedup_estimate))
  }
  
  cat("\n=== EXTRACTION COMPLETE ===\n\n")
  
  return(neighborhoods)
}