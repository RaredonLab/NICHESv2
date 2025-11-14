library(Seurat)
library(igraph)
library(dplyr)
library(RANN)

CreateEdgelist <- function(seurat_obj,
                           weighted_by_distance = TRUE, 
                           distance_threshold = NULL, 
                           num_neighbors = NULL, 
                           downsample_col = NULL, 
                           downsample_fraction = 0.1) {
  
  cells <- colnames(seurat_obj)
  coords <- seurat_obj@meta.data[, c("x", "y")]
  
  if (!is.null(coords) && ncol(coords) == 2) {
    has_spatial_coords <- TRUE
  } else {
    has_spatial_coords <- FALSE
    weighted_by_distance <- FALSE
  }
  
  if (weighted_by_distance && !is.null(distance_threshold)) {
    nn_result <- RANN::nn2(coords, coords, k = 100, searchtype = "radius", radius = distance_threshold)
    idx.result <- nn_result$nn.idx*(nn_result$nn.dists<=25)
    
    edges <- do.call(rbind, lapply(seq_len(nrow(idx.result)), function(i) {
      data.frame(from = rep(i, sum(idx.result[i, ] != 0)),
                 to = idx.result[i, ][idx.result[i, ] != 0])
    }))
    
    edges <- edges %>%
      mutate(from = cells[from],
             to = cells[to])
    
    edges <- as.matrix(edges)
  } else if (weighted_by_distance && !is.null(num_neighbors)) {
    # Calculate distances and filter by k nearest neighbors using KD-tree
    knn <- FNN::get.knn(as.matrix(coords), k = num_neighbors)
    idx.result <- knn$nn.index
    idx.result <- cbind(seq_len(nrow(coords)), idx.result)
    
    edges <- do.call(rbind, lapply(seq_len(nrow(idx.result)), function(i) {
      data.frame(from = rep(i, sum(idx.result[i, ] != 0)),
                 to = idx.result[i, ][idx.result[i, ] != 0])
    }))
    
    # Map numeric indices to cell barcodes
    edges <- edges %>%
      mutate(from = cells[from],
             to = cells[to])
    
    edges <- as.matrix(edges)
  } else if (!weighted_by_distance && !is.null(downsample_col)) {
    # Downsample cells based on the specified metadata column
    metadata <- seurat_obj@meta.data
    sampled_cells <- metadata %>%
      group_by(!!as.symbol(downsample_col)) %>%
      sample_frac(downsample_fraction) %>%
      pull(names)
    
    sampled_cells_idx <- match(sampled_cells, cells)
    edges <- expand.grid(sampled_cells_idx, sampled_cells_idx)
    edges <- edges[edges[, 1] != edges[, 2], ]
  } else {
    warning("Do you really want to do this? The graph size will be extremely and unnecessarily large. Recommend user filter it.")
    if (has_spatial_coords) {
      edges <- expand.grid(1:length(cells), 1:length(cells))
      edges <- edges[edges[, 1] != edges[, 2], ]
    } else {
      edges <- expand.grid(cells, cells)
      edges <- edges[edges[, 1] != edges[, 2], ]
    }
  }
  
  # Create the graph
  g <- graph_from_edgelist(as.matrix(edges),directed = TRUE)
  g <- igraph::as_data_frame(g)
  return(g)
}
  
  # if (weighted_by_distance && has_spatial_coords) {
  #   E(g)$weight <- sapply(1:nrow(edges), function(i) dist(rbind(coords[edges[i, 1], ], coords[edges[i, 2], ])))
  # }
  
test <- CreateEdgelist(sub,
                       weighted_by_distance = TRUE,
                       distance_threshold = 25)
