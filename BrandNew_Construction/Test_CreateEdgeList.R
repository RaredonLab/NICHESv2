library(Seurat)
library(igraph)
library(dplyr)
library(RANN)

setwd("~/Documents/GitHub/NICHESv2-NW/BrandNew_Construction")

load("/Volumes/T7/Pneum/Xenium/PPLR.annotated_Type.5_SEE_08-20-2025.Robj")

load("/Volumes/T7/Sophie/data/systems.combined_02-04-2025.Robj")

sub <- subset(PPLR.annotated, subset = TMA == 'PNX3_M1_Sp')

CreateEdgelist <- function(seurat_obj, 
                           radius = NULL, 
                           num_neighbors = NULL, 
                           downsample_col = NULL, 
                           downsample_fraction = NULL, 
                           downsample_num = NULL, 
                           check_coords = c("x", "y")) {
  
  cells <- colnames(seurat_obj)
  metadata <- seurat_obj@meta.data
  
  if (!is.null(downsample_fraction) && !is.null(downsample_num)) {
    stop("Conflict: Both downsample_fraction and downsample_num are provided. Please input only one of them.")
  }
  
  if (is.null(downsample_fraction) && is.null(downsample_num)) {
    downsample_fraction <- 1.0
  }
  
  if (!is.null(downsample_fraction) && (downsample_fraction <= 0 || downsample_fraction > 1)) {
    stop("Invalid downsample_fraction value. It must be >0 and <=1.")
  }
  
  has_spatial_coords <- all(check_coords %in% colnames(metadata))
  
  if (has_spatial_coords) {
    coords <- metadata[, check_coords]
    message("Detected spatial data. Proceeding with distance-based neighborhood thresholding.")
    if (!is.null(radius) && !is.null(num_neighbors)) {
      stop("Conflict: Both radius and num_neighbors are provided. Please use only one of them.")
    }
    if (is.null(radius) && is.null(num_neighbors)) {
      num_neighbors <- 4
      message("No radius or num_neighbors provided. Default value for K-nearest neighbors is set to k = ", num_neighbors)
    }
  } else {
    coords <- NULL
    message("No spatial coordinates '", paste(check_coords, collapse = "', '"), "' found. Please check if coordinates are present in your data under different names and consider renaming them.")
    if (!is.null(colnames(metadata))) {
      message("Available columns in the metadata are: ", paste(colnames(metadata), collapse = ", "))
    }
  }
  
  if (has_spatial_coords) {
    if (!is.null(radius)) {
      message("Thresholding neighborhoods based on radius with radius = ", radius)
      nn_result <- RANN::nn2(coords, coords, k = 100, searchtype = "radius", radius = radius)
      idx.result <- nn_result$nn.idx * (nn_result$nn.dists <= radius)
      
      edges <- do.call(rbind, lapply(seq_len(nrow(idx.result)), function(i) {
        neighbors <- idx.result[i, ][idx.result[i, ] != 0]
        data.frame(from = rep(cells[i], length(neighbors)), to = cells[neighbors])
      }))
      
      edges <- as.matrix(edges)
    } else if (!is.null(num_neighbors)) {
      message("Thresholding neighborhoods based on K-nearest neighbors with k = ", num_neighbors)
      knn <- FNN::get.knn(as.matrix(coords), k = num_neighbors - 1)
      idx.result <- knn$nn.index
      idx.result <- cbind(seq_len(nrow(coords)), idx.result)
      
      edges <- do.call(rbind, lapply(seq_len(nrow(idx.result)), function(i) {
        neighbors <- idx.result[i, ][idx.result[i, ] != 0]
        data.frame(from = rep(cells[i], length(neighbors)), to = cells[neighbors])
      }))
      
      edges <- as.matrix(edges)
    }
  } else {
    if (!is.null(downsample_col) && !downsample_col %in% colnames(metadata)) {
      warning("The provided downsample_col does not exist in the metadata. Treating the entire dataset as one group.")
      downsample_col <- NULL
    }
    
    if (is.null(downsample_col)) {
      if (!is.null(downsample_num)) {
        message("Downsampling the entire dataset by the number of cells.")
        sampled_cells <- metadata %>%
          sample_n(min(downsample_num, n()), replace = FALSE) %>%
          pull()
      } else if (downsample_fraction < 1 && downsample_fraction > 0) {
        message("Downsampling the entire dataset by fraction.")
        sampled_cells <- metadata %>%
          sample_frac(downsample_fraction) %>%
          pull()
      } else {
        message("No downsampling performed -- compute full graph.")
        sampled_cells <- cells
      }
    } else {
      if (!is.null(downsample_num)) {
        message("Downsampling non-spatial data based on provided metadata column: ", downsample_col)
        sampled_cells <- metadata %>%
          group_by(!!sym(downsample_col)) %>%
          sample_n(downsample_num, replace = TRUE) %>%  # ensure sampling required num
          ungroup() %>%
          pull()
        
        # Check and message for groups sampled more than available
        sampled_groups <- metadata %>%
          group_by(!!sym(downsample_col)) %>%
          summarize(group_size = n()) %>%
          filter(group_size < downsample_num)
        
        if (nrow(sampled_groups) > 0) {
          message("Warning: Sample number exceeds actual number in the following groups. Sampled all available cells.")
          for (i in 1:nrow(sampled_groups)) {
            message("Group: ", sampled_groups[[downsample_col]][i], ", Available: ", sampled_groups$group_size[i])
          }
        }
      } else if (downsample_fraction < 1 && downsample_fraction > 0) {
        message("Downsampling non-spatial data based on provided metadata column: ", downsample_col)
        sampled_cells <- metadata %>%
          group_by(!!sym(downsample_col)) %>%
          sample_frac(downsample_fraction) %>%
          ungroup() %>%
          pull()
      } else {
        message("No downsampling performed -- compute full graph.")
        sampled_cells <- cells
      }
    }
    
    if (length(sampled_cells) > 50000) {
      stop("Total number of cells exceeds 50,000 so full graph computation will not be performed. Please use downsampling options.")
    }
    
    edges <- expand.grid(from = sampled_cells, to = sampled_cells)
    edges <- as.matrix(edges)
  }
  
  # Create the graph
  g <- graph_from_edgelist(edges, directed = TRUE)
  g <- igraph::as_data_frame(g)
  
  return(g)
}


test_rad <- CreateEdgelist(sub, radius = 25)
test_kNN <- CreateEdgelist(sub, num_neighbors = 9)
test_spatial_noinput <- CreateEdgelist(sub)
test_spatial_conflict <- CreateEdgelist(sub, radius = 25, num_neighbors = 10)
test_spatial_nocoord <- CreateEdgelist(sub, radius = 25, check_coords = c('x_row','y_col'))

seurat_obj <- systems.combined
downsample_col <- 'CellType_Final'
downsample_fraction <- 0.1

test_downsample_conflict <- CreateEdgelist(systems.combined, downsample_col = "CellType_Final", downsample_fraction = 0.1, downsample_num = 1000)
test_downsample_noinput <- CreateEdgelist(systems.combined)
test_downsample_frac <- CreateEdgelist(systems.combined, radius = 25, downsample_fraction = 0.1)
test_downdownsample_num <- CreateEdgelist(systems.combined, downsample_num = 500)

test_downsample_frac_bycol <- CreateEdgelist(systems.combined, downsample_col = downsample_col, downsample_fraction = 0.1)
test_downdownsample_num_bycol <- CreateEdgelist(systems.combined, downsample_col = downsample_col, downsample_num = 500)


