library(Seurat)
library(ggplot2)
library(dplyr)
library(spNNGP)

# setwd('/Volumes/T7/Pneum/Xenium/Test_Composition')
# 
# load("/Volumes/T7/Pneum/Xenium/PPLR.annotated_Type.5_SEE_08-20-2025.Robj")

sub <- subset(PPLR.annotated, subset = TMA == 'PNX3_F1_Sp')

compute_edgelist <- function(filtered.obj,
                             position.x,
                             position.y,
                             k=4,
                             rad.set=NULL,
                             weighted = FALSE
){
  
  ### CREATE MAPPING ###
  
  # Create adjacency matrix
  # Adapted from :: https://stackoverflow.com/questions/16075232/how-to-create-adjacency-matrix-from-grid-coordinates-in-r
  # Setup numbering and labeling
  # jc: possible bug, change object to filtered.obj
  df <- data.frame(x = filtered.obj[[position.x]], y = filtered.obj[[position.y]])
  df$barcode <- rownames(df)
  df$x <- as.character(df$x)
  df$y <- as.character(df$y)
  df$x <- as.numeric(df$x)
  df$y <- as.numeric(df$y)
  df <- df[,c('x','y')]
  
  # Compute the euclidean distance matrix
  distance_mat <- apply(df, 1, function(pt)
    (sqrt(abs(pt["x"] - df$x)^2 + abs(pt["y"] - df$y)^2))
  )
  # keep the cell names
  rownames(distance_mat) <- colnames(distance_mat)
  
  if(!is.null(k)){
    message("Compute edgelist based on mutual nearest neighbors.")
    if(!is.null(rad.set))warning("'k' is not NULL. Parameter 'rad.set' will be ignored.")
    # neighbor_mat: the indices of k nearest neighbors  
    neighbor_mat <- apply(distance_mat,1,function(dis_vec){
      order(dis_vec)[1:(k+1)]
    })
    
    # k nearest neighbor adjacency matrix
    adj_mat <- matrix(data=0,nrow=nrow(distance_mat),ncol = ncol(distance_mat))
    rownames(adj_mat) <- rownames(distance_mat)
    colnames(adj_mat) <- colnames(distance_mat)
    for(cell in colnames(neighbor_mat)) adj_mat[cell,neighbor_mat[,cell]] <- 1
    # mutual nearest neighbor adjacency matrix
    adj_mat_final <- 1*(adj_mat & t(adj_mat))
    
  }
  else if(!is.null(rad.set)){
    message("\n Compute edgelist based on spatial radius threshold.")
    # Alternatively the adjacency matrix is defined by the absolute radius
    adj_mat_final <- 1*(distance_mat <= rad.set)
  }
  else{
    stop("Both k and rad.set are NULL.")
  }
  
  # Convert adj matrix to edgelist
  if(weighted == TRUE){
    edgelist <- igraph::graph_from_adjacency_matrix(distance_mat*adj_mat_final, weighted = TRUE)
    edgelist <- igraph::as_data_frame(edgelist)
  }
  else{
    edgelist <- igraph::graph_from_adjacency_matrix(adj_mat_final)
    edgelist <- igraph::as_data_frame(edgelist)
  }
  
  return(edgelist)
  
}

compute_Neighborhood_Composition <- function(obj, 
                                             position.x, 
                                             position.y, 
                                             k = 4, 
                                             rad.set = NULL,
                                             composition)
{
  edge.list <- compute_edgelist(obj, 'x', 'y', k, rad.set)
  
  base.df.from <- data.frame(from = colnames(obj),
                             type = obj[[composition]])
  
  base.df.to <- data.frame(to = colnames(obj),
                           type = obj[[composition]])
  
  
  edge.list <- dplyr::left_join(edge.list, base.df.from, by = 'from')
  names(edge.list)[names(edge.list) == composition] <- "SendingType"
  
  edge.list <- dplyr::left_join(edge.list, base.df.to, by = 'to')
  names(edge.list)[names(edge.list) == composition] <- "ReceivingType"
  
  # Receiving neighborhood calculation -- raw counts version
  receiving_counts_df <- edge.list %>%
    group_by(from, ReceivingType) %>%
    summarise(count = n(), .groups = 'drop') %>%
    group_by(from) %>%
    mutate(rawcounts = count) %>%
    ungroup() %>%
    select(from, ReceivingType, rawcounts) 
  
  
  receiving_wide_counts_df <- receiving_counts_df %>%
    tidyr::pivot_wider(names_from = ReceivingType,
                       values_from = rawcounts,
                       names_prefix = paste0(composition,"_"),
                       values_fill = list(rawcounts = 0))
  
  
  # Receiving neighborhood calculation -- percentage version
  receiving_fraction_df <- edge.list %>%
    group_by(from, ReceivingType) %>%
    summarise(count = n(), .groups = 'drop') %>%
    group_by(from) %>%
    mutate(total = sum(count),
           fraction = (count / total)) %>%
    ungroup() %>%
    select(from, ReceivingType, fraction)  # Select relevant columns
  
  # Step 2: Pivot data to wide format
  receiving_wide_fraction_df <- receiving_fraction_df %>%
    tidyr::pivot_wider(names_from = ReceivingType,
                       values_from = fraction,
                       names_prefix = paste0(composition,"_"),
                       values_fill = list(fraction = 0))
  
  
  # Sam found that "+" sign will cause problem in further analysis so let's delete them here.
  new_count_colnames <- gsub("\\+", "", colnames(receiving_wide_counts_df))
  # Assign new column names to the data frame
  colnames(receiving_wide_counts_df) <- new_count_colnames
  
  new_frac_colnames <- gsub("\\+", "", colnames(receiving_wide_fraction_df))
  # Assign new column names to the data frame
  colnames(receiving_wide_fraction_df) <- new_frac_colnames
  
  return(list(raw.df = receiving_wide_counts_df, pct.df = receiving_wide_fraction_df))
}


#### This is how I test writing my function ######

# edge.list <- compute_edgelist(favorite.data, 'x', 'y', k = NULL, rad.set = 25)
# 
# 
# 
# base.df.from <- data.frame(from = colnames(favorite.data),
#                            type = favorite.data[[decomposition]])
# 
# base.df.to <- data.frame(to = colnames(favorite.data),
#                          type = favorite.data[[decomposition]])
# 
# 
# edge.list <- dplyr::left_join(edge.list, base.df.from, by = 'from')
# names(edge.list)[names(edge.list) == decomposition] <- "SendingType"
# 
# edge.list <- dplyr::left_join(edge.list, base.df.to, by = 'to')
# names(edge.list)[names(edge.list) == decomposition] <- "ReceivingType"
# 
# # Receiving neighborhood calculation
# receiving_fraction_df <- edge.list %>%
#   group_by(from, ReceivingType) %>%
#   summarise(count = n(), .groups = 'drop') %>%
#   group_by(from) %>%
#   mutate(total = sum(count),
#          fraction = (count / total)) %>%
#   ungroup() %>%
#   select(from, ReceivingType, fraction)  # Select relevant columns
# 
# # Step 2: Pivot data to wide format
# receiving_wide_fraction_df <- receiving_fraction_df %>%
#   tidyr::pivot_wider(names_from = ReceivingType,
#                      values_from = fraction,
#                      names_prefix = "ReceivingType_",
#                      values_fill = list(fraction = 0))
# 
# 
# sending_fraction_df <- edge.list %>%
#   group_by(to, SendingType) %>%
#   summarise(count = n(), .groups = 'drop') %>%
#   group_by(to) %>%
#   mutate(total = sum(count),
#          fraction = (count / total)) %>%
#   ungroup() %>%
#   select(to, SendingType, fraction)  # Select relevant columns
# 
# # Step 2: Pivot data to wide format
# sending_wide_fraction_df <- sending_fraction_df %>%
#   tidyr::pivot_wider(names_from = SendingType,
#                      values_from = fraction,
#                      names_prefix = "SendingType_",
#                      values_fill = list(fraction = 0))
# 
# rownames(sending_wide_fraction_df) <- sending_wide_fraction_df$to
# sending_wide_fraction_df <- sending_wide_fraction_df[,-1]
# 
# Neighborhood.Sending.decomposition <- as.matrix(sending_wide_fraction_df)
# rownames(Neighborhood.Sending.decomposition) <- rownames(sending_wide_fraction_df)
# 
# dim(Neighborhood.Sending.decomposition)
# head(rownames(Neighborhood.Sending.decomposition))
