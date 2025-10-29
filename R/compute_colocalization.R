library(Seurat)
library(ggplot2)
library(dplyr)

# setwd('/Volumes/T7/Pneum/Xenium/Test_Composition')
# load("/Volumes/T7/Pneum/Xenium/PPLR.annotated_Type.5_SEE_08-20-2025.Robj")
# Neighborhood.Assay <- readRDS('NeighborhoodComposition.byRad25.Type5.rds')

compute_colocalization <- function(Neighborhood.Assay, annotated.obj, col.name){
  
  test <- as.data.frame(t(Neighborhood.Assay$rawcounts$counts))
  test$cell.name <- rownames(test)
  
  celltype.info <- data.frame(cell.name = rownames(annotated.obj@meta.data),
                              new.col = annotated.obj@meta.data[[col.name]])
  
  if(class(celltype.info$new.col) == 'factor'){
    celltype.info$new.col <- as.character(celltype.info$new.col)
  }
  
  celltype.joined <- test %>%
    left_join(celltype.info, by = c("cell.name" = "cell.name"))
  
  # List all feature columns except 'cell' and 'cell_type'
  feature_cols <- setdiff(names(celltype.joined), c("cell.name", "new.col"))
  
  # Method 1: Binarize the values (>0 becomes 1) and then aggregate by cell type
  df_binary_aggregated <- celltype.joined %>%
    mutate(across(all_of(feature_cols), ~ ifelse(. > 0, 1, 0))) %>%
    group_by(new.col) %>%
    summarise(across(all_of(feature_cols), ~ sum(.x, na.rm = TRUE)))
  
  # Method 2: Aggregate the original values by cell type
  df_aggregated <- celltype.joined %>%
    group_by(new.col) %>%
    summarise(across(all_of(feature_cols), ~ sum(.x, na.rm = TRUE)))
  
  # Some naming style change
  new.Name <- paste0(col.name, '.orig')
  df_binary_aggregated[[new.Name]] <- df_binary_aggregated[["new.col"]]
  df_aggregated[[new.Name]] <- df_aggregated[['new.col']]
  
  df_binary_aggregated[['new.col']] <- NULL
  df_aggregated[['new.col']] <- NULL
  
  # Print the results
  return(list(binary_aggregated = df_binary_aggregated, aggregated = df_aggregated))
}


# output <- compute_colocalization(Neighborhood.Assay, PPLR.annotated, col.name = 'Type.5')





