library(Seurat)
library(ggplot2)
library(dplyr)

# setwd('/Volumes/T7/Pneum/Xenium/Clean_Objects')
# 
# load("/Volumes/T7/Pneum/Xenium/PPLR.annotated_Type.5_SEE_08-20-2025.Robj")
# 
# PPLR.CTC <- readRDS("/Volumes/T7/Pneum/Xenium/NICHES_test/PNX0_F1_Ad.CTC.byRad.normalized.rds")

AddMetaData <- function(CTC.obj,
                        annotated.obj,
                        col.name){
  
  CTC.edges <- rownames(CTC.obj@meta.data)
  split_list <- strsplit(CTC.edges, "—")
  
  # Use 'sapply' to extract the sending and receiving parts into two separate vectors
  sending.cell <- sapply(split_list, `[`, 1)
  receiving.cell <- sapply(split_list, `[`, 2)
  
  # Create the CTC meta data so that we can join
  CTC.meta <- CTC.obj@meta.data
  CTC.meta$edge <- CTC.edges
  CTC.meta$sending.cell <- sending.cell
  CTC.meta$receiving.cell <- receiving.cell
  
  # Extract the new meta data column in the annotated obj and create the necessary data frame for joining.
  temp.df <- data.frame(cell.name = rownames(annotated.obj@meta.data),
                        new.col = annotated.obj@meta.data[[col.name]])
  
  # We do not like factor variables so let's change it to character if needed.
  if(class(temp.df$new.col) == 'factor'){
    temp.df$new.col <- as.character(temp.df$new.col)
  }
  
  # Add sending cell meta data
  sending.joined <- CTC.meta %>%
    left_join(temp.df, by = c("sending.cell" = "cell.name")) %>%
    rename(new.col.sending = new.col)
  
  # Add receiving cell meta data
  receiving.joined <- sending.joined %>%
    left_join(temp.df, by = c("receiving.cell" = "cell.name")) %>%
    rename(new.col.receiving = new.col)
  
  # Some naming change. Allow user-defined column names.
  new.sendingName <- paste0(col.name,'.Sending')
  new.receivingName <- paste0(col.name,'.Receiving')
  new.jointName <- paste0(col.name, '.Joint')
  
  receiving.joined[[new.sendingName]] <- receiving.joined[["new.col.sending"]]
  receiving.joined[[new.receivingName]] <- receiving.joined[['new.col.receiving']]
  receiving.joined[[new.jointName]] <- paste0(receiving.joined[[new.sendingName]],' - ',receiving.joined[[new.receivingName]])
  
  receiving.joined[["new.col.sending"]] <- NULL
  receiving.joined[["new.col.receiving"]] <- NULL
  
  # Add a joint column
  
  
  # return a new data frame which can be used to update the existing Seurat object.
  return(receiving.joined)
  
}

# An example on how to use the function.
# col.name <- 'Type.5'
# new.meta.data <- AddMetaData(PPLR.CTC, PPLR.annotated, col.name)
