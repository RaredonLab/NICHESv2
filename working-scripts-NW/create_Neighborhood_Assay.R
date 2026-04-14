library(Seurat)
library(ggplot2)
library(dplyr)

create_New_Assay <- function(obj, 
                             position.x = 'x',
                             position.y = 'y',
                             k = NULL,
                             rad.set = 25,
                             composition = composition,
                             TMA = NULL,
                             TMA.list = NULL){
  
  if(is.null(TMA) || is.null(TMA.list)){
    Neighborhood_composition <- compute_Neighborhood_Composition(obj, 
                                                                 position.x, 
                                                                 position.y, 
                                                                 k, 
                                                                 rad.set, 
                                                                 composition)
    
    rawcounts <- Neighborhood_composition$raw.df
    pct <- Neighborhood_composition$pct.df
    
    print('Working on whole data, did not split it by sample/TMA. Want to double check?')
    
    rownames(rawcounts) <- rawcounts$from
    rownames(pct) <- pct$from
    
    Neighborhood_composition.mat.rawcounts <- t(as.matrix(rawcounts)[,-1])
    Neighborhood_composition.mat.rawcounts[is.na(Neighborhood_composition.mat.rawcounts)] <- 0
    
    Neighborhood_composition.mat.pct <- t(as.matrix(pct)[,-1])
    Neighborhood_composition.mat.pct[is.na(Neighborhood_composition.mat.pct)] <- 0
    
    assay.rawcounts <- CreateAssay5Object(counts = Neighborhood_composition.mat.rawcounts)
    assay.pct <- CreateAssay5Object(counts = Neighborhood_composition.mat.pct)
    return(list(rawcounts = assay.rawcounts, pct = assay.pct))
  }
  
  else{
    
    for(i in 1:length(TMA.list)){
      cur.TMA <- TMA.list[i] # Take the current TMA
      print(paste0('Working on ', i, 'th sample, sample name: ', cur.TMA))
      df <- data.frame(TMA = obj[[TMA]], cell.id = rownames(obj[[TMA]]))
      rownames(df) <- NULL
      selected <- df[df[[TMA]] == cur.TMA,]$cell.id
      cur.data <- subset(obj, cells = selected) # Subset the current TMA
      
      cur.Neighborhood_Composition <- compute_Neighborhood_Composition(obj = cur.data,
                                                                       position.x,
                                                                       position.y,
                                                                       k,
                                                                       rad.set,
                                                                       composition)
      cur.rawcounts <- cur.Neighborhood_Composition$raw.df
      cur.pct <- cur.Neighborhood_Composition$pct.df
      
      if(i == 1){
        concatenated.Neighborhood_Composition.rawcounts <- cur.rawcounts
        concatenated.Neighborhood_Composition.pct <- cur.pct
      }
      else{
        concatenated.Neighborhood_Composition.rawcounts <- bind_rows(concatenated.Neighborhood_Composition.rawcounts,cur.rawcounts)
        concatenated.Neighborhood_Composition.pct <- bind_rows(concatenated.Neighborhood_Composition.pct,cur.pct )
      }
    }
    
    rownames(concatenated.Neighborhood_Composition.rawcounts) <- concatenated.Neighborhood_Composition.rawcounts$from
    rownames(concatenated.Neighborhood_Composition.pct) <- concatenated.Neighborhood_Composition.pct$from
    
    Neighborhood_Composition.mat.rawcounts <- t(as.matrix(concatenated.Neighborhood_Composition.rawcounts)[,-1])
    Neighborhood_Composition.mat.rawcounts[is.na(Neighborhood_Composition.mat.rawcounts)] <- 0
    
    Neighborhood_Composition.mat.pct <- t(as.matrix(concatenated.Neighborhood_Composition.pct)[,-1])
    Neighborhood_Composition.mat.pct[is.na(Neighborhood_Composition.mat.pct)] <- 0
    
    
    assay.rawcounts <- CreateAssay5Object(counts = Neighborhood_Composition.mat.rawcounts)
    assay.pct <- CreateAssay5Object(counts = Neighborhood_Composition.mat.pct)
    
    return(list(rawcounts = assay.rawcounts, pct = assay.pct))
  }
  
}


#### How to use this function #### 

# setwd('/Volumes/T7/Pneum/Xenium/Test_Composition')
# 
# load("/Volumes/T7/Pneum/Xenium/PPLR.annotated_Type.5_SEE_08-20-2025.Robj")
# 
# source('/Volumes/T7/Pneum/Xenium/Test_Composition/CellNeighborhood_Composition.R')
# 
# 
# PPLR.annotated$Type.5 <- as.character(PPLR.annotated$Type.5)
# TMA.list <- names(PPLR.annotated@images)[-1]
# 
# 
# NeighborhoodComposition.byRad25 <- create_New_Assay(obj = PPLR.annotated,
#                                                     position.x = 'x',
#                                                     position.y = 'y',
#                                                     k = NULL,
#                                                     rad.set = 25,
#                                                     composition = 'Type.5',
#                                                     TMA = 'TMA',
#                                                     TMA.list = TMA.list)
# 
# NeighborhoodComposition.byKNN9 <- create_New_Assay(obj = PPLR.annotated,
#                                                    position.x = 'x',
#                                                    position.y = 'y',
#                                                    k = 9,
#                                                    rad.set = NULL,
#                                                    composition = 'Type.5',
#                                                    TMA = 'TMA',
#                                                    TMA.list = TMA.list)
# 
# saveRDS(NeighborhoodComposition.byRad25, file = 'NeighborhoodComposition.byRad25.Type5.rds')
# saveRDS(NeighborhoodComposition.byKNN9, file = 'NeighborhoodComposition.byKNN9.Type5.rds')