CoLocalization <- function(object,
                       neighborhood.assay = "Neighborhood_Type5_Rad25",
                       group.by = 'Type.6',
                       split.by = 'TMA',
                       tag.by = c('Timepoint')){
  
  
  DefaultAssay(PPLR.master) <- neighborhood.assay
  binary <- PPLR.master@assays[[neighborhood.assay]]@layers$counts > 0
  
  coi <- PPLR.master$Type.6 == 'ATI' & PPLR.master$Timepoint == 'Day 0'
  temp <- binary[,coi] # pull just the neighborhoods of interest
  out <- data.frame(Neighborhood = 'ATI',
                    Probability = rowSums(temp)/ncol(temp),
                    Feature = rownames(PPLR.master),
                    Timepoint = 'Day 0')
  
  coi <- PPLR.master$Type.6 == 'ATI' & PPLR.master$Timepoint == 'Day 3'
  temp <- binary[,coi] # pull just the neighborhoods of interest
  out.2 <- data.frame(Neighborhood = 'ATI',
                      Probability = rowSums(temp)/ncol(temp),
                      Feature = rownames(PPLR.master),
                      Timepoint = 'Day 3')
  
  coi <- PPLR.master$Type.6 == 'ATI' & PPLR.master$Timepoint == 'Day 7'
  temp <- binary[,coi] # pull just the neighborhoods of interest
  out.3 <- data.frame(Neighborhood = 'ATI',
                      Probability = rowSums(temp)/ncol(temp),
                      Feature = rownames(PPLR.master),
                      Timepoint = 'Day 7')
  
  coi <- PPLR.master$Type.6 == 'ATI' & PPLR.master$Timepoint == 'Day 14'
  temp <- binary[,coi] # pull just the neighborhoods of interest
  out.4 <- data.frame(Neighborhood = 'ATI',
                      Probability = rowSums(temp)/ncol(temp),
                      Feature = rownames(PPLR.master),
                      Timepoint = 'Day 14')
  ati <- rbind(out,out.2,out.3,out.4)
}

