DefineReceptorClass <- function(feature.metadata){

  ###### RECEPTOR.CLASS    #####

  ## Plexin
  indices.temp <- grep(pattern = 'Plexin',
                       x = feature.metadata$RECEPTOR.FAMILY,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$RECEPTOR.CLASS <- 'Plexin'

  ## Activin
  indices.temp <- grep(pattern = 'Activin',
                       x = feature.metadata$RECEPTOR.FAMILY,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$RECEPTOR.CLASS <- 'Activin'

  ## Toll-Like Interactions
  indices.temp <- grep(pattern = 'Toll-Like Receptor',
                       x = feature.metadata$RECEPTOR.FAMILY,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$RECEPTOR.CLASS <- 'Toll-Like Receptor'

  ## Ephrin Interactions
  indices.temp <- grep(pattern = 'Ephrin',
                       x = feature.metadata$RECEPTOR.FAMILY,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$RECEPTOR.CLASS <- 'Ephrin Family Receptor'

  ## Integrin Interactions
  indices.temp <- grep(pattern = 'Integrin',
                       x = feature.metadata$RECEPTOR.FAMILY,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$RECEPTOR.CLASS <- 'Integrin'

  ## Growth Factor Receptor
  indices.temp <- grep(pattern = 'WNT|CSF|PDGF|VEGF|NGF|PGF|BMP|TGFB|FGF|IGF|EGF|GDF|ANGPT|BDNF|HGF|GDNF|NEGF|CNTF',
                       x = feature.metadata$RECEPTOR.FAMILY,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$RECEPTOR.CLASS <- 'Growth Factor Receptor'

  ## Cytokine Receptor
  indices.temp <- grep(pattern = '^CXC Family|^CC Family',
                       x = feature.metadata$RECEPTOR.FAMILY,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$RECEPTOR.CLASS <- 'Cytokine Receptor'
  ## Metabolic
  indices.temp <- grep(pattern = 'Lipid',
                       x = feature.metadata$RECEPTOR.FAMILY,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$RECEPTOR.CLASS <- 'Metabolic'


  ## Definition Level (for MSBR reference)
  temp <- feature.metadata[,c('LIGAND.CLASS',
                              'RECEPTOR.CLASS',
                              "LIGAND.FAMILY",
                              "RECEPTOR.FAMILY")]
  temp <- !is.na(temp)

  feature.metadata$DEFINITION.LEVEL <- as.character(rowSums(temp))

  return(feature.metadata)
}
