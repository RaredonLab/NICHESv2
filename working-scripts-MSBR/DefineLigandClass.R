DefineLigandClass <- function(feature.metadata){

  ###### LIGAND.CLASS    #####

  ## Neuroproteins
  indices.temp <- grep(pattern = 'Neuropeptide|Neuromedin',
                       x = feature.metadata$LIGAND.FAMILY,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$LIGAND.CLASS <- 'Neuroproteins'

  ## Catalysis
  indices.temp <- grep(pattern = 'Presenilin',
                       x = feature.metadata$LIGAND.FAMILY,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$LIGAND.CLASS <- 'Catalysis'

  ## Cytosolic (Non-Ligands)
  indices.temp <- grep(pattern = 'Cingulin|Intracellular|G Protein Complex|LRPAP|Heat-shock Protein|Myocilin|Intracellular Enzyme|Intracellular (uncategorized)',
                       x = feature.metadata$LIGAND.FAMILY,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$LIGAND.CLASS <- 'Cytosolic (Non-Ligands)'

  ## Protease Inhibitors
  indices.temp <- grep(pattern = 'Serpin|Alpha Globulin',
                       x = feature.metadata$LIGAND.FAMILY,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$LIGAND.CLASS <- 'Protease Inhibitors'

  ## Kinin-Kallikrein System
  # indices.temp <- grep(pattern = 'Kininogen',
  #                      x = feature.metadata$LIGAND.FAMILY,
  #                      ignore.case = FALSE)
  # feature.metadata[indices.temp,]$LIGAND.CLASS <- 'Kinin-Kallikrein System'

  ## Cell-Cell Contact
  indices.temp <- grep(pattern = 'Cadherin|Cell Adhesion|Neuroligin|Contactin|Neurexophilin|Amyloid Beta Precursor|Selectin Ligand',
                       x = feature.metadata$LIGAND.FAMILY,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$LIGAND.CLASS <- 'Cell-Cell Contact'

  ## Spatial Guidance
  indices.temp <- grep(pattern = 'Semaphorin|SLIT|Netrin|Ephrin|Reticulon|Repulsive guidance molecule',
                       x = feature.metadata$LIGAND.FAMILY,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$LIGAND.CLASS <- 'Spatial Guidance'

  ## Metabolic
  indices.temp <- grep(pattern = 'Growth Hormone|Ghrelin|Adiponectin|Apolipoprotein|Lipoprotein|Tryptophan hydroxylase|^PPY$',
                       x = feature.metadata$LIGAND.FAMILY,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$LIGAND.CLASS <- 'Metabolic'

  ## Calcium
  indices.temp <- grep(pattern = 'Calmodulin|Calreticulin|Calcitonin',
                       x = feature.metadata$LIGAND.FAMILY,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$LIGAND.CLASS <- 'Calcium'

  ## Blood Homeostasis
  indices.temp <- grep(pattern = 'Complement|VWF|Transferrin|Coagulation',
                       x = feature.metadata$LIGAND.FAMILY,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$LIGAND.CLASS <- 'Blood Homeostasis'

  ## Matrix Molecules
  indices.temp <- grep(pattern = 'Extracellular Enzyme|Plasmin|Reelin|ECM (misc)|Fibulin|Fibrillin|MFAP|Osteopontin|TIMP|Collagen|Laminin|Fibronectin|Vitronectin|BGN|Decorin|TNC|Agrin|CTGF|CYR61|Versican|Neuregulin|Nephronectin|Thrombospondin|CTHRC1|Nidogen',
                       x = feature.metadata$LIGAND.FAMILY,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$LIGAND.CLASS <- 'Matrix Molecules'

  ## Hormones
  indices.temp <- grep(pattern = 'Vasopressin|Natriuretic|Adrenomedullin|Endothelin|Agouti-related|Renin-Angiotensin|VIP|Urocortin|Relaxin|Sex Hormone Binding Protein',
                       x = feature.metadata$LIGAND.FAMILY,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$LIGAND.CLASS <- 'Hormones'

  ## Antigen Presentation Type
  indices.temp <- grep(pattern = 'Antigen Presentation',
                       x = feature.metadata$LIGAND.FAMILY,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$LIGAND.CLASS <- 'Antigen Presentation'

  ## GF Type
  indices.temp <- grep(pattern = 'KIT|Apelin|WNT|CSF|PDGF|VEGF|NGF|PGF|BMP|TGFB|FGF|IGF|EGF|GDF|ANGPT|BDNF|HGF|GDNF|NEGF|CNTF|NOTCH|Activin|Inhibin|Hedgehog',
                       x = feature.metadata$LIGAND.FAMILY,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$LIGAND.CLASS <- 'Growth Factor'

  ## Cytokine Category
  indices.temp <- grep(pattern = 'Interleukin|Interferon|C chemokine|CCL|CXC|TNF|Oncostatin M|Calprotectin|Cathelicidin',
                       x = feature.metadata$LIGAND.FAMILY,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$LIGAND.CLASS <- 'Cytokine'

  ## Difficult to categorize
  indices.temp <- grep(pattern = 'Prosaposin',
                       x = feature.metadata$LIGAND.FAMILY,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$LIGAND.CLASS <- 'Difficult to Categorize'

 return(feature.metadata)
}
