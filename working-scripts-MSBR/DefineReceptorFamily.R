DefineReceptorFamily <- function(feature.metadata){
  ###### RECEPTOR.FAMILY ######
  ## Growth Hormone Family Receptor
  indices.temp <- grep(pattern = '^GHR',
                       x = feature.metadata$RECEPTOR,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$RECEPTOR.FAMILY <- 'Growth Hormone Family Receptor'
  ## Glucagon-like peptide receptor
  indices.temp <- grep(pattern = '^GLP',
                       x = feature.metadata$RECEPTOR,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$RECEPTOR.FAMILY <- 'Glucagon-like peptide receptor'
  ## LPA Receptor
  indices.temp <- grep(pattern = '^LPA',
                       x = feature.metadata$RECEPTOR,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$RECEPTOR.FAMILY <- 'LPA receptor'
  ## Asialoglycoprotein receptor
  indices.temp <- grep(pattern = '^ASGR',
                       x = feature.metadata$RECEPTOR,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$RECEPTOR.FAMILY <- 'Asialoglycoprotein receptor'
  ## GDNF Family Receptor
  indices.temp <- grep(pattern = '^RET',
                       x = feature.metadata$RECEPTOR,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$RECEPTOR.FAMILY <- 'GDNF Family Receptor'

  ## Prokineticin family receptor
  indices.temp <- grep(pattern = '^PROK',
                       x = feature.metadata$RECEPTOR,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$RECEPTOR.FAMILY <- 'Prokineticin family receptor'
  ## Guanylate cyclase
  indices.temp <- grep(pattern = '^GUCY',
                       x = feature.metadata$RECEPTOR,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$RECEPTOR.FAMILY <- 'Guanylate cyclase'
  ## Prolactin Family Receptor
  indices.temp <- grep(pattern = '^PRL',
                       x = feature.metadata$RECEPTOR,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$RECEPTOR.FAMILY <- 'Prolactin Family Receptor'
  ## TAM Family Kinase
  indices.temp <- grep(pattern = '^AXL$|^MERTK$|^TYRO3$',
                       x = feature.metadata$RECEPTOR,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$RECEPTOR.FAMILY <- 'TAM Family Kinase'
  ## Formyl peptide receptor
  indices.temp <- grep(pattern = '^FPR',
                       x = feature.metadata$RECEPTOR,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$RECEPTOR.FAMILY <- 'Formyl peptide receptor'
  ## RAMP
  indices.temp <- grep(pattern = '^RAMP',
                       x = feature.metadata$RECEPTOR,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$RECEPTOR.FAMILY <- 'RAMP'
  ## Cell Adhesion Molecule
  indices.temp <- grep(pattern = 'CAM',
                       x = feature.metadata$RECEPTOR,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$RECEPTOR.FAMILY <- 'Cell Adhesion Molecule'
  ## Syndecan
  indices.temp <- grep(pattern = '^SDC',
                       x = feature.metadata$RECEPTOR,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$RECEPTOR.FAMILY <- 'Syndecan'
  ## ABC Transporter
  indices.temp <- grep(pattern = '^ABC',
                       x = feature.metadata$RECEPTOR,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$RECEPTOR.FAMILY <- 'ABC Transporter'
  ## Parathyroid Hormone Receptor
  indices.temp <- grep(pattern = '^PTH',
                       x = feature.metadata$RECEPTOR,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$RECEPTOR.FAMILY <- 'Parathyroid Hormone Receptor'
  ## Neuropeptide Receptor
  indices.temp <- grep(pattern = '^NPY|^GALR',
                       x = feature.metadata$RECEPTOR,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$RECEPTOR.FAMILY <- 'NPY'
  ## IGF
  indices.temp <- grep(pattern = '^IGF',
                       x = feature.metadata$RECEPTOR,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$RECEPTOR.FAMILY <- 'IGF'
  ## Solute Carrier Family
  indices.temp <- grep(pattern = '^SLC',
                       x = feature.metadata$RECEPTOR,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$RECEPTOR.FAMILY <- 'Solute Carrier Family'

  ## TNFSF Receptor Subunit
  indices.temp <- grep(pattern = '^TNFR|^CD27$',
                       x = feature.metadata$RECEPTOR,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$RECEPTOR.FAMILY <- 'TNFSF Receptor Subunit'
  ## Opiod Receptor Subunit
  indices.temp <- grep(pattern = '^OPR',
                       x = feature.metadata$RECEPTOR,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$RECEPTOR.FAMILY <- 'Opiod Receptor Subunit'

  ## Interleukin Receptor Subunit
  indices.temp <- grep(pattern = '^IL',
                       x = feature.metadata$RECEPTOR,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$RECEPTOR.FAMILY <- 'Interleukin Receptor Subunit'
  ## G-Protein Receptor
  indices.temp <- grep(pattern = '^GPR',
                       x = feature.metadata$RECEPTOR,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$RECEPTOR.FAMILY <- 'G-Protein Receptor'
  ## Cadherin
  indices.temp <- grep(pattern = '^CDH',
                       x = feature.metadata$RECEPTOR,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$RECEPTOR.FAMILY <- 'Cadherin'
  ## Neuropeptide Receptor
  indices.temp <- grep(pattern = '^NPFFR',
                       x = feature.metadata$RECEPTOR,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$RECEPTOR.FAMILY <- 'Neuropeptide Receptor'
  ## LRP Family
  indices.temp <- grep(pattern = '^LRP',
                       x = feature.metadata$RECEPTOR,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$RECEPTOR.FAMILY <- 'LRP Family'
  ## Insulin Receptor
  indices.temp <- grep(pattern = '^INSR',
                       x = feature.metadata$RECEPTOR,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$RECEPTOR.FAMILY <- 'Insulin Receptor'

  ## CXC Family Receptor
  indices.temp <- grep(pattern = '^CXC',
                       x = feature.metadata$RECEPTOR,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$RECEPTOR.FAMILY <- 'CXC Family Receptor'

  ## CC Family Receptor
  indices.temp <- grep(pattern = '^CCR',
                       x = feature.metadata$RECEPTOR,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$RECEPTOR.FAMILY <- 'CC Family Receptor'


  ## Neurexin
  indices.temp <- grep(pattern = '^NRX',
                       x = feature.metadata$RECEPTOR,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$RECEPTOR.FAMILY <- 'Neurexin'


  ## EGF Family Receptor
  indices.temp <- grep(pattern = '^ERBB|EGFR',
                       x = feature.metadata$RECEPTOR,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$RECEPTOR.FAMILY <- 'EGF Family Receptor'

  ## Integrins
  indices.temp <- grep(pattern = 'ITG',
                       x = feature.metadata$RECEPTOR,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$RECEPTOR.FAMILY <- 'Integrin'

  ## NOTCH
  indices.temp <- grep(pattern = 'NOTCH',
                       x = feature.metadata$RECEPTOR,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$RECEPTOR.FAMILY <- 'NOTCH Family Receptor'

  ## TLR
  indices.temp <- grep(pattern = 'TLR',
                       x = feature.metadata$RECEPTOR,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$RECEPTOR.FAMILY <- 'Toll-Like Receptor'

  ## Lipid Handling
  indices.temp <- grep(pattern = 'LDL',
                       x = feature.metadata$RECEPTOR,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$RECEPTOR.FAMILY <- 'Lipid Handling'


  ## Reticulon Family Receptor
  indices.temp <- grep(pattern = 'RTN',
                       x = feature.metadata$RECEPTOR,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$RECEPTOR.FAMILY <- 'Reticulon Family Receptor'

  ## Ephrin Family
  indices.temp <- grep(pattern = '^EPH',
                       x = feature.metadata$RECEPTOR,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$RECEPTOR.FAMILY <- 'Ephrin Family Receptor'

  ## WNT Family
  indices.temp <- grep(pattern = '^LGR|FZD|RYK',
                       x = feature.metadata$RECEPTOR,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$RECEPTOR.FAMILY <- 'WNT Family Receptor'

  ## BMP Family
  indices.temp <- grep(pattern = '^BMPR',
                       x = feature.metadata$RECEPTOR,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$RECEPTOR.FAMILY <- 'BMP Family Receptor'

  ## Collagen
  indices.temp <- grep(pattern = '^COL',
                       x = feature.metadata$RECEPTOR,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$RECEPTOR.FAMILY <- 'Collagen'

  ## FGF Family
  indices.temp <- grep(pattern = '^FGFR',
                       x = feature.metadata$RECEPTOR,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$RECEPTOR.FAMILY <- 'FGF Family Receptor'
  ## NGF Family
  indices.temp <- grep(pattern = '^NGF|^NTRK',
                       x = feature.metadata$RECEPTOR,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$RECEPTOR.FAMILY <- 'NGF Family Receptor'

  ## GDNF Family Receptor
  indices.temp <- grep(pattern = '^GFR', #GFRA = GDNF receptor https://en.wikipedia.org/wiki/GFRA1
                       x = feature.metadata$RECEPTOR,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$RECEPTOR.FAMILY <- 'GDNF Family Receptor'

  ## EGF Family
  indices.temp <- grep(pattern = '^EGF|^EREG|^AREG',
                       x = feature.metadata$RECEPTOR,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$RECEPTOR.FAMILY <- 'EGF Family Receptor'
  ## PDGF Family
  indices.temp <- grep(pattern = '^PDGF',
                       x = feature.metadata$RECEPTOR,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$RECEPTOR.FAMILY <- 'PDGF Family Receptor'
  ## TGFBR
  indices.temp <- grep(pattern = '^TGFBR',
                       x = feature.metadata$RECEPTOR,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$RECEPTOR.FAMILY <- 'TGFB Family Receptor'

  ## ANGPT
  indices.temp <- grep(pattern = '^TIE|^TEK',
                       x = feature.metadata$RECEPTOR,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$RECEPTOR.FAMILY <- 'ANGPT Family Receptor'
  ## Plexin
  indices.temp <- grep(pattern = '^PLXN',
                       x = feature.metadata$RECEPTOR,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$RECEPTOR.FAMILY <- 'Plexin'

  ## Selectin
  indices.temp <- grep(pattern = '^SEL',
                       x = feature.metadata$RECEPTOR,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$RECEPTOR.FAMILY <- 'Selectin'

  ## Activin
  indices.temp <- grep(pattern = '^ACVR',
                       x = feature.metadata$RECEPTOR,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$RECEPTOR.FAMILY <- 'Activin'

  return(feature.metadata)
}
