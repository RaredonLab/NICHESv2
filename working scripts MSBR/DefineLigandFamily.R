DefineLigandFamily <- function(feature.metadata){


  ###### LIGAND.FAMILY  ######

  ## Cytokines
  ### Chemokines

  #### CCL Family
  indices.temp <- grep(pattern = 'CCL',
                       x = feature.metadata$LIGAND,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'CCL Family'
  #### CXC Family
  indices.temp <- grep(pattern = '^CXC|^CX3|^PF4',
                       x = feature.metadata$LIGAND,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'CXC Family'
  #### C Chemokine Family (XCL Family)
  indices.temp <- grep(pattern = '^XCL1$|^XCL2$',
                       x = feature.metadata$LIGAND,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'C Chemokine Family (XCL Family)'
  #### Colony Stimulating Factors (CSF Family)
  indices.temp <- grep(pattern = '^CSF',
                       x = feature.metadata$LIGAND,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Colony Stimulating Factors (CSF Family)'

  ### Interferons
  #### IFN Family
  indices.temp <- grep(pattern = '^IFN',
                       x = feature.metadata$LIGAND,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'IFN Family'

  ### Interleukins
  #### CLCF Family
  indices.temp <- grep(pattern = '^CLCF1$',
                       x = feature.metadata$LIGAND,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'CLCF Family'
  #### IL Family
  indices.temp <- grep(pattern = '^IL',
                       x = feature.metadata$LIGAND,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'IL Family'
  #### Oncostatin M
  indices.temp <- grep(pattern = '^OSM',
                       x = feature.metadata$LIGAND,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Oncostatin M'
  #### Leukemia Inhibitory Factor

  ### Tumor Necrosis Factors
  #### TNF Family
  #### TNFSF Family

  ### Lymphoid Interaction Factors
  #### HLA Complex Molecules
  indices.temp <- grep(pattern = '^HLA-|^B2M',
                       x = feature.metadata$LIGAND,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'HLA Complex Molecules'
  #### CD40LG
  indices.temp <- grep(pattern = '^CD40LG$',
                       x = feature.metadata$LIGAND,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'CD40LG'
  #### BTLA
  indices.temp <- grep(pattern = '^BTLA$',
                       x = feature.metadata$LIGAND,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'BTLA'

  ### Miscellaneous
  #### "Apoptosis inhibitor expressed by macrophages"
  indices.temp <- grep(pattern = '^CD5L',
                       x = feature.metadata$LIGAND,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Apoptosis inhibitor expressed by macrophages'
  #### Calprotectin
  indices.temp <- grep(pattern = '^S100',
                       x = feature.metadata$LIGAND,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Calprotectin'
  #### Cathelicidin
  indices.temp <- grep(pattern = '^CAMP',
                       x = feature.metadata$LIGAND,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Cathelicidin'
  #### Lactoferrin
  indices.temp <- grep(pattern = '^LTF',
                       x = feature.metadata$LIGAND,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Lactoferrin'
  #### Annexins
  indices.temp <- grep(pattern = '^ANXA',
                       x = feature.metadata$LIGAND,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Annexins'

  ## Growth Factors

  ### Canonical Growth Factors
  #### Epidermal Growth Factors (EGF Family)
  indices.temp <- grep(pattern = '^EGF|^EREG|^AREG|^TGFA|^BTC',
                       x = feature.metadata$LIGAND,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Epidermal Growth Factors (EGF Family)'

  #### Fibroblast Growth Factors (FGF Family)
  indices.temp <- grep(pattern = 'FGF',
                       x = feature.metadata$LIGAND,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Fibroblast Growth Factors (FGF Family)'
  #### Heparin-binding EGF-like Growth factor (HBEGF Family)
  indices.temp <- grep(pattern = '^HBEGF',
                       x = feature.metadata$LIGAND,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Heparin-binding EGF-like Growth factor (HBEGF Family)'
  #### Hepatocyte Growth Factor (HGF Family)
  indices.temp <- grep(pattern = '^HGF',
                       x = feature.metadata$LIGAND,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Hepatocyte Growth Factor (HGF Family)'
  #### Insulin-like Growth Factor (IGF Family)
  indices.temp <- grep(pattern = '^IGF',
                       x = feature.metadata$LIGAND,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Insulin-like Growth Factor (IGF Family)'
  #### Neurite Promoting Growth Factor (NEGF Family)
  indices.temp <- grep(pattern = '^PTN|^MDK',
                       x = feature.metadata$LIGAND,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Neurite Promoting Growth Factor (NEGF Family)'
  #### Nerve Growth Factor (NGF Family)
  indices.temp <- grep(pattern = '^NGF|^NTF', #NTF is part of NGF family https://en.wikipedia.org/wiki/Neurotrophin-3',
                       x = feature.metadata$LIGAND,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Nerve Growth Factor (NGF Family)'
  #### Neuregulins (NRG Family)
  indices.temp <- grep(pattern = '^NRG',
                       x = feature.metadata$LIGAND,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Neuregulins (NRG Family)'
  #### Platelet Derived Growth Factor (PDGF Family)
  indices.temp <- grep(pattern = '^PDGF',
                       x = feature.metadata$LIGAND,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Platelet Derived Growth Factor (PDGF Family)'
  #### Placental Growth Factor (PGF Family)
  indices.temp <- grep(pattern = '^PGF',
                       x = feature.metadata$LIGAND,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Placental Growth Factor (PGF Family)'
  #### Transforming Growth Factor Beta (TGFB Family)
  indices.temp <- grep(pattern = 'TGFB',
                       x = feature.metadata$LIGAND,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Transforming Growth Factor Beta (TGFB Family)'
  #### Vascular Endothelial Growth Factor (VEGF Family)
  indices.temp <- grep(pattern = 'VEGF|FIGF', # FIGF = VEGFD!
                       x = feature.metadata$LIGAND,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Vascular Endothelial Growth Factor (VEGF Family)'
  #### Angiopoietin (ANGPT Family)
  indices.temp <- grep(pattern = '^ANGPT',
                       x = feature.metadata$LIGAND,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Angiopoietin (ANGPT Family)'
  #### Brain-derived Neurotrophic Factor (BDNF Family)
  indices.temp <- grep(pattern = '^BDNF',
                       x = feature.metadata$LIGAND,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Brain-derived Neurotrophic Factor (BDNF Family)'
  #### Bone Morphogenetic Protein (BMP Family)
  indices.temp <- grep(pattern = 'BMP',
                       x = feature.metadata$LIGAND,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Bone Morphogenetic Protein (BMP Family)'
  #### Apelin
  indices.temp <- grep(pattern = '^APLN',
                       x = feature.metadata$LIGAND,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Apelin'
  #### Ciliary Neurotrophic factor
  #### Growth Differentiation Factors (GDF Family)
  indices.temp <- grep(pattern = 'GDF',
                       x = feature.metadata$LIGAND,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Growth Differentiation Factors (GDF Family)'
  #### Stem Cell Factor (KITLG)
  indices.temp <- grep(pattern = '^KITLG',
                       x = feature.metadata$LIGAND,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Stem Cell Factor (KITLG)'
  #### Glial cell line-derived neurotrophic factors (GDNF Family)
  indices.temp <- grep(pattern = '^GDNF',
                       x = feature.metadata$LIGAND,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Glial cell line-derived neurotrophic factors (GDNF Family)'
  ### Canonical Morphogens
  #### NOTCH Family
  indices.temp <- grep(pattern = '^JAG|^DLL|^DLK|^MFNG',
                       x = feature.metadata$LIGAND,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'NOTCH Family'
  #### WNT Family
  indices.temp <- grep(pattern = '^WNT|^RSPO|^SFRP|^NDP',
                       x = feature.metadata$LIGAND,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'WNT Family'
  #### Dickkopf Family (WNT-inhibitor)
  indices.temp <- grep(pattern = '^DKK',
                       x = feature.metadata$LIGAND,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Dickkopf Family (WNT-inhibitor)'

  #### Hedgehog Family
  indices.temp <- grep(pattern = '^DHH|^SHH|^IHH',
                       x = feature.metadata$LIGAND,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Hedgehog Family'
  #### Activins
  #### Inhibins
  indices.temp <- grep(pattern = '^INH',
                       x = feature.metadata$LIGAND,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Inhibins'

  ## Endocrine Molecules

  ### Calcium Homeostasis
  #### Calcitonin
  indices.temp <- grep(pattern = '^CALCA|^CALCB',
                       x = feature.metadata$LIGAND,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Calcitonin'
  #### Calmodulins
  # indices.temp <- grep(pattern = '^CALM',
  #                          x = feature.metadata$LIGAND,
  #                          ignore.case = FALSE)
  # feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Calmodulins'
  #### Calreticulin
  indices.temp <- grep(pattern = '^CALR',
                       x = feature.metadata$LIGAND,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Calreticulin'
  ### Metabolic Homeostasis
  #### Growth Hormone
  indices.temp <- grep(pattern = '^GH1$|^GH2$|^GHRH$',
                       x = feature.metadata$LIGAND,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Growth Hormone'
  #### Ghrelin
  indices.temp <- grep(pattern = '^GHRL$',
                       x = feature.metadata$LIGAND,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Ghrelin'
  #### Pancreatic Polypeptide
  #### Insulin-like Family
  indices.temp <- grep(pattern = '^INSL',
                       x = feature.metadata$LIGAND,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Insulin-like Family'
  #### Somatostatin
  indices.temp <- grep(pattern = '^SST',
                       x = feature.metadata$LIGAND,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Somatostatin'
  #### Agouti-signaling protein
  indices.temp <- grep(pattern = '^ASIP',
                       x = feature.metadata$LIGAND,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Agouti-signaling protein'
  ### Lipid Homeostasis
  #### Adiponectin
  indices.temp <- grep(pattern = '^ADIPOQ',
                       x = feature.metadata$LIGAND,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Adiponectin'
  #### Apolipoproteins
  indices.temp <- grep(pattern = '^APO',
                       x = feature.metadata$LIGAND,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Apolipoproteins'

  #### Lipoproteins
  indices.temp <- grep(pattern = '^LPL',
                       x = feature.metadata$LIGAND,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Lipoproteins'

  ### Cardiovascular & Osmotic Homeostasis
  #### Adrenomedullin
  indices.temp <- grep(pattern = '^ADM',
                       x = feature.metadata$LIGAND,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Adrenomedullin'
  #### Endothelin
  indices.temp <- grep(pattern = '^EDN',
                       x = feature.metadata$LIGAND,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Endothelin'

  #### Natriuretic Protein Family
  indices.temp <- grep(pattern = '^NPP',
                       x = feature.metadata$LIGAND,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Natriuretic Protein Family'
  #### Relaxin
  indices.temp <- grep(pattern = '^RLN',
                       x = feature.metadata$LIGAND,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Relaxin'
  #### Renin-Angiotensin System
  indices.temp <- grep(pattern = '^ACE|^AGT|^REN',
                       x = feature.metadata$LIGAND,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Renin-Angiotensin System'
  #### Secretin
  indices.temp <- grep(pattern = '^SCT',
                       x = feature.metadata$LIGAND,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Secretin'
  #### Vasopressin
  indices.temp <- grep(pattern = '^AVP$',
                       x = feature.metadata$LIGAND,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Vasopressin'
  #### VIP
  indices.temp <- grep(pattern = '^VIP|^ADCYAP',
                       x = feature.metadata$LIGAND,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'VIP'
  ### Reproductive Homeostasis
  #### Oxytocin
  indices.temp <- grep(pattern = '^OXT',
                       x = feature.metadata$LIGAND,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Oxytocin'
  #### Sex Hormone Binding Protein
  indices.temp <- grep(pattern = '^SHBG',
                       x = feature.metadata$LIGAND,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Sex Hormone Binding Protein'
  #### Human chorionic gonadotropin (hCG)
  indices.temp <- grep(pattern = '^CGA|^CGB',
                       x = feature.metadata$LIGAND,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Human chorionic gonadotropin (hCG)'
  ### Stress Homeostasis
  #### Urocortin
  indices.temp <- grep(pattern = '^UCN',
                       x = feature.metadata$LIGAND,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Urocortin'
  #### CRH
  indices.temp <- grep(pattern = '^CRH',
                       x = feature.metadata$LIGAND,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'CRH'


  ## Matrix Molecules

  ### Structural
  #### Collagens
  indices.temp <- grep(pattern = '^COL',
                       x = feature.metadata$LIGAND,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Collagens'
  #### Laminins
  indices.temp <- grep(pattern = '^LAM',
                       x = feature.metadata$LIGAND,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Laminins'
  #### Elastins
  #### Fibronectins
  indices.temp <- grep(pattern = '^FN1',
                       x = feature.metadata$LIGAND,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Fibronectins'
  ### Matricellular
  #### Agrin (AGRN)
  indices.temp <- grep(pattern = '^AGRN',
                       x = feature.metadata$LIGAND,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Agrin (AGRN)'
  #### Biglycan (BGN)
  indices.temp <- grep(pattern = '^BGN',
                       x = feature.metadata$LIGAND,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Biglycan (BGN)'
  #### Connective Tissue Growth Factor (CTGF)
  indices.temp <- grep(pattern = '^CTGF',
                       x = feature.metadata$LIGAND,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Connective Tissue Growth Factor (CTGF)'
  #### CTHRC1
  indices.temp <- grep(pattern = '^CTHRC1',
                       x = feature.metadata$LIGAND,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'CTHRC1'
  #### CYR61
  indices.temp <- grep(pattern = '^CYR61',
                       x = feature.metadata$LIGAND,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'CYR61'
  #### Decorins
  indices.temp <- grep(pattern = '^DCN',
                       x = feature.metadata$LIGAND,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Decorins'
  #### ECM1
  indices.temp <- grep(pattern = '^ECM1$',
                       x = feature.metadata$LIGAND,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'ECM1'
  #### Fibulins
  indices.temp <- grep(pattern = '^EFEMP|^FBLN',
                       x = feature.metadata$LIGAND,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Fibulins'
  #### Fibrillins
  indices.temp <- grep(pattern = '^FBN',
                       x = feature.metadata$LIGAND,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Fibrillins'
  #### Lecticans
  indices.temp <- grep(pattern = '^VCAN|^BCAN|^ACAN|^NCAN',
                       x = feature.metadata$LIGAND,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Lecticans'
  #### "MFAP"
  indices.temp <- grep(pattern = '^MFAP',
                       x = feature.metadata$LIGAND,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'MFAP'
  #### Myocilin
  indices.temp <- grep(pattern = '^MYOC',
                       x = feature.metadata$LIGAND,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Myocilin'
  #### Nephronectins
  indices.temp <- grep(pattern = '^NPNT',
                       x = feature.metadata$LIGAND,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Nephronectins'
  #### Nidogens
  indices.temp <- grep(pattern = '^NID1|^NID2',
                       x = feature.metadata$LIGAND,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Nidogens'
  #### Osteopontins
  indices.temp <- grep(pattern = '^SPP1',
                       x = feature.metadata$LIGAND,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Osteopontins'
  #### Prosaposins
  indices.temp <- grep(pattern = '^PSAP',
                       x = feature.metadata$LIGAND,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Prosaposins'
  #### Reelin
  indices.temp <- grep(pattern = '^RELN',
                       x = feature.metadata$LIGAND,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Reelin'
  #### Tenascins
  indices.temp <- grep(pattern = '^TNC',
                       x = feature.metadata$LIGAND,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Tenascins'
  #### Thrombospondins
  indices.temp <- grep(pattern = '^THBS',
                       x = feature.metadata$LIGAND,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Thrombospondins'
  #### Vitronectins
  indices.temp <- grep(pattern = '^VTN',
                       x = feature.metadata$LIGAND,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Vitronectins'
  #### ZP3
  indices.temp <- grep(pattern = '^ZP3$',
                       x = feature.metadata$LIGAND,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'ZP3'
  ### Proteases
  #### Urokinase (PLAU)
  indices.temp <- grep(pattern = '^PLAU',
                       x = feature.metadata$LIGAND,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Urokinase (PLAU)'
  #### Metalloprotease Family
  indices.temp <- grep(pattern = '^MMP',
                       x = feature.metadata$LIGAND,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Metalloprotease Family'
  #### ADAM Family
  indices.temp <- grep(pattern = '^ADAM',
                       x = feature.metadata$LIGAND,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'ADAM Family'
  #### Tissue factor pathway inhibitor (TFPI)
  indices.temp <- grep(pattern = 'TFPI',
                       x = feature.metadata$LIGAND,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Tissue factor pathway inhibitor (TFPI)'
  #### Transglutaminase Family
  indices.temp <- grep(pattern = '^TGM',
                       x = feature.metadata$LIGAND,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Transglutaminase Family'
  #### Presenilin Family
  indices.temp <- grep(pattern = '^PSEN',
                       x = feature.metadata$LIGAND,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Presenilin Family'
  #### Granzymes
  indices.temp <- grep(pattern = '^GZMB',
                       x = feature.metadata$LIGAND,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Granzymes'
  #### Plasmin
  indices.temp <- grep(pattern = '^PLG',
                       x = feature.metadata$LIGAND,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Plasmin'
  ### Protease Inhibition
  #### Serpins
  indices.temp <- grep(pattern = '^SERP',
                       x = feature.metadata$LIGAND,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Serpins'
  #### Globulins
  indices.temp <- grep(pattern = '^A2M',
                       x = feature.metadata$LIGAND,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Globulins'
  #### Tissue inhibitors of metalloproteinases (TIMPs)
  indices.temp <- grep(pattern = '^TIMP',
                       x = feature.metadata$LIGAND,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Tissue inhibitors of metalloproteinases (TIMPs)'
  ### Circulatory Factors
  #### Kinin-Kallikrein System
  # indices.temp <- grep(pattern = '^KNG',
  #                      x = feature.metadata$LIGAND,
  #                      ignore.case = FALSE)
  # feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Kinin-Kallikrein System'
  #### Complement Molecules
  indices.temp <- grep(pattern = '^C1Q|^C3|^C4|^CD55$|^C5$|^CRP$',
                       x = feature.metadata$LIGAND,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Complement Molecules'
  #### Transferrins
  indices.temp <- grep(pattern = '^TF$',
                       x = feature.metadata$LIGAND,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Transferrins'
  #### Coagulation Factors
  indices.temp <- grep(pattern = '^VWF|^F2$|^F10$|^F12$|^F8$|^F7$|^F9$',
                       x = feature.metadata$LIGAND,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Coagulation Factors'
  ### Miscellaneous/Uncategorized
  #### LYPD3
  #### Secretory Proteins
  indices.temp <- grep(pattern = '^MUC|^SCGB|^SFTP',
                       x = feature.metadata$LIGAND,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Secretory Proteins'

  ## Spatial Guidance Molecules

  ### Cell Adhesion Molecules

  #### Amyloid Beta Precursor
  indices.temp <- grep(pattern = '^APP',
                       x = feature.metadata$LIGAND,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Amyloid Beta Precursor'
  #### Cadherins
  indices.temp <- grep(pattern = '^CDH',
                       x = feature.metadata$LIGAND,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Cadherins'

  #### CAM Family
  indices.temp <- grep(pattern = '^VCAM|^NCAM|^ICAM|^L1CAM',
                       x = feature.metadata$LIGAND,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'CAM Family'
  #### Contactins
  indices.temp <- grep(pattern = '^CNTN',
                       x = feature.metadata$LIGAND,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Contactins'
  #### Neuroligins
  indices.temp <- grep(pattern = '^NLGN',
                       x = feature.metadata$LIGAND,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Neuroligins'
  #### Neurexophilins
  indices.temp <- grep(pattern = '^NXP',
                       x = feature.metadata$LIGAND,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Neurexophilins'
  #### Lectins
  indices.temp <- grep(pattern = '^CLEC|^SELPL',
                       x = feature.metadata$LIGAND,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Lectins'
  ### Ephrins
  indices.temp <- grep(pattern = '^EFN',
                       x = feature.metadata$LIGAND,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Ephrins'
  ### Netrins
  indices.temp <- grep(pattern = '^NTN',
                       x = feature.metadata$LIGAND,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Netrins'
  ### Repulsive Guidance Molecule
  indices.temp <- grep(pattern = '^RGM',
                       x = feature.metadata$LIGAND,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Repulsive Guidance Molecule'
  ### Reticulons
  indices.temp <- grep(pattern = '^RTN',
                       x = feature.metadata$LIGAND,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Reticulons'
  ### Semaphorins
  indices.temp <- grep(pattern = '^SEMA',
                       x = feature.metadata$LIGAND,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Semaphorins'
  ### SLITs
  indices.temp <- grep(pattern = 'SLIT',
                       x = feature.metadata$LIGAND,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'SLITs'
  ## Neuroproteins


  ### Neuropeptides
  #### Agouti-Related Peptide
  indices.temp <- grep(pattern = '^AGRP',
                       x = feature.metadata$LIGAND,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Agouti-Related Peptide'
  #### CNTF
  indices.temp <- grep(pattern = '^CNTF',
                       x = feature.metadata$LIGAND,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'CNTF'
  #### Cortistatin
  indices.temp <- grep(pattern = '^CORT',
                       x = feature.metadata$LIGAND,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Cortistatin'
  #### Tachykinin
  indices.temp <- grep(pattern = '^TAC',
                       x = feature.metadata$LIGAND,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Tachykinin'
  #### Neuropeptide, Unsorted
  indices.temp <- grep(pattern = '^NPW$|^NPY$|^NPB$|^POMC$|^PMCH$|^HCRT$|^NTS$',
                       x = feature.metadata$LIGAND,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Neuropeptide, Unsorted'
  #### Neuromedins
  indices.temp <- grep(pattern = '^NMU$|^NMB$|^NMS$',
                       x = feature.metadata$LIGAND,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Neuromedin Family'

  ## Non-Ligands (Cytosolic Molecules)

  ### Intracellular (Structural)

  #### Afadin
  indices.temp <- grep(pattern = '^MLLT4$',
                       x = feature.metadata$LIGAND,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Afadin'
  #### Cingulin
  indices.temp <- grep(pattern = '^CGN',
                       x = feature.metadata$LIGAND,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Cingulin'
  #### Vimentin
  indices.temp <- grep(pattern = '^VIM$',
                       x = feature.metadata$LIGAND,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Vimentin'

  ### Intracellular (Enzyme)
  indices.temp <- grep(pattern = '^DUSP|^AANAT|^HDC$',
                       x = feature.metadata$LIGAND,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Intracellular Enzymes'

  ### Intracellular (Uncategorized)

  #### G Protein Complex Molecules
  indices.temp <- grep(pattern = '^GNA',
                       x = feature.metadata$LIGAND,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'G Protein Complex Molecules'
  #### LDL-receptor-related protein-associated protein (LRPAP)
  indices.temp <- grep(pattern = '^LRPAP',
                       x = feature.metadata$LIGAND,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'LDL-receptor-related protein-associated protein (LRPAP)'
  #### PIGF
  indices.temp <- grep(pattern = '^PIGF',
                       x = feature.metadata$LIGAND,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Phosphatidylinositol-glycan biosynthesis class F protein'
  #### Tryptophan Hydroxylase
  # indices.temp <- grep(pattern = '^TPH1',
  #                      x = feature.metadata$LIGAND,
  #                      ignore.case = FALSE)
  # feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Tryptophan hydroxylase'
  #### Heat-shock Proteins
  indices.temp <- grep(pattern = '^HSP',
                       x = feature.metadata$LIGAND,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Heat-shock Proteins'

  #### Intracellular (uncategorized)
  indices.temp <- grep(pattern = '^FARP|^HRAS$|^KRTAP4|^ARF1$|^SHANK1$|^C1orf|^C5orf|^C6orf',
                       x = feature.metadata$LIGAND,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$LIGAND.FAMILY <- 'Intracellular (uncategorized)'



################################
######## UNSORTED BELOW ########
################################


  return(feature.metadata)
}
