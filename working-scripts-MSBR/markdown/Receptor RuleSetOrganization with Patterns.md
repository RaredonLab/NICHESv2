Receptor Categorization Reference Markdown
Version: 1.0.0
Author: MSBR
Date Last Modified: 2024-12-31

# FANTOM5 Receptors

%%%%%%%%%%

## Cytokine Receptors

### Chemokine Family Receptors
#### CCL Family Receptor
##### '^CCR'
#### CXC Family Receptor
##### '^CXC'
#### C Chemokine (XCL) Family Receptors
##### '^XCR'
#### Colony Stimulating Factor (CSF) Family Receptors
##### '^CSF'
#### Prokineticin Family Receptors
##### '^PROKR'

### Interferon Family Receptors
#### IFN Family Receptors
##### '^INF'

### Interleukin Family Receptors
#### CLCF Family
##### '^CRLF'
#### IL Family Receptor Subunis
##### '^IL'
#### Oncostatin M
##### '^OSMR'
#### Leukemia Inhibitory Factor
##### '^LIFR'

### Tumor Necrosis Factor Family Receptors
#### TNF Family
##### '^LTBR$'
#### TNFSF Family Receptor Subunit
##### '^TNFR|^CD27$'

### Lymphoid Interaction Receptors
#### HLA Complex Molecules
##### '^LILRB|^KIR|^CD1A$|^CD1B$|^CD74$'
#### CD40LG
##### '^CD40'
#### BTLA
##### '^VTCN|^CD79A'
#### B-Cell Antigen Receptor Complex (BCR)
##### '^CD19$'

### Pattern Recognition Receptors
#### Toll-Like Receptors
##### 'TLR'

### Miscellaneous (Receptors)
#### "Apoptosis inhibitor expressed by macrophages"
##### '^CD5$'
#### Calprotectin Family Receptor
##### NA
#### Cathelicidin Family Receptor
##### '^FPR'
#### Lactotransferrin Family Receptor
##### NA
#### Annexin Family Receptor
##### '^DYSF'

%%%%%%%%%%

## Growth Factor Receptors

### Canonical Growth Factor Receptors
#### Epidermal Growth Factors (EGF) Family Receptor
##### '^ERBB|EGFR'
#### Fibroblast Growth Factors (FGF) Family Receptor
##### '^FGFR'
#### Heparin-binding EGF-like Growth factor (HBEGF Family)
##### '^CD9$'
#### Hepatocyte Growth Factor (HGF Family)
##### '^MET$'
#### Insulin-like Growth Factor (IGF) Family Receptor
##### '^IGF'
#### Neurite Promoting Growth Factor (NEGF Family)
##### '^PTPR'
#### Nerve Growth Factor (NGF) Family Receptor
##### '^NGF|^NTRK'
#### Neuregulins (NRG Family)
##### NA
#### Platelet Derived Growth Factor (PDGF) Family Receptor
##### '^PDGF'
#### Placental Growth Factor (PGF Family)
##### NA
#### Transforming Growth Factor Beta (TGFB) Family Receptor
##### '^TGFBR'
#### Vascular Endothelial Growth Factor (VEGF Family)
##### '^KDR$|^FLT'
#### Angiopoietin (ANGPT) Family Receptor
##### '^TIE|^TEK'
#### Brain-derived Neurotrophic Factor (BDNF Family)
##### NA
#### Bone Morphogenetic Protein (BMP) Family Receptor
##### '^BMPR'
#### Apelin
##### '^APLNR$'
#### Ciliary Neurotrophic factor
##### '^CNTFR'
#### Growth Differentiation Factors (GDF Family)
##### NA
#### Stem Cell Factor (KITLG) Receptor
##### '^KIT$'
#### Glial cell line-derived neurotrophic factors (GDNF) Family Receptor
##### '^RET|^GFR'

### Canonical Morphogen Receptors
#### NOTCH Family Receptors
##### 'NOTCH'
#### WNT Family Receptors
##### '^LGR|FZD|RYK'
#### Dickkopf Family (WNT-inhibitor) Family Receptors
##### '^KREMEN'
#### Hedgehog Family Receptors
##### '^PTCH|^HH|^CDON'
#### Inhibin Family Receptors
##### '^ACVR|^BAMBI'

%%%%%%%%%%

## Endocrine Molecule Receptors
  ## Parathyroid Hormone Receptor
  indices.temp <- grep(pattern = '^PTH',
                       x = feature.metadata$RECEPTOR,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$RECEPTOR.FAMILY <- 'Parathyroid Hormone Receptor'
  ## Prolactin Family Receptor
  indices.temp <- grep(pattern = '^PRL',
                       x = feature.metadata$RECEPTOR,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$RECEPTOR.FAMILY <- 'Prolactin Family Receptor'

  ## Glucagon-like peptide receptor
  indices.temp <- grep(pattern = '^GLP',
                       x = feature.metadata$RECEPTOR,
                       ignore.case = FALSE)
  feature.metadata[indices.temp,]$RECEPTOR.FAMILY <- 'Glucagon-like peptide receptor'
### Calcium Homeostasis
#### Calcitonin Family Receptor
##### '^CALCR|^RAMP'
#### Calreticulin
##### NA

### Metabolic Homeostasis
#### Growth Hormone Family Receptor
##### '^GHR'
#### Ghrelin
##### NA
#### Pancreatic Polypeptide
##### NA
#### "Insulin-like"
##### NA
#### Somatostatin
##### ^'SSTR'
#### Agouti-signaling protein
##### '^ATRN$'

### Lipid Homeostasis
#### Adiponectin Family Receptor
##### '^ADIPOR'

#### Lipoprotein Family Receptor
##### '^LDL|^LRP'

### Cardiovascular & Osmotic Homeostasis
#### Adrenomedullin Family Receptor
##### '^RAMP'
#### Endothelin Family Receptor
##### '^EDNR'
#### Natriuretic Family Receptor
##### '^NPR'
#### Relaxin Family Receptor
##### '^RXFP'
#### Renin-Angiotensin System Receptor
##### '^AGTR'
#### Secretin Family Receptor
##### '^SCTR$'
#### Vasopressin Family Receptor
##### '^AVPR'
#### VIP Family Receptor
##### '^VIPR'

### Reproductive Homeostasis
#### Oxytocin Receptor
##### '^OXTR'
#### Sex Hormone Binding Protein Receptors
##### NA
#### Gonadotropin Family Receptors
##### '^FSH|^LHCGR'

### Stress Homeostasis
#### Urocortin
##### '^CRHR'
#### CRH
##### NA

### Opiod System
#### Opiod System Receptors
##### '^OPR'

%%%%%%%%%%

## Matrix Molecule Receptors


### Heparan sulfate proteoglycans (HSPGs)
#### Glypicans
##### '^GPC'
#### Syndecans
##### '^SYN'
### Integrins
#### Integrins
##### '^ITG'

### Matricellular Family Receptors
#### Agrin (AGRN)
##### '^MUSK$'
#### Biglycan (BGN)
##### NA
#### Connective Tissue Growth Factor (CTGF)
##### NA
#### CTHRC1
##### NA
#### CYR61
##### NA
#### Decorins
##### NA
#### ECM1
##### NA
#### Fibulins
##### NA
#### Fibrillins
##### NA
#### Lecticans
##### NA
#### Microfibrillar-associated proteins
##### NA
#### Myocilin
##### NA
#### Nephronectins
##### NA
#### Nidogens 
##### NA
#### Osteopontins
##### NA
#### Prosaposins
##### '^GPR37|^CELSR1$'
#### Reelin
##### NA
#### Tenascins
##### 
#### Thrombospondins
##### '^CD36$|^CD47$'
#### Vitronectins
##### NA
#### ZP3
##### NA

### Protease Family Receptors
#### Urokinase (PLAU)
##### '^PLAUR$'
#### Metalloprotease Family
##### NA
#### ADAM Family
##### NA
#### Tissue factor pathway inhibitor (TFPI) 
##### NA
#### Transglutaminase Family
##### '^GPR56$'
#### Presenilin Family
##### '^NCSTN'
#### Granzymes
##### NA
#### Plasmin Family Receptor
##### '^PLGRKT$'

### Protease Inhibition
#### Serpins
##### NA
#### Globulins
##### NA
#### Tissue inhibitors of metalloproteinases (TIMPs)
##### '^CD63$'

### Circulatory Factors
#### Kinin-Kallikrein System
##### '^BDKRB'
#### Complement Molecules
##### '^CR1$|^C3AR|^C5AR|^CD97'
#### Transferrins
##### '^TFR'
#### Coagulation Factors
##### '^GP1BA$|^F3$|^F2|^ASGR2$'
### Miscellaneous/Uncategorized
#### LYPD3
##### NA
#### Secretory Proteins
##### NA

  
%%%%%%%%%%

## Spatial-Guidance Molecule Receptors

### Cell-Adhesion Molecules
#### Amyloid Beta Precursor
##### NA
#### Cadherin Family Receptors
##### '^CDH'
#### CAM Family Receptor Subunit
##### 'CAM'
#### Contactins
##### '^CNTN'
#### Neurexins
##### '^NRXN'
#### Lectin Family Receptors
##### '^SEL|^CD93'

### ECM-like Molecules
#### Netrin Family Receptors
##### '^UNC5'
#### SLIT Family Receptors
##### '^ROBO'

### Growth & Migration Inhibitors
#### Repulsive guidance molecules (RGM Family)
##### '^NEO'
#### Reticulon (Nogo) Family Receptors
##### '^RTN'

### Context Signaling Molecules
#### Ephrin Family Receptors
##### ^EPH'
#### Semaphorin Family Receptors
##### ^PLXN'

%%%%%%%%%%

## Neuropeptide Receptors

### Central Melanocortin System
#### Agouti-Related Peptide
##### NA
#### NPY
##### '^NPY'
#### POMC
##### '^MC1|^MC2|^MC3|^MC4|^MC5|'
### Neuropeptide, Unsorted
#### Neuropeptide Receptor
##### '^NPFFR|^GALR''
#### Pro-melanin stimulating hormone (PMCH)
##### '^MCHR'
#### Orexin (HCRT)
##### ^'HCRT'
#### Neurotensin (NTS)
##### '^NTSR'
#### Neuromedins
##### '^NM'
#### NPW
##### NA
#### NPB 
##### '^NPBW'
#### CNTF
##### '^CNTFR$'
#### Cortistatin
##### NA
#### Tachykinin
##### '^TACR'

## Pleiotrophic Receptors
### Transmembrane Glycoproteins
#### Neuropilins
##### '^NRP'

%%%%

## Non-Receptors (Cytosolic Molecules)

### Intracellular (Structural)
#### Afadin
##### NA
#### Vimentin
##### NA
#### Cingulin
##### NA
#### Scaffold Proteins
##### NA
#### Keratin-associated
##### NA 
### Intracellular (Signaling)
#### G Protein Complex Molecules
##### NA
#### LDL-receptor-related protein-associated protein (LRPAP)
##### NA 
#### Phosphatidylinositol-glycan biosynthesis class F protein (PIGF)
##### NA 
#### Heat-shock Proteins
##### NA
#### GTPase Signaling
##### NA 
### Intracellular (Enzymes)
#### Tryptophan Hydroxylase
##### NA 
#### Enzyme, Unsorted
##### NA 
### Intracellular (Uncategorized)
#### Intracellular (Uncategorized)
##### NA 

%%%%%%%%%%%%%%%%%%
