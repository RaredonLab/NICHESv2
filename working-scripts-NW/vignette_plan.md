# NICHESv2 Vignette Plan & Formatting Reference
# Pass 5 — carried forward from 2026-04-17 session

================================================================================
FORMATTING STANDARD (apply to all 3 vignettes identically)
================================================================================

YAML front matter template:
```yaml
---
title: "<Vignette Title>"
date: "`r Sys.Date()`"
output:
  prettydoc::html_pretty:
    theme: cayman
    highlight: github
vignette: >
  %\VignetteIndexEntry{<Vignette Title>}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---
```

knitr setup chunk (all vignettes):
```r
knitr::opts_chunk$set(
  collapse  = TRUE,
  comment   = "#>"
)
```

DESCRIPTION Suggests must include: knitr, rmarkdown, prettydoc
DESCRIPTION must include: VignetteBuilder: knitr

File location: vignettes/<name>.Rmd
Rendering: eval=TRUE for vignette 1 (bundled data only);
           eval=FALSE + pre-rendered HTML committed for vignettes 2 and 3


================================================================================
VIGNETTE STATUS
================================================================================

Vignette 1 — getting-started.Rmd         STATUS: COMPLETE
  Data:    Bundled NICHESv2_inputs + NICHESv2_example
  eval:    TRUE (no external deps, fast)
  Purpose: Orient users to NICHESObject slot structure
  Notes:   Column name confirmed: Sending.celltype


Vignette 2 — scrna-sampled.Rmd           STATUS: NOT STARTED
  Data:    ifnb dataset from SeuratData (~14k cells, 14 cell types,
           CTRL + STIM conditions)
  eval:    FALSE + pre-rendered HTML committed to repo
  Purpose: Sampled mode workflow for annotated scRNA-seq without spatial coords

  Content outline:
  1. Biological motivation: no spatial coords → no principled neighborhood
     definition; all-pairs is combinatorially intractable; sampled mode
     provides stratified-by-cell-type reproducible pairs
  2. extract_NICHESInputs_Seurat() — ingest from Seurat object
  3. Sample splitting — prominently commented code block showing
     sample.col = "stim"; explain consequences of omitting it
     (ifnb has CTRL + STIM: natural, non-contrived demonstration)
  4. create_NICHESObject(mode = "sampled", sample.col = "stim")
  5. Inspect result: print(), $edge.data, $edge.meta
  6. compute_PseudoBulk_Connectivity() — show 2-3 rows x 5-6 cols of
     output data.table, then ggplot2 heatmap (SendingType x ReceivingType,
     colored by aggregated score)
  7. Explicit note: no in-situ plotting without spatial coords; see
     spatial vignette

  Key packages needed in Suggests: SeuratData, SeuratObject


Vignette 3 — spatial.Rmd                 STATUS: NOT STARTED
  Data:    stxBrain anterior section from SeuratData (~2800 spots,
           mouse brain, strong spatial tissue organization)
  eval:    FALSE + pre-rendered HTML committed to repo
  Purpose: Spatial mode workflow with radius selection and in-situ plot

  Content outline:
  1. Biological motivation: physical proximity defines signaling
     neighborhoods; two parameter paths: kNN vs radius
  2. Explore_Radius_Neighborhood() — show diagnostic output, explain
     how to read it to pick radius (users think in avg neighbors per cell)
  3. create_NICHESObject(mode = "spatial", rad = <chosen value>)
  4. Sample splitting — code comment only (stxBrain is single section):
     "If your experiment has multiple tissue sections or capture areas,
      set sample.col here to avoid cross-section edges."
  5. Inspect result: print(), $edge.data with spatial columns, $edge.meta
  6. CellToCell_InSitu_Single_LRM() — one representative LRM pair
     chosen after exploratory run (see prerequisite below)
  7. No mention of batch plotting function

  PREREQUISITE before writing: run stxBrain through create_NICHESObject(),
  inspect $LRM.list and $edge.data, pick an LRM pair whose in-situ plot
  shows spatially structured connectivity. Lock in that pair, then write
  the vignette around it.

  Key packages needed in Suggests: SeuratData, SeuratObject


================================================================================
PENDING DECISIONS (resolve before or during Pass 5 writing)
================================================================================

- LRM pair for spatial in-situ plot: requires exploratory run on stxBrain
- Cell type column name in ifnb: confirm before writing vignette 2
  (convention: use whatever column name is in the Seurat object metadata)
- Column names to show in pseudobulk data.table head(): confirm against
  actual compute_PseudoBulk_Connectivity() output
