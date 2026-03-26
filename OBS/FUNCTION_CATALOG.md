# NICHES Package - Complete Function Catalog

## Table of Contents
1. [Main Entry Point](#main-entry-point)
2. [S3 Dispatch Methods](#s3-dispatch-methods)
3. [Core Organization Functions (Run* family)](#core-organization-functions)
4. [Ground Truth Loading (Load* family)](#ground-truth-loading-functions)
5. [Helper/Utility Functions](#helperutility-functions)
6. [Visualization Functions](#visualization-functions)
7. [Analysis Functions](#analysis-functions)
8. [Function Call Hierarchy](#function-call-hierarchy)

---

## MAIN ENTRY POINT

### RunNICHES (Generic)
**File:** `/Users/msbr/GitHub/NICHESv2-msbr/R/RunNICHES.R`

**Type:** Generic S3 method dispatcher

**Parameters:**
- `object`: Seurat object or normalized data matrix
- `...`: Arguments passed to specific methods

**Returns:** List of Seurat objects or raw matrices (depending on output_format)

**External Dependencies:** Seurat

**Description:** Main user-facing function that dispatches to appropriate method (Seurat, matrix, or Matrix) based on input type.

---

## S3 DISPATCH METHODS

### RunNICHES.Seurat
**File:** `/Users/msbr/GitHub/NICHESv2-msbr/R/RunNICHES.R` (lines 280-372)

**Type:** S3 method for Seurat objects

**Parameters:**
- `object`: Seurat object (required)
- `assay`: string, default "RNA" - which assay to use
- `LR.database`: string, default "fantom5" - "fantom5", "omnipath", or "custom"
- `species`: string (required) - "human", "mouse", "rat", "pig" (depends on database)
- `min.cells.per.ident`: integer, default NULL - minimum cells per cell type
- `min.cells.per.gene`: integer, default NULL - minimum cells expressing gene
- `meta.data.to.map`: character vector, default all columns - metadata to propagate
- `position.x`: string, default NULL - spatial x-coordinate metadata column
- `position.y`: string, default NULL - spatial y-coordinate metadata column
- `cell_types`: string, default NULL - metadata column containing cell type labels
- `custom_LR_database`: dataframe, default NULL - custom ligand-receptor database
- `k`: integer, default 4 - number of neighbors for k-NN
- `rad.set`: numeric, default NULL - radius threshold for spatial neighbors
- `blend`: string, default "mean" - "sum", "mean", or "mean.adj" for aggregation
- `CellToCell`: logical, default TRUE
- `CellToSystem`: logical, default FALSE
- `SystemToCell`: logical, default FALSE
- `CellToCellSpatial`: logical, default FALSE
- `CellToNeighborhood`: logical, default FALSE
- `NeighborhoodToCell`: logical, default FALSE
- `output_format`: string, default "seurat" - "seurat" or "raw"

**Returns:** List of Seurat objects or raw matrices (named by organization type)

**Internal Calls:**
- `RunNICHES.default()`

**External Dependencies:** Seurat

**Description:** Validates inputs, checks metadata columns, sets cell type identities, then calls default method.

---

### RunNICHES.matrix
**File:** `/Users/msbr/GitHub/NICHESv2-msbr/R/RunNICHES.R` (lines 412-508)

**Type:** S3 method for matrix objects

**Parameters:**
- `object`: matrix - expression matrix (rows=genes, columns=cells) (required)
- `LR.database`: string, default "fantom5"
- `species`: string (required)
- `min.cells.per.ident`: integer, default NULL
- `min.cells.per.gene`: integer, default NULL
- `meta.data.df`: dataframe, default NULL - cell metadata
- `position.x`: string, default NULL
- `position.y`: string, default NULL
- `cell_types`: string (required) - column name in meta.data.df for cell types
- `custom_LR_database`: dataframe, default NULL
- `k`: integer, default 4
- `rad.set`: numeric, default NULL
- `blend`: string, default "mean"
- `CellToCell`: logical, default TRUE
- `CellToSystem`: logical, default FALSE
- `SystemToCell`: logical, default FALSE
- `CellToCellSpatial`: logical, default FALSE
- `CellToNeighborhood`: logical, default FALSE
- `NeighborhoodToCell`: logical, default FALSE
- `output_format`: string, default "seurat"

**Returns:** List of Seurat objects or raw matrices

**Internal Calls:**
- `RunNICHES.default()`

**External Dependencies:** Seurat

**Description:** Converts matrix input to Seurat object, validates cell_types, then calls default method.

---

### RunNICHES.Matrix
**File:** `/Users/msbr/GitHub/NICHESv2-msbr/R/RunNICHES.R` (line 518)

**Type:** S3 method alias

**Parameters:** Same as RunNICHES.matrix

**Returns:** Same as RunNICHES.matrix

**Description:** Wrapper that points to RunNICHES.matrix for sparse Matrix objects.

---

### RunNICHES.default
**File:** `/Users/msbr/GitHub/NICHESv2-msbr/R/RunNICHES.R` (lines 62-247)

**Type:** Core implementation method

**Parameters:** (Same comprehensive list as RunNICHES.Seurat)

**Returns:** List of outputs named by enabled organization types

**Internal Calls:**
- `PrepNodeObject()` - preprocesses input
- `LoadGroundTruth()` - loads L-R database
- `ComputeEdgelist()` - (if spatial analysis requested)
- `RunCellToCell()` - (if enabled)
- `RunCellToSystem()` - (if enabled)
- `RunSystemToCell()` - (if enabled)
- `RunCellToCellSpatial()` - (if enabled)
- `RunCellToNeighborhood()` - (if enabled)
- `RunNeighborhoodToCell()` - (if enabled)

**External Dependencies:** Seurat

**Description:** Main preprocessing pipeline. Orchestrates validation, preprocessing, ground truth loading, edgelist computation (if spatial), and calls all enabled organization functions.

---

## CORE ORGANIZATION FUNCTIONS

### RunCellToCell
**File:** `/Users/msbr/GitHub/NICHESv2-msbr/R/RunCellToCell.R`

**Type:** Exported user-facing function

**Parameters:**
- `node.object`: Seurat object (preprocessed, required)
- `ground.truth`: list with source.subunits and target.subunits matrices (required)
- `assay`: string - assay to use (required)
- `meta.data.to.map`: character vector - metadata columns to propagate (required)
- `output_format`: string - "seurat" or "raw" (required)

**Returns:** 
- If `output_format == "seurat"`: Seurat object with assay "CellToCell"
- If `output_format == "raw"`: List with `$CellToCellMatrix` and `$metadata`

**Internal Calls:**
- `return_celltypes()` - extracts unique cell types
- `GetSeuratAssay()` - retrieves expression data (Seurat version-aware)
- `Seurat::subset()` - subset by cell type
- `Seurat::Idents()`
- `Seurat::CreateSeuratObject()`
- `Seurat::NormalizeData()` - (if Seurat v5)
- `Seurat::AddMetaData()`
- `abind::abind()` - for metadata array binding

**External Dependencies:** 
- Seurat
- dplyr
- abind

**Description:** 
Creates cell-to-cell (CellToCell) signaling matrix by:
1. Building ligand expression matrices per cell type (multiplying subunit expressions)
2. Building receptor expression matrices per cell type
3. For each sender type, downsampling to min(sender_size, receiver_size) edges per pair
4. Computing interaction strength as ligand * receptor expression
5. Creating Seurat object with mechanismsxpairs dimensions
6. Adding metadata (SendingType, ReceivingType, VectorType, and mapped metadata)

---

### RunCellToSystem
**File:** `/Users/msbr/GitHub/NICHESv2-msbr/R/RunCellToSystem.R`

**Type:** Exported user-facing function

**Parameters:**
- `node.object`: Seurat object (preprocessed, required)
- `ground.truth`: list with source.subunits and target.subunits matrices (required)
- `assay`: string - assay to use (required)
- `blend`: string, default "mean" - "sum", "mean", or "mean.adj"
- `meta.data.to.map`: character vector (required)
- `output_format`: string - "seurat" or "raw" (required)

**Returns:**
- If `output_format == "seurat"`: Seurat object with assay "CellToSystem"
- If `output_format == "raw"`: List with `$CellToSystemMatrix` and `$metadata`

**Internal Calls:**
- `GetSeuratAssay()` - retrieves expression data
- `Matrix::rowSums()`, `Matrix::rowMeans()` - aggregate across system
- `Seurat::CreateSeuratObject()`
- `Seurat::NormalizeData()` - (if Seurat v5)
- `Seurat::Idents()`
- `Seurat::AddMetaData()`

**External Dependencies:** 
- Seurat
- Matrix

**Description:**
Creates cell-to-system (outgoing signaling) matrix by:
1. Building ligand expression across ALL cells
2. Aggregating ligand expression (sum/mean/mean.adj across all receiving cells)
3. Building receptor expression across ALL cells
4. Computing interaction strength as aggregated_ligand * receptor_per_cell
5. Creating Seurat object with mechanismsxsending_cells dimensions
6. Columns represent each sending cell and its total output to system

---

### RunSystemToCell
**File:** `/Users/msbr/GitHub/NICHESv2-msbr/R/RunSystemToCell.R`

**Type:** Exported user-facing function

**Parameters:**
- `node.object`: Seurat object (preprocessed, required)
- `ground.truth`: list with source.subunits and target.subunits matrices (required)
- `assay`: string - assay to use (required)
- `blend`: string, default "mean" - "sum", "mean", or "mean.adj"
- `meta.data.to.map`: character vector (required)
- `output_format`: string - "seurat" or "raw" (required)

**Returns:**
- If `output_format == "seurat"`: Seurat object with assay "SystemToCell"
- If `output_format == "raw"`: List with `$SystemToCellMatrix` and `$metadata`

**Internal Calls:**
- `GetSeuratAssay()` - retrieves expression data
- `Matrix::rowSums()`, `Matrix::rowMeans()` - aggregate across system
- `Seurat::CreateSeuratObject()`
- `Seurat::NormalizeData()` - (if Seurat v5)
- `Seurat::Idents()`
- `Seurat::AddMetaData()`

**External Dependencies:**
- Seurat
- Matrix

**Description:**
Creates system-to-cell (cellular niche/microenvironment) matrix by:
1. Building ligand expression across ALL cells
2. Building receptor expression per cell
3. Aggregating receptor expression (sum/mean/mean.adj across all receiving cells)
4. Computing interaction strength as ligand_per_cell * aggregated_receptor
5. Creating Seurat object with mechanismsxreceiving_cells dimensions
6. Columns represent each receiving cell and its total input from system

---

### RunCellToCellSpatial
**File:** `/Users/msbr/GitHub/NICHESv2-msbr/R/RunCellToCellSpatial.R`

**Type:** Exported user-facing function

**Parameters:**
- `node.object`: Seurat object (preprocessed, required)
- `ground.truth`: list with source.subunits and target.subunits matrices (required)
- `assay`: string - assay to use (required)
- `meta.data.to.map`: character vector (required)
- `edgelist`: dataframe with "from" and "to" columns (required)
- `output_format`: string - "seurat" or "raw" (required)

**Returns:**
- If `output_format == "seurat"`: Seurat object with assay "CellToCellSpatial"
- If `output_format == "raw"`: List with `$CellToCellSpatialMatrix` and `$metadata`

**Internal Calls:**
- `GetSeuratAssay()` - retrieves expression from spatially-connected cells only
- `Seurat::CreateSeuratObject()`
- `Seurat::NormalizeData()` - (if Seurat v5)
- `Seurat::Idents()`
- `Seurat::AddMetaData()`
- `abind::abind()`

**External Dependencies:**
- Seurat
- abind

**Description:**
Creates spatially-restricted cell-to-cell matrix by:
1. Building ligand expression only for "from" (sending) cells in edgelist
2. Building receptor expression only for "to" (receiving) cells in edgelist
3. Computing interaction strength for each spatial edge pair
4. Creating Seurat object with mechanismsxspatial_pairs dimensions
5. Similar metadata structure to RunCellToCell but restricted to spatial neighbors

---

### RunCellToNeighborhood
**File:** `/Users/msbr/GitHub/NICHESv2-msbr/R/RunCellToNeighborhood.R`

**Type:** Exported user-facing function

**Parameters:**
- `node.object`: Seurat object (preprocessed, required)
- `ground.truth`: list with source.subunits and target.subunits matrices (required)
- `assay`: string - assay to use (required)
- `meta.data.to.map`: character vector (required)
- `blend`: string, default "mean" - "sum" or "mean"
- `edgelist`: dataframe with "from" and "to" columns (required)
- `output_format`: string - "seurat" or "raw" (required)

**Returns:**
- If `output_format == "seurat"`: Seurat object with assay "CellToNeighborhood"
- If `output_format == "raw"`: List with `$CellToNeighborhoodMatrix` and `$metadata`

**Internal Calls:**
- `GetSeuratAssay()` - retrieves expression from spatial neighbors
- `rowsum()`, `rowMeans()` - aggregate by sending cell
- `Seurat::CreateSeuratObject()`
- `Seurat::NormalizeData()` - (if Seurat v5)
- `Seurat::Idents()`
- `Seurat::AddMetaData()`

**External Dependencies:**
- Seurat

**Description:**
Creates spatially-restricted cell-to-neighborhood (outgoing to spatial neighbors) matrix by:
1. Building ligand expression for sending cells in spatial edges
2. Building receptor expression for receiving cells in spatial edges
3. Computing interactions per edge, then aggregating by sender (sum/mean)
4. Creating Seurat object with mechanismsxsending_cells dimensions
5. Represents signaling output from each cell to its spatial neighborhood only

---

### RunNeighborhoodToCell
**File:** `/Users/msbr/GitHub/NICHESv2-msbr/R/RunNeighborhoodToCell.R`

**Type:** Exported user-facing function

**Parameters:**
- `node.object`: Seurat object (preprocessed, required)
- `ground.truth`: list with source.subunits and target.subunits matrices (required)
- `assay`: string - assay to use (required)
- `meta.data.to.map`: character vector (required)
- `blend`: string, default "mean" - "sum" or "mean"
- `edgelist`: dataframe with "from" and "to" columns (required)
- `output_format`: string - "seurat" or "raw" (required)

**Returns:**
- If `output_format == "seurat"`: Seurat object with assay "NeighborhoodToCell"
- If `output_format == "raw"`: List with `$NeighborhoodToCellMatrix` and `$metadata`

**Internal Calls:**
- `GetSeuratAssay()` - retrieves expression from spatial neighbors
- `rowsum()`, `rowMeans()` - aggregate by receiving cell
- `Seurat::CreateSeuratObject()`
- `Seurat::NormalizeData()` - (if Seurat v5)
- `Seurat::Idents()`
- `Seurat::AddMetaData()`

**External Dependencies:**
- Seurat

**Description:**
Creates spatially-restricted neighborhood-to-cell (spatial niche) matrix by:
1. Building ligand expression for sending cells in spatial edges
2. Building receptor expression for receiving cells in spatial edges
3. Computing interactions per edge, then aggregating by receiver (sum/mean)
4. Creating Seurat object with mechanismsxreceiving_cells dimensions
5. Represents spatial microenvironment input to each cell from spatial neighbors only

---

## GROUND TRUTH LOADING FUNCTIONS

### LoadGroundTruth
**File:** `/Users/msbr/GitHub/NICHESv2-msbr/R/LoadGroundTruth.R`

**Type:** Dispatcher/wrapper

**Parameters:**
- `LR.database`: string - "fantom5", "omnipath", or "custom"
- `custom_LR_database`: dataframe, default NULL - required if LR.database=="custom"
- `species`: string - "human", "mouse", "rat", "pig"
- `input_rownames`: character vector - gene names from input data

**Returns:** List with:
- `$source.subunits`: matrix of ligand subunits (rows=mechanisms, cols=subunit_columns, NA for absent subunits)
- `$target.subunits`: matrix of receptor subunits (rows=mechanisms, cols=subunit_columns, NA for absent subunits)

**Internal Calls:**
- `LoadFantom5()` - if LR.database=="fantom5"
- `LoadOmniPath()` - if LR.database=="omnipath"
- `LoadCustom()` - if LR.database=="custom"
- `FilterGroundTruth()` - filters to genes in input_rownames

**External Dependencies:** None (wrappers are loaded separately)

**Description:** 
Entry point for loading ligand-receptor databases. Routes to appropriate loader and filters results to only include mechanisms with genes present in the input dataset.

---

### LoadFantom5
**File:** `/Users/msbr/GitHub/NICHESv2-msbr/R/LoadFantom5.R`

**Type:** Exported user-facing function

**Parameters:**
- `species`: string (required) - "human", "mouse", "rat", or "pig"

**Returns:** List with:
- `$source.subunits`: matrix (mechanisms x 1 column for single ligand subunits)
- `$target.subunits`: matrix (mechanisms x 1 column for single receptor subunits)

**Internal Calls:**
- `get(utils::data())` - loads packaged data objects

**External Dependencies:**
- utils (for data loading)

**Description:**
Loads FANTOM5 ligand-receptor database for specified species. FANTOM5 assumes single subunits per ligand/receptor. Data is pre-packaged in R package (data/*.rda files).

Database Objects:
- `ncomms8866_human` - Human L-R pairs
- `ncomms8866_mouse` - Mouse L-R pairs
- `ncomms8866_rat` - Rat L-R pairs
- `ncomms8866_pig` - Pig L-R pairs

---

### LoadOmniPath
**File:** `/Users/msbr/GitHub/NICHESv2-msbr/R/LoadOmniPath.R`

**Type:** Exported user-facing function

**Parameters:**
- `species`: string (required) - "human", "mouse", or "rat"

**Returns:** List with:
- `$source.subunits`: matrix (mechanisms x subunit_columns)
- `$target.subunits`: matrix (mechanisms x subunit_columns)

**Internal Calls:**
- `OmnipathR::import_ligrecextra_interactions()` - downloads OmniPath data
- `dplyr::select()`, `dplyr::distinct()` - data processing
- `tidyr::separate()` - splits multi-subunit proteins by "_" separator
- `sapply()`, `strsplit()` - determine max subunit counts

**External Dependencies:**
- OmnipathR (online data source)
- dplyr
- tidyr

**Description:**
Loads OmniPath ligand-receptor interactions. OmniPath supports multi-subunit complexes where subunits are separated by underscores (e.g., "GATA6_GATA3" for a 2-subunit ligand). Downloads data dynamically from OmnipathR package.

---

### LoadCustom
**File:** `/Users/msbr/GitHub/NICHESv2-msbr/R/LoadCustom.R`

**Type:** Exported user-facing function

**Parameters:**
- `custom_LR_database`: dataframe (required)
  - Column 1: ligand gene names (subunits separated by "_")
  - Column 2: receptor gene names (subunits separated by "_")

**Returns:** List with:
- `$source.subunits`: matrix (mechanisms x subunit_columns)
- `$target.subunits`: matrix (mechanisms x subunit_columns)

**Internal Calls:**
- `tidyr::separate()` - splits subunit strings by "_"

**External Dependencies:**
- tidyr

**Description:**
Loads user-provided ligand-receptor database. Expects same format as OmniPath (underscore-separated subunits). Converts input dataframe to standardized subunit matrix format.

---

### FilterGroundTruth
**File:** `/Users/msbr/GitHub/NICHESv2-msbr/R/FilterGroundTruth.R`

**Type:** Exported utility function

**Parameters:**
- `ground.truth`: list with source.subunits and target.subunits
- `input_rownames`: character vector of genes in input dataset

**Returns:** Filtered ground.truth list (same structure, fewer mechanisms)

**Internal Calls:**
- `apply()` with custom function checking if all subunits are in input

**External Dependencies:** Base R only

**Description:**
Filters ground truth mechanisms to only those where:
1. ALL ligand subunits are present in input_rownames (or NA)
2. ALL receptor subunits are present in input_rownames (or NA)
3. Removes subunit columns where all values are NA after filtering

Reduces computation by eliminating L-R pairs that cannot be analyzed.

---

## HELPER/UTILITY FUNCTIONS

### PrepNodeObject
**File:** `/Users/msbr/GitHub/NICHESv2-msbr/R/PrepNodeObject.R`

**Type:** Exported utility function

**Parameters:**
- `node.object`: Seurat object (required)
- `assay`: string (required) - which assay to filter
- `min.cells.per.ident`: integer, default NULL - minimum cells per cell type
- `min.cells.per.gene`: integer, default NULL - minimum cells expressing gene

**Returns:** Filtered Seurat object

**Internal Calls:**
- `Seurat::DefaultAssay()`
- `Seurat::Idents()`
- `Seurat::subset()`
- `GetSeuratAssay()` - retrieves counts
- `Matrix::rowSums()` - count non-zero cells per gene

**External Dependencies:**
- Seurat
- Matrix

**Description:**
Preprocessing function that:
1. Sets default assay
2. Filters to cell populations with > min.cells.per.ident cells
3. Filters to genes expressed in > min.cells.per.gene cells
4. Returns counts and prints summary

---

### ComputeEdgelist
**File:** `/Users/msbr/GitHub/NICHESv2-msbr/R/ComputeEdgelist.R`

**Type:** Exported spatial utility function

**Parameters:**
- `node.object`: Seurat object (required)
- `position.x`: string (required) - metadata column with x-coordinates
- `position.y`: string (required) - metadata column with y-coordinates
- `k`: integer, default 4 - number of nearest neighbors
- `rad.set`: numeric, default NULL - radius threshold (ignored if k provided)
- `nn.method`: string, default 'aoz' - "aoz" (fast) or NULL (legacy)

**Returns:** Dataframe with columns:
- `from`: sending cell barcode
- `to`: receiving cell barcode

**Internal Calls:**
- If nn.method==NULL (legacy):
  - `apply()` - compute distances
  - `order()` - find k-nearest
  - `igraph::graph.adjacency()` - convert to igraph format
  - `igraph::get.data.frame()` - extract edgelist
- If nn.method=='aoz':
  - `spNNGP:::mkNNIndxCB()` - fast k-NN with spatial ordering
  - `data.table::rbindlist()` - concatenate edges
  - `na.omit()` - remove NAs
  - `order()` - reindex to original coordinates

**External Dependencies:**
- igraph (legacy method)
- spNNGP (fast method)
- data.table

**Description:**
Computes spatial neighbor graph. Two methods:
1. **Legacy (nn.method=NULL)**: Full distance matrix, filter k-NN, convert to mutual NN
2. **Fast (nn.method='aoz')**: Uses spNNGP package for efficient mutual nearest neighbor computation

Output is directed edgelist where each row is an edge from "from" cell to "to" cell.

---

### ReturnCelltypes (also referenced as return_celltypes)
**File:** `/Users/msbr/GitHub/NICHESv2-msbr/R/ReturnCellTypes.R`

**Type:** Exported utility function

**Parameters:**
- `seurat_object`: Seurat object (required)

**Returns:** Character vector of unique cell types (from Idents())

**Internal Calls:**
- `Seurat::Idents()`
- `table()`
- `names()`

**External Dependencies:**
- Seurat

**Description:**
Extracts unique cell types from Seurat object's identity. Prints warning reminder that active Idents must correspond to cell types.

---

### GetSeuratAssay
**File:** `/Users/msbr/GitHub/NICHESv2-msbr/R/getSeuratAssay.R`

**Type:** Exported utility function

**Parameters:**
- `object`: Seurat object (required)
- `assay`: string (required) - assay name
- `slot`: string (required) - "data", "counts", or "scale.data"

**Returns:** Matrix/sparse matrix of assay data

**Internal Calls:**
- `SeuratObject::Version()` - check Seurat version
- `Seurat::GetAssayData()` - version-aware accessor

**External Dependencies:**
- Seurat
- SeuratObject

**Description:**
Wrapper for GetAssayData that handles Seurat v4 vs v5 API differences:
- v5: uses `layer` parameter
- v4: uses `slot` parameter

---

### FilterGroundTruth
**File:** `/Users/msbr/GitHub/NICHESv2-msbr/R/FilterGroundTruth.R`

**Type:** Exported utility function

**Parameters:**
- `ground.truth`: list with source.subunits and target.subunits
- `input_rownames`: character vector

**Returns:** Filtered ground.truth list

**Internal Calls:**
- `apply()` - check subunit presence
- `which()` - find matching indices
- `intersect()` - find mechanisms with both ligands and receptors
- `colSums()` - identify empty columns

**External Dependencies:** Base R only

**Description:**
Filters ground truth to mechanisms present in input data. Removes rows where any required subunit is missing, removes empty subunit columns.

---

## VISUALIZATION FUNCTIONS

### CircuitPlot
**File:** `/Users/msbr/GitHub/NICHESv2-msbr/R/CircuitPlot.R`

**Type:** Exported visualization function

**Parameters:**
- `transcr.obj`: Seurat object (node/transcriptomic data, required)
- `connect.obj`: Seurat object (edge/connectomic data, required)
- `feature`: string (required) - L-R mechanism with "|" separator (e.g., "TGFB1|TGFBR1")
- `plot.function`: string, default 'ggCircuit' - plotting function to use
- `group.by`: string, default 'CellClass' - metadata column for grouping
- `graph.angle`: numeric, default 45 - rotation angle
- `h`: numeric, default 0.2 - edge offset from node center
- `offset`: numeric, default 0.05 - spacing between paracrine edges
- `autocrine.offset`: numeric, default 0.03 - spacing for autocrine loops
- `edge.scale.factor`: numeric, default 20 - controls edge thickness
- `arrow.head.angle`: numeric, default 15 - arrowhead angle
- `arrow.head.length`: numeric, default 0.03 - arrowhead length
- `autocrine.arrow.curvature`: numeric, default 8 - curvature of autocrine loops
- `cols.use`: character vector - colors for cell types
- `edge.fixed.size`: logical, default TRUE - fixed or variable edge thickness
- `global.node.list`: character vector - node labels to include
- `min.edge.value`: numeric, default 0 - minimum connectivity for visualization
- `max.edge.value`: numeric, default NULL - maximum connectivity (auto if NULL)
- `unity.normalize`: logical, default FALSE - normalize connectivity values to 0-1

**Returns:** ggplot2 object

**Internal Calls:**
- `DefineNodeObject()` - extract node-level data
- `DefineEdgeObject()` - extract edge-level data
- `AggregateNodeData()` - aggregate by group.by
- `AggregateEdgeData()` - aggregate by metadata column
- `ggCircuit()` - create plot

**External Dependencies:**
- ggplot2
- scales
- RColorBrewer

**Description:**
Creates circular "circuit" plot of cell-cell signaling. Shows:
- Cell types as nodes around circle (sized by system fraction)
- Arrows between nodes (paracrine signaling)
- Curved loops (autocrine signaling)
- Color/size represent connectivity strength

---

### ggCircuit
**File:** `/Users/msbr/GitHub/NICHESv2-msbr/R/ggCircuit.R`

**Type:** Internal visualization helper

**Parameters:**
- `edge.aggregate`: aggregated edge dataframe
- `node.aggregate`: aggregated node dataframe
- `graph.angle`: numeric - rotation
- `h`: numeric - edge offset
- `offset`: numeric - paracrine spacing
- `autocrine.offset`: numeric - autocrine spacing
- `edge.scale.factor`: numeric
- `arrow.head.angle`: numeric
- `arrow.head.length`: numeric
- `autocrine.arrow.curvature`: numeric
- `cols.use`: character vector
- `edge.fixed.size`: logical
- `split.by`: logical, default FALSE
- `min.edge.value`: numeric, default NULL
- `max.edge.value`: numeric, default NULL
- `title`: string, default NULL

**Returns:** ggplot2 object

**Internal Calls:**
- `netCoin::layoutCircle()` - compute node positions
- `atan()`, `cos()`, `sin()` - trigonometry for edge placement
- `ggplot()`, `geom_segment()`, `geom_curve()`, `geom_point()` - plotting

**External Dependencies:**
- ggplot2
- netCoin
- grid (for arrows)
- Seurat

**Description:**
Low-level plotting function that:
1. Places nodes on circle
2. Computes edge start/end coordinates
3. Offsets paracrine edges to avoid overlap
4. Creates curved loops for autocrine signals
5. Colors and sizes based on connectivity

---

### NetworkPlot
**File:** `/Users/msbr/GitHub/NICHESv2-msbr/R/NetworkPlot.R`

**Type:** Exported visualization function

**Parameters:**
- `transcriptome.object`: Seurat object (required) - node data with UMAP reduction
- `connectome.object`: Seurat object (required) - edge data
- `mechanism.of.interest`: string (required) - L-R mechanism
- `legends.to.plot`: character vector, default NULL - metadata columns for legends
- `split.by`: string, default NULL - metadata column for splitting plots
- `transpose`: logical, default FALSE - plot orientation
- `reduction`: string, default 'umap.rpca' - dimensionality reduction to plot
- `min.connectivity.thresh`: numeric, default 0 - minimum connectivity to display
- `point.size`: numeric, default 0.1
- `alpha.min`: numeric, default 0
- `alpha.max`: numeric, default 0.85
- `line.thickness`: numeric, default 0.2
- `arrowhead.size`: numeric, default 0
- `connectivity.color.min`: numeric, default 0
- `connectivity.color.max`: numeric, default NULL (auto)
- `ligand.color.min`: numeric, default 0
- `receptor.color.min`: numeric, default 0
- `legend.palettes`: list, default NULL - named color palettes
- `legend.alpha`: numeric, default 0.25
- `black.points`: logical, default TRUE
- `node.assay`: string, default 'RNA'
- `edge.assay`: string, default 'CellToCell'

**Returns:** ggplot2 object (or patchwork if split.by or legends present)

**Internal Calls:**
- `strsplit()` - parse mechanism
- `ggplot()`, `geom_segment()`, `geom_point()` - plotting
- `ggnewscale::new_scale_color()` - multiple color scales
- `Seurat::DimPlot()` - legend plots
- `cowplot::plot_grid()` - arrange subplots
- `stringr::str_split_fixed()` - parse barcodes

**External Dependencies:**
- ggplot2
- stringr
- ggnewscale
- cowplot
- Seurat

**Description:**
Creates network visualization overlaid on UMAP/embedding:
- Background: dimensionality reduction plot with cell positions
- Edges: directed lines showing signaling (color/alpha = strength)
- Optional: Ligand/receptor expression colors on cells
- Optional: Split views by metadata groups
- Optional: Legend plots for metadata annotations

---

### DotPlotSCC
**File:** `/Users/msbr/GitHub/NICHESv2-msbr/R/DotPlotSCC.R`

**Type:** Exported visualization function

**Parameters:**
- `object`: Seurat object (required)
- `assay`: string, default NULL (uses DefaultAssay)
- `features`: character vector (required) - L-R mechanisms to plot
- `cols`: character vector, default c("lightgrey", "blue") - color scale
- `col.min`: numeric, default -2.5
- `col.max`: numeric, default 2.5
- `dot.min`: numeric, default 0
- `dot.scale`: numeric, default 6
- `idents`: character vector, default NULL
- `group.by`: string, default NULL
- `split.by`: string, default NULL
- `cluster.idents`: logical, default FALSE
- `scale`: logical, default TRUE
- `scale.by`: string, default "radius"
- `scale.min`: numeric, default NA
- `scale.max`: numeric, default NA

**Returns:** ggplot2 object

**Internal Calls:**
- `Seurat::DefaultAssay()`
- `Seurat::CellsByIdentities()`
- `Seurat::FetchData()`
- `Seurat::Idents()`
- `ggplot()`, `geom_point()` - plotting
- `stats::hclust()`, `stats::dist()` - hierarchical clustering
- `RColorBrewer::brewer.pal.info` - color checking

**External Dependencies:**
- ggplot2
- RColorBrewer
- cowplot
- Seurat

**Description:**
Modified dot plot for signaling (CellToCell) data:
- X-axis: L-R mechanisms
- Y-axis: cell types
- Dot size: percent of cells with interaction
- Dot color: average interaction strength
- Supports faceting by features, splitting, clustering

---

## ANALYSIS FUNCTIONS

### AggregateNodeData
**File:** `/Users/msbr/GitHub/NICHESv2-msbr/R/AggregateNodeData.R`

**Type:** Helper function for CircuitPlot

**Parameters:**
- `node.object`: dataframe (required) - node-level data
- `group.by`: string (required) - metadata column for grouping
- `global.node.list`: character vector (required) - node labels

**Returns:** Dataframe with:
- `node.label`: node name
- `system.fraction`: proportion in system

**Internal Calls:**
- `table()`, `prop.table()` - compute proportions

**External Dependencies:** Base R only

**Description:**
Aggregates node-level data by grouping variable. Ensures all nodes in global.node.list are represented (0 fraction if absent).

---

### AggregateEdgeData
**File:** `/Users/msbr/GitHub/NICHESv2-msbr/R/AggregateEdgeData.R`

**Type:** Helper function for CircuitPlot

**Parameters:**
- `edge.object`: dataframe (required) - edge-level data
- `group.by`: string (required) - metadata column for grouping

**Returns:** Dataframe with aggregated edge features and parsed sending/receiving labels

**Internal Calls:**
- `stats::aggregate()` - aggregate by group
- `stringr::str_split_fixed()` - parse edge labels

**External Dependencies:**
- stringr

**Description:**
Aggregates edge-level data by grouping variable, computing mean connectivity. Extracts sending and receiving labels from edge identifiers.

---

### DefineNodeObject
**File:** `/Users/msbr/GitHub/NICHESv2-msbr/R/DefineNodeObject.R`

**Type:** Internal helper

**Parameters:**
- `transcr.obj`: Seurat object (required)
- `feature`: string (required) - gene name
- `assay`: string, default 'RNA'

**Returns:** Dataframe with:
- `feature.value`: expression values
- `barcode`: cell barcodes
- Columns from transcr.obj@meta.data

**Internal Calls:**
- Accesses `transcr.obj@assays[[assay]]@layers$data` (Seurat v5)

**External Dependencies:**
- Seurat (v5)

**Description:**
Extracts expression data for a single gene and combines with metadata.

---

### DefineEdgeObject
**File:** `/Users/msbr/GitHub/NICHESv2-msbr/R/DefineEdgeObject.R`

**Type:** Internal helper

**Parameters:**
- `connect.obj`: Seurat object (required)
- `feature`: string (required) - L-R mechanism
- `assay`: string, default 'CellToCell'

**Returns:** Dataframe with:
- `feature.value`: interaction strength
- `barcode.pair`: cell pair identifiers
- Columns from connect.obj@meta.data

**Internal Calls:**
- Accesses `connect.obj@assays[[assay]]@layers$data` (Seurat v5)

**External Dependencies:**
- Seurat (v5)

**Description:**
Extracts interaction data for a single mechanism and combines with metadata.

---

### AddMetaData
**File:** `/Users/msbr/GitHub/NICHESv2-msbr/R/AddMetaData.R`

**Type:** Internal/utility function (incomplete/experimental)

**Parameters:**
- `CTC.obj`: Seurat object with edge data (required)
- `annotated.obj`: Seurat object with new metadata (required)
- `col.name`: string - metadata column to add (required)

**Returns:** Updated metadata dataframe with new columns

**Internal Calls:**
- `strsplit()` - parse barcodes from edge identifiers
- `dplyr::left_join()` - join metadata
- `dplyr::rename()` - rename columns

**External Dependencies:**
- dplyr

**Description:**
Utility for adding metadata from annotated object to edge object. Creates separate Sending/Receiving/Joint versions of metadata.

**Note:** Contains inline example code at bottom; appears to be under development.

---

### compute_Neighborhood_Composition
**File:** `/Users/msbr/GitHub/NICHESv2-msbr/R/CellNeighborhood_Composition.R` (referenced, not fully reviewed)

**Type:** Spatial analysis helper

**Description:** Computes cell type composition within spatial neighborhoods.

---

### create_Neighborhood_Assay
**File:** `/Users/msbr/GitHub/NICHESv2-msbr/R/create_Neighborhood_Assay.R`

**Type:** Spatial analysis function

**Parameters:**
- `obj`: Seurat object (required)
- `position.x`: string, default 'x'
- `position.y`: string, default 'y'
- `k`: integer, default NULL
- `rad.set`: numeric, default 25
- `composition`: string (required) - metadata column for cell type composition
- `TMA`: string, default NULL - sample identifier column
- `TMA.list`: character vector, default NULL - list of samples

**Returns:** List with:
- `$rawcounts`: Seurat Assay5 object
- `$pct`: Seurat Assay5 object

**Internal Calls:**
- `compute_Neighborhood_Composition()` - main computation
- `Seurat::CreateAssay5Object()` - create assays
- `dplyr::bind_rows()` - concatenate results

**External Dependencies:**
- Seurat
- dplyr

**Description:**
Creates neighborhood composition assay for spatial data. Can work on whole dataset or split by sample/TMA. Returns both raw counts and percentages of cell types in spatial neighborhoods.

---

### CoLocalization
**File:** `/Users/msbr/GitHub/NICHESv2-msbr/R/CoLocalization.R` (incomplete)

**Type:** Spatial analysis function (experimental/incomplete)

**Description:** Appears to compute co-localization probabilities but is hardcoded for specific dataset. Not generalized.

---

## FUNCTION CALL HIERARCHY

### Top-Level Entry Point
```
RunNICHES (generic)
  ├─ RunNICHES.Seurat
  │   └─ RunNICHES.default
  ├─ RunNICHES.matrix
  │   └─ RunNICHES.default
  └─ RunNICHES.Matrix (alias to RunNICHES.matrix)
```

### RunNICHES.default Workflow
```
RunNICHES.default
  ├─ [Validation & parameter checking]
  ├─ PrepNodeObject
  ├─ LoadGroundTruth
  │   ├─ LoadFantom5 OR LoadOmniPath OR LoadCustom
  │   └─ FilterGroundTruth
  ├─ ComputeEdgelist [if spatial modes enabled]
  │   ├─ (legacy: igraph, distance matrix)
  │   └─ (fast: spNNGP)
  └─ [Call enabled organization functions]
      ├─ RunCellToCell
      │   ├─ return_celltypes
      │   ├─ GetSeuratAssay (for each celltype)
      │   └─ [Metadata assembly with abind, AddMetaData]
      ├─ RunCellToSystem
      │   ├─ GetSeuratAssay
      │   ├─ Matrix::rowSums/rowMeans
      │   └─ [Metadata assembly]
      ├─ RunSystemToCell
      │   ├─ GetSeuratAssay
      │   ├─ Matrix::rowSums/rowMeans
      │   └─ [Metadata assembly]
      ├─ RunCellToCellSpatial [if spatial & CellToCellSpatial]
      │   ├─ GetSeuratAssay (subset to edgelist cells)
      │   └─ [Metadata assembly with abind]
      ├─ RunCellToNeighborhood [if spatial & CellToNeighborhood]
      │   ├─ GetSeuratAssay (subset to edgelist cells)
      │   └─ [Aggregate by sender]
      └─ RunNeighborhoodToCell [if spatial & NeighborhoodToCell]
          ├─ GetSeuratAssay (subset to edgelist cells)
          └─ [Aggregate by receiver]
```

### Visualization Workflow
```
CircuitPlot
  ├─ DefineNodeObject
  ├─ DefineEdgeObject
  ├─ AggregateNodeData
  │   └─ global.node.list provided by user
  ├─ AggregateEdgeData
  └─ ggCircuit
      ├─ netCoin::layoutCircle (node positions)
      ├─ [Trigonometry for edge placement]
      ├─ ggplot, geom_segment, geom_curve, geom_point
      └─ Seurat::NoAxes, gg_color_hue

NetworkPlot
  ├─ [Data extraction & validation]
  ├─ [UMAP coordinate extraction]
  ├─ [Connectivity data preparation]
  ├─ [Edge thresholding & log transform]
  ├─ ggplot, geom_segment, geom_point (base layer)
  ├─ ggnewscale::new_scale_color (multiple scales)
  ├─ cowplot::plot_grid (if split.by or legends)
  └─ Seurat::DimPlot (for legends)

DotPlotSCC
  ├─ Seurat::CellsByIdentities
  ├─ Seurat::FetchData
  ├─ [Scaling & normalization]
  ├─ [Hierarchical clustering if requested]
  └─ ggplot, geom_point, facet_grid
```

### Spatial Analysis Workflow
```
create_Neighborhood_Assay
  ├─ [Sample iteration if TMA provided]
  ├─ compute_Neighborhood_Composition (for each sample)
  │   └─ [Spatial neighborhood computation]
  ├─ dplyr::bind_rows (if multi-sample)
  ├─ Seurat::CreateAssay5Object
  └─ Return list(rawcounts, pct)
```

---

## EXTERNAL PACKAGE DEPENDENCIES

### Core Analysis
- **Seurat** (>= 4.0) - Single-cell object handling, core functions
- **SeuratObject** - Seurat v5 compatibility
- **Matrix** - Sparse matrix operations
- **dplyr** - Data manipulation (pipes, joins, select)
- **tidyr** - Data tidying (separate)

### Spatial Analysis
- **igraph** - Graph operations (legacy edgelist method)
- **spNNGP** - Fast spatial nearest neighbor computation (fast edgelist method)
- **data.table** - Fast data operations
- **netCoin** - Circular layout for circuit plots

### Visualization
- **ggplot2** - Core visualization
- **ggnewscale** - Multiple color scales in ggplot
- **cowplot** - Plot arrangement and theming
- **RColorBrewer** - Color palettes

### Spatial Data Integration
- **OmnipathR** - Download OmniPath database

### Utility
- **abind** - Array binding for metadata
- **stringr** - String manipulation
- **scales** - Scale utilities (rescale)
- **utils** - Data loading utilities
- **grid** - Arrow graphics
- **stats** - Statistical functions (aggregate, hclust, dist)

---

## GLOBAL VARIABLES

**File:** `/Users/msbr/GitHub/NICHESv2-msbr/R/globals.R`

Global variables registered for R CMD check compliance:
- `ncomms8866_human` - FANTOM5 human database
- `ncomms8866_mouse` - FANTOM5 mouse database
- `ncomms8866_rat` - FANTOM5 rat database
- `ncomms8866_pig` - FANTOM5 pig database
- `source_genesymbol` - Column reference (OmniPath)
- `target_genesymbol` - Column reference (OmniPath)
- `non.zero.cells` - Data frame column reference

---

## DATA OBJECTS

**Location:** `data/*.rda` files

Pre-packaged ligand-receptor databases:
- `ncomms8866_human.rda` - Human L-R pairs (FANTOM5)
- `ncomms8866_mouse.rda` - Mouse L-R pairs (FANTOM5)
- `ncomms8866_rat.rda` - Rat L-R pairs (FANTOM5)
- `ncomms8866_pig.rda` - Pig L-R pairs (FANTOM5)

---

## EXAMPLE WORKFLOW

### Non-spatial analysis:
```r
library(NICHES)

# Create NICHES objects
niches.output <- RunNICHES(
  object = seurat_object,
  species = "human",
  LR.database = "fantom5",
  CellToCell = TRUE,
  CellToSystem = TRUE,
  SystemToCell = TRUE,
  output_format = "seurat"
)

# Extract outputs
ctc <- niches.output$CellToCell  # Cell-to-cell interactions
cts <- niches.output$CellToSystem # Outgoing signaling
stc <- niches.output$SystemToCell # Incoming signaling
```

### Spatial analysis:
```r
# Spatial NICHES analysis with k-NN
niches.spatial <- RunNICHES(
  object = spatial_seurat,
  species = "human",
  LR.database = "omnipath",
  position.x = "x",
  position.y = "y",
  k = 5,
  CellToCellSpatial = TRUE,
  CellToNeighborhood = TRUE,
  NeighborhoodToCell = TRUE,
  output_format = "seurat"
)

# Extract spatial outputs
ctc.spatial <- niches.spatial$CellToCellSpatial
ctn <- niches.spatial$CellToNeighborhood
ntc <- niches.spatial$NeighborhoodToCell
```

### Visualization:
```r
# Circuit plot for a specific mechanism
circuit <- CircuitPlot(
  transcr.obj = seurat_object,
  connect.obj = niches.output$CellToCell,
  feature = "TGFB1|TGFBR1",
  group.by = "cell_type"
)

# Network plot
network <- NetworkPlot(
  transcriptome.object = seurat_object,
  connectome.object = niches.output$CellToCell,
  mechanism.of.interest = "IL6|IL6R",
  reduction = "umap"
)
```

