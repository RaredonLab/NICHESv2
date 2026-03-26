# NICHES Package - Function Architecture Diagram

## Executive Summary

The NICHES package follows a **modular pipeline architecture** with:
- **1 generic entry point** (RunNICHES) using S3 dispatch
- **6 core transformation functions** (Run* family) that convert gene expression → signaling space
- **3 database loaders** (Load* family) for ligand-receptor ground truth
- **4 visualization functions** with helper utilities
- **Clear separation** between spatial and non-spatial analyses

---

## 1. MAIN DATA FLOW ARCHITECTURE

```
┌─────────────────────────────────────────────────────────────────────────┐
│                         USER INPUT                                      │
│         (Seurat object, matrix, or sparse Matrix)                       │
└────────────────────────────┬────────────────────────────────────────────┘
                             │
                             ▼
┌─────────────────────────────────────────────────────────────────────────┐
│                        RunNICHES (Generic)                              │
│                       S3 Method Dispatcher                              │
└──────────┬──────────────────┬────────────────┬─────────────────────────┘
           │                  │                │
           ▼                  ▼                ▼
    .Seurat()          .matrix()        .Matrix()
           │                  │                │
           └──────────────────┴────────────────┘
                             │
                             ▼
┌─────────────────────────────────────────────────────────────────────────┐
│                     RunNICHES.default()                                 │
│                    CORE ORCHESTRATION                                   │
│                                                                         │
│  ┌──────────────────────────────────────────────────────────────────┐  │
│  │ Step 1: PrepNodeObject()                                         │  │
│  │   • Filter cells by min.cells.per.ident                          │  │
│  │   • Filter genes by min.cells.per.gene                           │  │
│  │   • Returns cleaned Seurat object                                │  │
│  └──────────────────────────────────────────────────────────────────┘  │
│                             │                                           │
│  ┌──────────────────────────▼───────────────────────────────────────┐  │
│  │ Step 2: LoadGroundTruth()                                        │  │
│  │   ├─ LoadFantom5() OR LoadOmniPath() OR LoadCustom()            │  │
│  │   └─ FilterGroundTruth() (restrict to genes in dataset)         │  │
│  │   • Returns list(source.subunits, target.subunits)               │  │
│  └──────────────────────────────────────────────────────────────────┘  │
│                             │                                           │
│  ┌──────────────────────────▼───────────────────────────────────────┐  │
│  │ Step 3: [IF SPATIAL] ComputeEdgelist()                           │  │
│  │   • Compute k-NN or radius-based spatial neighbor graph          │  │
│  │   • Fast method: spNNGP (default)                                │  │
│  │   • Legacy method: igraph with full distance matrix              │  │
│  │   • Returns edgelist: dataframe(from, to)                        │  │
│  └──────────────────────────────────────────────────────────────────┘  │
│                             │                                           │
│  ┌──────────────────────────▼───────────────────────────────────────┐  │
│  │ Step 4: BRANCHING - Call enabled organization functions          │  │
│  │                                                                   │  │
│  │   NON-SPATIAL MODES (all-to-all connectivity):                   │  │
│  │   ├─ [IF CellToCell] RunCellToCell()                            │  │
│  │   ├─ [IF CellToSystem] RunCellToSystem()                        │  │
│  │   └─ [IF SystemToCell] RunSystemToCell()                        │  │
│  │                                                                   │  │
│  │   SPATIAL MODES (neighbor-constrained):                          │  │
│  │   ├─ [IF CellToCellSpatial] RunCellToCellSpatial()              │  │
│  │   ├─ [IF CellToNeighborhood] RunCellToNeighborhood()            │  │
│  │   └─ [IF NeighborhoodToCell] RunNeighborhoodToCell()            │  │
│  └──────────────────────────────────────────────────────────────────┘  │
└────────────────────────────┬────────────────────────────────────────────┘
                             │
                             ▼
┌─────────────────────────────────────────────────────────────────────────┐
│                       OUTPUT                                            │
│   Named list of Seurat objects (or raw matrices)                       │
│   - $CellToCell, $CellToSystem, $SystemToCell                          │
│   - $CellToCellSpatial, $CellToNeighborhood, $NeighborhoodToCell       │
└─────────────────────────────────────────────────────────────────────────┘
```

---

## 2. FUNCTION DEPENDENCY GRAPH

### Core Entry Point Dependencies

```
RunNICHES (generic)
├── RunNICHES.Seurat ───────┐
│   └── [validates metadata] │
├── RunNICHES.matrix ────────┼──→ RunNICHES.default ─┐
│   └── [creates Seurat]     │                       │
└── RunNICHES.Matrix ────────┘                       │
    └── [alias to .matrix]                           │
                                                     │
                                                     ▼
                              ┌────────────────────────────────────────┐
                              │     RunNICHES.default                  │
                              │  (MAIN ORCHESTRATION PIPELINE)         │
                              └────────┬───────────────────────────────┘
                                       │
           ┌───────────────────────────┼───────────────────────────┐
           │                           │                           │
           ▼                           ▼                           ▼
    PrepNodeObject            LoadGroundTruth            ComputeEdgelist
           │                           │                    [if spatial]
           │                           │                           │
           │          ┌────────────────┴───────────┐              │
           │          │                            │              │
           │          ▼                            ▼              │
           │    LoadFantom5()              LoadOmniPath()         │
           │    LoadCustom()                     │               │
           │          │                           │               │
           │          └────────────┬──────────────┘               │
           │                       ▼                              │
           │              FilterGroundTruth()                     │
           │                       │                              │
           └───────────────────────┴──────────────────────────────┘
                                   │
                ┌──────────────────┴──────────────────┐
                │         BRANCHING LOGIC             │
                │   (Based on boolean flags)          │
                └──────────────────┬──────────────────┘
                                   │
      ┌────────────┬───────────────┼───────────────┬────────────┬────────────┐
      │            │               │               │            │            │
      ▼            ▼               ▼               ▼            ▼            ▼
RunCellToCell  RunCellToSystem  RunSystemToCell  RunCellTo   RunCellTo  RunNeighborhood
                                                  CellSpatial  Neighborhood  ToCell
      │            │               │               │            │            │
      └────────────┴───────────────┴───────────────┴────────────┴────────────┘
                                   │
                                   ▼
                          Named list output
```

### Organization Function Dependencies

```
┌───────────────────────────────────────────────────────────────────┐
│                 ALL Run* Functions Share:                         │
│                                                                   │
│  ┌─────────────────────────────────────────────────────────────┐ │
│  │ GetSeuratAssay()                                            │ │
│  │   • Version-aware Seurat accessor (v4/v5 compatibility)    │ │
│  │   • Called multiple times per Run* function                │ │
│  │   • Inputs: object, assay, slot                            │ │
│  │   • Returns: expression matrix                             │ │
│  └─────────────────────────────────────────────────────────────┘ │
│                                                                   │
│  ┌─────────────────────────────────────────────────────────────┐ │
│  │ Seurat object creation pipeline                            │ │
│  │   1. CreateSeuratObject(counts = interaction_matrix)       │ │
│  │   2. NormalizeData()                                       │ │
│  │   3. [Seurat v5 ONLY] Manually create data layer          │ │
│  │   4. Set Idents()                                          │ │
│  │   5. AddMetaData()                                         │ │
│  └─────────────────────────────────────────────────────────────┘ │
│                                                                   │
│  ┌─────────────────────────────────────────────────────────────┐ │
│  │ Metadata propagation                                       │ │
│  │   • For each column in meta.data.to.map:                   │ │
│  │     - .Sending (from sender cell)                          │ │
│  │     - .Receiving (from receiver cell)                      │ │
│  │     - .Joint (concatenated "sender - receiver")            │ │
│  └─────────────────────────────────────────────────────────────┘ │
└───────────────────────────────────────────────────────────────────┘
```

---

## 3. DETAILED FUNCTION SPECIFICATIONS

### 3A. Core Organization Functions

| Function | Input Dimensions | Output Dimensions | Key Operation | Spatial? |
|----------|-----------------|-------------------|---------------|----------|
| **RunCellToCell** | genes × cells | mechanisms × cell_pairs | Downsample to min(sender,receiver) per type pair | No |
| **RunCellToSystem** | genes × cells | mechanisms × sending_cells | Aggregate all receptors (sum/mean/mean.adj) | No |
| **RunSystemToCell** | genes × cells | mechanisms × receiving_cells | Aggregate all ligands (sum/mean/mean.adj) | No |
| **RunCellToCellSpatial** | genes × cells | mechanisms × spatial_pairs | Filter to edgelist neighbors only | Yes |
| **RunCellToNeighborhood** | genes × cells | mechanisms × sending_cells | Aggregate receptors of neighbors per sender | Yes |
| **RunNeighborhoodToCell** | genes × cells | mechanisms × receiving_cells | Aggregate ligands of neighbors per receiver | Yes |

### 3B. Organization Function Details

#### RunCellToCell
```
INPUTS:
  • node.object: Seurat (preprocessed)
  • ground.truth: list(source.subunits, target.subunits)
  • assay: string
  • meta.data.to.map: character vector
  • output_format: "seurat" or "raw"

ALGORITHM:
  1. Extract unique cell types → return_celltypes()
  2. FOR each sending celltype:
       a. Build ligand expression matrix (multiply subunits)
          → GetSeuratAssay() × num_mechanisms
       b. FOR each receiving celltype:
            i. Build receptor expression matrix
            ii. Downsample edges: min(n_sender, n_receiver)
            iii. Compute interaction: ligand[sender] × receptor[receiver]
  3. Concatenate all pairs into single matrix
  4. Create Seurat object with metadata (SendingType, ReceivingType, VectorType)
  5. Propagate meta.data.to.map columns (.Sending, .Receiving, .Joint)

OUTPUT:
  • Seurat object: rows=mechanisms, cols=barcode_pairs
  • OR raw matrix + metadata dataframe

DEPENDENCIES:
  • return_celltypes()
  • GetSeuratAssay() [MANY CALLS]
  • Seurat: subset(), CreateSeuratObject(), NormalizeData(), AddMetaData()
  • abind::abind()
```

#### RunCellToSystem
```
INPUTS: (same parameters as RunCellToCell except no downsampling needed)

ALGORITHM:
  1. Build ligand expression matrix across ALL cells
  2. Build receptor expression matrix across ALL cells
  3. Aggregate receptor expression:
       • "sum": rowSums()
       • "mean": rowMeans()
       • "mean.adj": rowMeans() × sqrt(ncol)
  4. FOR each sending cell:
       Compute interaction: ligand[sender] × aggregated_receptor
  5. Create Seurat object (mechanisms × sending_cells)
  6. Propagate metadata from sending cells only

OUTPUT:
  • Seurat object: rows=mechanisms, cols=sending_cells
  • Represents total signaling OUTPUT from each cell

DEPENDENCIES:
  • GetSeuratAssay() [2 CALLS]
  • Matrix::rowSums(), Matrix::rowMeans()
  • Seurat: CreateSeuratObject(), NormalizeData(), AddMetaData()
```

#### RunSystemToCell
```
INPUTS: (same parameters as RunCellToSystem)

ALGORITHM:
  1. Build ligand expression matrix across ALL cells
  2. Build receptor expression matrix across ALL cells
  3. Aggregate ligand expression:
       • "sum": rowSums()
       • "mean": rowMeans()
       • "mean.adj": rowMeans() × sqrt(ncol)
  4. FOR each receiving cell:
       Compute interaction: aggregated_ligand × receptor[receiver]
  5. Create Seurat object (mechanisms × receiving_cells)
  6. Propagate metadata from receiving cells only

OUTPUT:
  • Seurat object: rows=mechanisms, cols=receiving_cells
  • Represents total signaling INPUT to each cell (microenvironment/niche)

DEPENDENCIES:
  • GetSeuratAssay() [2 CALLS]
  • Matrix::rowSums(), Matrix::rowMeans()
  • Seurat: CreateSeuratObject(), NormalizeData(), AddMetaData()
```

#### RunCellToCellSpatial
```
INPUTS:
  • node.object, ground.truth, assay, meta.data.to.map, output_format
  • edgelist: dataframe(from, to) ← OUTPUT FROM ComputeEdgelist()

ALGORITHM:
  1. Build ligand expression for cells in edgelist$from
  2. Build receptor expression for cells in edgelist$to
  3. FOR each row in edgelist:
       Compute interaction: ligand[from] × receptor[to]
  4. Create Seurat object (mechanisms × spatial_pairs)
  5. Propagate metadata (SendingType, ReceivingType, VectorType)

OUTPUT:
  • Seurat object: rows=mechanisms, cols=spatial_barcode_pairs
  • Similar to RunCellToCell but restricted to spatial neighbors

DEPENDENCIES:
  • GetSeuratAssay() [EDGELIST-DRIVEN]
  • Seurat: CreateSeuratObject(), NormalizeData(), AddMetaData()
  • abind::abind()
```

#### RunCellToNeighborhood
```
INPUTS: (adds blend parameter, requires edgelist)

ALGORITHM:
  1. Build ligand/receptor matrices for edgelist cells
  2. Compute interactions for each edge
  3. Aggregate by sending cell (group by edgelist$from):
       • "sum": sum all outgoing edges per sender
       • "mean": mean all outgoing edges per sender
  4. Create Seurat object (mechanisms × sending_cells)

OUTPUT:
  • Seurat object: rows=mechanisms, cols=sending_cells
  • Represents signaling OUTPUT to spatial neighborhood

DEPENDENCIES:
  • GetSeuratAssay()
  • rowsum(), rowMeans() [base R aggregation]
  • Seurat: CreateSeuratObject(), NormalizeData(), AddMetaData()
```

#### RunNeighborhoodToCell
```
INPUTS: (same as RunCellToNeighborhood)

ALGORITHM:
  1. Build ligand/receptor matrices for edgelist cells
  2. Compute interactions for each edge
  3. Aggregate by receiving cell (group by edgelist$to):
       • "sum": sum all incoming edges per receiver
       • "mean": mean all incoming edges per receiver
  4. Create Seurat object (mechanisms × receiving_cells)

OUTPUT:
  • Seurat object: rows=mechanisms, cols=receiving_cells
  • Represents signaling INPUT from spatial neighborhood (spatial niche)

DEPENDENCIES:
  • GetSeuratAssay()
  • rowsum(), rowMeans()
  • Seurat: CreateSeuratObject(), NormalizeData(), AddMetaData()
```

---

## 4. GROUND TRUTH LOADING SYSTEM

```
┌────────────────────────────────────────────────────────────┐
│                   LoadGroundTruth()                        │
│                    (Dispatcher)                            │
│                                                            │
│  Inputs:                                                   │
│    • LR.database: "fantom5", "omnipath", or "custom"      │
│    • species: "human", "mouse", "rat", "pig"              │
│    • custom_LR_database: dataframe (if custom)            │
│    • input_rownames: genes in dataset                     │
│                                                            │
│  Routing Logic:                                            │
│    IF LR.database == "fantom5" → LoadFantom5()            │
│    IF LR.database == "omnipath" → LoadOmniPath()          │
│    IF LR.database == "custom" → LoadCustom()              │
│                                                            │
│  Post-processing:                                          │
│    → FilterGroundTruth() (filter to input genes)          │
│                                                            │
│  Output:                                                   │
│    list(source.subunits, target.subunits)                 │
│      • Matrices: rows=mechanisms, cols=subunits           │
│      • NA values for unused subunit columns               │
└────────────────────────────────────────────────────────────┘
                             │
         ┌───────────────────┼───────────────────┐
         │                   │                   │
         ▼                   ▼                   ▼
┌─────────────────┐  ┌──────────────────┐  ┌──────────────┐
│ LoadFantom5()   │  │ LoadOmniPath()   │  │ LoadCustom() │
│                 │  │                  │  │              │
│ • Single subunit│  │ • Multi-subunit  │  │ • User-      │
│   per L/R       │  │   complexes      │  │   provided   │
│ • Pre-packaged  │  │ • Downloaded     │  │ • Flexible   │
│   in data/*.rda │  │   from OmnipathR │  │   format     │
│ • Species:      │  │ • Species:       │  │              │
│   h/m/r/p       │  │   h/m/r only     │  │              │
│                 │  │ • Uses tidyr to  │  │              │
│ Data objects:   │  │   separate "_"   │  │              │
│  ncomms8866_*   │  │   subunits       │  │              │
└─────────────────┘  └──────────────────┘  └──────────────┘
         │                   │                   │
         └───────────────────┼───────────────────┘
                             │
                             ▼
                   ┌──────────────────────┐
                   │ FilterGroundTruth()  │
                   │                      │
                   │ For each mechanism:  │
                   │  • Check if ALL      │
                   │    ligand subunits   │
                   │    exist in input    │
                   │  • Check if ALL      │
                   │    receptor subunits │
                   │    exist in input    │
                   │  • Keep only valid   │
                   │    mechanisms        │
                   │  • Remove empty      │
                   │    subunit columns   │
                   └──────────────────────┘
```

### Ground Truth Data Format

```
ground.truth = list(
  source.subunits = matrix(
    # Rows: L-R mechanisms (e.g., "TGFB1—TGFBR1", "IL6—IL6R_IL6ST")
    # Cols: subunit_1, subunit_2, ..., subunit_n
    # Values: gene names or NA

    Example row for single-subunit:
      subunit_1 = "TGFB1", subunit_2 = NA, ...

    Example row for multi-subunit:
      subunit_1 = "IL6R", subunit_2 = "IL6ST", subunit_3 = NA, ...
  ),

  target.subunits = matrix(
    # Same structure as source.subunits but for receptors
  )
)
```

---

## 5. SPATIAL ANALYSIS SYSTEM

```
┌─────────────────────────────────────────────────────────────┐
│                    ComputeEdgelist()                        │
│                                                             │
│  Required for:                                              │
│    • CellToCellSpatial, CellToNeighborhood,                 │
│      NeighborhoodToCell                                     │
│                                                             │
│  Inputs:                                                    │
│    • node.object: Seurat with spatial coordinates          │
│    • position.x, position.y: metadata columns              │
│    • k: number of neighbors (default 4)                    │
│    • rad.set: radius threshold (alternative to k)          │
│    • nn.method: "aoz" (fast) or NULL (legacy)              │
│                                                             │
│  Methods:                                                   │
│  ┌────────────────────────────────────────────────────┐    │
│  │ FAST (nn.method='aoz', DEFAULT)                    │    │
│  │  • Uses spNNGP::mkNNIndxCB()                       │    │
│  │  • Mutual nearest neighbors via spatial ordering   │    │
│  │  • Efficient for large datasets                    │    │
│  │  • Returns directed edges                          │    │
│  └────────────────────────────────────────────────────┘    │
│                                                             │
│  ┌────────────────────────────────────────────────────┐    │
│  │ LEGACY (nn.method=NULL)                            │    │
│  │  • Computes full distance matrix                   │    │
│  │  • Filters to k-NN or radius                       │    │
│  │  • Converts to igraph object                       │    │
│  │  • Extracts mutual edges                           │    │
│  └────────────────────────────────────────────────────┘    │
│                                                             │
│  Output:                                                    │
│    dataframe(from=cell_barcode, to=cell_barcode)           │
│    • Each row is a directed spatial edge                   │
│    • Mutual nearest neighbors only                         │
└─────────────────────────────────────────────────────────────┘
```

---

## 6. VISUALIZATION SYSTEM

### 6A. CircuitPlot Workflow

```
CircuitPlot()
  │
  ├─ Inputs:
  │    • transcr.obj: Seurat (node/gene expression)
  │    • connect.obj: Seurat (edge/interaction data)
  │    • feature: "LIGAND|RECEPTOR" (mechanism to plot)
  │    • group.by: metadata column for aggregation
  │    • [styling parameters: colors, sizes, offsets]
  │
  ├─ Step 1: Extract node data
  │    └─ DefineNodeObject(transcr.obj, feature, assay='RNA')
  │         • Returns: dataframe with expression + metadata
  │
  ├─ Step 2: Extract edge data
  │    └─ DefineEdgeObject(connect.obj, feature, assay='CellToCell')
  │         • Returns: dataframe with connectivity + metadata
  │
  ├─ Step 3: Aggregate node data
  │    └─ AggregateNodeData(node.object, group.by, global.node.list)
  │         • Groups by cell type
  │         • Computes system.fraction (cell type proportions)
  │         • Returns: node.aggregate
  │
  ├─ Step 4: Aggregate edge data
  │    └─ AggregateEdgeData(edge.object, group.by)
  │         • Groups by sending/receiving types
  │         • Computes mean connectivity
  │         • Parses sending/receiving labels
  │         • Returns: edge.aggregate
  │
  └─ Step 5: Create plot
       └─ ggCircuit(edge.aggregate, node.aggregate, [styling params])
            │
            ├─ netCoin::layoutCircle() → node positions on circle
            │
            ├─ FOR each edge:
            │    • Compute start/end coordinates (trigonometry)
            │    • Offset paracrine edges to avoid overlap
            │    • Create curved loops for autocrine signals
            │    • Color/thickness based on connectivity
            │
            └─ ggplot() + geom_segment() + geom_curve() + geom_point()
                 • Nodes: circles sized by system.fraction
                 • Edges: arrows colored by connectivity strength
                 • Autocrine: curved loops
```

### 6B. NetworkPlot Workflow

```
NetworkPlot()
  │
  ├─ Inputs:
  │    • transcriptome.object: Seurat with UMAP/reduction
  │    • connectome.object: Seurat (CellToCell or spatial)
  │    • mechanism.of.interest: "LIGAND|RECEPTOR"
  │    • reduction: UMAP, tSNE, etc.
  │    • [filtering, styling parameters]
  │
  ├─ Step 1: Parse mechanism
  │    └─ strsplit(mechanism, "|") → ligand, receptor
  │
  ├─ Step 2: Extract coordinates
  │    └─ transcriptome.object@reductions[[reduction]]@cell.embeddings
  │
  ├─ Step 3: Extract connectivity data
  │    └─ connectome.object assay for mechanism
  │         • Filter by min.connectivity.thresh
  │         • Log-transform if needed
  │
  ├─ Step 4: Extract ligand/receptor expression
  │    └─ transcriptome.object assay for ligand and receptor genes
  │
  ├─ Step 5: Parse barcode pairs
  │    └─ stringr::str_split_fixed() → sender, receiver barcodes
  │
  ├─ Step 6: Map coordinates
  │    └─ Match sender/receiver to UMAP positions
  │         • Create segment dataframe (x, y, xend, yend)
  │
  └─ Step 7: Create plot
       └─ ggplot() +
            geom_segment(connectivity, alpha/color by strength) +  # Edges
            geom_point(UMAP positions) +                           # Cells
            [OPTIONAL] ggnewscale::new_scale_color() +             # Multiple scales
            geom_point(ligand expression, color scale) +           # Ligand overlay
            geom_point(receptor expression, color scale) +         # Receptor overlay
            [OPTIONAL] facet_wrap(split.by) +                      # Split views
            [OPTIONAL] cowplot::plot_grid(main + legends)          # Add legends
```

### 6C. DotPlotSCC Workflow

```
DotPlotSCC()
  │
  ├─ Inputs:
  │    • object: Seurat (edge/CellToCell data)
  │    • features: character vector of mechanisms
  │    • group.by: metadata column (e.g., VectorType)
  │    • [styling: colors, dot size scaling]
  │
  ├─ Step 1: Fetch data
  │    └─ Seurat::FetchData(object, features)
  │
  ├─ Step 2: Group cells by identity
  │    └─ Seurat::CellsByIdentities() or group.by
  │
  ├─ Step 3: Compute statistics per group
  │    └─ FOR each group × feature:
  │         • percent.expressed: % cells > 0
  │         • avg.expression: mean expression
  │         • avg.expression.scaled: z-score
  │
  ├─ Step 4: [OPTIONAL] Hierarchical clustering
  │    └─ IF cluster.idents:
  │         • stats::dist() → distance matrix
  │         • stats::hclust() → dendrogram
  │         • Reorder identities by cluster
  │
  └─ Step 5: Create dot plot
       └─ ggplot() +
            geom_point(
              x = feature,
              y = identity,
              size = percent.expressed,
              color = avg.expression.scaled
            ) +
            scale_color_gradient(low=cols[1], high=cols[2]) +
            [OPTIONAL] facet_grid(split.by)
```

---

## 7. HELPER UTILITIES

### 7A. Version Compatibility Layer

```
┌──────────────────────────────────────────────────────────────┐
│                     GetSeuratAssay()                         │
│                                                              │
│  Purpose: Abstract Seurat v4 vs v5 API differences          │
│                                                              │
│  Inputs:                                                     │
│    • object: Seurat object                                  │
│    • assay: assay name                                      │
│    • slot: "data", "counts", or "scale.data"                │
│                                                              │
│  Logic:                                                      │
│    IF SeuratObject::Version() >= 5.0:                       │
│      → Seurat::GetAssayData(object, assay, layer=slot)      │
│    ELSE:                                                     │
│      → Seurat::GetAssayData(object, assay, slot=slot)       │
│                                                              │
│  Returns: Matrix or sparse Matrix                           │
│                                                              │
│  Called by: ALL Run* functions (multiple times each)        │
└──────────────────────────────────────────────────────────────┘
```

### 7B. Preprocessing & Filtering

```
┌──────────────────────────────────────────────────────────────┐
│                    PrepNodeObject()                          │
│                                                              │
│  Purpose: Filter input before transformation                 │
│                                                              │
│  Inputs:                                                     │
│    • node.object: Seurat object                             │
│    • assay: assay to filter                                 │
│    • min.cells.per.ident: filter rare cell types            │
│    • min.cells.per.gene: filter lowly expressed genes       │
│                                                              │
│  Operations:                                                 │
│    1. Set default assay                                     │
│    2. Count cells per identity → Seurat::Idents()           │
│    3. Filter identities with < min.cells.per.ident          │
│    4. Count non-zero cells per gene → Matrix::rowSums()     │
│    5. Filter genes with < min.cells.per.gene                │
│    6. Subset object → Seurat::subset()                      │
│                                                              │
│  Returns: Filtered Seurat object                            │
│                                                              │
│  Called by: RunNICHES.default (once at beginning)           │
└──────────────────────────────────────────────────────────────┘

┌──────────────────────────────────────────────────────────────┐
│                   FilterGroundTruth()                        │
│                                                              │
│  Purpose: Remove L-R mechanisms not in dataset               │
│                                                              │
│  Inputs:                                                     │
│    • ground.truth: list(source.subunits, target.subunits)   │
│    • input_rownames: genes in dataset                       │
│                                                              │
│  Operations:                                                 │
│    1. FOR each mechanism (row):                             │
│         • Check if ALL ligand subunits in input_rownames    │
│         • Check if ALL receptor subunits in input_rownames  │
│         • Keep only if both conditions met                  │
│    2. Remove subunit columns that are all NA                │
│                                                              │
│  Returns: Filtered ground.truth (same structure)            │
│                                                              │
│  Called by: LoadGroundTruth (after loading database)        │
└──────────────────────────────────────────────────────────────┘
```

### 7C. Cell Type & Metadata Utilities

```
┌──────────────────────────────────────────────────────────────┐
│                   ReturnCelltypes()                          │
│  (also lowercase: return_celltypes)                          │
│                                                              │
│  Purpose: Extract unique cell types from Seurat object       │
│                                                              │
│  Inputs:                                                     │
│    • seurat_object: Seurat object                           │
│                                                              │
│  Operations:                                                 │
│    • Seurat::Idents(seurat_object) → factor                 │
│    • table() → count per identity                           │
│    • names() → unique identities                            │
│    • Print warning about active Idents requirement          │
│                                                              │
│  Returns: Character vector of unique cell types              │
│                                                              │
│  Called by: RunCellToCell, RunCellToCellSpatial             │
└──────────────────────────────────────────────────────────────┘

┌──────────────────────────────────────────────────────────────┐
│                     AddMetaData()                            │
│                 (Experimental/Incomplete)                    │
│                                                              │
│  Purpose: Add metadata from annotated object to edge object  │
│                                                              │
│  Inputs:                                                     │
│    • CTC.obj: Seurat (edge data)                            │
│    • annotated.obj: Seurat (node data with new metadata)    │
│    • col.name: metadata column to add                       │
│                                                              │
│  Operations:                                                 │
│    1. Parse barcode pairs from CTC.obj                      │
│    2. Join annotated.obj metadata by sender barcode         │
│    3. Join annotated.obj metadata by receiver barcode       │
│    4. Create .Sending, .Receiving, .Joint versions          │
│                                                              │
│  Returns: Updated metadata dataframe                         │
│                                                              │
│  Note: Appears to be development version; contains inline   │
│        example code at bottom of file                       │
└──────────────────────────────────────────────────────────────┘
```

---

## 8. PARAMETER STANDARDIZATION ANALYSIS

### 8A. Input Object Parameters (INCONSISTENT)

| Function | Object Parameter | Type Expected |
|----------|-----------------|---------------|
| RunNICHES.Seurat | `object` | Seurat |
| RunNICHES.matrix | `object` | matrix |
| RunNICHES.default | `object` | Seurat |
| PrepNodeObject | `node.object` | Seurat |
| **RunCellToCell** | **`node.object`** | **Seurat** |
| **RunCellToSystem** | **`node.object`** | **Seurat** |
| **RunSystemToCell** | **`node.object`** | **Seurat** |
| **RunCellToCellSpatial** | **`node.object`** | **Seurat** |
| **RunCellToNeighborhood** | **`node.object`** | **Seurat** |
| **RunNeighborhoodToCell** | **`node.object`** | **Seurat** |
| ComputeEdgelist | `node.object` | Seurat |
| CircuitPlot | `transcr.obj` | Seurat |
| CircuitPlot | `connect.obj` | Seurat |
| NetworkPlot | `transcriptome.object` | Seurat |
| NetworkPlot | `connectome.object` | Seurat |
| DotPlotSCC | `object` | Seurat |
| GetSeuratAssay | `object` | Seurat |
| ReturnCelltypes | `seurat_object` | Seurat |

**RECOMMENDATION:** Standardize to `object` for primary input across all functions.

---

### 8B. Spatial Parameters (CONSISTENT)

| Parameter | Used In | Default |
|-----------|---------|---------|
| `position.x` | RunNICHES, ComputeEdgelist, create_Neighborhood_Assay | NULL or 'x' |
| `position.y` | RunNICHES, ComputeEdgelist, create_Neighborhood_Assay | NULL or 'y' |
| `k` | RunNICHES, ComputeEdgelist, create_Neighborhood_Assay | 4 or NULL |
| `rad.set` | RunNICHES, ComputeEdgelist, create_Neighborhood_Assay | NULL or 25 |

**RECOMMENDATION:** Good consistency. Consider consolidating defaults (use k=4 everywhere).

---

### 8C. Assay Parameters (MOSTLY CONSISTENT)

| Parameter | Used In | Default |
|-----------|---------|---------|
| `assay` | All Run* functions, PrepNodeObject, GetSeuratAssay | "RNA" (in RunNICHES) |
| `node.assay` | NetworkPlot | 'RNA' |
| `edge.assay` | NetworkPlot | 'CellToCell' |

**RECOMMENDATION:** `assay` is good standard. NetworkPlot should align naming.

---

### 8D. Ground Truth Parameters (CONSISTENT)

| Parameter | Used In | Type |
|-----------|---------|------|
| `ground.truth` | All Run* functions | list(source.subunits, target.subunits) |
| `LR.database` | RunNICHES, LoadGroundTruth | string: "fantom5", "omnipath", "custom" |
| `species` | RunNICHES, LoadGroundTruth, Load* functions | string: "human", "mouse", "rat", "pig" |
| `custom_LR_database` | RunNICHES, LoadGroundTruth, LoadCustom | dataframe |

**RECOMMENDATION:** Excellent consistency.

---

### 8E. Aggregation Parameters (INCONSISTENT)

| Parameter | Used In | Options | Default |
|-----------|---------|---------|---------|
| `blend` | RunCellToSystem, RunSystemToCell | "sum", "mean", "mean.adj" | "mean" |
| `blend` | RunCellToNeighborhood, RunNeighborhoodToCell | "sum", "mean" | "mean" |

**RECOMMENDATION:** Harmonize options. Should spatial functions also support "mean.adj"?

---

### 8F. Metadata Parameters (CONSISTENT)

| Parameter | Used In | Type |
|-----------|---------|------|
| `meta.data.to.map` | All Run* functions | character vector |
| `group.by` | CircuitPlot, DotPlotSCC | string (single column) |
| `split.by` | NetworkPlot, DotPlotSCC | string (single column) |

**RECOMMENDATION:** Good naming. `meta.data.to.map` is clear.

---

### 8G. Output Parameters (CONSISTENT)

| Parameter | Used In | Options | Default |
|-----------|---------|---------|---------|
| `output_format` | All Run* functions | "seurat", "raw" | "seurat" |

**RECOMMENDATION:** Perfect consistency.

---

## 9. CODE DUPLICATION ANALYSIS

### 9A. HIGH DUPLICATION: Seurat Object Creation

**Duplicated ~6 times** in RunCellToCell, RunCellToSystem, RunSystemToCell, RunCellToCellSpatial, RunCellToNeighborhood, RunNeighborhoodToCell

```r
# Pattern repeated in every Run* function:
output.seurat <- CreateSeuratObject(counts = interaction_matrix,
                                     project = 'NICHES',
                                     assay = assay_name)
output.seurat <- NormalizeData(output.seurat)

# Seurat v5 compatibility (manual data layer creation)
if (SeuratObject::Version() >= '5.0.0') {
  output.seurat@assays[[assay_name]]@layers$data <-
    output.seurat@assays[[assay_name]]@layers$counts
}

Idents(output.seurat) <- output.seurat$some_metadata_column
output.seurat <- AddMetaData(output.seurat, metadata_to_add)
```

**RECOMMENDATION:** Create helper function:
```r
CreateNichesSeurat <- function(counts, assay_name, metadata, ident_column) {
  # Consolidate Seurat object creation + v5 compatibility + metadata
}
```

---

### 9B. HIGH DUPLICATION: Ligand/Receptor Matrix Construction

**Duplicated ~6 times** - Core algorithm for multiplying subunit expressions

```r
# Pattern for ligand matrix:
ligand.matrix <- matrix(1, nrow=nrow(source.subunits), ncol=ncol(expression))
for (subunit_col in 1:ncol(source.subunits)) {
  subunit_genes <- source.subunits[, subunit_col]
  subunit_genes <- subunit_genes[!is.na(subunit_genes)]
  ligand.matrix <- ligand.matrix * expression[subunit_genes, , drop=FALSE]
}
rownames(ligand.matrix) <- mechanism_names

# Similar for receptor matrix
```

**RECOMMENDATION:** Create helper function:
```r
ComputeComplexExpression <- function(subunit.matrix, expression.matrix) {
  # Returns mechanism × cell matrix
}
```

---

### 9C. MEDIUM DUPLICATION: Metadata Propagation

**Duplicated ~6 times** - Logic for creating .Sending, .Receiving, .Joint metadata columns

```r
# Pattern:
for (meta_col in meta.data.to.map) {
  metadata_sending <- node.object@meta.data[sending_cells, meta_col]
  metadata_receiving <- node.object@meta.data[receiving_cells, meta_col]
  metadata_joint <- paste(metadata_sending, metadata_receiving, sep=" - ")

  output.metadata[[paste0(meta_col, ".Sending")]] <- metadata_sending
  output.metadata[[paste0(meta_col, ".Receiving")]] <- metadata_receiving
  output.metadata[[paste0(meta_col, ".Joint")]] <- metadata_joint
}
```

**RECOMMENDATION:** Create helper function:
```r
PropagateMetadata <- function(node.object, sending_cells, receiving_cells,
                               meta.data.to.map, direction = "both") {
  # direction: "both", "sending", "receiving"
  # Returns dataframe with propagated columns
}
```

---

### 9D. MEDIUM DUPLICATION: Aggregation Logic

**Duplicated ~4 times** in RunCellToSystem, RunSystemToCell, RunCellToNeighborhood, RunNeighborhoodToCell

```r
# Pattern:
if (blend == "sum") {
  aggregated <- Matrix::rowSums(expression_matrix)
} else if (blend == "mean") {
  aggregated <- Matrix::rowMeans(expression_matrix)
} else if (blend == "mean.adj") {
  aggregated <- Matrix::rowMeans(expression_matrix) * sqrt(ncol(expression_matrix))
}
```

**RECOMMENDATION:** Create helper function:
```r
AggregateExpression <- function(expression.matrix, method = "mean") {
  # method: "sum", "mean", "mean.adj"
  # Returns aggregated vector
}
```

---

### 9E. LOW DUPLICATION: Edge Parsing

**Duplicated ~2 times** in NetworkPlot, AddMetaData, visualization helpers

```r
# Pattern:
barcode_pairs <- stringr::str_split_fixed(colnames(edge.object), "—", 2)
sender_barcodes <- barcode_pairs[, 1]
receiver_barcodes <- barcode_pairs[, 2]
```

**RECOMMENDATION:** Create helper function:
```r
ParseBarcodePairs <- function(barcode.pairs, separator = "—") {
  # Returns list(sender, receiver)
}
```

---

## 10. REFACTORING ROADMAP

### Phase 1: Extract Core Utilities (NO API CHANGES)

**Priority: High (reduce duplication immediately)**

1. **CreateNichesSeurat()** - Consolidate Seurat object creation (6 places)
2. **ComputeComplexExpression()** - Ligand/receptor matrix construction (6 places)
3. **PropagateMetadata()** - Metadata propagation (6 places)
4. **AggregateExpression()** - Sum/mean/mean.adj logic (4 places)
5. **ParseBarcodePairs()** - Edge barcode parsing (2-3 places)

**Impact:** Reduce code by ~30-40%, improve maintainability, fix bugs in one place

---

### Phase 2: Standardize Parameter Names (BREAKING CHANGES)

**Priority: Medium (improves user experience)**

1. **Rename all object parameters to `object`**
   - `node.object` → `object`
   - `transcr.obj` → `node.object` or `transcriptome.object`
   - `connect.obj` → `edge.object` or `connectome.object`
   - `seurat_object` → `object`

2. **Harmonize blend options**
   - Add "mean.adj" to spatial functions
   - Document when each method is appropriate

3. **Consolidate spatial defaults**
   - Use k=4 everywhere (or make explicit in docs)
   - Clarify rad.set vs k priority

**Impact:** Cleaner API, easier to remember, better docs

---

### Phase 3: Modularize Transformation Pipeline (BREAKING CHANGES)

**Priority: Low (architectural improvement, major refactor)**

1. **Extract pre-processing stage**
   ```r
   PrepareNichesInput(object, assay, min.cells.per.ident, min.cells.per.gene)
   ```

2. **Extract transformation stage**
   ```r
   ComputeInteractions(ligand.expr, receptor.expr, ground.truth)
   ```

3. **Extract post-processing stage**
   ```r
   BuildNichesOutput(interaction.matrix, metadata, output_format)
   ```

4. **Refactor Run* functions to use pipeline**
   - Each Run* function becomes a thin wrapper
   - Core logic extracted to shared functions

**Impact:** Maximum flexibility, easier testing, enable custom pipelines

---

## 11. QUICK REFERENCE: FUNCTION INDEX

### By Category

**Entry Points (User-facing):**
- RunNICHES (generic + 4 S3 methods)

**Core Transformations:**
- RunCellToCell
- RunCellToSystem
- RunSystemToCell
- RunCellToCellSpatial
- RunCellToNeighborhood
- RunNeighborhoodToCell

**Ground Truth:**
- LoadGroundTruth (dispatcher)
- LoadFantom5
- LoadOmniPath
- LoadCustom
- FilterGroundTruth

**Preprocessing:**
- PrepNodeObject
- ComputeEdgelist

**Visualization:**
- CircuitPlot
- NetworkPlot
- DotPlotSCC
- ggCircuit (internal)

**Data Helpers:**
- DefineNodeObject
- DefineEdgeObject
- AggregateNodeData
- AggregateEdgeData
- create_Neighborhood_Assay

**Utilities:**
- GetSeuratAssay
- ReturnCelltypes
- AddMetaData
- gg_color_hue

---

## 12. EXTERNAL DEPENDENCIES SUMMARY

**Critical (cannot run without):**
- Seurat (>= 4.0, < 6.0)
- SeuratObject
- Matrix
- dplyr
- tidyr

**Important (some features unavailable):**
- spNNGP (fast spatial neighbors)
- igraph (legacy spatial neighbors)
- OmnipathR (OmniPath database)

**Visualization:**
- ggplot2, cowplot, ggnewscale
- RColorBrewer, netCoin

**Utility:**
- abind, stringr, data.table, scales

---

## END OF ARCHITECTURE DIAGRAM
