# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

NICHES (Niche Interactions and Communication Heterogeneity in Extracellular Signaling) is an R package for analyzing cell-cell signaling at single-cell resolution. It transforms single-cell atlases into single-cell-signaling atlases by creating unique one-to-one cell pairings and characterizing cellular interactions across ligand-receptor signaling mechanisms.

**Key Concept**: NICHES transforms cell-level expression data into edge-level signaling data, where each "edge" represents a sender-receiver cell pair. This enables dimensional reduction and analysis of cell-cell interactions in "signal-space" rather than gene-expression space.

## Development Commands

### Building and Testing
```r
# Build package documentation
devtools::document()

# Install the package locally
devtools::install()

# Build package from source
R CMD build .

# Check package (runs R CMD check)
devtools::check()

# Load package for interactive development
devtools::load_all()
```

### Installation
Users install via:
```r
library(devtools)
install_github('msraredon/NICHES', ref = 'master')
```

## Architecture

### Core Transformation Pipeline

The package implements a modular transformation pipeline through `RunNICHES()`:

1. **Input Processing** (`PrepNodeObject.R`): Filters Seurat objects by minimum cells per identity and minimum gene expression thresholds
2. **Ground Truth Loading** (`LoadGroundTruth.R`): Loads ligand-receptor databases (FANTOM5, OmniPath, or custom)
3. **Edgelist Computation** (`ComputeEdgelist.R`): For spatial data, computes neighbor relationships via k-NN or radius-based methods
4. **Organization Functions** (`RunCellToCell.R`, `RunSystemToCell.R`, etc.): Generate different views of cell signaling

### Transformation Modes

NICHES supports six analysis modes, controlled by boolean flags in `RunNICHES()`:

**Non-spatial modes** (all-to-all connectivity):
- `CellToCell`: Pairwise cell-cell interactions with downsampling to manage combinatorial explosion
- `CellToSystem`: Summed signaling output from each cell to entire system
- `SystemToCell`: Summed signaling input to each cell from entire system (microenvironment)

**Spatial modes** (neighbor-constrained):
- `CellToCellSpatial`: Pairwise interactions between spatial neighbors only
- `CellToNeighborhood`: Summed signaling output to spatial neighbors
- `NeighborhoodToCell`: Summed signaling input from spatial neighbors (spatial niche)

### S3 Method Dispatch

`RunNICHES()` uses S3 dispatch with methods for:
- `RunNICHES.Seurat`: Takes Seurat objects (most common)
- `RunNICHES.matrix`: Takes expression matrices with metadata dataframes
- `RunNICHES.Matrix`: Handles sparse Matrix objects
- `RunNICHES.default`: Core implementation called by all methods

The dispatch pattern: Seurat/matrix methods validate inputs → call `RunNICHES.default()` → default calls organization functions (`RunCellToCell`, etc.) → returns list of Seurat objects or raw matrices

### Ligand-Receptor Databases

Three database options:
- **FANTOM5** (`LoadFantom5.R`): Simple ligand-receptor pairs stored in `data/*.rda`. Supports human, mouse, rat, pig
- **OmniPath** (`LoadOmniPath.R`): Multi-subunit complexes (ligands/receptors with multiple protein subunits separated by `_`). Downloads from OmnipathR. Supports human, mouse, rat
- **Custom** (`LoadCustom.R`): User-provided dataframe with source and target gene columns

Ground truth format: List with `source.subunits` and `target.subunits` matrices where rows are mechanisms and columns are subunits.

### Cell-Cell Interaction Matrix Construction

Core algorithm in `RunCellToCell()`:

1. For each celltype, compute ligand expression matrix (multiply expression of all ligand subunits per mechanism)
2. For each celltype, compute receptor expression matrix (multiply expression of all receptor subunits per mechanism)
3. For each sending celltype, sample edges to all receiving celltypes (downsampling to min of sender/receiver population size)
4. Create sender-receiver pairs: `ligand_expression * receptor_expression` per mechanism
5. Build Seurat object where:
   - Rows = ligand-receptor mechanisms (e.g., `TGFB1—TGFBR1`)
   - Columns = barcode pairs (e.g., `AAACCTGAGAAGGCCT-1—AAACCTGAGCAGATCA-1`)
   - Values = interaction strength (product of ligand and receptor expression)
   - Metadata includes `SendingType`, `ReceivingType`, `VectorType` (concatenation of sender-receiver types)

### Spatial Edge Computation

`ComputeEdgelist()` supports two methods:
- **k-NN based** (default, `nn.method='aoz'`): Uses spNNGP package for fast mutual nearest neighbor computation. Efficient for large datasets
- **Legacy** (`nn.method=NULL`): Computes full distance matrix then filters. Slower but more straightforward

Spatial coordinates must be in metadata columns specified by `position.x` and `position.y`.

## Important Implementation Details

### Seurat Version Compatibility

The package handles both Seurat v4 and v5. Key difference: v5 uses layered assays and requires manual creation of the `data` slot after `NormalizeData()` (see `RunCellToCell.R`:123-126).

### Cell Type Handling

- For Seurat objects: `cell_types` parameter specifies which metadata column contains cell type labels. If NULL, uses `Idents(object)`
- The active identity (`Idents`) is used for subsetting and must contain cell types for `CellToCell` mode
- Cell types are validated and set as identities before processing

### Metadata Propagation

`meta.data.to.map` parameter: Vector of metadata column names to carry forward into NICHES objects. For each metadata field, three versions are created:
- `.Sending`: Value from sending cell
- `.Receiving`: Value from receiving cell
- `.Joint`: Concatenation as "sending - receiving"

### Output Formats

`output_format` parameter controls return type:
- `"seurat"` (default): Returns list of Seurat objects, one per enabled organization
- `"raw"`: Returns list of lists with `$Matrix` (interaction matrix) and `$metadata` (compiled metadata dataframe)

## File Organization

- `R/RunNICHES.R`: Main entry point and S3 dispatch
- `R/RunCellToCell.R`, `R/RunSystemToCell.R`, etc.: Individual organization implementations
- `R/ComputeEdgelist.R`: Spatial neighbor graph construction
- `R/PrepNodeObject.R`: Input filtering and validation
- `R/Load*.R`: Database loading functions
- `data/`: Pre-packaged FANTOM5 ligand-receptor databases
- `man/`: Roxygen-generated documentation
- `vignettes/`: Usage examples organized by use case (spatial, differential analysis, pseudotime)
- `tests/`: Ad-hoc test scripts (not formal unit tests)

## Notes on Code Quality

- Extensive TODOs in `RunNICHES.R` indicate ongoing refactoring
- Parameter validation is partially implemented, with plans for more rigorous checking
- The package uses older Seurat patterns; v5 compatibility is being added incrementally
- Global variables are managed through `globals.R` for R CMD check compliance
