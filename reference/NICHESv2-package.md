# NICHESv2: Spatial Cell-to-Cell Ligand-Receptor Communication Analysis

NICHESv2 computes ligand-receptor (LR) signaling scores between pairs of
cells in spatial and non-spatial single-cell transcriptomics data. It
improves on the original NICHES package with:

- A memory-efficient long-format output object (`NICHESObject`, 15
  slots).

- Support for multi-gene **ligand and receptor complexes**
  (underscore-delimited, e.g. `"COL1A1_COL1A2"`).

- Built-in loaders for the **ConnectomeDB2025** (default), **FANTOM5**,
  and **OmniPath** LR databases.

- Per-edge and per-cell **normalized scores** alongside raw counts.

- Parallel LR computation via socket clusters on all platforms.

- Flexible aggregation (`aggregate_NICHESObject`).

## Typical workflow


    library(NICHESv2)

    # 1. Load / prepare inputs
    #    count.mtx : gene x cell sparse matrix  (rownames = genes, colnames = cells)
    #    meta.data : data.frame                  (rownames = cells; spatial mode
    #                                             needs "x" and "y" columns)
    #    LRM.db    : "connectomedb2025" | "fantom5" | "omnipath" | user data.frame

    # 2. Build the NICHESObject (one call does everything)
    niches.obj <- create_NICHESObject(
      count.mtx = count.mtx,
      meta.data = meta.data,
      LRM.db    = "connectomedb2025",
      k         = 6L,           # 6-nearest-neighbor graph
      method    = "product"
    )

    # 3. Inspect
    print(niches.obj)
    summary(niches.obj)
    dim(niches.obj)            # cells / genes / edges / LRMs

    # 4. Aggregate neighborhood edge signals
    niches.obj <- aggregate_NICHESObject(niches.obj, cell.type.col = "celltype")

    # 5. Subset to a cell population (aggregation rows carry over unchanged)
    b.cells <- niches.obj$cell.meta[celltype == "B_cell", cell]
    sub.obj  <- niches.obj[b.cells]

## NICHESObject slot reference

|  |  |  |
|----|----|----|
| **Slot** | **Type** | **Contents** |
| `$cell.list` | character | Cell barcodes |
| `$cell.data` | data.table | (cell, gene, score, score.norm) – sparse |
| `$cell.meta` | data.table | Per-cell metadata (x, y, celltype, ...) |
| `$edge.list` | data.table | (edge, SendingCell, ReceivingCell) |
| `$edge.data` | data.table | (edge, LRM, score, score.norm) – sparse |
| `$edge.meta` | data.table | Per-edge metadata with Sending/Receiving/Joint cols |
| `$neighborhood.cell.list` | list | Per-cell neighbor barcodes |
| `$neighborhood.edge.list` | list | Per-cell neighborhood edge IDs |
| `$neighborhood.meta` | data.table | (cell, n.neighbors, n.edges) |
| `$gene.list` | character | Gene names |
| `$gene.meta` | data.table | Per-gene annotation |
| `$LRM.list` | character | "ligand\|receptor" strings |
| `$LRM.meta` | data.table | Per-LRM annotation (pathway, class, ...) |
| `$segmentation` | sf or NULL | Optional cell polygons |
| `$aggregations` | list | Aggregation tables computed on demand |

## Key separators

- Edge IDs: `"SendingCell|ReceivingCell"`

- LRM IDs: `"ligand|receptor"`

- Gene complexes in database: `"gene1_gene2"` (underscore)

## Key naming conventions

- Function names use underscores:
  [`create_NICHESObject()`](https://raredonlab.github.io/NICHESv2/reference/create_NICHESObject.md)

- Variable / argument names use dots: `n.cores`, `bg.color`

- Score columns: `score` (raw), `score.norm` (normalized)

## References

Ramilowski, J.A. et al. (2015). A draft network of
ligand-receptor-mediated multicellular signalling in human. *Nature
Communications*, 6, 7866.

Turei, D. et al. (2021). Integrated intra- and intercellular signaling
knowledge for multicellular omics analysis. *Molecular Systems Biology*,
17(3), e9923.

## See also

[`create_NICHESObject`](https://raredonlab.github.io/NICHESv2/reference/create_NICHESObject.md),
[`aggregate_NICHESObject`](https://raredonlab.github.io/NICHESv2/reference/aggregate_NICHESObject.md),
[`load_LRM_database`](https://raredonlab.github.io/NICHESv2/reference/load_LRM_database.md)

## Author

NICHESv2 Authors
