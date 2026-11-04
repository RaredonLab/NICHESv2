# Package index

## Build an object

- [`create_NICHESObject()`](https://raredonlab.github.io/NICHESv2/reference/create_NICHESObject.md)
  : Create a NICHESObject from Single-Cell Data
- [`extract_NICHESInputs_Seurat()`](https://raredonlab.github.io/NICHESv2/reference/extract_NICHESInputs_Seurat.md)
  : Extract inputs for create_NICHESObject from a Seurat object
- [`extract_NICHESInputs_AnnData()`](https://raredonlab.github.io/NICHESv2/reference/extract_NICHESInputs_AnnData.md)
  : Extract inputs for create_NICHESObject from an AnnData (.h5ad) file
- [`load_LRM_database()`](https://raredonlab.github.io/NICHESv2/reference/load_LRM_database.md)
  : Load a Ligand-Receptor Pair Database
- [`merge_NICHESObjects()`](https://raredonlab.github.io/NICHESv2/reference/merge_NICHESObjects.md)
  : Merge multiple NICHESObjects into one

## Neighborhoods and aggregation

- [`Explore_Radius_Neighborhood()`](https://raredonlab.github.io/NICHESv2/reference/Explore_Radius_Neighborhood.md)
  : Explore the relationship between search radius and neighborhood size
- [`add_neighborhoods()`](https://raredonlab.github.io/NICHESv2/reference/add_neighborhoods.md)
  : Add Neighborhood Slots to a NICHESObject
- [`extract_neighborhoods()`](https://raredonlab.github.io/NICHESv2/reference/extract_neighborhoods.md)
  : Extract Neighborhoods from an Edge List
- [`aggregate_NICHESObject()`](https://raredonlab.github.io/NICHESv2/reference/aggregate_NICHESObject.md)
  : Run the Full Aggregation Pipeline on a NICHESObject
- [`add_neighborhood_edge_agg()`](https://raredonlab.github.io/NICHESv2/reference/add_neighborhood_edge_agg.md)
  : Add Unified Neighborhood Edge Aggregation to a NICHESObject
- [`add_neighborhood_composition()`](https://raredonlab.github.io/NICHESv2/reference/add_neighborhood_composition.md)
  : Add Neighborhood Cell-Type Composition to a NICHESObject

## Work with an object

- [`` `[`( ``*`<NICHESObject>`*`)`](https://raredonlab.github.io/NICHESv2/reference/sub-.NICHESObject.md)
  : Subset a NICHESObject to a Set of Cells

- [`print(`*`<NICHESObject>`*`)`](https://raredonlab.github.io/NICHESv2/reference/print.NICHESObject.md)
  : Print a NICHESObject

- [`summary(`*`<NICHESObject>`*`)`](https://raredonlab.github.io/NICHESv2/reference/summary.NICHESObject.md)
  : Summarize a NICHESObject

- [`dim(`*`<NICHESObject>`*`)`](https://raredonlab.github.io/NICHESv2/reference/dim.NICHESObject.md)
  : Dimensions of a NICHESObject

- [`length(`*`<NICHESObject>`*`)`](https://raredonlab.github.io/NICHESv2/reference/length.NICHESObject.md)
  : Number of Cells in a NICHESObject

- [`is_NICHESObject()`](https://raredonlab.github.io/NICHESv2/reference/is_NICHESObject.md)
  : Test Whether an Object Is a NICHESObject

- [`add_meta_column()`](https://raredonlab.github.io/NICHESv2/reference/add_meta_column.md)
  :

  Append or Overwrite a Single Column in `$cell.meta`

- [`add_cell_meta()`](https://raredonlab.github.io/NICHESv2/reference/add_cell_meta.md)
  :

  Replace `$cell.meta` on a NICHESObject

- [`add_cell_rawcounts()`](https://raredonlab.github.io/NICHESv2/reference/add_cell_rawcounts.md)
  :

  Add Raw UMI Counts to `$cell.data`

- [`add_segmentation()`](https://raredonlab.github.io/NICHESv2/reference/add_segmentation.md)
  : Attach Cell Segmentation Polygons to a NICHESObject

## Export

- [`export_to_TissuePlex()`](https://raredonlab.github.io/NICHESv2/reference/export_to_TissuePlex.md)
  : Export a NICHESObject to TissuePlex edges.parquet

## Lower-level steps

- [`compute_CellToCell()`](https://raredonlab.github.io/NICHESv2/reference/compute_CellToCell.md)
  : Compute Cell-to-Cell Ligand-Receptor Signals
- [`create_Edgelist_Spatial()`](https://raredonlab.github.io/NICHESv2/reference/create_Edgelist_Spatial.md)
  : Build a Spatial Neighborhood Edge List from 2-D Cell Coordinates
- [`create_Edgelist_SamplingByType()`](https://raredonlab.github.io/NICHESv2/reference/create_Edgelist_SamplingByType.md)
  : Build a Type-Stratified Sampled Edge List for Non-Spatial scRNA Data

## Example data

- [`download_NICHESv2_data()`](https://raredonlab.github.io/NICHESv2/reference/download_NICHESv2_data.md)
  : Download a Public Tutorial Dataset
- [`NICHESv2_inputs`](https://raredonlab.github.io/NICHESv2/reference/NICHESv2_inputs.md)
  : Example NICHESv2 inputs (raw triple)
- [`NICHESv2_example`](https://raredonlab.github.io/NICHESv2/reference/NICHESv2_example.md)
  : Example prebuilt NICHESObject

## LR databases

- [`connectomedb2025_human`](https://raredonlab.github.io/NICHESv2/reference/connectomedb2025_human.md)
  : ConnectomeDB2025 ligand-receptor database – human
- [`connectomedb2025_mouse`](https://raredonlab.github.io/NICHESv2/reference/connectomedb2025_mouse.md)
  : ConnectomeDB2025 ligand-receptor database – mouse
- [`connectomedb2025_pig`](https://raredonlab.github.io/NICHESv2/reference/connectomedb2025_pig.md)
  : ConnectomeDB2025 ligand-receptor database – pig
- [`connectomedb2025_rat`](https://raredonlab.github.io/NICHESv2/reference/connectomedb2025_rat.md)
  : ConnectomeDB2025 ligand-receptor database – rat
- [`ncomms8866`](https://raredonlab.github.io/NICHESv2/reference/ncomms8866.md)
  : FANTOM5 ligand-receptor database – original unfiltered table
- [`ncomms8866_human`](https://raredonlab.github.io/NICHESv2/reference/ncomms8866_human.md)
  : FANTOM5 ligand-receptor database – human
- [`ncomms8866_mouse`](https://raredonlab.github.io/NICHESv2/reference/ncomms8866_mouse.md)
  : FANTOM5 ligand-receptor database – mouse
- [`ncomms8866_pig`](https://raredonlab.github.io/NICHESv2/reference/ncomms8866_pig.md)
  : FANTOM5 ligand-receptor database – pig
- [`ncomms8866_rat`](https://raredonlab.github.io/NICHESv2/reference/ncomms8866_rat.md)
  : FANTOM5 ligand-receptor database – rat
