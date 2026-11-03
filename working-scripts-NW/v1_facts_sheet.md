# NICHESv2 1.0.0 facts sheet

Internal note. Every document (vignette, README, NEWS, website, manuscript)
copies these facts, so they cannot drift apart. Each fact was checked against
the code or a real run on 2026-11-02. If the code changes, update this file
first.

## Package

- Name: NICHESv2. Free on CRAN (including the archive) and Bioconductor.
- Version 1.0.0. R >= 4.1.0.
- License: GNU AGPL v3, version 3 only (`License: AGPL-3`), plus a commercial
  license from the Raredon Lab, the same as TissuePlex. Copyright holders:
  Nuoya Wang and Micha Sam Brickman Raredon; Yale is not listed, as for
  TissuePlex. Versions up to 0.1.0 were released under MIT and stay MIT.
- Authors: Nuoya Wang (maintainer), Micha Sam Brickman Raredon.
- Repository: https://github.com/RaredonLab/NICHESv2
- Install today: `remotes::install_github("RaredonLab/NICHESv2", ref = "main")`,
  which works once dev is merged into main.
  Change to `install.packages("NICHESv2")` after CRAN acceptance.
- Imports: data.table, Matrix, methods, RANN, dbscan, ggplot2, parallel, stats,
  utils.
- Suggests: SeuratObject, rhdf5, OmnipathR, arrow, testthat, withr, sf, knitr,
  rmarkdown.
- 20 exported functions. S3 methods: `[`, `dim`, `length`, `print`, `summary`.

## Not known yet

- Website URL, the Zenodo code archive DOI, paper title and citation.
- Do not write these into any document until they exist.

## Analysis modes

- `mode = "spatial"` (default): supply exactly one of `k` (kNN, mutual edges
  only) or `rad` (fixed radius). `meta.data` needs `x` and `y`. Every cell also
  gets an autocrine self-loop edge.
- `mode = "nonspatial"`: type-stratified random sampling. Requires
  `cell.type.col`. `seed` defaults to 42. `max.cells` defaults to `NULL`.
- There is no `"sampled"` mode. That name was removed.

## Scoring

- `method`: `"product"` (default) or `"minimum"`.
- `normalize.method`: `"prop"` (default) or `"none"`.
- `create_NICHESObject()` has no default for `LRM.db`.
- `n.cores` defaults to 1 and uses socket clusters on every platform.
- Gene complexes use underscores, such as `COL1A1_COL1A2`.
- Key separators: edge `"SendingCell|ReceivingCell"`, LRM `"ligand|receptor"`.

## LR databases, rows in the bundled tables

| Database | human | mouse | rat | pig |
|---|---|---|---|---|
| `connectomedb2025` (default) | 3,550 | 3,252 | 3,347 | 3,358 |
| `fantom5` | 2,557 | 2,557 | 2,557 | 2,552 |

`omnipath` needs internet access and supports human, mouse and rat, not pig.
Custom data frames need ligand and receptor columns.

## NICHESObject

- S3 list with 15 slots at four levels, plus two extras:
  - cell: `$cell.list`, `$cell.data`, `$cell.meta`
  - edge: `$edge.list`, `$edge.data`, `$edge.meta`
  - neighborhood: `$neighborhood.cell.list`, `$neighborhood.edge.list`,
    `$neighborhood.meta`
  - gene and LRM: `$gene.list`, `$gene.meta`, `$LRM.list`, `$LRM.meta`
  - extras: `$segmentation`, `$aggregations`
- `$cell.data` and `$edge.data` are keyed long-format `data.table`s holding
  non-zero values only.
- Neighborhood slots are filled at construction with
  `edge.filter.mode = "both_in_neighborhood"` and `remove.autocrine = FALSE`.
  `add_neighborhoods(overwrite = TRUE)` changes them.
- Attributes: `mode`, `sample.col`, `edge.filter.mode`, `remove.autocrine`.
- `add_cell_meta()` replaces all metadata. `add_meta_column()` changes one column.

## Aggregation

- `aggregate_NICHESObject()` runs `add_neighborhood_edge_agg()` and, when
  `cell.type.col` is given, `add_neighborhood_composition()`.
- Edge aggregation: 27 columns per (cell, LRM). Five categories: `in`, `out`,
  `cross`, `self.autocrine`, `other.autocrine`. `chunk.size` defaults to 50.
- Composition needs `mode = "spatial"`, so it isn't available for `ifnb`.

## Subsetting

- `obj[cells]` keeps each retained cell's neighborhood entries and aggregation
  rows unchanged, so they still describe the full tissue.
- Edge slots keep only edges between retained cells.
- Aggregate the full object first, then subset.
- Recomputing on a subset stops with an error, unless
  `add_neighborhoods(obj, overwrite = TRUE)` runs first. That treats the subset
  as a separate tissue and clears the carried-over aggregations.

## Threads

- `n.cores` is the only setting the package asks users about.
- NICHESv2 doesn't change `data.table`'s thread setting. Thread count made no
  measurable difference on `stxBrain` or `ifnb` (2, 8 and 16 threads).

## Bundled data

- `data(NICHESv2_inputs)`: list with `count.mtx` (50 genes x 80 cells),
  `meta.data` (`x`, `y`, `nCount_RNA`, `celltype`) and `LRM.db` (15 pairs).
  Cell types: TypeA 30, TypeB 25, TypeC 25.
- `data(NICHESv2_example)`: spatial NICHESObject with 80 cells, 50 genes, 354
  edges and 15 LRMs. `$aggregations` is empty.
- `system.file("extdata", "NICHESv2_example.h5ad", package = "NICHESv2")`: small
  AnnData file with the same synthetic data.
- Both are synthetic. Use `data()` to load them, because `LazyData` is false.

## Public test data

Runs below used raw counts, `method = "product"`, the default database and
`n.cores = 8` on a 16-core Apple M3 Max with 64 GB of memory.

| | `stxBrain` (anterior1) | `ifnb` |
|---|---|---|
| Size | 31,053 genes x 2,696 spots | 14,053 genes x 13,999 cells |
| Download | 143 MB | 413 MB |
| Mode | spatial, `k = 6`, mouse | nonspatial, human |
| Other settings | assay `"Spatial"` | `sample.col = "stim"`, `cell.type.col = "seurat_annotations"` |
| Scored rows in `$edge.data` | 5,191,520 | 605,973 |
| Create | about 77 s | about 27 s |
| Aggregate | about 12 s | not run |
| TissuePlex export | about 4 s | not run |

- Both datasets come from SeuratData, which isn't on CRAN. NICHESv2 hosts no
  data. Install with `remotes::install_github("satijalab/seurat-data")`, then
  `SeuratData::InstallData("stxBrain")` and `SeuratData::InstallData("ifnb")`.
  SeuratData imports Seurat, so installing it also installs Seurat. NICHESv2
  itself only needs SeuratObject to read Seurat objects. SeuratData stays out of
  `DESCRIPTION`.
- `UpdateSeuratObject()` is not needed. `extract_NICHESInputs_Seurat()` gives
  identical results with or without it on both datasets (SeuratObject 5.0.2 and
  5.1.0). SeuratData objects keep their own Seurat version: `ifnb` is a V4-style
  `Assay` and `anterior1` is `Assay5`.
- `stxBrain` needs `assay = "Spatial"`. Add coordinates to the metadata with
  `SeuratObject::GetTissueCoordinates()`, which returns `x`, `y` and `cell`.
- `stxBrain` has no cell-type labels, so composition needs labels from the user.
- `ifnb` has 13 cell types with no missing labels. Conditions: CTRL 6,548 cells,
  STIM 7,451 cells.
- Earlier dry run: `stxBrain` took 7.5 min on 1 core and 3.9 min on 2 cores with
  identical scores. Its object was 684 MB, its parquet file 42 MB. `ifnb` on 2
  cores took 62 s, with a 421 MB object and 3.7 GB peak memory.
- Never use the PPLR data. It is private.

### Licenses, checked 2026-11-02

- `stxBrain` anterior1 is the 10x Genomics dataset "Mouse Brain Serial Section 1
  (Sagittal-Anterior)", Space Ranger 1.0.0, under CC BY 4.0. Any copy must credit
  10x Genomics, link to the license and state what was changed.
- `ifnb` comes from Kang et al., Nat Biotechnol 2018;36(1):89-94 (GSE96583), with
  cell-type labels added by the Satija Lab. The `ifnb.SeuratData` package is
  labeled CC BY 4.0, but no license from the original authors was found.

## TissuePlex

- NICHESv2 stops at the exported `edges.parquet` file.
- Viewing is covered by the TissuePlex manual, so link to it rather than repeat
  it: https://raredonlab.github.io/TissuePlex/
- Worked examples: https://raredonlab.github.io/TissuePlex/vignettes.html

## Wording

- US spelling in all prose. Keep published titles as printed.
- Say "ligand-receptor (LR)", "LRM", "NICHESObject", "nonspatial".
- Don't describe NICHESv2 as a Seurat or Bioconductor package.
