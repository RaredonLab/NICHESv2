# For reviewers

This page takes you from a clean R session to real results on public
data. You need R 4.1.0 or later, and no Python or sign-up.

The two public datasets here let you test the software. They don’t
reproduce the paper’s figures, which used a separate dataset.

## 1. Install

``` r

install.packages(c("remotes", "SeuratObject", "arrow"))
remotes::install_github("RaredonLab/NICHESv2", ref = "main")
```

## 2. Confirm it works

This takes a few seconds and uses data that ship with the package.

``` r

library(NICHESv2)

data(NICHESv2_inputs)
obj <- create_NICHESObject(
  count.mtx = NICHESv2_inputs$count.mtx,
  meta.data = NICHESv2_inputs$meta.data,
  LRM.db    = NICHESv2_inputs$LRM.db,
  mode      = "spatial",
  k         = 5L
)
dim(obj)
```

You should see 80 cells, 50 genes, 354 edges and 15 LRMs.

## 3. Get the test data

Both datasets come from
[SeuratData](https://github.com/satijalab/seurat-data). They download
once, and together take about 560 MB. SeuratData also installs the
Seurat package, so this step can take a while. NICHESv2 itself doesn’t
need Seurat.

``` r

remotes::install_github("satijalab/seurat-data")
SeuratData::InstallData("stxBrain")
SeuratData::InstallData("ifnb")
```

## 4. Run the tutorials

| Tutorial | Data | Time with `n.cores = 2` |
|----|----|----|
| [Spatial: mouse brain](https://raredonlab.github.io/NICHESv2/articles/spatial-tutorial.md) | 2,696 Visium spots | About 4 minutes |
| [Non-spatial: stimulated PBMCs](https://raredonlab.github.io/NICHESv2/articles/nonspatial-tutorial.md) | 13,999 cells | About 1 minute, 4 GB of memory |

Each tutorial’s code runs from top to bottom, so you can copy it
directly into R. Results should match the printed output, because
scoring is deterministic and nonspatial sampling uses a fixed seed.

For a quick check, these are the expected object sizes:

| Tutorial    | Cells  | Genes  | Edges  | LRMs  |
|-------------|--------|--------|--------|-------|
| Spatial     | 2,696  | 31,053 | 18,500 | 3,097 |
| Non-spatial | 13,999 | 14,053 | 87,957 | 1,047 |

The spatial tutorial stores 5,191,520 non-zero scores, and the
non-spatial tutorial stores 605,973.

## 5. Explore further

- [Input
  formats](https://raredonlab.github.io/NICHESv2/articles/input-formats.md)
  shows how to use your own data.
- [The
  NICHESObject](https://raredonlab.github.io/NICHESv2/articles/niches-object.md)
  explains each slot.
- [Neighborhood
  aggregation](https://raredonlab.github.io/NICHESv2/articles/aggregation.md)
  explains the 27 aggregation columns and how subsets work.
- [Exporting to
  TissuePlex](https://raredonlab.github.io/NICHESv2/articles/tissueplex-export.md)
  links to the interactive viewer.

## If something goes wrong

- **A download fails:** the data are served by the Satija Lab. Try again
  later.
- **A step runs out of memory:** close other programs. The non-spatial
  tutorial needs about 4 GB.
- **A step feels slow:** raise `n.cores` if your computer has more
  cores.
- **Something else:** please open an
  [issue](https://github.com/RaredonLab/NICHESv2/issues).
