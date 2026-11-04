# Starting from a Seurat object

Most single-cell data in R lives in Seurat objects. This page shows how
[`extract_NICHESInputs_Seurat()`](https://raredonlab.github.io/NICHESv2/reference/extract_NICHESInputs_Seurat.md)
turns two real Seurat objects into NICHESv2 inputs. They are the same
two datasets the tutorials use, and these are the steps that made the
tutorial inputs.

The objects come from
[SeuratData](https://github.com/satijalab/seurat-data), which the Satija
Lab maintains. The tutorials don’t need it, because they download the
finished inputs with
[`download_NICHESv2_data()`](https://raredonlab.github.io/NICHESv2/reference/download_NICHESv2_data.md).

## Install SeuratData

SeuratData isn’t on CRAN. Run this once; the downloads are about 143 MB
for `stxBrain` and 413 MB for `ifnb`. SeuratData also installs the
Seurat package, so the first install can take a while.

``` r

install.packages(c("remotes", "SeuratObject"))
remotes::install_github("satijalab/seurat-data")
SeuratData::InstallData(c("stxBrain", "ifnb"))
```

``` r

library(NICHESv2)
library(SeuratObject)
library(SeuratData)
```

SeuratData objects keep the Seurat version they were saved with.
[`extract_NICHESInputs_Seurat()`](https://raredonlab.github.io/NICHESv2/reference/extract_NICHESInputs_Seurat.md)
reads both V4 and V5 objects, so you don’t need to update them first.
The spatial object below is V5 and the non-spatial one is V4.

## A spatial object: Visium mouse brain

`stxBrain` holds four sections, and the spatial tutorial uses
`anterior1`.

``` r

data("anterior1")

spatial <- extract_NICHESInputs_Seurat(
  seurat.obj = anterior1,
  assay      = "Spatial",
  layer      = "counts",
  verbose    = FALSE
)
dim(spatial$count.mtx)
#> [1] 31053  2696
```

Spatial coordinates sit in the image part of the Seurat object, so add
them to the metadata as `x` and `y`. Matching by barcode keeps every
spot aligned.

``` r

coords <- GetTissueCoordinates(anterior1)
spatial$meta.data$x <- coords[rownames(spatial$meta.data), "x"]
spatial$meta.data$y <- coords[rownames(spatial$meta.data), "y"]
head(spatial$meta.data[, c("x", "y")], 3)
```

|                    |     x |     y |
|:-------------------|------:|------:|
| AAACAAGTATCTCCCA-1 | 7,475 | 8,501 |
| AAACACCAATAACTGC-1 | 8,553 | 2,788 |
| AAACAGAGCGACTCCT-1 | 3,164 | 7,950 |

## A non-spatial object: stimulated PBMCs

``` r

data("ifnb")

nonspatial <- extract_NICHESInputs_Seurat(
  seurat.obj = ifnb,
  assay      = "RNA",
  layer      = "counts",
  verbose    = FALSE
)
dim(nonspatial$count.mtx)
#> [1] 14053 13999
table(nonspatial$meta.data$stim)
#> 
#> CTRL STIM 
#> 6548 7451
```

No coordinates are needed here, because non-spatial mode pairs cells by
sampling within cell types.

## Next steps

`spatial` and `nonspatial` hold the same inputs that
[`download_NICHESv2_data()`](https://raredonlab.github.io/NICHESv2/reference/download_NICHESv2_data.md)
returns. Pass either to
[`create_NICHESObject()`](https://raredonlab.github.io/NICHESv2/reference/create_NICHESObject.md)
and continue with the [spatial
tutorial](https://raredonlab.github.io/NICHESv2/articles/spatial-tutorial.md)
or the [non-spatial
tutorial](https://raredonlab.github.io/NICHESv2/articles/nonspatial-tutorial.md).
For other Seurat options, such as choosing a layer, see [Input
formats](https://raredonlab.github.io/NICHESv2/articles/input-formats.md).
