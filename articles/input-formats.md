# Input formats

[`create_NICHESObject()`](https://raredonlab.github.io/NICHESv2/reference/create_NICHESObject.md)
needs an expression matrix, cell metadata and an LR database. This guide
shows how to prepare each one.

``` r

library(NICHESv2)
data(NICHESv2_inputs)
```

## Expression matrix

`count.mtx` is a gene-by-cell matrix, sparse or dense. Row names are
gene symbols and column names are cell barcodes.

``` r

class(NICHESv2_inputs$count.mtx)
#> [1] "dgCMatrix"
#> attr(,"package")
#> [1] "Matrix"
dim(NICHESv2_inputs$count.mtx)
#> [1] 50 80
```

The values you pass set what the scores mean. Raw counts give absolute
magnitudes, and normalized expression gives relative scores. NICHESv2
uses the matrix as it is, so this choice affects which edges score
highest.

## Cell metadata

`meta.data` is a data frame with cell barcodes as row names, matching
the matrix columns. What else it needs depends on the mode:

| Mode           | Required columns                    |
|----------------|-------------------------------------|
| `"spatial"`    | `x` and `y`                         |
| `"nonspatial"` | the column named by `cell.type.col` |

Other columns carry through to `$cell.meta` and appear on each edge in
`$edge.meta`.

``` r

head(NICHESv2_inputs$meta.data, 3)
```

|        |     x |     y | nCount_RNA | celltype |
|:-------|------:|------:|-----------:|:---------|
| Cell_1 | 0.136 | 0.298 |         66 | TypeA    |
| Cell_2 | 0.238 | 0.172 |         53 | TypeA    |
| Cell_3 | 0.225 | 0.240 |         55 | TypeA    |

### Several samples

When your data has several samples, set `sample.col` to the column that
identifies them. Each sample is processed separately and the results are
merged, so no edges cross samples.

Barcodes must be unique across all samples, and a duplicate stops the
merge with an error. Add a prefix first if needed, such as
`paste0(sample.id, "_", barcode)`.

## LR databases

Pass one of three database names to `LRM.db`, or your own data frame.

| Name                 | Species                | Needs internet      |
|----------------------|------------------------|---------------------|
| `"connectomedb2025"` | human, mouse, rat, pig | No                  |
| `"fantom5"`          | human, mouse, rat, pig | No                  |
| `"omnipath"`         | human, mouse, rat      | Yes, plus OmnipathR |

[`create_NICHESObject()`](https://raredonlab.github.io/NICHESv2/reference/create_NICHESObject.md)
has no default database, so always set `LRM.db`. Set `species` to match
your data. To look at a database first:

``` r

db <- load_LRM_database("connectomedb2025", species = "mouse", verbose = FALSE)
nrow(db)
#> [1] 3252
head(db[, c("ligand", "receptor")], 3)
```

| ligand | receptor |
|:-------|:---------|
| A2m    | Lrp1     |
| Ace    | Bdkrb2   |
| Adam10 | Epha3    |

### Your own database

A custom data frame needs ligand and receptor columns. If they have
other names, set `ligand.col` and `receptor.col`. Multi-gene complexes
are written with underscores, and a pair is dropped if any of its genes
is missing from the matrix.

``` r

my.db <- data.frame(
  L = c("Gene1", "Gene3"),
  R = c("Gene2", "Gene4_Gene5")
)
obj <- create_NICHESObject(
  count.mtx    = NICHESv2_inputs$count.mtx,
  meta.data    = NICHESv2_inputs$meta.data,
  LRM.db       = my.db,
  ligand.col   = "L",
  receptor.col = "R",
  mode         = "spatial",
  k            = 5L,
  verbose      = FALSE
)
obj$LRM.list
#> [1] "Gene1|Gene2"       "Gene3|Gene4_Gene5"
```

## Seurat objects

[`extract_NICHESInputs_Seurat()`](https://raredonlab.github.io/NICHESv2/reference/extract_NICHESInputs_Seurat.md)
needs the SeuratObject package and reads both V4 and V5 objects, so you
don’t need to update an older object first. It returns `count.mtx` and
`meta.data`.

``` r

seu <- SeuratObject::CreateSeuratObject(
  counts    = NICHESv2_inputs$count.mtx,
  meta.data = NICHESv2_inputs$meta.data
)
inputs <- extract_NICHESInputs_Seurat(
  seu, assay = "RNA", layer = "counts", verbose = FALSE
)
dim(inputs$count.mtx)
#> [1] 50 80
```

Spatial coordinates aren’t copied automatically, because Seurat keeps
them with the image. For spatial objects, add them yourself, as shown in
the [spatial
tutorial](https://raredonlab.github.io/NICHESv2/articles/spatial-tutorial.md):

``` r

coords <- SeuratObject::GetTissueCoordinates(seurat.obj)
inputs$meta.data$x <- coords[rownames(inputs$meta.data), "x"]
inputs$meta.data$y <- coords[rownames(inputs$meta.data), "y"]
```

## AnnData files

[`extract_NICHESInputs_AnnData()`](https://raredonlab.github.io/NICHESv2/reference/extract_NICHESInputs_AnnData.md)
reads `.h5ad` files through the rhdf5 package, so Python isn’t needed.
`layer = "counts"` reads the counts layer, and `obs.cols` limits which
metadata columns are kept.

``` r

h5ad <- system.file("extdata", "NICHESv2_example.h5ad", package = "NICHESv2")
inputs <- extract_NICHESInputs_AnnData(h5ad, layer = "counts", verbose = FALSE)
dim(inputs$count.mtx)
#> [1] 50 80
colnames(inputs$meta.data)
#> [1] "celltype"   "nCount_RNA" "x"          "y"
```

Install rhdf5 from Bioconductor:

``` r

install.packages("BiocManager")
BiocManager::install("rhdf5")
```
