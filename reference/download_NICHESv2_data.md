# Download a Public Tutorial Dataset

Downloads one of the two public datasets used in the website tutorials,
checks its MD5 sum, and returns a count matrix and metadata ready for
[`create_NICHESObject()`](https://raredonlab.github.io/NICHESv2/reference/create_NICHESObject.md).
The files are attached to a release of the NICHESv2 GitHub repository,
so no other package is needed. Each dataset is downloaded once and then
read from `cache.dir`.

## Usage

``` r
download_NICHESv2_data(
  dataset,
  cache.dir = tools::R_user_dir("NICHESv2", "cache"),
  force = FALSE,
  verbose = TRUE
)
```

## Arguments

- dataset:

  `"stxBrain"` or `"ifnb"`.

- cache.dir:

  Folder to keep downloaded files in. Default: the NICHESv2 user cache
  folder, `tools::R_user_dir("NICHESv2", "cache")`. Delete the folder to
  free the space.

- force:

  Download again even if the dataset is already cached. Default `FALSE`.

- verbose:

  Print progress messages. Default `TRUE`.

## Value

A list with `count.mtx` (a sparse genes x cells matrix of raw counts)
and `meta.data` (a data frame with one row per cell), the same shape
that
[`extract_NICHESInputs_Seurat()`](https://raredonlab.github.io/NICHESv2/reference/extract_NICHESInputs_Seurat.md)
returns. If the download fails, or the file doesn't match its MD5 sum, a
message explains why and `NULL` is returned invisibly.

## Datasets

- `"stxBrain"`:

  10x Genomics Visium mouse brain, sagittal-anterior section 1: 2,696
  spots and 31,053 genes, with spot coordinates in `x` and `y`. About 14
  MB. CC BY 4.0, credit 10x Genomics.

- `"ifnb"`:

  Human PBMCs, untreated or stimulated with interferon-beta (Kang et
  al., 2018): 13,999 cells and 14,053 genes, with the condition in
  `stim` and cell-type labels from the Satija Lab in
  `seurat_annotations`. About 11 MB. CC BY 4.0.

Each dataset's folder in `cache.dir` also holds a `PROVENANCE.md` file
stating its source, license and the changes made.

## References

Kang HM, Subramaniam M, Targ S, et al. (2018). Multiplexed droplet
single-cell RNA-sequencing using natural genetic variation. *Nature
Biotechnology*, 36(1), 89–94.
[doi:10.1038/nbt.4042](https://doi.org/10.1038/nbt.4042)

## Examples

``` r
# \donttest{
# Downloads about 11 MB, kept in a temporary folder here.
inputs <- download_NICHESv2_data("ifnb", cache.dir = tempdir())
#> [download_NICHESv2_data] Downloading 'ifnb' (11 MB) ...
if (!is.null(inputs)) dim(inputs$count.mtx)
#> [1] 14053 13999
# }
```
