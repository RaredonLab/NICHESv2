# Attach Cell Segmentation Polygons to a NICHESObject

Stores a pre-built `sf` object of cell segmentation polygons in the
`$segmentation` slot of a `NICHESObject`.

NICHESv2 does not extract segmentation data from Seurat objects – that
step is platform-specific and belongs in the user's own preprocessing
script. Once you have built an `sf` object (using sf), call this
function to attach it so that downstream plotting utilities can overlay
cell boundaries.

The `sf` object must have a `cell` column whose values are cell
barcodes, and every cell in `niches.obj$cell.list` must have a
corresponding polygon – since both the count matrix and the segmentation
are extracted from the same source object, missing polygons indicate a
barcode mismatch that must be resolved upstream. An error is raised if
any cell is absent. Extra polygons for barcodes not in the object are
permitted with a warning.

## Usage

``` r
add_segmentation(niches.obj, seg.sf)
```

## Arguments

- niches.obj:

  A `NICHESObject`.

- seg.sf:

  An `sf` object with at least a `cell` column containing cell barcodes
  and an `sf` geometry column of cell polygons. Barcodes must match
  those in `niches.obj$cell.list`.

## Value

The input `NICHESObject` with `$segmentation` populated.

## See also

[`NICHESv2_example`](https://raredonlab.github.io/NICHESv2/reference/NICHESv2_example.md)
for an example object whose `$segmentation` slot is `NULL` and ready to
receive an attached `sf` object;
[`merge_NICHESObjects`](https://raredonlab.github.io/NICHESv2/reference/merge_NICHESObjects.md)
which uses `sf::rbind()` to combine `$segmentation` slots across merged
objects.

## Examples

``` r
# NICHESv2 expects 'seg.sf' to be an sf data.frame with:
#   * a 'cell' column matching niches.obj$cell.list
#   * an sf geometry column of polygons (one per cell)
# The exact construction is platform-specific (Visium HD, Xenium, MERFISH,
# CosMx, etc.). The skeleton below shows the required shape, using a unit
# square for every cell:
if (requireNamespace("sf", quietly = TRUE)) {
  data(NICHESv2_example)

  # Replace 'polygon.list' with polygons from your segmentation pipeline:
  polygon.list <- lapply(NICHESv2_example$cell.list, function(cell) {
    sf::st_polygon(list(matrix(c(0, 0, 1, 0, 1, 1, 0, 1, 0, 0),
                               ncol = 2, byrow = TRUE)))
  })

  seg.sf <- sf::st_sf(
    cell     = NICHESv2_example$cell.list,
    geometry = sf::st_sfc(polygon.list)
  )

  niches.obj <- add_segmentation(NICHESv2_example, seg.sf)
}
#> [add_segmentation] Segmentation attached: 80 polygon(s) for 80 cell(s).
```
