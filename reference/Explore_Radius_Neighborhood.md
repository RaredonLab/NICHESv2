# Explore the relationship between search radius and neighborhood size

Estimates the mean number of neighboring cells as a function of search
radius, using a randomly sampled subset of cells for efficiency.
Designed to help users choose a radius parameter before calling
[`create_NICHESObject()`](https://raredonlab.github.io/NICHESv2/reference/create_NICHESObject.md).

## Usage

``` r
Explore_Radius_Neighborhood(
  coord,
  r.min = NULL,
  r.max = NULL,
  r.step = NULL,
  n.sample = 2000L,
  seed = 42L,
  n.cores = 1L,
  show.plot = TRUE
)
```

## Arguments

- coord:

  data.frame. First two columns must be numeric x and y coordinates. Row
  names must be unique cell barcodes.

- r.min:

  numeric(1). Minimum radius to evaluate. Default: 1% of the
  bounding-box diagonal.

- r.max:

  numeric(1). Maximum radius to evaluate. Default: 20% of the
  bounding-box diagonal.

- r.step:

  numeric(1). Step size between evaluated radii. Default:
  `(r.max - r.min) / 50` (50 evenly spaced evaluation points).

- n.sample:

  integer(1) or NULL. Number of cells to randomly sample for the
  neighborhood estimate. Default 2000. Set to `NULL` to use all cells –
  not recommended for large datasets.

- seed:

  integer(1) or NULL. Random seed for reproducible sampling. Default 42.
  The caller's random number state is restored on exit.

- n.cores:

  integer(1). Number of parallel workers (socket cluster) for
  distributing the radius evaluation. Default 1.

- show.plot:

  logical. Whether to print the summary ggplot before returning. Default
  `TRUE`.

## Value

A `data.table` with columns:

- radius:

  evaluated radius value

- mean.neighbors:

  mean neighbor count across sampled cells

- sd.neighbors:

  standard deviation of neighbor count

- median.neighbors:

  median neighbor count (robust to density outliers)

Returned invisibly when `show.plot = TRUE`.

## Details

A single kd-tree query via
[`dbscan::frNN()`](https://michael.hahsler.net/dbscan/reference/frNN.html)
is run at `r.max`; all smaller radii are evaluated by thresholding the
already-retrieved distance lists, so the kd-tree is built only once
regardless of how many radius steps are requested.

## Examples

``` r
# \donttest{
data(NICHESv2_inputs)

# Build a coord data.frame from the example metadata
coord <- NICHESv2_inputs$meta.data[, c("x", "y")]
rownames(coord) <- rownames(NICHESv2_inputs$meta.data)

# Sweep radius across the bounding box; small dataset -> use all 80 cells
radius.curve <- Explore_Radius_Neighborhood(
  coord     = coord,
  n.sample  = 80L,
  n.cores   = 1L,
  show.plot = FALSE
)
head(radius.curve)
#>        radius mean.neighbors sd.neighbors median.neighbors
#>         <num>          <num>        <num>            <num>
#> 1: 0.01317801          0.100    0.3765264                0
#> 2: 0.01828785          0.200    0.5372221                0
#> 3: 0.02339769          0.350    0.7647288                0
#> 4: 0.02850753          0.525    0.8711261                0
#> 5: 0.03361737          0.675    0.9246963                0
#> 6: 0.03872721          0.775    0.9933322                1
# }
```
