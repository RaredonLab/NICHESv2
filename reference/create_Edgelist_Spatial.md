# Build a Spatial Neighborhood Edge List from 2-D Cell Coordinates

Constructs a directed edge list representing spatial cell neighborhoods.
Two computation modes are available for inter-cell edges:

- **kNN mode** – finds `k` nearest neighbors per cell via a kd-tree
  ([`RANN::nn2`](https://jefferislab.github.io/RANN/reference/nn2.html)).
  Only *mutual* edges are retained (A-\>B kept only if B-\>A also
  exists). Non-mutual edges are discarded.

- **Radius mode** – finds all cells within Euclidean distance `rad` via
  a kd-tree
  ([`dbscan::frNN`](https://michael.hahsler.net/dbscan/reference/frNN.html)).
  Edges are inherently symmetric.

Autocrine self-loop edges (A-\>A, `distance = 0`) are always appended
for every cell in both modes.

## Usage

``` r
create_Edgelist_Spatial(coord, k = NULL, rad = NULL, n.cores = 1L)
```

## Arguments

- coord:

  `data.frame`. First column = x, second column = y. **Row names must be
  unique cell identifiers.**

- k:

  `integer` or `NULL`. Number of nearest neighbors (kNN mode). Positive
  integer.

- rad:

  `numeric` or `NULL`. Neighborhood radius (radius mode). Positive
  finite number in the same units as `coord`.

- n.cores:

  `integer`. CPU cores. Default 1. Uses socket-based clusters
  ([`parallel::makeCluster`](https://rdrr.io/r/parallel/makeCluster.html))
  which work on all platforms (Windows, macOS, Linux).

## Value

A `data.table` keyed on `(from, to)` with columns:

- `from`:

  Sending cell barcode (character).

- `to`:

  Receiving cell barcode (character).

- `distance`:

  Euclidean distance; 0 for self-loops (numeric).

- `x.from`, `y.from`:

  Sending cell coordinates (numeric).

- `x.to`, `y.to` :

  Receiving cell coordinates (numeric).

Carries metadata attributes: `mode`, `param`, `n.cells`, `cell.names`.

## Details

This function is a low-level primitive: it builds the spatial edge list
that
[`create_NICHESObject`](https://raredonlab.github.io/NICHESv2/reference/create_NICHESObject.md)
consumes internally. Most users should call
`create_NICHESObject(mode = "spatial", ...)`, which invokes this
function and assembles the result into a full `NICHESObject` in one
step. `create_Edgelist_Spatial()` is exposed for power users who want to
inspect or post-process the edge list, benchmark alternative `k` or
`rad` values without re-running the scoring step, or feed the edges into
a non-NICHESv2 scoring pipeline.

**Parameter resolution (priority):**

1.  Both NULL / invalid -\> default `k = 4`.

2.  Both valid -\> radius mode; user notified.

3.  Only `rad` valid -\> radius mode.

4.  Only `k` valid -\> kNN mode (mutual edges only).

## Examples

``` r
set.seed(1)
co <- data.frame(x = runif(50), y = runif(50),
                 row.names = paste0("cell_", 1:50))
el <- create_Edgelist_Spatial(co, k = 5)
#> [create_Edgelist_Spatial] Mode: kNN    | Param: 5 | Cells: 50 | Inter-cell: 198 | Autocrine: 50
el.rad <- create_Edgelist_Spatial(co, rad = 0.2)
#> [create_Edgelist_Spatial] Mode: radius | Param: 0.2 | Cells: 50 | Inter-cell: 286 | Autocrine: 50
```
