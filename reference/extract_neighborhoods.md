# Extract Neighborhoods from an Edge List

For each cell in the edge list, identifies its spatial neighbors and the
edges that fall within its neighborhood. Returns three structures
suitable for populating the neighborhood slots of a `NICHESObject`.

## Usage

``` r
extract_neighborhoods(
  edge.list,
  edge.filter.mode = "both_in_neighborhood",
  remove.autocrine = FALSE,
  n.cores = 1,
  verbose = TRUE
)
```

## Arguments

- edge.list:

  A `data.table` with columns `edge`, `SendingCell`, and
  `ReceivingCell`. Typically `niches.obj$edge.list` or output of
  [`create_Edgelist_Spatial()`](https://raredonlab.github.io/NICHESv2/reference/create_Edgelist_Spatial.md).

- edge.filter.mode:

  Character. One of:

  `"both_in_neighborhood"`

  :   (default) Both endpoints of an edge must be in the cell's
      neighborhood (cell itself + its neighbors).

  `"any_in_neighborhood"`

  :   At least one endpoint must be in the neighborhood.

  `"central_cell_only"`

  :   The central cell must be an endpoint.

- remove.autocrine:

  Logical. If `TRUE`, removes self-loops from the adjacency (neighbor)
  list. Autocrine edges may still appear in `neighborhood.edge.list` if
  they pass the edge filter. Default `FALSE`.

- n.cores:

  Integer. Number of cores for parallel processing. Default 1 (serial).

- verbose:

  Logical. Print progress messages. Default `TRUE`.

## Value

A list with three elements:

- `neighborhood.cell.list`:

  Named list of character vectors. Each element is keyed by a cell
  barcode and contains the barcodes of its neighbors.

- `neighborhood.edge.list`:

  Named list of character vectors. Each element is keyed by a cell
  barcode and contains the edge IDs (format
  `"SendingCell|ReceivingCell"`) belonging to that cell's neighborhood.

- `neighborhood.meta`:

  A `data.table` with columns `cell`, `n.neighbors`, `n.edges`, keyed on
  `cell`.

Attributes: `"edge.filter.mode"`, `"remove.autocrine"`.

## Details

Neighborhoods are defined topologically from the edge list: cell B is a
neighbor of cell A if any edge connects them (in either direction). The
`edge.filter.mode` parameter controls which edges are assigned to each
cell's neighborhood.

This function is a low-level primitive: it operates on a raw edge list
(a `data.table` of `(edge, SendingCell, ReceivingCell)` rows) rather
than on a full `NICHESObject`. Most users should call
[`add_neighborhoods`](https://raredonlab.github.io/NICHESv2/reference/add_neighborhoods.md)
instead, which extracts neighborhoods from `niches.obj$edge.list` and
writes them back into the three neighborhood slots in one step.
`extract_neighborhoods()` is exposed for power users and for testing
pipelines that build edge lists outside the standard
[`create_NICHESObject()`](https://raredonlab.github.io/NICHESv2/reference/create_NICHESObject.md)
workflow.

## Examples

``` r
data(NICHESv2_example)
nbhd <- extract_neighborhoods(NICHESv2_example$edge.list, verbose = FALSE)
head(nbhd$neighborhood.meta)
#> Key: <cell>
#>       cell n.neighbors n.edges
#>     <char>       <int>   <int>
#> 1:  Cell_1           5      17
#> 2: Cell_10           6      30
#> 3: Cell_11           6      28
#> 4: Cell_12           4      10
#> 5: Cell_13           2       4
#> 6: Cell_14           6      30
```
