# Add Neighborhood Slots to a NICHESObject

Extracts neighborhoods from the edge list of a `NICHESObject` and
populates the three neighborhood slots: `$neighborhood.cell.list`,
`$neighborhood.edge.list`, and `$neighborhood.meta`.

## Usage

``` r
add_neighborhoods(
  niches.obj,
  edge.filter.mode = "both_in_neighborhood",
  remove.autocrine = FALSE,
  overwrite = FALSE,
  n.cores = 1L,
  verbose = TRUE
)
```

## Arguments

- niches.obj:

  A `NICHESObject` with a non-empty `$edge.list`.

- edge.filter.mode:

  Character. Passed to
  [`extract_neighborhoods`](https://raredonlab.github.io/NICHESv2/reference/extract_neighborhoods.md).
  Default `"both_in_neighborhood"`.

- remove.autocrine:

  Logical. Passed to
  [`extract_neighborhoods`](https://raredonlab.github.io/NICHESv2/reference/extract_neighborhoods.md).
  Default `FALSE`.

- overwrite:

  Logical. If `FALSE` (default), raises a hard error when the
  neighborhood slots are already populated. Set to `TRUE` to recompute;
  stale neighborhood-dependent aggregation slots are cleared
  automatically.

- n.cores:

  Integer. Number of cores for parallel processing. Default `1L`.

- verbose:

  Logical. Print progress messages. Default `TRUE`.

## Value

The input `NICHESObject` with the three neighborhood slots Attributes
`"edge.filter.mode"` and `"remove.autocrine"` are set on the object.

## Details

This function calls
[`extract_neighborhoods`](https://raredonlab.github.io/NICHESv2/reference/extract_neighborhoods.md)
internally. If the neighborhood slots are already populated and
`overwrite = FALSE` (the default), a hard error is raised. Pass
`overwrite = TRUE` to recompute; any stale neighborhood-dependent
aggregation slots (`$aggregations$neighborhood.edge.agg` and
`$aggregations$neighborhood.composition`) are cleared automatically.

## Examples

``` r
# \donttest{
data(NICHESv2_example)

# Default neighborhood: edges where BOTH endpoints lie in the same
# neighborhood as the central cell
obj.nb <- add_neighborhoods(NICHESv2_example, overwrite = TRUE)
#> [add_neighborhoods] Recomputing neighborhood slots.
#> Step 1/4: Building integer index mapping...
#>   80 cells, 354 edges
#> Step 2/4: Building adjacency lists...
#>   Adjacency lists built
#> Step 3/4: Building edge lookup indices...
#>   Edge lookup indices built
#> Step 4/4: Extracting neighborhoods [mode: both_in_neighborhood]...
#>   Running in serial mode...
#>   Complete: 0.0 seconds (44384 cells/sec)
#>   Avg neighbors/cell: 4.4
#>   Avg edges/cell: 17.1

length(obj.nb$neighborhood.cell.list)   # one entry per cell
#> [1] 80
head(obj.nb$neighborhood.meta)          # n.neighbors and n.edges per cell
#> Key: <cell>
#>       cell n.neighbors n.edges
#>     <char>       <int>   <int>
#> 1:  Cell_1           5      17
#> 2: Cell_10           6      30
#> 3: Cell_11           6      28
#> 4: Cell_12           4      10
#> 5: Cell_13           2       4
#> 6: Cell_14           6      30

# Inspect one cell's neighborhood
first.cell <- names(obj.nb$neighborhood.cell.list)[1]
obj.nb$neighborhood.cell.list[[first.cell]]
#> [1] "Cell_1"  "Cell_10" "Cell_14" "Cell_25" "Cell_4"  "Cell_6" 
obj.nb$neighborhood.edge.list[[first.cell]]
#>  [1] "Cell_10|Cell_1"  "Cell_10|Cell_10" "Cell_10|Cell_14" "Cell_10|Cell_25"
#>  [5] "Cell_10|Cell_4"  "Cell_10|Cell_6"  "Cell_1|Cell_1"   "Cell_1|Cell_10" 
#>  [9] "Cell_1|Cell_14"  "Cell_14|Cell_1"  "Cell_14|Cell_10" "Cell_14|Cell_14"
#> [13] "Cell_14|Cell_25" "Cell_14|Cell_4"  "Cell_14|Cell_6"  "Cell_25|Cell_10"
#> [17] "Cell_25|Cell_14" "Cell_25|Cell_25" "Cell_25|Cell_4"  "Cell_25|Cell_6" 
#> [21] "Cell_4|Cell_10"  "Cell_4|Cell_14"  "Cell_4|Cell_25"  "Cell_4|Cell_4"  
#> [25] "Cell_4|Cell_6"   "Cell_6|Cell_10"  "Cell_6|Cell_14"  "Cell_6|Cell_25" 
#> [29] "Cell_6|Cell_4"   "Cell_6|Cell_6"  
# }
```
