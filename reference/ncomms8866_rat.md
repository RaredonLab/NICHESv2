# FANTOM5 ligand-receptor database – rat

Rat ortholog-mapped ligand-receptor pairs derived from the human FANTOM5
table (Ramilowski et al., 2015) using the methodology of Raredon et al.
(2019). This is the table loaded when
`load_LRM_database("fantom5", species = "rat")` is called.

## Format

A `data.frame` with one row per ligand-receptor pair and the following
key columns:

- `Ligand.ApprovedSymbol`:

  Rat gene symbol for the ligand (character).

- `Receptor.ApprovedSymbol`:

  Rat gene symbol for the receptor (character).

Additional annotation columns are retained where ortholog mapping was
possible.

## Source

Ramilowski, J.A. et al. (2015). A draft network of
ligand-receptor-mediated multicellular signalling in human. *Nature
Communications*, **6**, 7866.
[doi:10.1038/ncomms8866](https://doi.org/10.1038/ncomms8866)

## See also

[`load_LRM_database`](https://raredonlab.github.io/NICHESv2/reference/load_LRM_database.md)
for the recommended user interface;
[`ncomms8866`](https://raredonlab.github.io/NICHESv2/reference/ncomms8866.md),
[`ncomms8866_human`](https://raredonlab.github.io/NICHESv2/reference/ncomms8866_human.md),
[`ncomms8866_mouse`](https://raredonlab.github.io/NICHESv2/reference/ncomms8866_mouse.md),
[`ncomms8866_pig`](https://raredonlab.github.io/NICHESv2/reference/ncomms8866_pig.md).

## Examples

``` r
# \donttest{
db <- load_LRM_database("fantom5", species = "rat")
#> [load_LRM_database] Loading FANTOM5 (rat) from bundled package data ...
#> [load_LRM_database] FANTOM5 (rat): 2557 unique LR pairs loaded.
head(db)
#>   ligand receptor
#> 1    A2m     Lrp1
#> 2  Aanat   Mtnr1a
#> 3  Aanat   Mtnr1b
#> 4    Ace    Agtr2
#> 5    Ace   Bdkrb2
#> 6 Adam10      Axl
nrow(db)
#> [1] 2557
# }
```
