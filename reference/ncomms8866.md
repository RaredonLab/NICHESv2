# FANTOM5 ligand-receptor database – original unfiltered table

The complete, unfiltered set of human ligand-receptor pairs from the
FANTOM5 consortium as published by Ramilowski et al. (2015). This is the
source table from which the four species-specific variants
([`ncomms8866_human`](https://raredonlab.github.io/NICHESv2/reference/ncomms8866_human.md),
[`ncomms8866_mouse`](https://raredonlab.github.io/NICHESv2/reference/ncomms8866_mouse.md),
[`ncomms8866_rat`](https://raredonlab.github.io/NICHESv2/reference/ncomms8866_rat.md),
[`ncomms8866_pig`](https://raredonlab.github.io/NICHESv2/reference/ncomms8866_pig.md))
are derived.

## Format

A `data.frame` with one row per ligand-receptor pair and the following
key columns:

- `Ligand.ApprovedSymbol`:

  HGNC-approved gene symbol for the ligand (character).

- `Receptor.ApprovedSymbol`:

  HGNC-approved gene symbol for the receptor (character).

Additional annotation columns from the original publication are
retained.

## Source

Ramilowski, J.A. et al. (2015). A draft network of
ligand-receptor-mediated multicellular signalling in human. *Nature
Communications*, **6**, 7866.
[doi:10.1038/ncomms8866](https://doi.org/10.1038/ncomms8866)

## See also

[`load_LRM_database`](https://raredonlab.github.io/NICHESv2/reference/load_LRM_database.md)
for the recommended user interface;
[`ncomms8866_human`](https://raredonlab.github.io/NICHESv2/reference/ncomms8866_human.md),
[`ncomms8866_mouse`](https://raredonlab.github.io/NICHESv2/reference/ncomms8866_mouse.md),
[`ncomms8866_rat`](https://raredonlab.github.io/NICHESv2/reference/ncomms8866_rat.md),
[`ncomms8866_pig`](https://raredonlab.github.io/NICHESv2/reference/ncomms8866_pig.md)
for the species-specific filtered tables.

## Examples

``` r
# \donttest{
db <- load_LRM_database("fantom5", species = "human")
#> [load_LRM_database] Loading FANTOM5 (human) from bundled package data ...
#> [load_LRM_database] FANTOM5 (human): 2557 unique LR pairs loaded.
head(db)
#>   ligand receptor
#> 1    A2M     LRP1
#> 2  AANAT   MTNR1A
#> 3  AANAT   MTNR1B
#> 4    ACE    AGTR2
#> 5    ACE   BDKRB2
#> 6 ADAM10      AXL
nrow(db)
#> [1] 2557
# }
```
