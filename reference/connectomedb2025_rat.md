# ConnectomeDB2025 ligand-receptor database – rat

Rat ortholog-mapped ligand-receptor pairs from ConnectomeDB2025, bundled
within the package. This is the table loaded when
`load_LRM_database("connectomedb2025", species = "rat")` is called.

## Format

A `data.frame` with one row per ligand-receptor pair. Key columns:

- `ligand`:

  Rat gene symbol for the ligand (character).

- `receptor`:

  Rat gene symbol for the receptor (character).

Additional annotation columns are retained where ortholog mapping was
possible.

## Source

ConnectomeDB2025. <https://connectomedb.org>.

## See also

[`load_LRM_database`](https://raredonlab.github.io/NICHESv2/reference/load_LRM_database.md)
for the recommended user interface;
[`connectomedb2025_human`](https://raredonlab.github.io/NICHESv2/reference/connectomedb2025_human.md),
[`connectomedb2025_mouse`](https://raredonlab.github.io/NICHESv2/reference/connectomedb2025_mouse.md),
[`connectomedb2025_pig`](https://raredonlab.github.io/NICHESv2/reference/connectomedb2025_pig.md).

## Examples

``` r
# \donttest{
db <- load_LRM_database("connectomedb2025", species = "rat")
#> [load_LRM_database] Loading ConnectomeDB2025 (rat) from bundled package data ...
#> [load_LRM_database] ConnectomeDB2025 (rat): 3347 unique LR pairs loaded.
head(db)
#>   ligand receptor interaction.id      lr.pair evidence         ligand.location
#> 1    A2m     Lrp1  CDB15:0000001     A2m Lrp1 Inferred                secreted
#> 2    Ace   Bdkrb2  CDB15:0000002   Ace Bdkrb2 Inferred secreted, cell membrane
#> 3 Adam10    Epha3  CDB15:0000003 Adam10 Epha3 Inferred           cell membrane
#> 4 Adam12    Itga9  CDB15:0000004 Adam12 Itga9 Inferred secreted, cell membrane
#> 5 Adam12    Itgb1  CDB15:0000005 Adam12 Itgb1 Inferred secreted, cell membrane
#> 6 Adam12     Sdc4  CDB15:0000006  Adam12 Sdc4 Inferred secreted, cell membrane
#>   receptor.location
#> 1     cell membrane
#> 2     cell membrane
#> 3     cell membrane
#> 4     cell membrane
#> 5     cell membrane
#> 6     cell membrane
nrow(db)
#> [1] 3347
# }
```
