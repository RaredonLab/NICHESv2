# ConnectomeDB2025 ligand-receptor database – human

Human ligand-receptor pairs from ConnectomeDB2025, bundled within the
package. This is the table loaded when
`load_LRM_database("connectomedb2025", species = "human")` is called,
and also when
[`load_LRM_database()`](https://raredonlab.github.io/NICHESv2/reference/load_LRM_database.md)
is called with no arguments (the package default).

## Format

A `data.frame` with one row per ligand-receptor pair. Key columns:

- `ligand`:

  Gene symbol for the ligand (character).

- `receptor`:

  Gene symbol for the receptor (character).

Additional annotation columns (evidence category, subcellular location,
literature support, etc.) are retained from the source database.

## Source

ConnectomeDB2025. <https://connectomedb.org>.

## See also

[`load_LRM_database`](https://raredonlab.github.io/NICHESv2/reference/load_LRM_database.md)
for the recommended user interface;
[`connectomedb2025_mouse`](https://raredonlab.github.io/NICHESv2/reference/connectomedb2025_mouse.md),
[`connectomedb2025_rat`](https://raredonlab.github.io/NICHESv2/reference/connectomedb2025_rat.md),
[`connectomedb2025_pig`](https://raredonlab.github.io/NICHESv2/reference/connectomedb2025_pig.md).

## Examples

``` r
# \donttest{
db <- load_LRM_database("connectomedb2025", species = "human")
#> [load_LRM_database] Loading ConnectomeDB2025 (human) from bundled package data ...
#> [load_LRM_database] ConnectomeDB2025 (human): 3550 unique LR pairs loaded.
head(db)
#>   ligand receptor interaction.id      lr.pair evidence         ligand.location
#> 1    A2M     LRP1  CDB15:0000001     A2M LRP1   Direct                secreted
#> 2    ACE   BDKRB2  CDB15:0000002   ACE BDKRB2   Direct secreted, cell membrane
#> 3 ADAM10    EPHA3  CDB15:0000003 ADAM10 EPHA3   Direct           cell membrane
#> 4 ADAM12    ITGA9  CDB15:0000004 ADAM12 ITGA9   Direct secreted, cell membrane
#> 5 ADAM12    ITGB1  CDB15:0000005 ADAM12 ITGB1   Direct secreted, cell membrane
#> 6 ADAM12     SDC4  CDB15:0000006  ADAM12 SDC4 Inferred secreted, cell membrane
#>   receptor.location
#> 1     cell membrane
#> 2     cell membrane
#> 3     cell membrane
#> 4     cell membrane
#> 5     cell membrane
#> 6     cell membrane
nrow(db)
#> [1] 3550
# }
```
