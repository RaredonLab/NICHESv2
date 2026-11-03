# Load a Ligand-Receptor Pair Database

Returns a standardized ligand-receptor pair table for use as the
`LRM.db` argument of
[`create_NICHESObject()`](https://raredonlab.github.io/NICHESv2/reference/create_NICHESObject.md)
and
[`compute_CellToCell()`](https://raredonlab.github.io/NICHESv2/reference/compute_CellToCell.md).
Two curated databases can be loaded by name; alternatively, any
user-supplied `data.frame` is accepted and column names are remapped to
the canonical `"ligand"` / `"receptor"` names.

## Usage

``` r
load_LRM_database(
  db = "connectomedb2025",
  species = "human",
  ligand.col = "ligand",
  receptor.col = "receptor",
  verbose = TRUE
)
```

## Arguments

- db:

  Character string (`"connectomedb2025"` (default), `"fantom5"`, or
  `"omnipath"`) OR a `data.frame` / `data.table` with at least two
  columns containing ligand and receptor gene symbols.

- species:

  Species for built-in databases. One of `"human"` (default), `"mouse"`,
  `"rat"`, or `"pig"` (ConnectomeDB2025 and FANTOM5 only; OmniPath does
  not support pig). Ignored silently when `db` is a `data.frame`.

- ligand.col:

  Column name for ligand symbols when `db` is a `data.frame`. Default:
  `"ligand"`.

- receptor.col:

  Column name for receptor symbols when `db` is a `data.frame`. Default:
  `"receptor"`.

- verbose:

  Print progress messages. Default `TRUE`.

## Value

A `data.frame` with columns `"ligand"`, `"receptor"`, and any additional
annotation columns present in the source. Rows with `NA` in ligand or
receptor are silently removed.

## Built-in databases

- `"connectomedb2025"`:

  ConnectomeDB2025 ligand-receptor pairs (default), bundled within the
  package for human, mouse, rat, and pig. No internet access or
  additional packages required. Includes experimentally verified and
  orthology-inferred interactions with evidence labels and subcellular
  location annotations.

- `"fantom5"`:

  Ramilowski et al. (2015) FANTOM5 ligand-receptor pairs, bundled within
  the package for human, mouse, rat, and pig. No internet access or
  additional packages required. Species-specific datasets were derived
  using the methodology in Raredon et al. 2019 (DOI:
  10.1126/sciadv.aaw3851).

- `"omnipath"`:

  OmniPath ligand-receptor extra interactions (Turei et al. 2021).
  Requires internet access and the OmnipathR package. Install with:
  `BiocManager::install("OmnipathR")`. Supported species: `"human"`,
  `"mouse"`, `"rat"`.

## Complex ligands/receptors

Both built-in databases and user-supplied tables may contain gene
complexes. Complexes are represented as underscore-joined gene symbols,
e.g. `"COL1A1_COL1A2"` for the collagen I heterodimer. The scoring
functions in NICHESv2 handle complexes automatically.

## Examples

``` r
# ConnectomeDB2025 is the default -- no db argument needed.
# \donttest{
lrm.human <- load_LRM_database(species = "human")
#> [load_LRM_database] Loading ConnectomeDB2025 (human) from bundled package data ...
#> [load_LRM_database] ConnectomeDB2025 (human): 3550 unique LR pairs loaded.
head(lrm.human)
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

# FANTOM5 -- bundled, always available, no internet.
lrm.f5  <- load_LRM_database("fantom5", species = "human")
#> [load_LRM_database] Loading FANTOM5 (human) from bundled package data ...
#> [load_LRM_database] FANTOM5 (human): 2557 unique LR pairs loaded.
lrm.rat <- load_LRM_database("fantom5", species = "rat")
#> [load_LRM_database] Loading FANTOM5 (rat) from bundled package data ...
#> [load_LRM_database] FANTOM5 (rat): 2557 unique LR pairs loaded.

# User-supplied table -- species argument is ignored. Column names are
# remapped to canonical "ligand" / "receptor".
my.db <- data.frame(
  Ligand   = c("MIF", "VEGFA", "Cxcl12_Pf4"),
  Receptor = c("CD74", "KDR",  "Ackr3")
)
lrm.user <- load_LRM_database(my.db,
                              ligand.col   = "Ligand",
                              receptor.col = "Receptor")
#> [load_LRM_database] Custom table: 3 LR pairs loaded.
# }

# OmniPath needs internet access, so this runs in interactive sessions only.
if (interactive()) {
  lrm.op <- load_LRM_database("omnipath", species = "mouse")
}
```
