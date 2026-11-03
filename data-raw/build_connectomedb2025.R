## build_connectomedb2025.R ----
# Developer script -- run once to build bundled ConnectomeDB2025 .rda objects.
# Output: data/connectomedb2025_{species}.rda for human, mouse, rat, pig.
# Re-run whenever ConnectomeDB2025 releases a new version.
#
# Usage (from package root):
#   source("data-raw/build_connectomedb2025.R")
#
# Requires: nothing beyond base R. No internet-connected user dependency.

## -- 1. Configuration --------------------------------------------------------

BASE.URL <- "https://connectomedb.org/downloads/Current-Release/CSV"

SPECIES.MAP <- list(
  human = "ConnectomeDB2025_human.csv",
  mouse = "ConnectomeDB2025_mouse.csv",
  rat   = "ConnectomeDB2025_rat.csv",
  pig   = "ConnectomeDB2025_pig.csv"
)

## -- 2. Helper: extract primary gene symbol ----------------------------------
#
# Input format: "GENE1 (ALIAS1, ALIAS2, ...)"  or just "GENE1"
# Output: "GENE1"
# All commas in the raw data are inside the parenthesis block (aliases only).

.extract_primary <- function(x) {
  trimws(sub(" \\(.*", "", x))
}

## -- 3. Helper: process one species CSV --------------------------------------

.process_species_csv <- function(species.name, file.name, base.url, verbose = TRUE) {

  url <- paste0(base.url, "/", file.name)

  if (verbose)
    message(sprintf("[ConnectomeDB2025] Downloading %s ...", file.name))

  raw <- tryCatch(
    read.csv(url(url), check.names = TRUE, stringsAsFactors = FALSE),
    error = function(e)
      stop(sprintf("[build_connectomedb2025] Failed to download '%s': %s",
                   url, conditionMessage(e)))
  )

  if (verbose)
    message(sprintf("[ConnectomeDB2025] %d rows downloaded.", nrow(raw)))

  ## -- Required columns (names after check.names = TRUE) --------------------
  ## Raw CSV headers and their R-safe equivalents:
  ##   "Interaction ID"        -> Interaction.ID
  ##   "LR Pair"               -> LR.Pair
  ##   "Evidence"              -> Evidence
  ##   "Ligand Symbols"        -> Ligand.Symbols
  ##   "Receptor Symbols"      -> Receptor.Symbols
  ##   "Ligand Location"       -> Ligand.Location
  ##   "Receptor Location"     -> Receptor.Location

  required <- c("Interaction.ID", "LR.Pair", "Evidence",
                 "Ligand.Symbols", "Receptor.Symbols",
                 "Ligand.Location", "Receptor.Location")

  missing.cols <- setdiff(required, names(raw))
  if (length(missing.cols) > 0L)
    stop(sprintf(
      "[build_connectomedb2025] Species '%s': missing expected columns: %s",
      species.name, paste(missing.cols, collapse = ", ")))

  ## -- Extract primary symbols -----------------------------------------------
  ligand.raw   <- raw[["Ligand.Symbols"]]
  receptor.raw <- raw[["Receptor.Symbols"]]

  ligand.primary   <- .extract_primary(ligand.raw)
  receptor.primary <- .extract_primary(receptor.raw)

  ## -- Sanity check: no empty symbols after extraction -----------------------
  empty.lig <- which(nchar(ligand.primary) == 0L)
  empty.rec <- which(nchar(receptor.primary) == 0L)
  if (length(empty.lig) > 0L)
    warning(sprintf(
      "[build_connectomedb2025] Species '%s': %d empty ligand symbol(s) after extraction (rows: %s). These rows will be dropped.",
      species.name, length(empty.lig), paste(head(empty.lig, 5L), collapse = ",")))
  if (length(empty.rec) > 0L)
    warning(sprintf(
      "[build_connectomedb2025] Species '%s': %d empty receptor symbol(s) after extraction (rows: %s). These rows will be dropped.",
      species.name, length(empty.rec), paste(head(empty.rec, 5L), collapse = ",")))

  keep <- which(nchar(ligand.primary) > 0L & nchar(receptor.primary) > 0L)

  ## -- Assemble output data.frame -------------------------------------------
  ## Standard NICHESv2 columns: ligand | receptor | <extras>
  ## Extras follow dot-naming convention.
  out <- data.frame(
    ligand            = ligand.primary[keep],
    receptor          = receptor.primary[keep],
    interaction.id    = raw[["Interaction.ID"]][keep],
    lr.pair           = raw[["LR.Pair"]][keep],
    evidence          = raw[["Evidence"]][keep],
    ligand.location   = raw[["Ligand.Location"]][keep],
    receptor.location = raw[["Receptor.Location"]][keep],
    stringsAsFactors  = FALSE
  )

  if (verbose)
    message(sprintf(
      "[ConnectomeDB2025] %s: %d valid LR pairs retained.",
      species.name, nrow(out)))

  out
}

## -- 4. Build and save each species ------------------------------------------

data.dir <- file.path("data")
if (!dir.exists(data.dir))
  dir.create(data.dir, recursive = TRUE)

for (sp in names(SPECIES.MAP)) {

  db.obj  <- .process_species_csv(
    species.name = sp,
    file.name    = SPECIES.MAP[[sp]],
    base.url     = BASE.URL,
    verbose      = TRUE
  )

  ## Object name follows FANTOM5 convention: connectomedb2025_{species}
  obj.name <- paste0("connectomedb2025_", sp)
  assign(obj.name, db.obj)

  out.path <- file.path(data.dir, paste0(obj.name, ".rda"))
  save(list = obj.name, file = out.path, compress = "xz")

  message(sprintf("[ConnectomeDB2025] Saved: %s  (%d KB)",
                  out.path,
                  round(file.size(out.path) / 1024L)))

  rm(list = obj.name)
}

message("\n[ConnectomeDB2025] All species built successfully.")
message("Re-run devtools::document() twice after adding new data objects.")
