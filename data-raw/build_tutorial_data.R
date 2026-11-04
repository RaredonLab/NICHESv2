# build_tutorial_data.R
#
# PURPOSE: Builds the two public datasets the website tutorials download with
#          download_NICHESv2_data(). Each dataset becomes one zip file holding
#          inputs.rds (a sparse count matrix and its metadata, ready for
#          create_NICHESObject()) and PROVENANCE.md (source, license, changes).
#          The zips and a SHA256SUMS file are uploaded by hand as assets of the
#          GitHub release named in .niches_data_url(). Both datasets are
#          CC BY 4.0, so this copy credits the sources and states the changes.
#
# NEEDS:   SeuratObject, SeuratData, and its datasets installed once with
#          SeuratData::InstallData(c("stxBrain", "ifnb")).
#
# USAGE:   Run from the package root:
#          Rscript data-raw/build_tutorial_data.R
#          Then copy the printed MD5 sums into .niches_datasets() in
#          R/download_NICHESv2_data.R. Output goes to data-raw/tutorial-data/,
#          which git ignores.

suppressPackageStartupMessages({
  library(SeuratObject)
  library(SeuratData)
})
devtools::load_all(quiet = TRUE)

out.dir <- "data-raw/tutorial-data"
unlink(out.dir, recursive = TRUE)
dir.create(out.dir, recursive = TRUE)

.write_dataset <- function(name, inputs, provenance) {
  folder <- file.path(out.dir, name)
  dir.create(folder)
  saveRDS(inputs, file.path(folder, "inputs.rds"), compress = "xz")
  writeLines(provenance, file.path(folder, "PROVENANCE.md"))
  zip.file <- paste0(name, ".zip")
  old.wd <- setwd(out.dir)
  on.exit(setwd(old.wd))
  utils::zip(zip.file, files = name, flags = "-r9Xq")
  unlink(name, recursive = TRUE)
  invisible(file.path(out.dir, zip.file))
}


## stxBrain: 10x Visium mouse brain, anterior section 1 ----

data("anterior1")
spatial <- extract_NICHESInputs_Seurat(
  seurat.obj = anterior1,
  assay      = "Spatial",
  layer      = "counts",
  verbose    = FALSE
)
coords <- GetTissueCoordinates(anterior1)
spatial$meta.data$x <- coords[rownames(spatial$meta.data), "x"]
spatial$meta.data$y <- coords[rownames(spatial$meta.data), "y"]

.write_dataset("stxBrain", spatial, c(
  "# stxBrain: 10x Visium mouse brain, sagittal-anterior section 1",
  "",
  "## Source",
  "",
  "10x Genomics, \"Mouse Brain Serial Section 1 (Sagittal-Anterior)\", Space",
  "Ranger 1.0.0. Taken from the `anterior1` object in the SeuratData package",
  sprintf("`stxBrain.SeuratData` %s (Satija Lab).",
          utils::packageVersion("stxBrain.SeuratData")),
  "",
  "## License",
  "",
  "CC BY 4.0, https://creativecommons.org/licenses/by/4.0/",
  "Credit: 10x Genomics.",
  "",
  "## Changes",
  "",
  "- Kept the raw counts of the `Spatial` assay (genes x spots) and the spot",
  "  metadata.",
  "- Added spot coordinates `x` and `y` from `GetTissueCoordinates()`.",
  "- Saved them as a list with `count.mtx` (sparse matrix) and `meta.data`",
  "  (data frame) for NICHESv2. No values were changed."
))


## ifnb: interferon-beta stimulated and control PBMCs ----

data("ifnb")
nonspatial <- extract_NICHESInputs_Seurat(
  seurat.obj = ifnb,
  assay      = "RNA",
  layer      = "counts",
  verbose    = FALSE
)

.write_dataset("ifnb", nonspatial, c(
  "# ifnb: interferon-beta stimulated and control human PBMCs",
  "",
  "## Source",
  "",
  "Kang HM, Subramaniam M, Targ S, et al. Multiplexed droplet single-cell",
  "RNA-sequencing using natural genetic variation. Nat Biotechnol.",
  "2018;36(1):89-94. doi:10.1038/nbt.4042. GEO: GSE96583.",
  "Taken from the `ifnb` object in the SeuratData package",
  sprintf("`ifnb.SeuratData` %s (Satija Lab), which adds the cell-type labels",
          utils::packageVersion("ifnb.SeuratData")),
  "in `seurat_annotations`.",
  "",
  "## License",
  "",
  "CC BY 4.0, https://creativecommons.org/licenses/by/4.0/, as stated by",
  "`ifnb.SeuratData`. Credit: Kang et al. (2018) and the Satija Lab.",
  "",
  "## Changes",
  "",
  "- Kept the raw counts of the `RNA` assay (genes x cells) and the cell",
  "  metadata, including `stim` and `seurat_annotations`.",
  "- Saved them as a list with `count.mtx` (sparse matrix) and `meta.data`",
  "  (data frame) for NICHESv2. No values were changed."
))


## Checksums ----

zips <- file.path(out.dir, c("stxBrain.zip", "ifnb.zip"))
sha <- vapply(zips, function(f) {
  strsplit(system2("shasum", c("-a", "256", shQuote(f)), stdout = TRUE), " ")[[1L]][1L]
}, character(1L))
writeLines(paste0(sha, "  ", basename(zips)), file.path(out.dir, "SHA256SUMS"))

for (f in zips) {
  cat(sprintf("%-12s %6.1f MB  md5 %s\n", basename(f), file.size(f) / 1e6,
              unname(tools::md5sum(f))))
}
