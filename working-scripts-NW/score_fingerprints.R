# score_fingerprints.R
#
# PURPOSE: Fingerprint NICHESv2 outputs so code changes can be shown not to
#          alter results. Each table is sorted, rounded to 10 significant
#          digits, written to a temp file, and hashed with MD5. All runs use
#          raw counts + product scoring. Temporary release tool for 1.0.0.
#
# USAGE:   Rscript working-scripts-NW/score_fingerprints.R <label>
#          Writes working-scripts-NW/score_fingerprints_<label>.txt
#          Compare with: diff score_fingerprints_baseline.txt score_fingerprints_<label>.txt

args  <- commandArgs(trailingOnly = TRUE)
label <- if (length(args)) args[[1L]] else "baseline"
repo  <- "/Users/nw426/Documents/GitHub/NICHESv2-NW"
out   <- file.path(repo, "working-scripts-NW",
                   sprintf("score_fingerprints_%s.txt", label))

suppressMessages(devtools::load_all(repo, quiet = TRUE))
suppressPackageStartupMessages({
  library(data.table)
  library(SeuratData)
  library(SeuratObject)
})

fp <- function(dt) {
  dt <- as.data.table(copy(dt))
  num <- names(dt)[vapply(dt, is.numeric, logical(1L))]
  for (col in num) set(dt, j = col, value = sprintf("%.10g", dt[[col]]))
  setcolorder(dt, sort(names(dt)))
  setorderv(dt, names(dt))
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp), add = TRUE)
  fwrite(dt, tmp)
  unname(tools::md5sum(tmp))
}

lines <- character()
rec <- function(dataset, table, dt) {
  lines <<- c(lines, sprintf("%-10s %-24s %10d %s", dataset, table, nrow(dt), fp(dt)))
}

obj_tables <- function(dataset, obj) {
  rec(dataset, "edge.list", obj$edge.list)
  rec(dataset, "edge.data", obj$edge.data)
  rec(dataset, "cell.data", obj$cell.data)
}

## Bundled synthetic data ------------------------------------------------------
data("NICHESv2_inputs", package = "NICHESv2", envir = environment())
b <- NICHESv2_inputs
b.sp <- suppressMessages(create_NICHESObject(
  count.mtx = b$count.mtx, meta.data = b$meta.data, LRM.db = b$LRM.db,
  mode = "spatial", k = 5L, method = "product", n.cores = 2L, verbose = FALSE
))
obj_tables("bundled", b.sp)
b.ag <- suppressMessages(aggregate_NICHESObject(b.sp, cell.type.col = "celltype", verbose = FALSE))
rec("bundled", "neighborhood.edge.agg", b.ag$aggregations$neighborhood.edge.agg)
rec("bundled", "neighborhood.composition", b.ag$aggregations$neighborhood.composition)

b.ns <- suppressMessages(create_NICHESObject(
  count.mtx = b$count.mtx, meta.data = b$meta.data, LRM.db = b$LRM.db,
  mode = "nonspatial", cell.type.col = "celltype", method = "product",
  n.cores = 2L, verbose = FALSE
))
rec("bundled", "nonspatial.edge.list", b.ns$edge.list)
rec("bundled", "nonspatial.edge.data", b.ns$edge.data)

b.ex <- Explore_Radius_Neighborhood(
  data.frame(x = b$meta.data$x, y = b$meta.data$y), show.plot = FALSE
)
rec("bundled", "explore.radius", as.data.table(b.ex))

## ifnb ------------------------------------------------------------------------
data("ifnb", envir = environment())
ifnb <- suppressWarnings(UpdateSeuratObject(ifnb))
i <- suppressMessages(extract_NICHESInputs_Seurat(ifnb, assay = "RNA", layer = "counts", verbose = FALSE))
i.ns <- suppressMessages(create_NICHESObject(
  count.mtx = i$count.mtx, meta.data = i$meta.data, LRM.db = "connectomedb2025",
  species = "human", mode = "nonspatial", sample.col = "stim",
  cell.type.col = "seurat_annotations", method = "product",
  n.cores = 2L, verbose = FALSE
))
obj_tables("ifnb", i.ns)
i.ag <- suppressMessages(aggregate_NICHESObject(i.ns, verbose = FALSE))
rec("ifnb", "neighborhood.edge.agg", i.ag$aggregations$neighborhood.edge.agg)
rm(ifnb, i, i.ns, i.ag); invisible(gc())

## stxBrain --------------------------------------------------------------------
data("anterior1", envir = environment())
a <- suppressWarnings(UpdateSeuratObject(anterior1))
s <- suppressMessages(extract_NICHESInputs_Seurat(a, assay = "Spatial", layer = "counts", verbose = FALSE))
cc <- a@images$anterior1@boundaries$centroids@coords
rownames(cc) <- a@images$anterior1@boundaries$centroids@cells
s$meta.data$x <- cc[rownames(s$meta.data), "x"]
s$meta.data$y <- cc[rownames(s$meta.data), "y"]
s.ex <- Explore_Radius_Neighborhood(
  data.frame(x = s$meta.data$x, y = s$meta.data$y), r.max = 800, show.plot = FALSE
)
rec("stxBrain", "explore.radius", as.data.table(s.ex))
s.sp <- suppressMessages(create_NICHESObject(
  count.mtx = s$count.mtx, meta.data = s$meta.data, LRM.db = "connectomedb2025",
  species = "mouse", mode = "spatial", k = 6L, method = "product",
  n.cores = 2L, verbose = FALSE
))
obj_tables("stxBrain", s.sp)
s.ag <- suppressMessages(aggregate_NICHESObject(s.sp, verbose = FALSE))
rec("stxBrain", "neighborhood.edge.agg", s.ag$aggregations$neighborhood.edge.agg)

writeLines(lines, out)
cat(paste(lines, collapse = "\n"), "\n")
