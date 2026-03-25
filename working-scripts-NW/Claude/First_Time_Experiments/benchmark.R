# ==============================================================================
# FILE:    benchmark.R
# PURPOSE: Benchmark computeCellToCellSpatial() + summarizeCellSignals()
#          against NICHES::RunNICHES() across datasets of varying size.
#
# STRUCTURE:
#   Section 0  — configuration
#   Section 1  — helper utilities (timing, error reporting, table printing)
#   Section 2  — main benchmark loop (NICHES vs. ours, per dataset)
#   Section 3  — parallelization scaling test (our function only)
#   Section 4  — final consolidated report
#
# TIMING BREAKDOWN — our pipeline:
#   T1  extract   : pull count matrix + coords from Seurat object
#   T2  edgelist  : createEdgelistSpatial()
#   T3  compute   : computeCellToCellSpatial()
#   T4  summarize : summarizeCellSignals()
#   T_total = T1 + T2 + T3 + T4
#
# TIMING — NICHES:
#   Single wall-clock time for RunNICHES() (it does everything internally)
#
# ERROR HANDLING:
#   NICHES errors (including disk/memory) are caught, logged, and skipped.
#   The script always continues to the next dataset.
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
  library(NICHES)
})

source("createEdgelistSpatial.R")
source("computeCellToCellSpatial.R")
source("summarizeCellSignals.R")

# ==============================================================================
# SECTION 0 ── Configuration
# ==============================================================================

PSEUDO_DIR  <- "/Volumes/T5 EVO/NICHESv2/pseudo_data"     # where generate_pseudo_data.R saved .rds files
LRM_PATH    <- "~/Library/CloudStorage/OneDrive-YaleUniversity/Lung_Tissue_CK/Connectome2025.LRMs.human.csv"

# Spatial neighborhood parameters (same for NICHES and our functions)
BENCH_RAD   <- 50     # radius mode  — used in rad.set / rad benchmarks
BENCH_K     <- 5      # kNN mode     — used in k / kNN benchmarks
BENCH_MODE  <- "rad"  # "rad" or "knn" — which mode to use for main benchmark

# Parallelization test: cores to test
PARALLEL_CORES <- c(1L, 2L, 4L)
# Datasets large enough to show parallelization benefit
# (run parallel test only if n_genes * n_cells >= this threshold)
PARALLEL_THRESHOLD <- 5000L * 1000L   # 5M cells*genes

# Report output file — saved in the same directory as this script.
# A timestamp is appended so repeated runs never overwrite each other.
# Set to NULL to skip saving.
REPORT_FILE <- sprintf("benchmark_report_%s.txt",
                       format(Sys.time(), "%Y%m%d_%H%M%S"))

# ==============================================================================
# SECTION 1 ── Helper utilities
# ==============================================================================

# ── Precise wall-clock timer ──────────────────────────────────────────────────
.tik  <- function() proc.time()[["elapsed"]]
.tok  <- function(t0) round(proc.time()[["elapsed"]] - t0, 3)

# ── Disk space check (in GB) ─────────────────────────────────────────────────
.free_disk_gb <- function() {
  info <- tryCatch(
    as.numeric(system("df -Pk . | tail -1 | awk '{print $4}'",
                      intern = TRUE)) / (1024^2),
    error = function(e) NA_real_
  )
  info
}

# ── Format seconds nicely ─────────────────────────────────────────────────────
.fmt_sec <- function(s) {
  if (is.na(s)) return("  ERROR ")
  if (s < 60)   return(sprintf("%6.2fs", s))
  return(sprintf("%5.1fmin", s / 60))
}

# ── Right-pad string ──────────────────────────────────────────────────────────
.rpad <- function(s, w) formatC(as.character(s), width = -w, flag = "-")

# ── Print section banner ──────────────────────────────────────────────────────
.banner <- function(txt) {
  cat("\n", strrep("=", 70), "\n", sep = "")
  cat(" ", txt, "\n")
  cat(strrep("=", 70), "\n\n", sep = "")
}

# ── Structured error reporter ─────────────────────────────────────────────────
# Builds a list with enough context to reconstruct what went wrong later.
.make_error_report <- function(tag, step, err, disk_before, disk_after) {
  list(
    dataset      = tag,
    step         = step,
    message      = conditionMessage(err),
    class        = paste(class(err), collapse = ", "),
    disk_before  = disk_before,
    disk_after   = disk_after,
    disk_delta   = if (!any(is.na(c(disk_before, disk_after))))
                     round(disk_before - disk_after, 2) else NA,
    timestamp    = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  )
}

# ── Print a single error report entry ────────────────────────────────────────
.print_error <- function(er) {
  cat(sprintf("  [ERROR] Dataset : %s\n",  er$dataset))
  cat(sprintf("          Step    : %s\n",  er$step))
  cat(sprintf("          Message : %s\n",  er$message))
  cat(sprintf("          Class   : %s\n",  er$class))
  if (!is.na(er$disk_delta))
    cat(sprintf("          Disk    : %.2f GB used during step\n", er$disk_delta))
  cat("\n")
}

# ==============================================================================
# SECTION 2 ── Load shared resources
# ==============================================================================

.banner("Loading shared resources")

# Ligand-receptor database
message("Loading LRM database from: ", LRM_PATH)
LRM.database <- read.csv(LRM_PATH)
message(sprintf("  %d LR pairs loaded.", nrow(LRM.database)))

# Discover datasets
rds_files <- sort(list.files(PSEUDO_DIR, pattern = "\\.rds$", full.names = TRUE))
if (length(rds_files) == 0)
  stop("No .rds files found in '", PSEUDO_DIR, "'. Run generate_pseudo_data.R first.")

# Parse n_genes / n_cells from filenames (format: pseudo_gXXXXX_cXXXXX.rds)
parsed <- regmatches(basename(rds_files),
                     regexpr("g(\\d+)_c(\\d+)", basename(rds_files)))
dataset_info <- data.frame(
  path   = rds_files,
  tag    = sub("\\.rds$", "", basename(rds_files)),
  n_genes = as.integer(sub("g(\\d+)_c.*",   "\\1", parsed)),
  n_cells = as.integer(sub(".*_c(\\d+)",    "\\1", parsed)),
  stringsAsFactors = FALSE
)
dataset_info <- dataset_info[order(dataset_info$n_genes * dataset_info$n_cells), ]
rownames(dataset_info) <- NULL

message(sprintf("  Found %d datasets to benchmark.", nrow(dataset_info)))
print(dataset_info[, c("tag", "n_genes", "n_cells")])

# ==============================================================================
# SECTION 3 ── Main benchmark loop
# ==============================================================================

.banner("Main Benchmark: NICHES vs. Our Pipeline")
cat(sprintf("  Mode: %s | rad = %d | k = %d\n\n", BENCH_MODE, BENCH_RAD, BENCH_K))

# Results storage
results      <- vector("list", nrow(dataset_info))
error_log    <- list()

for (i in seq_len(nrow(dataset_info))) {

  di  <- dataset_info[i, ]
  tag <- di$tag
  cat(sprintf("── [%d/%d] %s  (%d genes × %d cells)\n",
              i, nrow(dataset_info), tag, di$n_genes, di$n_cells))

  # ── Load Seurat object ──────────────────────────────────────────────────────
  obj <- readRDS(di$path)

  # Row to accumulate results for this dataset
  row <- list(
    tag       = tag,
    n_genes   = di$n_genes,
    n_cells   = di$n_cells,
    niches_t  = NA_real_,
    niches_ok = FALSE,
    our_t1_extract   = NA_real_,
    our_t2_edgelist  = NA_real_,
    our_t3_compute   = NA_real_,
    our_t4_summarize = NA_real_,
    our_t_total      = NA_real_,
    our_ok           = FALSE,
    n_edges   = NA_integer_,
    n_nonzero = NA_integer_
  )

  # ════════════════════════════════════════════════════════════════════════════
  # 3a. NICHES benchmark
  # ════════════════════════════════════════════════════════════════════════════
  cat("    [NICHES]  ... ")
  disk_before <- .free_disk_gb()

  niches_result <- tryCatch({

    t0 <- .tik()
    niches_out <- RunNICHES(
      object            = obj,
      assay             = "RNA",
      species           = "human",
      LR.database       = "custom",
      custom_LR_database = LRM.database,
      cell_types        = NULL,
      position.x        = "x",
      position.y        = "y",
      k                 = if (BENCH_MODE == "knn") BENCH_K   else NULL,
      rad.set           = if (BENCH_MODE == "rad") BENCH_RAD else NULL,
      meta.data.to.map  = names(obj@meta.data),
      CellToCell           = FALSE,
      CellToSystem         = FALSE,
      SystemToCell         = FALSE,
      CellToCellSpatial    = TRUE,
      CellToNeighborhood   = TRUE,
      NeighborhoodToCell   = TRUE
    )
    elapsed <- .tok(t0)
    list(ok = TRUE, time = elapsed, obj = niches_out)

  }, error = function(e) {
    list(ok = FALSE, time = NA_real_, error = e)
  })

  disk_after <- .free_disk_gb()

  if (niches_result$ok) {
    row$niches_t  <- niches_result$time
    row$niches_ok <- TRUE
    cat(sprintf("done (%s)\n", .fmt_sec(niches_result$time)))
    rm(niches_result)
  } else {
    er <- .make_error_report(tag, "NICHES::RunNICHES",
                             niches_result$error, disk_before, disk_after)
    error_log[[length(error_log) + 1]] <- er
    cat(sprintf("FAILED (%s)\n", er$class))
  }
  gc()

  # ════════════════════════════════════════════════════════════════════════════
  # 3b. Our pipeline benchmark
  # ════════════════════════════════════════════════════════════════════════════
  cat("    [Ours]    ")

  our_result <- tryCatch({

    # T1: Extract data from Seurat object
    cat("extract...")
    t1 <- .tik()
    count_mtx <- GetAssayData(obj, assay = "RNA", layer = "counts")
    coord_df  <- data.frame(
      x = obj@meta.data$x,
      y = obj@meta.data$y,
      row.names = colnames(obj)
    )
    t1_elapsed <- .tok(t1)

    # T2: Build edge list
    cat(" edgelist...")
    t2 <- .tik()
    edge_list <- createEdgelistSpatial(
      coord   = coord_df,
      k       = if (BENCH_MODE == "knn") BENCH_K   else NULL,
      rad     = if (BENCH_MODE == "rad") BENCH_RAD else NULL,
      n_cores = 1L
    )
    t2_elapsed <- .tok(t2)
    n_edges    <- nrow(edge_list)

    # T3: Compute CTC signals (n_cores = 1 for main benchmark)
    cat(" compute...")
    t3 <- .tik()
    ctc <- computeCellToCellSpatial(
      count.mtx    = count_mtx,
      edge.list    = edge_list,
      LRM.list     = LRM.database,
      ligand.col   = "Ligand",
      receptor.col = "Receptor",
      method       = "product",
      n_cores      = 1L
    )
    t3_elapsed <- .tok(t3)

    # T4: Summarize cell signals
    cat(" summarize...")
    t4 <- .tik()
    smry <- summarizeCellSignals(ctc)
    t4_elapsed <- .tok(t4)

    cat(sprintf(" done (%.2fs + %.2fs + %.2fs + %.2fs)\n",
                t1_elapsed, t2_elapsed, t3_elapsed, t4_elapsed))

    list(ok        = TRUE,
         t1        = t1_elapsed,
         t2        = t2_elapsed,
         t3        = t3_elapsed,
         t4        = t4_elapsed,
         n_edges   = n_edges,
         n_nonzero = ctc$params$n_nonzero)

  }, error = function(e) {
    cat(sprintf(" FAILED\n"))
    list(ok = FALSE, error = e)
  })

  if (our_result$ok) {
    row$our_t1_extract   <- our_result$t1
    row$our_t2_edgelist  <- our_result$t2
    row$our_t3_compute   <- our_result$t3
    row$our_t4_summarize <- our_result$t4
    row$our_t_total      <- our_result$t1 + our_result$t2 +
                            our_result$t3 + our_result$t4
    row$our_ok    <- TRUE
    row$n_edges   <- our_result$n_edges
    row$n_nonzero <- our_result$n_nonzero
  } else {
    er <- .make_error_report(tag, "our pipeline",
                             .free_disk_gb(), .free_disk_gb(), our_result$error)
    error_log[[length(error_log) + 1]] <- er
  }

  results[[i]] <- row
  rm(obj, coord_df, count_mtx)
  gc()
  cat("\n")
}

# Flatten results to data.frame
results_df <- rbindlist(lapply(results, as.data.table), fill = TRUE)

# ==============================================================================
# SECTION 4 ── Parallelization scaling test
# ==============================================================================

.banner("Parallelization Scaling Test")

# Select datasets that meet the complexity threshold
para_mask <- (dataset_info$n_genes * dataset_info$n_cells) >= PARALLEL_THRESHOLD
para_sets <- dataset_info[para_mask, ]

if (nrow(para_sets) == 0) {
  cat("  No datasets meet the parallelization threshold — skipping.\n")
} else {
  cat(sprintf("  Testing %d datasets with cores: %s\n\n",
              nrow(para_sets), paste(PARALLEL_CORES, collapse = ", ")))

  para_results <- list()

  for (i in seq_len(nrow(para_sets))) {
    di  <- para_sets[i, ]
    tag <- di$tag
    cat(sprintf("  Dataset: %s (%d genes × %d cells)\n", tag, di$n_genes, di$n_cells))

    obj       <- readRDS(di$path)
    count_mtx <- GetAssayData(obj, assay = "RNA", layer = "counts")
    coord_df  <- data.frame(x = obj@meta.data$x, y = obj@meta.data$y,
                            row.names = colnames(obj))

    edge_list <- createEdgelistSpatial(
      coord   = coord_df,
      k       = if (BENCH_MODE == "knn") BENCH_K   else NULL,
      rad     = if (BENCH_MODE == "rad") BENCH_RAD else NULL,
      n_cores = 1L
    )

    core_times <- numeric(length(PARALLEL_CORES))
    names(core_times) <- paste0("cores_", PARALLEL_CORES)

    for (j in seq_along(PARALLEL_CORES)) {
      nc <- PARALLEL_CORES[j]
      cat(sprintf("    n_cores = %d ... ", nc))

      t_par <- tryCatch({
        t0 <- .tik()
        computeCellToCellSpatial(
          count.mtx    = count_mtx,
          edge.list    = edge_list,
          LRM.list     = LRM.database,
          ligand.col   = "Ligand",
          receptor.col = "Receptor",
          method       = "product",
          n_cores      = nc
        )
        .tok(t0)
      }, error = function(e) {
        cat(sprintf("FAILED (%s)\n", conditionMessage(e)))
        NA_real_
      })

      core_times[j] <- t_par
      if (!is.na(t_par)) cat(sprintf("%s\n", .fmt_sec(t_par)))
    }

    # Speedup relative to serial (1 core)
    speedup <- core_times[1] / core_times
    para_results[[tag]] <- list(
      tag       = tag,
      n_genes   = di$n_genes,
      n_cells   = di$n_cells,
      times     = core_times,
      speedup   = speedup
    )

    rm(obj, count_mtx, coord_df, edge_list)
    gc()
    cat("\n")
  }
}

# ==============================================================================
# SECTION 5 ── Final consolidated report
# ==============================================================================

# Open sink: captures all cat() output below to a text file.
# split = TRUE means output goes to BOTH the file and the console,
# so you can still watch progress in real time.
# message() calls (progress during the loop above) are NOT captured —
# the saved file contains only the clean final report.
if (!is.null(REPORT_FILE)) {
  sink(REPORT_FILE, split = TRUE)
  cat(sprintf("Benchmark Report  |  generated: %s\n",
              format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
  cat(sprintf("Mode: %s  |  rad = %d  |  k = %d\n",
              BENCH_MODE, BENCH_RAD, BENCH_K))
  cat(strrep("=", 70), "\n\n")
}

.banner("BENCHMARK REPORT")

cat(sprintf("  Mode: %s | rad = %d | k = %d\n\n",
            BENCH_MODE, BENCH_RAD, BENCH_K))

# ── 5a. Main timing comparison table ─────────────────────────────────────────
cat("── Timing Comparison (seconds) ──────────────────────────────────────────\n\n")

# Header
hdr <- sprintf("%-26s  %8s  %8s  %8s  %8s  %8s  %10s  %9s",
               "Dataset", "NICHES", "Extract", "Edgelist", "Compute", "Summarize",
               "Our Total", "Speedup")
cat(hdr, "\n")
cat(strrep("-", nchar(hdr)), "\n")

for (i in seq_len(nrow(results_df))) {
  r <- results_df[i, ]

  niches_str  <- if (r$niches_ok) .fmt_sec(r$niches_t)  else "  FAILED"
  our_str     <- if (r$our_ok)    .fmt_sec(r$our_t_total) else "  FAILED"
  t1_str      <- if (r$our_ok)    .fmt_sec(r$our_t1_extract)   else "       -"
  t2_str      <- if (r$our_ok)    .fmt_sec(r$our_t2_edgelist)  else "       -"
  t3_str      <- if (r$our_ok)    .fmt_sec(r$our_t3_compute)   else "       -"
  t4_str      <- if (r$our_ok)    .fmt_sec(r$our_t4_summarize) else "       -"

  speedup_str <- if (r$niches_ok && r$our_ok) {
    ratio <- r$niches_t / r$our_t_total
    sprintf("  %5.2fx", ratio)
  } else "        -"

  cat(sprintf("%-26s  %8s  %8s  %8s  %8s  %8s  %10s  %9s\n",
              r$tag, niches_str, t1_str, t2_str, t3_str, t4_str,
              our_str, speedup_str))
}
cat("\n")

# ── 5b. Edge / signal counts ──────────────────────────────────────────────────
cat("── Dataset Statistics ────────────────────────────────────────────────────\n\n")
cat(sprintf("%-26s  %8s  %8s  %12s\n", "Dataset", "Edges", "Non-zero", "Sparsity %"))
cat(strrep("-", 62), "\n")

for (i in seq_len(nrow(results_df))) {
  r <- results_df[i, ]
  if (!r$our_ok) next
  sparsity <- round(100 * (1 - r$n_nonzero / (r$n_edges * nrow(LRM.database))), 1)
  cat(sprintf("%-26s  %8d  %8d  %11.1f%%\n",
              r$tag, r$n_edges, r$n_nonzero, sparsity))
}
cat("\n")

# ── 5c. Parallelization report ────────────────────────────────────────────────
if (exists("para_results") && length(para_results) > 0) {
  cat("── Parallelization Scaling ───────────────────────────────────────────────\n\n")

  core_header <- paste(sprintf("%9s", paste0(PARALLEL_CORES, "c")), collapse = "")
  spdup_header <- paste(sprintf("%9s", paste0(PARALLEL_CORES, "c (x)")), collapse = "")
  cat(sprintf("%-26s  %s\n", "Dataset", core_header))
  cat(strrep("-", 26 + 2 + length(PARALLEL_CORES) * 9), "\n")

  for (tag in names(para_results)) {
    pr  <- para_results[[tag]]
    ts  <- sprintf("%9s", sapply(pr$times, .fmt_sec))
    cat(sprintf("%-26s  %s\n", tag, paste(ts, collapse = "")))

    sp  <- sprintf("%8.2fx", pr$speedup)
    cat(sprintf("%-26s  %s   (speedup vs 1 core)\n", "", paste(sp, collapse = "")))
    cat("\n")
  }
}

# ── 5d. Error report ──────────────────────────────────────────────────────────
if (length(error_log) > 0) {
  cat(sprintf("── Error Report (%d failures) ────────────────────────────────────────\n\n",
              length(error_log)))
  for (er in error_log) .print_error(er)

  # Disk-related errors flagged separately
  disk_errors <- Filter(function(er) {
    grepl("disk|space|quota|write|cannot open|no space",
          er$message, ignore.case = TRUE) ||
    (!is.na(er$disk_delta) && er$disk_delta > 0.5)
  }, error_log)

  if (length(disk_errors) > 0) {
    cat(sprintf("  ⚠  %d error(s) appear disk/memory related:\n", length(disk_errors)))
    for (er in disk_errors) {
      cat(sprintf("     - %s : disk consumed %.2f GB during step '%s'\n",
                  er$dataset,
                  ifelse(is.na(er$disk_delta), NA, er$disk_delta),
                  er$step))
    }
    cat("\n  Recommendation: NICHES writes large intermediate files during\n")
    cat("  CellToCellSpatial computation. Consider:\n")
    cat("    - Freeing disk space before running\n")
    cat("    - Running NICHES on smaller datasets first\n")
    cat("    - Using a machine with more /tmp space\n\n")
  }
} else {
  cat("── No errors encountered ─────────────────────────────────────────────────\n\n")
}

cat(strrep("=", 70), "\n")
cat("  Benchmark complete.\n")
cat(strrep("=", 70), "\n\n")

# Close sink and confirm
if (!is.null(REPORT_FILE)) {
  sink()
  message(sprintf("[benchmark] Report saved to: %s", REPORT_FILE))
}

# Optionally save results for further analysis
# saveRDS(list(results = results_df, errors = error_log,
#              parallel = para_results), "benchmark_results.rds")
