# WJCD

# Core dependencies
library(Seurat)
library(data.table)
library(Matrix)
library(parallel)
library(stats)

# Modeling package
library(iDESC)      # (will auto-install via devtools if missing)

# Optional helper (for auto-install)
library(devtools)   # only needed if iDESC isn't already installed

#### Function ####

# ------------------------------------------------------------------------------
# run_idesc_xen()
# ------------------------------------------------------------------------------
# Batch runner for iDESC differential expression across one or more metadata
# "levels" (columns) and their "labels" (unique values) with robust logging,
# resumability, file cataloging, and automatic API-compatibility with different
# iDESC versions. Designed for Seurat v5/v4 count access (layer -> slot fallback).
#
# ---- Arguments (detailed) -----------------------------------------------------
# obj                 : Seurat object. Must have obj@meta.data and counts in
#                       `assay_name` (layer="counts" preferred; slot fallback).
#
# subset_col          : (character(1) or NULL) Optional metadata column used to
#                       pre-filter cells (e.g., "TimePoint"). If NULL, no subset.
#
# subset_values       : (vector or NULL) If subset_col is provided, the allowed
#                       values for that column. Ignored when subset_col is NULL.
#
# group_col           : (character(1)) Metadata column with exactly two groups:
#                       one `control_group` and one `experimental_group`.
#                       Internally coerced to factor with levels
#                       (control, experimental) to stabilize coefficient signs.
#
# control_group       : (character(1)) The reference group name in `group_col`.
#
# experimental_group  : (character(1)) The non-reference group name in `group_col`.
#                       Coefficients/log2FC are interpreted as experimental vs control.
#
# subject_col         : (character(1); default "TMA") Per-cell subject ID column
#                       (biological replicate). Used to prune subjects with too
#                       few cells and to enforce a minimum subject count.
#
# ctrl_ids            : (character vector) Subject IDs considered controls.
# case_ids            : (character vector) Subject IDs considered cases.
#                       Only subjects in ctrl_ids ∪ case_ids are included.
#
# assay_name          : (character(1); default "Xenium") Seurat assay to pull counts from.
#
# levels              : (character vector; default c("All")) Names of metadata
#                       columns to stratify by. "All" means one aggregated run on
#                       the filtered set. Any other entry must be a column in meta.
#
# min_cells_per_subject : (integer; default 5) Minimum cells per subject retained
#                         for the current comparison. Subjects with fewer are dropped.
#
# min_subjects_total    : (integer; default 2) Minimum number of distinct subjects
#                         required after pruning; else the comparison is skipped.
#
# out_root            : (character(1); default "iDESC_results") Root output folder.
#                       A date-prefixed subset directory is created inside.
#
# run_tag             : (character(1); default "TMArandom") Free-text tag appended
#                       to filenames to differentiate runs.
#
# skip_on_error       : (logical(1); default TRUE) If TRUE, traps iDESC errors,
#                       logs them, and continues; if FALSE, errors propagate.
#
# resume_after        : (named list or NULL) Per-level label to resume *after*.
#                       Example: list(Gentype = "Lymphoid") → skip all labels up
#                       to and including "Lymphoid", then continue with the next.
#
# completed_labels    : (named list or NULL) Labels already completed per level.
#                       Example: list(Gentype = c("Cap","Epi_Dist")). Those labels
#                       will be skipped with a log note.
#
# skip_labels         : (named list or vector) Labels to always skip:
#                         If a *vector* is supplied, it applies to *all* levels
#                         (including "All").
#                         If a *named list*, names must be level names
#                         (including "All"); each value is a character vector.
#
# ---- Notable Behaviors & Outputs ---------------------------------------------
# Filtering: Cells are restricted to (optional subset) ∩ (group ∈ {control, experimental})
#              ∩ (subject ∈ ctrl_ids ∪ case_ids). Skips if <100 cells remain.
# Group factor ordering: Enforced as (control, experimental) to stabilize the
#   sign of the experimental effect for log2FC extraction.
# iDESC API compatibility: Detects whether iDESC::iDESC expects
#   (expr/meta/subject_col/group_col) or (mat/meta/subject_var/group_var) and
#   calls accordingly, with sensible defaults (no normalization, minimal filters).
# Gene name inference: Tries common gene column names; falls back to rownames
#   or to the first n genes from the counts matrix if needed.
# log2FC inference: Searches coefficient columns for the experimental effect
#   and converts from natural-log scale to log2 by dividing by log(2). Records
#   the source column name in `log2FC_from`.
# FDR computation: If any p-value-like column exists, adjusts with BH ("fdr").
# Output structure (date-prefixed):
#     - RUNLOG: timestamped progress/errors
#     - MANIFEST.csv: one row per attempted (level, label) with status & paths
#     - Results per comparison:
#         * CSV: tidy results with gene, log2FC, FDR, level, label, assay, contrast
#         * RDS: same as CSV
#         * *_qc.rds: list with counts, subjects, per-subject cell counts, etc.
# Return value: Invisibly returns a small “catalog” data.frame for successful
#   comparisons (paths + sample sizes). Skips/errors still appear in the MANIFEST.
#
# ---- Returns -----------------------------------------------------------------
# (Invisibly) a data.frame “catalog” of successful comparisons with file paths
# and sample sizes. All attempts are also logged and summarized in MANIFEST.csv.
#
# ---- Side Effects -------------------------------------------------------------
# Creates date-prefixed folders/files under `out_root`. Writes logs, manifests,
# CSV/RDS results and QC RDS files. May install iDESC from GitHub if missing.
#
# ---- Dependencies -------------------------------------------------------------
# Seurat (GetAssayData, DefaultAssay)
# data.table (fwrite)
# iDESC (from CRAN or GitHub "yl883/iDESC")
# Matrix (dgCMatrix)
#
# ---- Example (pseudo) ---------------------------------------------------------
# cat <- run_idesc_xen(
#   obj = tmp,
#   subset_col = "TimePoint", subset_values = "E21",
#   group_col = "CDHstate", control_group = "N", experimental_group = "CDH",
#   subject_col = "TMA",
#   ctrl_ids = c("E21_N_039","E21_N_041","E21_N_042"),
#   case_ids = c("E21_CDH_048","E21_CDH_049","E21_CDH_051","E21_CDH_052"),
#   assay_name = "Xenium",
#   levels = c("All","Type","Subtype"),
#   run_tag = "myrun",
#   skip_labels = list(Type = c("Unknown"))
# )
# ------------------------------------------------------------------------------

run_idesc_xen <- function(
    obj,
    subset_col      = NULL,
    subset_values   = NULL,
    group_col,
    control_group,
    experimental_group,
    subject_col  = "TMA",
    ctrl_ids,
    case_ids,
    assay_name   = "Xenium",
    levels       = c("All"),     # user may supply any metadata column names here
    min_cells_per_subject = 5,
    min_subjects_total    = 2,
    out_root     = "iDESC_results",
    run_tag      = "TMArandom",
    skip_on_error = TRUE,
    resume_after     = NULL,   # e.g., list(Gentype = "Lymphoid")
    completed_labels = NULL,   # e.g., list(Gentype = c("Cap","Epi_Dist",...))
    skip_labels      = list(   # always skip these labels per-level (keys must match entries in `levels`, or "All")
      All = character(0)
    )
) {
  # ---- helpers ----
  `%||%` <- function(a, b) if (!is.null(a)) a else b   # simple "or" for NULLs
  date_prefix <- format(Sys.Date(), "%Y%m%d")          # YYYYMMDD used for all emitted files/dirs
  
  # ---- argument checks ----
  # Enforce scalar character types and meta presence for required columns
  stopifnot(is.character(group_col), length(group_col) == 1)
  stopifnot(is.character(control_group), length(control_group) == 1)
  stopifnot(is.character(experimental_group), length(experimental_group) == 1)
  stopifnot(is.character(subject_col), length(subject_col) == 1)
  stopifnot(is.character(assay_name), length(assay_name) == 1)
  stopifnot(all(c(group_col, subject_col) %in% colnames(obj@meta.data)))
  # Optional subset: ensure column exists and values provided
  if (!is.null(subset_col)) {
    stopifnot(is.character(subset_col), length(subset_col) == 1)
    stopifnot(subset_col %in% colnames(obj@meta.data))
    stopifnot(length(subset_values) >= 1)
  }
  # Ensure all requested levels are valid (or "All")
  stopifnot(all(levels %in% c("All", colnames(obj@meta.data))))
  
  # ---- strip images/tools for iDESC compatibility (preprocessing) ----
  # iDESC only needs the counts matrix and metadata.
  # Large image payloads (Xenium/Visium) can cause serialization issues or
  # unnecessary memory usage, so drop them upfront in a safe/defensive way.
  suppressWarnings({
    # Drop any spatial images (Seurat stores them in @images).
    if (isS4(obj) && "images" %in% slotNames(obj)) {
      if (length(obj@images)) obj@images <- list()
    }
    # Drop Xenium tool payload (often stored at @tools$Xenium).
    if (!is.null(obj@tools) && !is.null(obj@tools[["Xenium"]])) {
      obj@tools[["Xenium"]] <- NULL
    }
  })
  
  # ---- subset to requested subset_col ∈ subset_values and group ∈ {control, experimental} ----
  DefaultAssay(obj) <- assay_name
  meta0 <- obj@meta.data
  
  # Boolean masks for subset and group membership
  in_subset <- if (is.null(subset_col)) rep(TRUE, nrow(meta0)) else meta0[[subset_col]] %in% subset_values
  in_groups <- meta0[[group_col]] %in% c(control_group, experimental_group)
  
  # Hard stop if no cells survive filters
  if (!any(in_subset & in_groups)) {
    if (is.null(subset_col)) {
      stop(sprintf("No cells after filtering: %s in {%s,%s}.",
                   group_col, control_group, experimental_group))
    } else {
      stop(sprintf("No cells after filtering: %s in {%s} and %s in {%s,%s}.",
                   subset_col, paste(subset_values, collapse = ","),
                   group_col, control_group, experimental_group))
    }
  }
  
  # Materialize subset
  keep_cells_ids <- rownames(meta0)[in_subset & in_groups]
  obj_sub <- subset(obj, cells = keep_cells_ids)
  
  # ---- defensive strip on subsetted object (in case future changes reintroduce payloads) ----
  suppressWarnings({
    if (isS4(obj_sub) && "images" %in% slotNames(obj_sub)) {
      if (length(obj_sub@images)) obj_sub@images <- list()
    }
    if (!is.null(obj_sub@tools) && !is.null(obj_sub@tools[["Xenium"]])) {
      obj_sub@tools[["Xenium"]] <- NULL
    }
  })
  
  # ---- meta and IDs ----
  meta  <- obj_sub@meta.data
  # Enforce ordered factor so modeled coefficient aligns with experimental - control
  meta[[group_col]] <- factor(as.character(meta[[group_col]]),
                              levels = c(control_group, experimental_group))
  # Subjects we will allow into comparisons
  all_ids  <- c(ctrl_ids, case_ids)
  
  # ---- output dirs for this subset (date-prefixed everywhere) ----
  subset_tag <- if (is.null(subset_col)) "ALL" else paste0(subset_col, "=", paste(subset_values, collapse = ","))
  dated_subset_dir <- paste0(date_prefix, "_", subset_tag)
  subset_root <- file.path(out_root, dated_subset_dir)
  dir.create(subset_root, recursive = TRUE, showWarnings = FALSE)
  
  # Set up run log and manifest paths
  log_path <- file.path(subset_root, sprintf("%s_RUNLOG_%s.txt", date_prefix, subset_tag))
  manifest_path <- file.path(subset_root, sprintf("%s_MANIFEST_%s.csv", date_prefix, subset_tag))
  .log <- function(...) {
    # Timestamped console + file log
    msg <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", sprintf(...))
    cat(msg, "\n"); write(msg, file = log_path, append = TRUE)
  }
  write("", file = log_path) # reset log at start
  
  # ---- skip_lists, completed, and label utilities ----
  # Normalize `skip_labels` so that a simple vector applies to all levels.
  normalize_skip <- function(x) {
    if (is.null(x)) return(list())
    if (is.list(x)) return(x)
    # if a vector was provided, apply to all requested levels including "All"
    out <- list(All = x)
    for (lv in setdiff(levels, "All")) out[[lv]] <- x
    out
  }
  skip_labels <- normalize_skip(skip_labels)
  
  # Helpers to read completed or always-skip sets for a given level
  get_completed_for_level <- function(level) {
    if (is.null(completed_labels)) return(character(0))
    if (is.null(completed_labels[[level]])) return(character(0))
    unique(as.character(stats::na.omit(completed_labels[[level]])))
  }
  get_skip_for_level <- function(level) {
    if (is.null(skip_labels[[level]])) return(character(0))
    unique(as.character(stats::na.omit(skip_labels[[level]])))
  }
  
  # Enumerate labels present for a level among *current subjects and groups*
  get_labels_for_level <- function(level) {
    if (level == "All") return("All")
    col <- level
    if (!col %in% colnames(meta)) return(character(0))
    idx <- (meta[[subject_col]] %in% all_ids) & (meta[[group_col]] %in% c(control_group, experimental_group))
    vals <- meta[idx, col, drop = TRUE]
    unique(stats::na.omit(as.character(vals)))
  }
  
  # ---- storage for manifest rows ----
  mani <- list()
  
  # ---- core runner for a single (level, label) ----
  run_one <- function(level, label) {
    # Select cells for this comparison (subject/group filters + optional label)
    idx <- (meta[[subject_col]] %in% all_ids) & (meta[[group_col]] %in% c(control_group, experimental_group))
    if (level != "All") {
      stopifnot(level %in% colnames(meta))
      idx <- idx & (meta[[level]] == label)
    }
    keep_cells <- rownames(meta)[idx]
    
    # Guard: need a minimum total cell count for iDESC to be meaningful
    if (length(keep_cells) < 100) {
      reason <- sprintf("Too few cells after subject/group filtering: %d < 100.", length(keep_cells))
      .log("[%s|%s] SKIP: %s", level, label, reason)
      return(list(status = "skipped", reason = reason, n_cells = length(keep_cells)))
    }
    
    # counts extraction (Seurat v5 layer → v4 slot fallback)
    DefaultAssay(obj_sub) <- assay_name
    counts_sub <- tryCatch({
      GetAssayData(obj_sub, assay = assay_name, layer = "counts")[, keep_cells, drop = FALSE]
    }, error = function(e) {
      .log("[%s|%s] Falling back to slot='counts': %s", level, label, conditionMessage(e))
      GetAssayData(obj_sub, assay = assay_name, slot = "counts")[, keep_cells, drop = FALSE]
    })
    if (!inherits(counts_sub, "dgCMatrix")) counts_sub <- as(counts_sub, "dgCMatrix")
    
    # Restrict meta to selected cells and re-enforce group factor order
    meta_sub <- meta[keep_cells, , drop = FALSE]
    meta_sub[[group_col]] <- factor(as.character(meta_sub[[group_col]]),
                                    levels = c(control_group, experimental_group))
    
    # Build iDESC meta (row-aligned to count columns)
    idesc_meta <- data.frame(
      subject   = meta_sub[[subject_col]],
      group     = meta_sub[[group_col]],
      cell_type = if (level == "All") "All" else paste(level, label, sep = "="),
      stringsAsFactors = FALSE,
      row.names = rownames(meta_sub)
    )
    
    # Ensure column/row alignment (iDESC requires colnames(mat) == rownames(meta))
    counts_cells <- colnames(counts_sub)
    idesc_meta   <- idesc_meta[counts_cells, , drop = FALSE]
    if (!identical(colnames(counts_sub), rownames(idesc_meta))) {
      reason <- "Alignment failed: colnames(counts_sub) != rownames(idesc_meta)."
      .log("[%s|%s] SKIP: %s", level, label, reason)
      return(list(status = "skipped", reason = reason, n_cells = ncol(counts_sub)))
    }
    
    # Prune subjects with insufficient cells for this comparison
    tab_subj  <- table(idesc_meta$subject)
    keep_subj <- names(tab_subj)[tab_subj >= min_cells_per_subject]
    keep_idx  <- idesc_meta$subject %in% keep_subj
    counts_sub <- counts_sub[, keep_idx, drop = FALSE]
    idesc_meta <- idesc_meta[keep_idx, , drop = FALSE]
    
    # QC sanity checks: need both groups and enough subjects
    if (length(unique(idesc_meta$group)) < 2) {
      reason <- sprintf("Need both %s and %s after filtering.", control_group, experimental_group)
      .log("[%s|%s] SKIP: %s", level, label, reason)
      return(list(status = "skipped", reason = reason, n_cells = ncol(counts_sub)))
    }
    if (length(unique(idesc_meta$subject)) < min_subjects_total) {
      reason <- sprintf("Need at least %d subjects after filtering.", min_subjects_total)
      .log("[%s|%s] SKIP: %s", level, label, reason)
      return(list(status = "skipped", reason = reason, n_cells = ncol(counts_sub)))
    }
    
    # ---- iDESC call + annotation wrapped to skip on error ----
    attempt <- function() {
      # Lazily install iDESC if not found (GitHub source)
      if (!requireNamespace("iDESC", quietly = TRUE)) {
        if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
        devtools::install_github("yl883/iDESC")
      }
      library(iDESC)
      
      # Output directory structure:
      # {out_root}/{YYYYMMDD_subsetTag}/{YYYYMMDD_level=Level}/{YYYYMMDD_label=Label or same-as-level-for-All}
      lev_dir_name <- paste0(date_prefix, "_level=", level)
      lev_dir <- file.path(subset_root, lev_dir_name)
      dir.create(lev_dir, recursive = TRUE, showWarnings = FALSE)
      
      lab_dir_name <- if (level == "All") lev_dir_name else paste0(date_prefix, "_label=", gsub("[/\\\\]", "_", label))
      lab_dir <- if (level == "All") lev_dir else file.path(lev_dir, lab_dir_name)
      if (level != "All") dir.create(lab_dir, recursive = TRUE, showWarnings = FALSE)
      
      # Determine cores: prefer SLURM_CPUS_PER_TASK, else all cores
      ncores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = parallel::detectCores()))
      ncores <- max(1, ncores)
      
      # Detect iDESC API signature and call appropriately
      fargs <- names(formals(iDESC::iDESC))
      res <- NULL
      
      if (all(c("expr","meta","subject_col","group_col") %in% fargs)) {
        # Newer signature: sparse counts OK
        res <- iDESC::iDESC(
          expr = counts_sub,
          meta = idesc_meta[, c("subject","group","cell_type")],
          subject_col = "subject",
          group_col   = "group",
          min_cell_per_sub = min_cells_per_subject,
          n.cores = ncores
        )
      } else if (all(c("mat","meta","subject_var","group_var") %in% fargs)) {
        # Older signature: expects dense matrix
        mat_use <- if (inherits(counts_sub, "dgCMatrix")) as.matrix(counts_sub) else counts_sub
        res <- iDESC::iDESC(
          mat                 = mat_use,
          meta                = idesc_meta,  # rownames(meta) == colnames(mat)
          subject_var         = "subject",
          group_var           = "group",
          norm_opt            = "None",
          sub_cell_filtering  = min_cells_per_subject,
          gene_sub_filtering  = 0,
          gene_cell_filtering = 0,
          ncell_filtering     = 0
        )
      } else {
        stop("Unrecognized iDESC::iDESC signature. Please run names(formals(iDESC::iDESC)).")
      }
      
      # ---- Safe annotation of genes + log2FC (+ FDR) ----
      res_df <- as.data.frame(res)
      
      # Gene column inference: pick a common name, else rownames fallback
      gcol <- intersect(c("gene","Gene","GENE","Feature","feature","symbol","Symbol"), names(res_df))
      if (length(gcol)) {
        res_df$gene <- as.character(res_df[[gcol[1]]])
      } else if (!is.null(rownames(res_df)) && all(nzchar(rownames(res_df)))) {
        res_df$gene <- rownames(res_df)
      } else {
        # If all else fails, synthesize from counts (keeps row counts aligned)
        res_df$gene <- rownames(counts_sub)[seq_len(nrow(res_df))]
      }
      
      # log2FC: derive from coefficient column matching the experimental effect
      bcols <- grep("^(Beta_|beta|Estimate|coef|Coef)", names(res_df), value = TRUE, ignore.case = TRUE)
      log2fc_from <- NA_character_
      if (length(bcols)) {
        # Prefer columns that contain the experimental group name
        exp_pat <- paste0("(", experimental_group, ")")
        pick <- grep(exp_pat, bcols, ignore.case = TRUE, value = TRUE)
        if (length(pick) == 0) {
          pick <- grep(paste0("(group|condition|treatment|", experimental_group, ")"),
                       bcols, ignore.case = TRUE, value = TRUE)
        }
        if (length(pick) == 0) pick <- bcols
        pick <- pick[1]
        log2fc_from <- pick
        suppressWarnings(res_df$log2FC <- res_df[[pick]] / log(2))  # convert ln to log2
      }
      res_df$log2FC_from <- log2fc_from
      
      # FDR: BH-adjust any available p-like column; prefer one matching the effect
      pcols <- grep("(Pval_|pval|^p$|Pr\\(>|<=?)", names(res_df), value = TRUE, ignore.case = TRUE)
      if (length(pcols)) {
        pickp <- grep(experimental_group, pcols, ignore.case = TRUE, value = TRUE)
        if (length(pickp) == 0) pickp <- pcols
        pickp <- pickp[1]
        suppressWarnings(res_df$FDR <- p.adjust(res_df[[pickp]], method = "fdr"))
      }
      
      # Provenance fields
      res_df <- res_df[, c("gene", setdiff(names(res_df), "gene"))]
      if (!is.null(subset_col)) res_df[[subset_col]] <- paste(subset_values, collapse = ",")
      res_df$level     <- level
      res_df$label     <- label
      res_df$assay     <- assay_name
      res_df$contrast  <- sprintf("%s_vs_%s", control_group, experimental_group)
      
      # Filenames (date-prefixed)
      label_tag <- if (level == "All") "All" else paste(level, label, sep = "=")
      prefix <- sprintf("%s_%s_%s_%s", subset_tag, label_tag, res_df$contrast[1], run_tag)
      
      rds_path <- file.path(lab_dir, paste0(date_prefix, "_iDESC_", prefix, ".rds"))
      csv_path <- file.path(lab_dir, paste0(date_prefix, "_iDESC_", prefix, ".csv"))
      qc_path  <- file.path(lab_dir, paste0(date_prefix, "_iDESC_", prefix, "_qc.rds"))
      
      # Persist results + QC
      saveRDS(res_df, rds_path)
      data.table::fwrite(res_df, csv_path)
      qc <- list(
        subset_col         = subset_col,
        subset_values      = subset_values,
        level              = level,
        label              = label,
        class_label        = if (level == "All") "All" else label_tag,
        assay_name         = assay_name,
        group_col          = group_col,
        control_group      = control_group,
        experimental_group = experimental_group,
        n_cells            = ncol(counts_sub),
        n_genes_attempted  = nrow(counts_sub),
        n_subjects         = length(unique(idesc_meta$subject)),
        subjects           = sort(unique(idesc_meta$subject)),
        group_tab          = table(idesc_meta$group),
        cells_per_subject  = sort(table(idesc_meta$subject), decreasing = TRUE),
        log2FC_from        = log2fc_from
      )
      saveRDS(qc, qc_path)
      
      .log("[%s|%s] OK: cells=%d, subjects=%d → %s",
           level, label, qc$n_cells, qc$n_subjects, basename(csv_path))
      
      # Summarize for manifest/catalog
      list(
        status  = "ok",
        reason  = NA_character_,
        n_cells = qc$n_cells,
        n_subj  = qc$n_subjects,
        paths   = list(csv = csv_path, rds = rds_path, qc = qc_path)
      )
    }
    
    # Error policy: skip vs. strict
    if (isTRUE(skip_on_error)) {
      out <- tryCatch(
        attempt(),
        error = function(e) {
          reason <- paste0("iDESC error: ", conditionMessage(e))
          .log("[%s|%s] SKIP (error): %s", level, label, reason)
          return(list(status = "skipped", reason = reason,
                      n_cells = ncol(counts_sub),
                      n_subj = length(unique(idesc_meta$subject)),
                      paths = NULL))
        }
      )
      return(out)
    } else {
      return(attempt())  # strict mode: propagate errors
    }
  }
  
  # ---- iterate levels and labels with resume/completion/skip controls ----
  catalog <- list()  # successful runs only
  for (lev in levels) {
    # Ensure per-level directory exists (date-prefixed)
    lev_dir <- file.path(subset_root, paste0(date_prefix, "_level=", lev))
    dir.create(lev_dir, showWarnings = FALSE, recursive = TRUE)
    
    # Determine labels to run for this level
    labs <- get_labels_for_level(lev)
    if (!length(labs)) labs <- if (lev == "All") "All" else character(0)
    if (!length(labs)) { .log("[Level=%s] No labels present; skipping level.", lev); next }
    
    # Resume logic: if resume_after[[lev]] is set, skip up to and including it
    resume_target <- if (!is.null(resume_after) && !is.null(resume_after[[lev]])) as.character(resume_after[[lev]]) else NULL
    resume_passed <- is.null(resume_target)
    
    # Completed/skip lists for this level
    already_done <- get_completed_for_level(lev)
    always_skip  <- get_skip_for_level(lev)
    
    for (lab in labs) {
      .log("Starting [%s|%s] ...", lev, lab)
      
      # RESUME gate
      if (!resume_passed) {
        if (identical(lab, resume_target)) {
          .log("[%s|%s] RESUME: target reached; will start with the next label.", lev, lab)
          resume_passed <- TRUE
          next
        } else {
          .log("[%s|%s] RESUME: before target (%s); skipping.", lev, lab, resume_target)
          next
        }
      }
      
      # Completed or always-skip labels
      if (lab %in% already_done) {
        .log("[%s|%s] SKIP: already completed.", lev, lab); next
      }
      if (lab %in% always_skip || is.na(lab) || identical(lab, "")) {
        .log("[%s|%s] SKIP: label in skip list or NA/empty.", lev, lab); next
      }
      
      # Run one comparison and record in manifest
      res <- run_one(lev, lab)
      
      mani[[length(mani) + 1]] <- data.frame(
        subset_tag   = subset_tag,
        level        = lev,
        label        = lab,
        status       = res$status,
        skipped_reason = ifelse(is.null(res$reason), NA, res$reason),
        n_cells      = res$n_cells %||% NA_integer_,
        n_subjects   = res$n_subj %||% NA_integer_,
        results_csv  = if (!is.null(res$paths)) res$paths$csv else NA_character_,
        results_rds  = if (!is.null(res$paths)) res$paths$rds else NA_character_,
        qc_rds       = if (!is.null(res$paths)) res$paths$qc  else NA_character_,
        stringsAsFactors = FALSE
      )
      
      # Build the success catalog (used as function return)
      if (identical(res$status, "ok")) {
        catalog[[length(catalog) + 1]] <- data.frame(
          subset_tag  = subset_tag,
          level       = lev,
          label       = lab,
          results_csv = res$paths$csv,
          results_rds = res$paths$rds,
          qc_rds      = res$paths$qc,
          n_cells     = res$n_cells,
          n_subjects  = res$n_subj,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  
  # ---- write manifest and return catalog (date-prefixed files) ----
  mani_df <- if (length(mani)) data.table::rbindlist(mani, fill = TRUE) else
    data.frame(subset_tag=character(),level=character(),label=character(),
               status=character(),skipped_reason=character(),n_cells=integer(),
               n_subjects=integer(),results_csv=character(),results_rds=character(),
               qc_rds=character(), stringsAsFactors = FALSE)
  data.table::fwrite(mani_df, manifest_path)
  .log("Wrote MANIFEST: %s", manifest_path)
  
  cat_df <- if (length(catalog)) data.table::rbindlist(catalog, fill = TRUE) else
    data.frame(subset_tag=character(),level=character(),label=character(),
               results_csv=character(),results_rds=character(),qc_rds=character(),
               n_cells=integer(),n_subjects=integer(), stringsAsFactors = FALSE)
  
  invisible(cat_df)  # return only successes; consult MANIFEST for full attempt log
}

#### Usage Examples (CDH Dataset) ####

# ----------------------------- iDESC E21: Neuron Class -----------------------------

# subset allows only Class == Neuron cells to be compared
# selection of FOVs (ctrl_ids / case_ids) uses only E21 pucks

idesc_E21_neuron <- run_idesc_xen(
  obj                = tmp,                     # your Seurat object
  subset_col         = "Class",                 # subset by Class
  subset_values      = "Neuron",                # only include Class == "Neuron"
  group_col          = "CDHstate",
  control_group      = "Norm",
  experimental_group = "CDH",
  subject_col        = "TMA",                   # default
  ctrl_ids           = c("E21_N_039", "E21_N_041", "E21_N_042"),
  case_ids           = c("E21_CDH_048", "E21_CDH_049", "E21_CDH_051", "E21_CDH_052"),
  assay_name         = "Xenium",                # default
  levels             = c("All"),                # compare across ctrl/exp groups for all cells in subset
  min_cells_per_subject = 5,                    # default
  min_subjects_total    = 2,                    # default
  out_root             = "iDESC_results",       # default output directory
  run_tag              = "E21_Neuron",          # tag to identify this run
  skip_on_error        = TRUE                   # continue past errors; details logged in RUNLOG
)


# ----------------------------- iDESC E19: All, Class, Type -----------------------------
# subset limits analysis to E19 cells only
# level selections will result in comparisons at tissue scale, for every class,
#   and for every type

idesc_E19_all <- run_idesc_grid(
  obj                = tmp,                     # your Seurat object
  subset_col         = "TimePoint",             # subset by embryonic timepoint
  subset_values      = "E19",                   # only include E19 cells
  group_col          = "CDHstate",
  control_group      = "Norm",                  
  experimental_group = "CDH",
  subject_col        = "TMA",                   # FOV / subject identity
  ctrl_ids           = c("E19_N_030","E19_N_031","E19_N_032","E19_N_034"),
  case_ids           = c("E19_CDH_053","E19_CDH_054","E19_CDH_056","E19_CDH_057"),
  assay_name         = "Xenium",                # default
  levels             = c("All","Class","Type"), # compare across All cells, then per-Class, then per-Type
  min_cells_per_subject = 5,                    # default
  min_subjects_total    = 2,                    # default
  out_root             = "iDESC_results",       # output directory root
  run_tag              = "E19_All_Class_Type",  # tag to identify this run
  skip_on_error        = TRUE                   # continue past errors; details logged in RUNLOG
)

