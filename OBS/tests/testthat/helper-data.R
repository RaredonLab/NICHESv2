# Helper functions for creating test data
# These functions are automatically loaded by testthat before running tests

#' Create a minimal Seurat object for testing
#'
#' @param n_cells Number of cells (default: 50)
#' @param n_genes Number of genes (default: 100)
#' @param n_celltypes Number of cell types (default: 3)
#' @return Seurat object
#' @export
create_test_seurat <- function(n_cells = 50, n_genes = 100, n_celltypes = 3) {
  # Generate fake expression data with some realistic structure
  set.seed(123)  # For reproducibility

  # Create gene names using real-looking gene symbols
  gene_names <- c(
    paste0("GENE", 1:n_genes)
  )

  # Create cell names
  cell_names <- paste0("Cell_", 1:n_cells)

  # Generate counts matrix with some cell type structure
  counts <- matrix(
    rpois(n_genes * n_cells, lambda = 3),
    nrow = n_genes,
    dimnames = list(gene_names, cell_names)
  )

  # Add some differential expression between cell types
  cells_per_type <- n_cells %/% n_celltypes
  for (i in 1:n_celltypes) {
    type_cells <- ((i-1) * cells_per_type + 1):min(i * cells_per_type, n_cells)
    type_genes <- ((i-1) * 10 + 1):min(i * 10, n_genes)
    counts[type_genes, type_cells] <- counts[type_genes, type_cells] * 3
  }

  # Create metadata
  metadata <- data.frame(
    cell_type = rep(
      paste0("CellType", 1:n_celltypes),
      length.out = n_cells
    ),
    sample = "test_sample",
    nCount_RNA = colSums(counts),
    nFeature_RNA = colSums(counts > 0),
    row.names = cell_names
  )

  # Create Seurat object
  obj <- Seurat::CreateSeuratObject(
    counts = counts,
    meta.data = metadata,
    project = "NICHES_test"
  )

  # Set identities
  Seurat::Idents(obj) <- obj$cell_type

  # Normalize
  obj <- Seurat::NormalizeData(obj, verbose = FALSE)

  return(obj)
}


#' Create a minimal spatial Seurat object for testing
#'
#' @param n_cells Number of cells (default: 50)
#' @param n_genes Number of genes (default: 100)
#' @param n_celltypes Number of cell types (default: 3)
#' @return Seurat object with spatial coordinates
#' @export
create_test_spatial_seurat <- function(n_cells = 50, n_genes = 100, n_celltypes = 3) {
  obj <- create_test_seurat(n_cells, n_genes, n_celltypes)

  # Add spatial coordinates
  # Create clustered spatial structure (cell types cluster spatially)
  set.seed(123)
  cells_per_type <- n_cells %/% n_celltypes

  x_coords <- numeric(n_cells)
  y_coords <- numeric(n_cells)

  for (i in 1:n_celltypes) {
    type_cells <- ((i-1) * cells_per_type + 1):min(i * cells_per_type, n_cells)
    # Each cell type clusters around a different center
    center_x <- (i - 1) * 50
    center_y <- (i - 1) * 50
    x_coords[type_cells] <- rnorm(length(type_cells), mean = center_x, sd = 10)
    y_coords[type_cells] <- rnorm(length(type_cells), mean = center_y, sd = 10)
  }

  obj$x <- x_coords
  obj$y <- y_coords

  return(obj)
}


#' Create minimal ground truth (ligand-receptor database)
#'
#' @param n_mechanisms Number of L-R mechanisms (default: 5)
#' @param genes Gene names to use (default: GENE1-GENE20)
#' @param multisubunit Include multi-subunit complexes (default: FALSE)
#' @return List with source.subunits and target.subunits matrices
#' @export
create_test_groundtruth <- function(n_mechanisms = 5, genes = NULL, multisubunit = FALSE) {
  if (is.null(genes)) {
    genes <- paste0("GENE", 1:20)
  }

  mechanism_names <- paste0("Mechanism_", 1:n_mechanisms)

  if (!multisubunit) {
    # Simple single-subunit mechanisms
    source.subunits <- matrix(
      genes[seq(1, n_mechanisms * 2, by = 2)],
      ncol = 1,
      dimnames = list(mechanism_names, "subunit_1")
    )

    target.subunits <- matrix(
      genes[seq(2, n_mechanisms * 2, by = 2)],
      ncol = 1,
      dimnames = list(mechanism_names, "subunit_1")
    )
  } else {
    # Multi-subunit mechanisms
    max_subunits <- 3

    source.subunits <- matrix(
      NA,
      nrow = n_mechanisms,
      ncol = max_subunits,
      dimnames = list(mechanism_names, paste0("subunit_", 1:max_subunits))
    )

    target.subunits <- matrix(
      NA,
      nrow = n_mechanisms,
      ncol = max_subunits,
      dimnames = list(mechanism_names, paste0("subunit_", 1:max_subunits))
    )

    gene_idx <- 1
    for (i in 1:n_mechanisms) {
      # Randomly assign 1-3 subunits per mechanism
      n_source_subunits <- sample(1:3, 1)
      n_target_subunits <- sample(1:3, 1)

      for (j in 1:n_source_subunits) {
        source.subunits[i, j] <- genes[gene_idx]
        gene_idx <- gene_idx + 1
        if (gene_idx > length(genes)) gene_idx <- 1
      }

      for (j in 1:n_target_subunits) {
        target.subunits[i, j] <- genes[gene_idx]
        gene_idx <- gene_idx + 1
        if (gene_idx > length(genes)) gene_idx <- 1
      }
    }
  }

  return(list(
    source.subunits = source.subunits,
    target.subunits = target.subunits
  ))
}


#' Create test expression matrix (no Seurat)
#'
#' @param n_cells Number of cells
#' @param n_genes Number of genes
#' @return Matrix of expression values
#' @export
create_test_matrix <- function(n_cells = 50, n_genes = 100) {
  set.seed(123)
  matrix(
    rpois(n_genes * n_cells, lambda = 3),
    nrow = n_genes,
    dimnames = list(
      paste0("GENE", 1:n_genes),
      paste0("Cell_", 1:n_cells)
    )
  )
}


#' Create test metadata dataframe
#'
#' @param n_cells Number of cells
#' @param n_celltypes Number of cell types
#' @return Data frame with metadata
#' @export
create_test_metadata <- function(n_cells = 50, n_celltypes = 3) {
  data.frame(
    cell_type = rep(
      paste0("CellType", 1:n_celltypes),
      length.out = n_cells
    ),
    sample = "test_sample",
    condition = rep(c("Control", "Treatment"), length.out = n_cells),
    row.names = paste0("Cell_", 1:n_cells)
  )
}


#' Create test edgelist for spatial analysis
#'
#' @param n_cells Number of cells
#' @param k Average neighbors per cell
#' @return Data frame with from/to columns
#' @export
create_test_edgelist <- function(n_cells = 50, k = 5) {
  set.seed(123)
  cell_names <- paste0("Cell_", 1:n_cells)

  # Create k-NN edgelist
  edges <- data.frame(
    from = character(),
    to = character(),
    stringsAsFactors = FALSE
  )

  for (i in 1:n_cells) {
    # Connect to k nearest neighbors (simulated)
    neighbors <- sample(setdiff(1:n_cells, i), min(k, n_cells - 1))
    edges <- rbind(
      edges,
      data.frame(
        from = cell_names[i],
        to = cell_names[neighbors],
        stringsAsFactors = FALSE
      )
    )
  }

  return(edges)
}


#' Load real FANTOM5 subset for integration tests
#'
#' @param species Species (human, mouse, rat, pig)
#' @param n_mechanisms Number of mechanisms to include
#' @return Ground truth list
#' @export
load_real_fantom5_subset <- function(species = "human", n_mechanisms = 10) {
  # Load real FANTOM5 data (if package is installed)
  tryCatch({
    if (species == "human") {
      data("ncomms8866_human", package = "NICHES", envir = environment())
      gt <- ncomms8866_human
    } else if (species == "mouse") {
      data("ncomms8866_mouse", package = "NICHES", envir = environment())
      gt <- ncomms8866_mouse
    } else if (species == "rat") {
      data("ncomms8866_rat", package = "NICHES", envir = environment())
      gt <- ncomms8866_rat
    } else if (species == "pig") {
      data("ncomms8866_pig", package = "NICHES", envir = environment())
      gt <- ncomms8866_pig
    } else {
      stop("Unknown species")
    }

    # Subset to n_mechanisms
    if (nrow(gt$source.subunits) > n_mechanisms) {
      rows_to_keep <- 1:n_mechanisms
      gt$source.subunits <- gt$source.subunits[rows_to_keep, , drop = FALSE]
      gt$target.subunits <- gt$target.subunits[rows_to_keep, , drop = FALSE]
    }

    return(gt)
  }, error = function(e) {
    warning("Could not load real FANTOM5 data. Using mock data instead.")
    return(create_test_groundtruth(n_mechanisms = n_mechanisms))
  })
}


#' Skip test if Seurat is not available
#' @export
skip_if_no_seurat <- function() {
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    testthat::skip("Seurat not available")
  }
}


#' Skip test if running on CI
#' @export
skip_if_ci <- function() {
  if (identical(Sys.getenv("CI"), "true")) {
    testthat::skip("Skipping on CI")
  }
}
