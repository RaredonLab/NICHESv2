# Helper utilities shared across test files.
# testthat automatically sources any helper-*.R file in tests/testthat/
# before running test blocks, so functions defined here are available
# everywhere without explicit source() calls.

# Returns an n×n regular grid data.frame with unit spacing.
# Row names are unique cell barcodes ("cell_1", "cell_2", ...).
make_grid_coord <- function(n, spacing = 1.0) {
  g <- expand.grid(
    x = seq(0, (n - 1L) * spacing, by = spacing),
    y = seq(0, (n - 1L) * spacing, by = spacing)
  )
  rownames(g) <- paste0("cell_", seq_len(nrow(g)))
  g
}
