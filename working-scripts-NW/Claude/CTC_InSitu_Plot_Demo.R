# ==============================================================================
# SCRIPT:  CTC_InSitu_Plot_Demo.R
# PURPOSE: Demonstrate usage of CellToCell_InSitu() for batch plotting and
#          CellToCell_InSitu_Single_LRM() for single-pair interactive plotting.
# ==============================================================================

library(NICHESv2)

# ── Assumed inputs ─────────────────────────────────────────────────────────────
# niches.obj  — a NICHESObject built by create_NICHESObject()
# PLOT.DIR    — character path to output folder

# ==============================================================================
# DEMO 1 — Batch plotting with CellToCell_InSitu()
#
# Plots a selected set of LRM pairs in parallel and saves one PNG per pair.
# Adjust n.cores to the number of physical cores available on your machine.
# ==============================================================================

my.LRM.list <- c("MIF|CD74", "VEGFA|KDR", "CXCL12|CXCR4")

result <- CellToCell_InSitu(
  niches.obj     = niches.obj,
  LRM.list       = my.LRM.list,   # NULL to plot all pairs in $LRM.list
  folder         = PLOT.DIR,
  bg.color       = "grey40",
  bg.linewidth   = 0.05,
  bg.alpha       = 0.1,
  sig.linewidth  = 0.1,
  sig.alpha      = 0.8,
  viridis.option = "C",
  width          = 35,
  height         = 25,
  dpi            = 800,
  n.cores        = 4L
)

# result is a data.frame with columns: LRM | status | path
print(result)

# ==============================================================================
# DEMO 2 — Single-pair interactive plotting with CellToCell_InSitu_Single_LRM()
#
# Returns a ggplot object for immediate inspection or manual saving.
# No parallelism — intended for exploratory use of one pair at a time.
# ==============================================================================

p <- CellToCell_InSitu_Single_LRM(
  niches.obj     = niches.obj,
  LRM            = "MIF|CD74",
  value.col      = "score",       # "score" (raw) or "score.norm" (per-edge)
  scale.upper    = NULL,          # NULL = auto 90th percentile of non-zero values
  bg.color       = "grey40",
  bg.linewidth   = 0.05,
  bg.alpha       = 0.1,
  sig.linewidth  = 0.1,
  sig.alpha      = 0.8,
  viridis.option = "C",
  title          = "MIF → CD74 Connectivity"
)

# View interactively
print(p)

# Save manually if desired
ggplot2::ggsave(
  filename = file.path(PLOT.DIR, "MIF_CD74.png"),
  plot     = p,
  width    = 35,
  height   = 25,
  dpi      = 800,
  bg       = "black"
)
