# ── Layer 1: all double-quoted "sampled" → "nonspatial" ──────────────────────
sed -i '' 's/"sampled"/"nonspatial"/g' R/*.R tests/testthat/*.R

# ── Layer 2a: common prose patterns ─────────────────────────────────────────
sed -i '' 's/sampled mode/nonspatial mode/g'   R/*.R tests/testthat/*.R
sed -i '' 's/Sampled mode/Nonspatial mode/g'   R/*.R tests/testthat/*.R
sed -i '' "s/mode = 'sampled'/mode = 'nonspatial'/g" R/*.R tests/testthat/*.R
sed -i '' 's/sampled-mode/nonspatial-mode/g'   R/*.R tests/testthat/*.R

# ── Layer 2b: one remaining single-quoted case in test description ───────────
# test-create_Edgelist_SamplingByType.R:119 has: mode attribute is 'sampled'
sed -i '' "s/is 'sampled'/is 'nonspatial'/g" tests/testthat/*.R

# ── Layer 2c: three targeted lines in source files ───────────────────────────
# compute_CellToCell.R:60 -- "spatial and sampled modes"
sed -i '' 's/and sampled modes/and nonspatial modes/g' R/compute_CellToCell.R

# create_NICHESObject_single.R:165 -- "Spatial and sampled (smooth = FALSE)"
sed -i '' 's/Spatial and sampled (smooth/Spatial and nonspatial (smooth/g' R/create_NICHESObject_single.R

# create_NICHESObject_single.R:166 -- "Sampled (smooth = TRUE): sampled barcode columns"
# Only the first word on this line is a mode reference; "sampled barcode columns" stays
sed -i '' 's/Sampled (smooth = TRUE)/Nonspatial (smooth = TRUE)/g' R/create_NICHESObject_single.R