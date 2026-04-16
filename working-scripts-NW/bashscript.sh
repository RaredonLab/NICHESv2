# 1. Each file has exactly one @keywords internal line
grep -c "^#' @keywords internal" R/extract_neighborhoods.R          # expect: 1
grep -c "^#' @keywords internal" R/create_Edgelist_Spatial.R        # expect: 1
grep -c "^#' @keywords internal" R/create_Edgelist_PseudoBulk.R     # expect: 1
grep -c "^#' @keywords internal" R/create_Edgelist_SamplingByType.R # expect: 1

# 2. Each file has the rationale paragraph (anchor on a unique phrase)
grep -c "low-level primitive" R/extract_neighborhoods.R             # expect: 1
grep -c "low-level primitive" R/create_Edgelist_Spatial.R           # expect: 1
grep -c "low-level primitive" R/create_Edgelist_PseudoBulk.R        # expect: 1
grep -c "low-level primitive" R/create_Edgelist_SamplingByType.R    # expect: 1

# 3. Each file still has @export (function stays callable)
grep -c "^#' @export" R/extract_neighborhoods.R                     # expect: 1
grep -c "^#' @export" R/create_Edgelist_Spatial.R                   # expect: 1
grep -c "^#' @export" R/create_Edgelist_PseudoBulk.R                # expect: 1
grep -c "^#' @export" R/create_Edgelist_SamplingByType.R            # expect: 1

# 4. @keywords internal is positioned immediately before @export in all 4
# (i.e., no @* tag between them)
for f in extract_neighborhoods create_Edgelist_Spatial create_Edgelist_PseudoBulk create_Edgelist_SamplingByType; do
  echo "=== $f ==="
  awk '/^#'"'"' @keywords internal/,/^#'"'"' @export/' R/${f}.R
done
# expect for each: exactly two lines printed -- "#' @keywords internal" then "#' @export"