# NICHESv2 — Test Suite Requirements

This document specifies the requirements for a clean rewrite of the
`tests/testthat/` test suite. It is intended as the primary specification
for a dedicated test-rewrite session and does NOT duplicate the main
package design in `NICHESv2_package_design.md` — refer to that document
for API surface, conventions, and NICHESObject structure.


# BACKGROUND

The existing test suite was written during early function development,
before the API stabilized. It works under `devtools::test()` but fails
under `R CMD check` due to a fixture-loading issue in
`tests/testthat/helper-testdata.R` that has resisted targeted fixes.
Rather than continue patching, the plan is to rewrite all test files
from scratch using the conventions below.


# GOALS

1. Clean `R CMD check --as-cran`: zero errors, zero warnings, zero test
   failures.
2. Each `test_that()` block is self-contained — no cross-block state,
   no dependence on test execution order.
3. Fixtures loaded from the bundled Pass-2 datasets
   (`data(NICHESv2_inputs)` / `data(NICHESv2_example)`) rather than
   from `inst/testdata/*.rds`.
4. Eliminate `tests/testthat/helper-testdata.R` and the
   `inst/testdata/` directory entirely once the rewrite is complete.
5. Preserve coverage of the current 373 passing tests; drop obsolete
   tests only with explicit justification.


# FIXTURE STRATEGY

## Primary: bundled data

Every test_that block that needs real NICHESv2 data should load it
inline:

```r
test_that("something about X", {
  data(NICHESv2_example)
  # ... assertions on NICHESv2_example ...
})
```

`data()` is cheap (lazy-loaded, cached by R after first call), so
calling it at the top of every block costs essentially nothing after
the first test in the session.

Use `NICHESv2_inputs` (list: count.mtx, meta.data, LRM.db) for tests
of upstream constructors: `create_NICHESObject`, `compute_CellToCell`,
the three `create_Edgelist_*` functions, `load_LRM_database`,
`compute_PseudoBulk_Connectivity`.

Use `NICHESv2_example` (prebuilt NICHESObject) for tests of downstream
functions: `add_neighborhoods`, `add_segmentation`, `add_cell_meta`,
`Aggregate_CellToCell_Edges`, `CellToCell_InSitu`,
`CellToCell_InSitu_Single_LRM`, `merge_NICHESObjects`, and all S3
methods of NICHESObject.

## Secondary: inline synthesis for edge cases

When a test needs inputs that differ from the bundled data (empty
matrices, all-zero counts, single-cell, barcode collisions, etc.),
synthesise them inline inside the `test_that` block. Keep synthesis
minimal — 6 genes x 8 cells is typically enough. Use `set.seed(1L)`
for reproducibility.

## What to DELETE

- `tests/testthat/helper-testdata.R` — replaced by direct `data()` calls
- `inst/testdata/` directory and all `.rds` files within it
- `inst/testdata/generate_testdata.R` (if no longer referenced)
- The `.setup_fixtures()`, `.load_fixture()`, `.td_path()` helpers
- Any top-of-file or between-block references to `count.mtx`,
  `meta.data`, `LRM.db`, `edge.list`, `niches.obj` that rely on
  helper-injected state

## What to KEEP

- `tests/testthat/helper-testutils.R` (contains `make_grid_coord()` —
  genuinely useful for spatial edge list tests and doesn't touch
  fixtures)


# STRUCTURAL REQUIREMENTS

## Test file organization

One test file per exported function, named
`test-<function_name>.R`. Keep the existing one-file-per-function
structure; don't consolidate.

## Per-test_that block

- Every block starts fresh. Load fixtures with `data(...)` inside the
  block, or synthesise inputs inline.
- No block may reference a variable that was created in a previous
  block.
- No block may rely on side effects of a previous block (e.g.
  environment variables, tempdir state, options settings) without
  explicitly re-creating them.

## Messages and warnings

Most NICHESv2 functions emit `message()` calls during normal
operation. Tests should not suppress these by default — `R CMD check`
tolerates messages. Use `expect_message()` when a specific message is
part of the contract being tested.

When a function legitimately emits a warning as part of correct
behaviour (e.g. `Explore_Radius_Neighborhood()` when `n.sample = NULL`
on a small dataset), suppress with `suppressWarnings()` narrowly —
only around the offending call, not around the whole block.

## Parallelization tests

All tests MUST use `n.cores = 1L`. CRAN check machines have unreliable
parallel support and tests that spin up socket clusters often fail on
their infrastructure. If you need to exercise the parallel code path,
do it in a separate test file wrapped in `skip_on_cran()`.

## Tempdir discipline

Any test that writes to disk (`CellToCell_InSitu`, file-creating
examples) MUST use `file.path(tempdir(), ...)` and clean up with
`unlink(..., recursive = TRUE)` at the end of the block. Do not write
anywhere else under any circumstances.


# SUGGESTS-PACKAGE GUARDING

The package uses four Suggests packages: `sf`, `SeuratObject`,
`rhdf5`, `OmnipathR`. Tests that exercise code paths requiring these
packages MUST guard with `skip_if_not_installed()`:

```r
test_that("add_segmentation attaches sf polygons", {
  skip_if_not_installed("sf")
  # ... test body ...
})
```

Without the guard, CRAN machines that don't have the Suggest installed
will report an ERROR instead of a clean skip.


# COVERAGE FLOOR

The rewrite must preserve coverage of at least the following
behaviours (reverse-engineered from the 373 currently passing tests —
cross-check against the actual files during the rewrite):

## Per exported function
- Happy-path call with minimal valid inputs produces the correct
  return type
- Required arguments missing -> clear error
- Column-schema invariants on returned `data.table`s hold (required
  columns present, keyed as documented)
- Attributes set on the return object match the @return section of the
  function's roxygen

## NICHESObject-specific
- 14 slots populated exactly as specified in the design doc
- `attr(obj, "mode")` and `attr(obj, "sample.col")` behaviour
- S3 methods (`print`, `summary`, `dim`, `length`, `[`) return the
  types documented in their @return sections
- Subsetting preserves the both-endpoints-in-i invariant for edges

## Sample-aware path
- `sample.col` correctly splits and merges
- Cross-sample barcode collision -> hard error
- Mode mismatch during merge -> hard error


# DELIVERY PROTOCOL

- One test file per delivery. Mode 2 (full file replacement).
- After each file lands, confirm `devtools::test(filter = "...")` and
  `devtools::check()` both pass cleanly.
- Final delivery: a single commit that deletes
  `tests/testthat/helper-testdata.R`, `inst/testdata/`, and any
  helper infrastructure no longer referenced.


# OUT OF SCOPE

- Adding new test behaviour not previously covered (e.g. fuzzing,
  property-based testing, performance benchmarks) — if desired, those
  are Pass 4 work.
- Touching function source code in `R/`. If a test genuinely cannot
  pass because the function has a bug, that's a finding to escalate,
  not a license to modify source.
- `DESCRIPTION`, `NEWS.md`, vignette, `tools::resaveRdaFiles()` — all
  deferred to Pass 3.
