# NICHESv2 — Session Focus: Remove EdgeLRResult

## Context
The `EdgeLRResult` S3 class in `R/compute_CellToCell.R` is vestigial design
debt predating the NICHESObject architecture. It was originally the return
type of `compute_CellToCell()` when that function was conceived as a
user-facing endpoint. Under the current architecture it has been reduced
to an internal intermediate format: `.create_NICHESObject_single()` calls
`compute_CellToCell()`, immediately unpacks the EdgeLRResult fields into
`$edge.list` / `$edge.data` / `$edge.meta` of a new NICHESObject, and
discards the wrapper. The class serves no purpose beyond adding a
construction-then-disassembly hop.

## Confirmed facts from the previous session's diagnostic greps
- `compute_CellToCell` is NOT currently exported (NAMESPACE grep returned
                                                  integer(0) for both `compute_CellToCell` and `EdgeLRResult`).
- Single internal caller: `R/create_NICHESObject_single.R` calls
`compute_CellToCell()` once and assigns to a local variable `elr`.
- Test references to EdgeLRResult exist in:
  tests/testthat/test-compute_CellToCell.R (lines 5, 12, 117)
tests/testthat/test-utilities.R (line 129)

## Scope of cleanup
(a) Change `compute_CellToCell()` to return a plain data.table with the
relevant columns (to be determined by reading the current
                  EdgeLRResult construction).
(b) Update the single caller in `create_NICHESObject_single.R` to read
the data.table directly.
(c) Update `test-compute_CellToCell.R` and `test-utilities.R` to assert
the new return shape.
(d) Delete the EdgeLRResult constructor, `print.EdgeLRResult`, `as_matrix`
generic, and `as_matrix.EdgeLRResult` method from `compute_CellToCell.R`.

## Open design question to resolve first
Should `compute_CellToCell()` remain a public exported function, or be
renamed to an internal helper?
  - Prior knowledge docs described it as public, but NAMESPACE shows it
currently isn't exported.
- If it should be public: re-export, write full roxygen, document as API.
- If it should be internal: rename to dot-prefix, tag @noRd, stop
  documenting as API.
This is a real API decision and should not be made mechanically.

## Package design reference
See NICHESv2_package_design.md (uploaded alongside this file) — noting
that its description of `compute_CellToCell` returning an EdgeLRResult
is what this cleanup is undoing.