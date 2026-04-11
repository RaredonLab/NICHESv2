# NICHESv2 — Engineering Aspirations
# Updated: 2026-04-10

This file describes what a mature R package looks like AT NICHESV2'S STAGE,
where NICHESv2 currently measures up, where it falls short, and — equally
importantly — which improvements would be premature to make right now.

The gap between this file and current practice is expected and okay.
Closing the gap is incremental. This is a reference for deliberate
improvement, not a checklist of shame.

================================================================================
A HONEST FRAMING NOTE
================================================================================

There is no single "right way" a senior engineer would build this package.
The right approach depends on what NICHESv2 is FOR:

  Methods-paper package:     optimize for clarity of the algorithmic
                             narrative, reproducibility of figures, ease
                             of someone else running the pipeline on their
                             data. Performance and broad API surface are
                             secondary.

  Community-adoption tool:   optimize for stable API, comprehensive docs,
                             breadth of input format support, Bioconductor
                             compliance, vignettes.

  Production infrastructure: optimize for rigorous testing, backward
                             compatibility, CI/CD, telemetry, performance
                             benchmarks.

Most R packages in computational biology are some mix of the first two.
NICHESv2 is currently mostly the first, with aspirations toward the
second. That is appropriate for the stage of work. A senior engineer who
came in and imposed CRAN/Bioconductor practices on a package still being
algorithmically defined would actually make it worse, because premature
standardization locks in design decisions before they're ready.

So the question is not "is my design perfect," it is "is my design
appropriate for what NICHESv2 is right now." Mostly yes, with specific
things that could improve as the package matures.


================================================================================
WHERE NICHESV2 CURRENTLY DOES WELL
================================================================================

**The NICHESObject 14-slot design.**
This is the heart of the package and it is well-thought-out. The
separation of cell.list / cell.data / cell.meta (and the analogous edge
trio) reflects a real conceptual distinction: "what cells exist" vs
"what we know about each cell quantitatively" vs "what we know about
each cell descriptively." Many packages collapse these into one bag and
regret it later. Neighborhood slots being NULL until populated is also
right — lazy population avoids paying memory costs you don't need.

**Naming convention discipline.**
Underscore for functions, dot for variables, distinct separators for
edges and LRMs. The specific choices matter less than the fact that the
conventions have been held rigidly across the entire codebase.
Consistency is the win.

**Mode dispatch architecture.**
`spatial` vs `sampled` as a mode attribute, with separate constructors
sharing machinery, is the right pattern. Many packages do this with
conditionals scattered through every function and it becomes a nightmare
to maintain.

**Asymmetric return types where appropriate.**
compute_PseudoBulk_Connectivity returns a plain data.table instead of
forcing a NICHESObject. This is a subtle but excellent decision —
pseudobulk is conceptually a different kind of output that doesn't need
the full object machinery. The instinct to force "consistency" here
would have made things worse.

**Platform-universal parallelism.**
Socket clusters everywhere instead of fork-based mclapply. Costs some
verbosity but works on all platforms and is reproducible. Many R
packages punt on this and then break on Windows or macOS with ggplot2
in the workers.


================================================================================
WHAT "GOOD" LOOKS LIKE AT NICHESV2'S STAGE
================================================================================

These are concrete criteria, not platitudes.

**File organization**
  - One exported function per file, named R/<function_name>.R
  - Internal helpers either (a) in utils.R if they are small universal
    utilities, or (b) in their own R/.<helper_name>.R file if they are
    substantial branch implementations of a public function
  - Dot-prefixed, @noRd tagged for internal helpers

**Test coverage**
  - Every exported function has at least one dedicated test file
  - Every error path (every stop() call) has at least one test that
    triggers it
  - Every conditional branch in the main path is exercised by at least
    one test
  - Use covr::package_coverage() to identify gaps; target is not 100%,
    target is "know what's not covered and have a reason for each gap"

**Documentation**
  - Every exported function has roxygen with @param for every argument,
    @return describing output structure, at least one @example (or
    \dontrun{} if running requires data), and @seealso linking to
    related functions
  - A NEWS.md file documenting every API change, one entry per release
  - At least one vignette demonstrating the intended user workflow
    end-to-end on a small example dataset

**Dependencies**
  - Imports: only what is used in every code path
  - Suggests: everything used in optional code paths (smoothing,
    segmentation, specific file formats)
  - No package used only in a single function's optional branch should
    be in Imports

**Reproducibility**
  - All stochastic operations accept a seed argument with a documented
    default
  - Default seed should be a specific integer (e.g., 42L), not
    sample.int(.Machine$integer.max, 1L)

**Performance**
  - Do NOT optimize before profiling
  - Run profvis::profvis() on a representative workload before making
    any performance claims or optimizations
  - Benchmark suites (e.g., via bench::mark()) for any function that
    users will call on large data


================================================================================
SPECIFIC GAPS IN NICHESV2 RIGHT NOW (actionable)
================================================================================

Ranked by impact-to-effort ratio.

**1. Zero test coverage on smooth = TRUE in create_NICHESObject().**
The smoothing path is entirely untested. A bug in C4 of the Pass 1
refactor (line 208, `.smooth_expression` → `smooth_expression`) sat
undetected because no test exercises it. Fix: write one minimal test
that calls create_NICHESObject(mode="sampled", smooth=TRUE) on a small
fixture and verifies the output NICHESObject is well-formed.

**2. No test coverage measurement.**
Set up covr::package_coverage() as a one-line check. Run once a month.
Current coverage is probably in the 60-80% range but nobody knows
exactly. Install: install.packages("covr"). Run: covr::package_coverage().

**3. No CI.**
Set up a minimal GitHub Actions workflow running R CMD check --as-cran
on Ubuntu R-release only. Not multiple platforms yet — start minimal.
This catches breakage that isn't caught by local devtools::test() (e.g.,
NAMESPACE drift, missing @importFrom, documentation mismatches).

**4. add_cell_meta lives in the wrong file.**
Currently in R/create_NICHESObject_single.R. Should be in
R/add_cell_meta.R. (NOTE: this was fixed in Pass 3 of the 2026-04-10
session — verify it has landed before marking this item done.)

**5. is_NICHESObject vs .is_NICHESObject inconsistency.**
Both names exist. The real definition is in R/NICHESObject_class.R as
is_NICHESObject (no dot). The dot-prefix variant is an alias or
duplicate. Call sites inconsistently use one or the other. Pick one
convention, apply everywhere, delete the other. Cosmetic but catches
attention during code review.

**6. No NEWS.md file.**
Every API change in this session (the 4 new exported functions from
Pass 1) should be documented somewhere users can find it. Create
NEWS.md and add a "## NICHESv2 (development version)" section with the
new exports listed.

**7. No vignette.**
A vignette forces you to think about user experience in ways that
writing function code doesn't. Start with one tiny vignette:
"Getting Started with NICHESv2 — Spatial Workflow on a Toy Dataset."
Use the existing test fixture as the dataset. The act of writing the
vignette will expose every awkwardness in the current API and let you
fix them before users complain.

**8. No lint configuration.**
lintr::lint_package() would catch style inconsistencies and latent
bugs like unused imports, undefined global variables, etc. Add a
.lintr file even if you don't enforce all rules — running it once
surfaces real issues.


================================================================================
THINGS THAT WOULD BE PREMATURE RIGHT NOW
================================================================================

Equally important: not everything a "mature" package has is appropriate
for NICHESv2 today. Don't do these yet.

**CRAN submission.**
The API is still being defined. Submit to CRAN only after the methods
paper is published and you're confident the public API is stable.
CRAN submission creates strong backward-compatibility expectations
that are harmful during active algorithm development.

**Multi-platform CI.**
Ubuntu R-release is enough. Don't add Windows + macOS + R-devel +
R-oldrel until you have a reason to support them (e.g., a user filing
a Windows-specific issue).

**Optimizing performance beyond profvis-identified hotspots.**
You have notes about "performance optimization for large datasets" as
a pending item. The right move is: profvis::profvis() on a
representative workload, then optimize whatever lights up red. Do NOT
optimize from intuition. "Pre-extract gene expression vectors once"
may or may not actually help — measure before changing.

**Aggressive refactoring toward senior-engineer aesthetics.**
The helper extraction done in Pass 1 of this session was arguably
already too aggressive. Exporting build_cell_data, build_edge_meta,
smooth_expression, and compute_LRM_Score is a commitment. If nobody
outside NICHESv2 actually needs to call those directly, they should
probably have stayed internal. Future refactoring should wait for
concrete external use cases.

**Comprehensive unit tests for every function.**
Integration tests that exercise the full create_NICHESObject pipeline
give you more coverage per line of test code than isolated unit tests.
Prioritize integration over unit until coverage measurement shows real
gaps.

**Adding a user-friendly wrapper layer.**
Your current API is decomposed by data flow (compute_CellToCell,
build_edge_meta), not by user goal (analyze_spatial_communication). A
wrapper layer would make the package easier for new users, but it's
premature — wait until the underlying API is stable AND until you have
evidence that new users are struggling with the current surface.


================================================================================
IF YOU EVER GIVE THIS PROJECT TO CLAUDE CODE FOR AN "ARTIST PASS"
================================================================================

Use a scoped, specific prompt. Do NOT ask for vague "improve quality" or
"make it more elegant" — those produce bad results. Here is a prompt that
would do meaningful work without breaking anything:

    I'm the author of the NICHESv2 R package, a methods-paper-stage
    package for spatial cell-cell communication analysis. I want you to
    do a deliberate quality-improvement pass without changing the public
    API or any algorithmic behavior.

    Focus on these specific areas, in priority order:

    1. Set up covr for test coverage and produce a coverage report.
       Identify the three lowest-coverage exported functions and propose
       minimal additional tests to cover them.

    2. Run lintr::lint_package() and goodpractice::gp(). Surface the top
       10 issues by impact, ignore stylistic noise.

    3. Audit the file organization: every exported function must live in
       R/<function_name>.R. Flag any violations.

    4. Audit roxygen coverage: every exported function must have @param
       for every argument, a @return describing the output structure, at
       least one @example (or \dontrun{} if running requires data), and
       a @seealso linking to related functions. Flag gaps.

    5. Set up a minimal GitHub Actions workflow running R CMD check
       --as-cran on Ubuntu R-release. Don't try multiple platforms yet.

    Do NOT touch algorithm code, do NOT rename anything, do NOT change
    function signatures, do NOT add new exports. The goal is
    infrastructure quality, not redesign. Present findings as a
    prioritized report before making any changes, and wait for my
    approval on each item.

Notice what is NOT in that prompt: "make it more elegant," "improve the
architecture," "make it more senior-engineer-ish." Vague prompts get
vague, usually harmful, results. Specific scoped prompts get useful work.


================================================================================
CLOSING NOTE
================================================================================

You are not behind schedule. The package is in good shape for a
methods-paper-stage package. The gaps above are normal for this stage
and closing them is incremental work that happens in parallel with
continued research, not a remediation project that blocks everything
else.

The most valuable single action from this list is probably #1 (zero
test coverage on smoothing) followed by #2 (set up covr). Everything
else can wait until there's a specific trigger.
