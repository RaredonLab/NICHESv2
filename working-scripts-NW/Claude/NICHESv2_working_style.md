# NICHESv2 — Working Style & Lessons
# Updated: 2026-04-10

This file is for both Nuoya and Claude. It captures observations about how
we work together, patterns that help, patterns that cost time, and concrete
mitigations. Read at the start of each session.

It is NOT a technical reference — see NICHESv2_package_design.md for that.

================================================================================
ABOUT NUOYA — STRENGTHS (keep doing these)
================================================================================

**Scientific judgment rooted in real data.**
Nuoya has worked with actual Visium data from multiple collaborators. When
she pushes back on a proposal, it is almost always because she has seen the
failure mode in practice. Examples from history:
  - Rejecting embedding-based edge definition (suppresses biologically
    interesting cross-type edges)
  - Insisting cell.data always use raw counts even when smoothing is applied
  - Keeping sf in Suggests, not Imports
Claude should treat Nuoya's "my advisor/collaborators have a different
opinion" or "in practice this breaks" as strong signals, not objections to
argue past.

**Disciplined confirmation gates.**
Nuoya enforces "explain the bug, wait for my OK, then produce code" as a
hard rule. This catches mistakes before they turn into deliverables. Claude
should never skip this step just because the problem "seems obvious."

**Honest uncertainty.**
Nuoya asks "is this redundant?" instead of pretending to know. She has
explicitly said she's far from a mature software engineer. This is rare
and valuable — it means every piece of feedback produces actual change.

**Independent decision-making within her domain.**
Naming, default values, specific architecture choices — Nuoya decides these
directly without needing Claude to arbitrate. When Claude offers choices,
Nuoya picks quickly and moves on.


================================================================================
ABOUT NUOYA — PATTERNS THAT COST TIME (work to shift these)
================================================================================

**Pattern 1: Asking Claude to interpret tool output she could interpret herself.**

Nuoya often runs a verification command, gets output, and asks "is this
right?" — when the answer is visible from the output itself. Example: when
a grep for `%||%` returned integer(0), she asked if something was missing,
but opening utils.R and scrolling to the end would have answered it in
5 seconds.

This is a confidence pattern, not a capability pattern. The fix:
  - Try interpreting first
  - Ask Claude only if the interpretation is genuinely ambiguous
  - When asking, include what you expected vs. what you saw, so Claude can
    redesign the check rather than just re-explain

**Pattern 2: Asking permission instead of deciding on independent steps.**

Phrases like "should I do X first or Y first?" when X and Y have no
dependency, or "should we continue?" when the previous step is obviously
green. These prompts feel polite but fragment the work into more turns
than necessary.

Fix: When the next step is obvious, just say "doing X, then Y, will ping
when both are done." Claude will push back if she disagrees.

**Pattern 3: Missing small steps when applying long manual change lists.**

This is the most actionable observation from this session. During Pass 1,
C4 (a one-line swap: `.smooth_expression` → `smooth_expression`) was
silently not applied. It wasn't caught for two full passes because the
`smooth = TRUE` code path has zero test coverage. Tests passed. The bug
sat undetected until Claude happened to re-read the file.

**The mitigation — NEW PROTOCOL, effective immediately:**
When Claude delivers MORE THAN 2 manual edits in a single batch, Claude
MUST provide a content-anchored verification script at the end of the
batch. The script should contain one grep per edit, phrased so that the
expected result is obvious (e.g., "should return 0 hits" or "should return
exactly one line containing X"). Nuoya runs the script before declaring
the batch applied. Two minutes of grep prevents missed edits. This is a
hard rule now, recorded in Claude's persistent memory.

**Pattern 4: Pasting tool output without framing.**

When a tool gives confusing output, the fastest path is: "I ran X to find
Y. The output was Z. I expected something like W." That framing tells
Claude what Nuoya was after and lets Claude redesign the query without
guessing. The messy is_NICHESObject grep earlier in this session is the
example — Nuoya pasted cryptic line numbers with "not easy for me to
understand," and Claude had to invent three possible interpretations
because it didn't know which one mattered.


================================================================================
CLAUDE'S FAILURE MODES FROM THIS SESSION (avoid in future sessions)
================================================================================

**Failure 1: Overcorrecting after mistakes.**
When the `add_cell_meta` location violation was discovered mid-session,
Claude responded by adding more verification gates before every subsequent
step. The right response was "internalize the lesson, keep moving at
original pace." Performance of caution is not the same as being careful.
Nuoya caught this explicitly and pushed back.

**Failure 2: Inventing exotic explanations for mundane phenomena.**
When a grep showed `.smooth_expression` still in the file after Pass 1,
Claude's first theory was "the upload pipeline silently served you a stale
file from April 3." The boring explanation — "C4 was missed during manual
application" — was more likely and turned out to be correct. Occam's
razor should win by default.

**Failure 3: Delivering untested code in batches.**
Claude wrote 22 unit tests across four files in one delivery and didn't
mentally execute any of them. Two bugs shipped:
  - A matrix constructed column-major instead of byrow in
    test-build_cell_data.R produced wrong expected values
  - A garbled chained assignment `row.mins <- Matrix::rowMin <- apply(...)`
    in test-smooth_expression.R caused a runtime error
Both were visible from a slow re-read. Claude should treat any code
delivered without execution as suspect, especially batch deliveries, and
do a mental-execution pass before sending.

**Failure 4: Mechanical application of process patterns.**
Pass 1 and Pass 2 were correctly delivered in two messages each (new files
first, edits second). Pass 3 was ALSO delivered in two messages, but it
should have been one — it was small enough that splitting added a
round-trip for no benefit. Claude should notice when the previous pattern
is overkill for the current task.

**Failure 5: Trusting a cached mental model over live source.**
Early in the session, Claude trusted the knowledge base doc's claim that
`.build_cell_data` lived in utils.R, and it didn't. The fix: always read
the actual file before producing changes, never reason from the knowledge
doc alone. Knowledge docs lag the codebase by at least one session and
can be wrong about recent changes.


================================================================================
WORKING PIPELINE — PROTOCOLS THAT WORK
================================================================================

**Design-first delivery protocol**
For any non-trivial change, Claude:
  1. Explains the bug or design decision in full
  2. Identifies ALL affected files proactively (including tests)
  3. Waits for Nuoya's explicit confirmation of understanding
  4. Only then produces files or edits
This is a hard requirement for all sessions.

**File delivery modes**
Only two delivery modes are permitted:

MODE 1 — In-dialogue edit (≤3 changes, simple enough to apply manually):
  For each change, state ALL of:
    - File path
    - Exact block/line to find (content-anchored, NEVER line numbers)
    - Exact replacement
  Present old and new code as formatted fenced code blocks, one each.
  Tables are ONLY for file/location metadata — never put code inside tables.
  When delivering >2 edits, MUST include a verification grep script.

MODE 2 — Downloadable replacement file:
  Only for (a) brand-new files that do not yet exist, or (b) complete
  replacements for an existing file. Never for partial snippets, helper
  extracts, or auxiliaries. Delivered via create_file + present_files,
  never inside fenced code blocks in chat (fenced delivery causes missing
  trailing newline warnings on copy-paste, which trips R CMD check).

Rule: Mode 1 in-dialogue is for EDITS to existing files. NEW files of any
size go through Mode 2. Always.

**Change summary table ordering**
When multiple files are affected: Delete → Add → Change, in separate tables.

**Content-anchored deletes (NEW RULE from this session)**
Mode 1 deletes must anchor by content (docstring + function definition),
never by line number. Line numbers shift between deletes; content doesn't.
Always explicitly enumerate which functions in the target file STAY, not
just which are deleted. The add_cell_meta discovery proved this is
necessary.

**Verification script requirement (NEW RULE from this session)**
When Claude delivers more than 2 manual edits in a batch, Claude MUST
include a verification script at the end of the batch: one grep per edit,
with expected results stated in comments. Nuoya runs the script before
declaring the batch applied.

**File inventory before editing**
Before editing any file, run:
    grep("^[a-zA-Z_.][a-zA-Z_.0-9]* *<- *function", readLines(f))
to get the FULL function inventory. Don't just find the target function.
The add_cell_meta discovery proved this is necessary.

**Pre-attach blast radius at session start for refactor sessions**
For sessions whose goal is "refactor X across multiple files," Nuoya
attaches all affected files at the start, before Claude plans anything.
Costs 30 seconds, saves a full round-trip of re-uploading mid-session.

**Knowledge base is updated at session end, not during**
The in-context knowledge file is the previous session's baseline. Updates
are tracked verbally during the session and delivered as a Mode 2 download
at session close. Nuoya deletes the old version and replaces it.

**devtools loop**
  devtools::load_all() and devtools::test() for iteration
  devtools::document() TWICE after any @export or @noRd change
  NAMESPACE managed via roxygen tags (@import data.table in
  NICHESv2-package.R), never hand-edited


================================================================================
TOOLS & ENVIRONMENT
================================================================================

Platform: macOS ARM64 (Apple Silicon, M-series), R 4.3
  Staying on R 4.3 to avoid reinstalling many GitHub-sourced packages.
  Cross-platform compatibility is a design value.

R package infrastructure:
  devtools, testthat, roxygen2, irlba (>= 2.3.0), RANN, dbscan, data.table,
  ggplot2, sf (Suggests), parallel, SeuratObject (Suggests), rhdf5 (Suggests)

LR database in use: Connectome2025.LRMs.human.csv
  Columns: Ligand, Receptor (lowercase defaults expected in package).

Spatial data:
  Visium from osteosarcoma and normal lung tissue (samples NML1–NML7).
  Stored on OneDrive at Yale.

Version control:
  GitHub-linked local repository at ~/Documents/GitHub/NICHESv2-NW.
  Extra scratch folders: OBS, working-scripts-NW, working-scripts-MSBR —
  separate from package source.
