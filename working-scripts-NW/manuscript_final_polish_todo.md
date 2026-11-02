# Manuscript final-polish list

Temporary working note. Do not edit the manuscript for these items now; apply
them together in one final polishing round, after the package is final.
Delete this file once that round is done.

Manuscript checked: `NICHESv2 paper draft Genome Biology Software 2026-10-19.docx`

## To fix

1. Fill the `[version]` placeholder in Methods, "Software implementation and
   availability".
2. Word the version statement to match the evidence. PPLR results came from an
   earlier version, so only say 1.0.0 if the score fingerprints support it.
3. Update "available on GitHub [38]" to name CRAN, the website and the Zenodo
   archive, whichever are live at that point.
4. Add the data and code availability text. The draft has none yet. Include the
   public test datasets (`stxBrain`, `ifnb`) and their sources.
5. Check how PPLR access works for reviewers against the journal's data policy.
   Public test data lets reviewers try the software, not reproduce PPLR figures.
6. Check Genome Biology's live guidelines for required Software article sections
   and any archive requirement. Confirmed so far: public repository and OSI
   license. Not confirmed: Zenodo or reviewer test-data rules.
7. Align the TissuePlex platform list. The manuscript names Xenium, MERSCOPE and
   CosMx; the TissuePlex README also lists Visium, Visium HD and seqFISH.
8. Update references 23 and 38 (titles, URLs, access dates).

## To re-verify against the final package

9. Requirements: R >= 4.1.0, R >= 4.2 for the AnnData reader, and the dependency
   list.
10. Functions, slots and defaults the manuscript cites, including 15 slots,
    27 aggregation columns, `n.sample` 2,000, seed 42 and `chunk.size` 50.
11. The Discussion sentence that the extraction functions "report the layer they
    use", so it still matches the message in both extract functions.

## After the title is final (package side)

12. Update the README citation and add `inst/CITATION`. The README currently uses
    an older title.
