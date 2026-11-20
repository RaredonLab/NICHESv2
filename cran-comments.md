## Resubmission

This is a resubmission. In response to the review:

- Software names in the Description field of DESCRIPTION are now in single
  quotes: 'NICHES', 'ConnectomeDB2025', 'FANTOM5' and 'OmniPath'. The Title
  names no software.

## Submission

This is the first submission of NICHESv2, version 1.0.0.

## Test environments

- Local: macOS 14.7 (aarch64), R 4.3.3
- GitHub Actions: macOS 26.6 (aarch64), Windows Server 2022 (x64) and
  Ubuntu 24.04 (x64), each with R 4.6.2
- GitHub Actions: Ubuntu 24.04, R-devel (4.7.0)
- win-builder: Windows Server 2022 (x64), R 4.6.2 and R-devel (2026-11-17 r90926)

## R CMD check results

0 errors | 0 warnings | 1 note

* checking CRAN incoming feasibility ... NOTE

  New submission

  Possibly misspelled words in DESCRIPTION: Raredon, et, al, Ligand, ligand,
  ligands, transcriptomics

  These are spelled correctly. Raredon et al. is the cited method paper, and
  the others are standard biology terms.

On some win-builder runs (R-devel, and one of two R-release runs of the same
file), a second NOTE appeared:

* checking top-level files ... NOTE

  Non-standard file/directory found at top level: 'omnipathr-log'

  This folder is not in the submitted package. It is created on the check
  machine by the suggested Bioconductor package OmnipathR, which writes an
  `omnipathr-log` folder into the working directory whenever it is loaded.
  NICHESv2 itself loads OmnipathR only inside `load_LRM_database("omnipath")`,
  which no example or test runs during the check.

## Notes for the reviewers

- `download_NICHESv2_data()` downloads one of two public CC BY 4.0 datasets
  (11 MB and 14 MB) from a release of the package's GitHub repository. Its
  example is in `\donttest{}` and writes to `tempdir()`. If the download fails
  or the file doesn't match its MD5 sum, it returns `NULL` with a message
  instead of an error. Its default cache folder is
  `tools::R_user_dir("NICHESv2", "cache")`. Its tests use a local `file://`
  URL and need no internet.
- `load_LRM_database("omnipath")` needs internet access. Its example runs only
  in interactive sessions, and no test contacts the service.
- Apart from the download example, no example is wrapped in `\donttest{}`;
  each runs in under 1.5 seconds. All examples use one core. The two tests that
  start two worker processes are skipped on CRAN.
- The suggested packages rhdf5 and OmnipathR are on Bioconductor. Code that
  uses any suggested package checks for it first, and the tests skip when it
  is missing.
