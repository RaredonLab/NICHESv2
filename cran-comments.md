## Submission

This is the first submission of NICHESv2, version 1.0.0.

## Test environments

- Local: macOS 14.7 (aarch64), R 4.3.3
- GitHub Actions: macOS 26.6 (aarch64), Windows Server 2022 (x64) and
  Ubuntu 24.04 (x64), each with R 4.6.2
- GitHub Actions: Ubuntu 24.04, R-devel (4.7.0)
- win-builder: R-devel (result to be added)

## R CMD check results

0 errors | 0 warnings | 1 note

* checking CRAN incoming feasibility ... NOTE

  New submission

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
