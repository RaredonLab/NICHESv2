# ==============================================================================
# FILE:    tests/testthat/test-download_NICHESv2_data.R
# PURPOSE: Unit tests for download_NICHESv2_data().
#          No internet needed: a small zip built from NICHESv2_inputs is served
#          from a local folder through a file:// URL.
# ==============================================================================

# Builds <src.dir>/tiny.zip holding tiny/inputs.rds and returns the dataset
# entry and the file:// URL of src.dir.
.make_test_release <- function(src.dir) {
  data("NICHESv2_inputs", package = "NICHESv2", envir = environment())
  inputs <- NICHESv2_inputs[c("count.mtx", "meta.data")]
  dir.create(file.path(src.dir, "tiny"), recursive = TRUE)
  saveRDS(inputs, file.path(src.dir, "tiny", "inputs.rds"))
  writeLines("# tiny", file.path(src.dir, "tiny", "PROVENANCE.md"))
  old.wd <- setwd(src.dir)
  on.exit(setwd(old.wd))
  utils::zip("tiny.zip", files = "tiny", flags = "-r9Xq")
  unlink(file.path(src.dir, "tiny"), recursive = TRUE)
  path <- normalizePath(src.dir, winslash = "/")
  list(
    inputs = inputs,
    entry  = list(file = "tiny.zip",
                  md5  = unname(tools::md5sum(file.path(src.dir, "tiny.zip"))),
                  size = "1 KB"),
    url    = paste0("file://", if (startsWith(path, "/")) "" else "/", path)
  )
}

# ── Argument checks ────────────────────────────────────────────────────────────

test_that("download_NICHESv2_data: unknown or missing dataset stops", {
  expect_error(download_NICHESv2_data("pbmc3k", cache.dir = tempfile()),
               "must be one of")
  expect_error(download_NICHESv2_data(cache.dir = tempfile()), "must be one of")
  expect_error(download_NICHESv2_data(c("stxBrain", "ifnb"),
                                      cache.dir = tempfile()),
               "must be one of")
})

test_that("download_NICHESv2_data: both datasets are registered", {
  datasets <- NICHESv2:::.niches_datasets()
  expect_setequal(names(datasets), c("stxBrain", "ifnb"))
  for (d in datasets) expect_match(d$md5, "^[0-9a-f]{32}$")
})

# ── Download, cache and failure paths ──────────────────────────────────────────

test_that("download_NICHESv2_data: downloads, unpacks and then reads the cache", {
  skip_if(Sys.which("zip") == "", "zip program not available")
  src.dir   <- withr::local_tempdir()
  cache.dir <- withr::local_tempdir()
  release   <- .make_test_release(src.dir)
  withr::local_options(NICHESv2.data.url = release$url)
  local_mocked_bindings(.niches_datasets = function() list(tiny = release$entry))

  result <- download_NICHESv2_data("tiny", cache.dir = cache.dir, verbose = FALSE)
  expect_equal(result, release$inputs)
  expect_true(file.exists(file.path(cache.dir, "tiny", "PROVENANCE.md")))
  expect_false(file.exists(file.path(cache.dir, "tiny.zip")))

  # The second call reads the cache, even with the source gone.
  unlink(file.path(src.dir, "tiny.zip"))
  expect_message(
    cached <- download_NICHESv2_data("tiny", cache.dir = cache.dir),
    "Reading cached"
  )
  expect_equal(cached, release$inputs)
})

test_that("download_NICHESv2_data: a wrong MD5 sum returns NULL with a message", {
  skip_if(Sys.which("zip") == "", "zip program not available")
  src.dir   <- withr::local_tempdir()
  cache.dir <- withr::local_tempdir()
  release   <- .make_test_release(src.dir)
  release$entry$md5 <- strrep("0", 32L)
  withr::local_options(NICHESv2.data.url = release$url)
  local_mocked_bindings(.niches_datasets = function() list(tiny = release$entry))

  expect_message(
    result <- download_NICHESv2_data("tiny", cache.dir = cache.dir, verbose = FALSE),
    "doesn't match its expected MD5 sum"
  )
  expect_null(result)
  expect_length(list.files(cache.dir), 0L)
})

test_that("download_NICHESv2_data: a failed download returns NULL with a message", {
  cache.dir <- withr::local_tempdir()
  missing.dir <- normalizePath(withr::local_tempdir(), winslash = "/")
  withr::local_options(NICHESv2.data.url = paste0(
    "file://", if (startsWith(missing.dir, "/")) "" else "/", missing.dir, "/none"
  ))

  expect_message(
    result <- download_NICHESv2_data("ifnb", cache.dir = cache.dir, verbose = FALSE),
    "Could not download 'ifnb'"
  )
  expect_null(result)
  expect_length(list.files(cache.dir), 0L)
})
