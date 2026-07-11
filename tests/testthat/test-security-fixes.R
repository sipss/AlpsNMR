# Regression tests for security fixes identified in a prior critical code
# review:
#
#  1. Zip-slip guard added inside download_MTBLS242() (R/download_MTBLS242.R)
#     before zip::unzip() is called on the (unauthenticated, plain FTP)
#     downloaded archive.
#  2. Class validation added to nmr_dataset_load() (R/nmr_dataset_load_save.R)
#     so that loading an .rds file that is not an AlpsNMR dataset object now
#     errors clearly instead of silently returning garbage.
#
# None of these tests perform any network access: download_MTBLS242()'s only
# network-touching call, the internal curl_download_retry(), is mocked out
# via testthat::local_mocked_bindings() so that the *real* zip-slip check
# inside download_MTBLS242() still runs against a locally-crafted zip file.

test_that("download_MTBLS242() zip-slip guard rejects a path-traversal entry", {
  skip_if_not_installed("zip")
  skip_if_not_installed("fs")

  # --- Build a malicious zip archive containing an entry that, once
  # combined with the extraction root, resolves outside of it. ---
  #
  # The relevant sample-selection logic in download_MTBLS242() only keeps
  # zip entries that start with "<filename_base>/3/" (the CPMG subfolder),
  # via plain startsWith() on the string -- so an entry name such as
  # "Obs0_0001s/3/../../../evil_target.txt" passes that filter while still
  # escaping the destination root on extraction. That's exactly the
  # scenario the zip-slip guard (fs::path_abs()/fs::path_norm() +
  # startsWith() check, aborting via cli::cli_abort()) is meant to catch.
  work_root <- withr::local_tempdir()
  dest_dir <- file.path(work_root, "mtbls_test")
  dst_rootdir <- file.path(dest_dir, "samples")
  dir.create(dst_rootdir, recursive = TRUE)

  # Create the (fake) legitimate subfolder structure so that the relative
  # path below is resolvable on disk when zip::zip() reads the file it
  # points to.
  withr::with_dir(dst_rootdir, {
    dir.create(file.path("Obs0_0001s", "3"), recursive = TRUE)
  })
  # The "escape hatch" file lives one directory above dst_rootdir, i.e.
  # exactly where 3 levels of ".." from "<dst_rootdir>/Obs0_0001s/3" lands.
  writeLines("evil-content", file.path(dest_dir, "evil_target.txt"))

  malicious_zip <- file.path(work_root, "malicious.zip")
  withr::with_dir(dst_rootdir, {
    entry <- file.path("Obs0_0001s", "3", "..", "..", "..", "evil_target.txt")
    stopifnot(file.exists(entry)) # sanity check on the traversal path itself
    zip::zip(zipfile = malicious_zip, files = entry)
  })
  # Confirm the archive really stores the literal ".."-containing entry
  # name (zip::zip() does not sanitize it away).
  entry_names <- zip::zip_list(malicious_zip)[["filename"]]
  expect_true(any(grepl("\\.\\.", entry_names)))

  # --- Mock the only network-touching call so no real download happens ---
  # The first call fetches the annotations metadata file, the second the
  # (here: malicious) per-sample zip archive.
  mock_curl_download_retry <- function(url, destfile, ...) {
    if (grepl("s_mtbls242\\.txt$", url)) {
      writeLines(
        c("Sample Name\tFactor Value[time point]", "0-0001-1\tpreop"),
        destfile
      )
    } else if (grepl("\\.zip$", url)) {
      file.copy(malicious_zip, destfile, overwrite = TRUE)
    } else {
      stop("Unexpected curl_download_retry() call in test mock: ", url)
    }
    invisible(destfile)
  }
  testthat::local_mocked_bindings(
    curl_download_retry = mock_curl_download_retry,
    .package = "AlpsNMR"
  )

  # The real download_MTBLS242() (with its real, inline zip-slip check)
  # must refuse to extract the malicious archive.
  expect_error(
    download_MTBLS242(
      dest_dir = dest_dir,
      force = TRUE,
      keep_only_CPMG_1r = TRUE,
      keep_only_preop_and_3months = TRUE,
      keep_only_complete_time_points = TRUE
    ),
    regexp = "zip-slip"
  )

  # Extraction must have been aborted before zip::unzip() ever ran. The
  # intermediate zip is legitimately copied into dst_rootdir before the
  # zip-slip check runs, but nothing must have been *extracted* from it:
  # the "Obs0_0001s/3" directory (created empty above) must still be empty,
  # and the pre-existing sentinel file at the escape target must be
  # untouched.
  expect_equal(list.files(file.path(dst_rootdir, "Obs0_0001s", "3")), character(0))
  expect_equal(readLines(file.path(dest_dir, "evil_target.txt")), "evil-content")
})

test_that("download_MTBLS242() still extracts a benign zip archive normally", {
  skip_if_not_installed("zip")
  skip_if_not_installed("fs")

  # Companion "no false positives" test: a well-formed archive (entries
  # confined to "<filename_base>/3/...") must still extract successfully,
  # i.e. the zip-slip guard does not break the legitimate code path.
  work_root <- withr::local_tempdir()
  dest_dir <- file.path(work_root, "mtbls_test")
  dst_rootdir <- file.path(dest_dir, "samples")
  dir.create(dst_rootdir, recursive = TRUE)

  benign_src <- file.path(work_root, "benign_src")
  dir.create(file.path(benign_src, "Obs0_0001s", "3"), recursive = TRUE)
  writeLines("spectrum-data", file.path(benign_src, "Obs0_0001s", "3", "1r"))
  benign_zip <- file.path(work_root, "benign.zip")
  withr::with_dir(benign_src, {
    zip::zip(zipfile = benign_zip, files = file.path("Obs0_0001s", "3", "1r"))
  })

  mock_curl_download_retry <- function(url, destfile, ...) {
    if (grepl("s_mtbls242\\.txt$", url)) {
      writeLines(
        c("Sample Name\tFactor Value[time point]", "0-0001-1\tpreop"),
        destfile
      )
    } else if (grepl("\\.zip$", url)) {
      file.copy(benign_zip, destfile, overwrite = TRUE)
    } else {
      stop("Unexpected curl_download_retry() call in test mock: ", url)
    }
    invisible(destfile)
  }
  testthat::local_mocked_bindings(
    curl_download_retry = mock_curl_download_retry,
    .package = "AlpsNMR"
  )

  expect_no_error(
    download_MTBLS242(
      dest_dir = dest_dir,
      force = TRUE,
      keep_only_CPMG_1r = TRUE,
      keep_only_preop_and_3months = TRUE,
      keep_only_complete_time_points = TRUE
    )
  )
  expect_true(file.exists(file.path(dst_rootdir, "Obs0_0001s.zip")))
})

test_that("nmr_dataset_load() rejects an .rds file that is not a valid AlpsNMR dataset", {
  bad_file_list <- withr::local_tempfile(fileext = ".rds")
  saveRDS(list(foo = "bar"), bad_file_list)
  expect_error(
    nmr_dataset_load(bad_file_list),
    regexp = "does not contain a valid AlpsNMR dataset"
  )

  bad_file_scalar <- withr::local_tempfile(fileext = ".rds")
  saveRDS(42, bad_file_scalar)
  expect_error(
    nmr_dataset_load(bad_file_scalar),
    regexp = "does not contain a valid AlpsNMR dataset"
  )
})

test_that("nmr_dataset_load() accepts a valid AlpsNMR dataset object", {
  fixture <- system.file("extdata", "nmr_dataset.rds", package = "AlpsNMR")
  skip_if(identical(fixture, ""), "nmr_dataset.rds fixture not found")

  valid_dataset <- readRDS(fixture)
  good_file <- withr::local_tempfile(fileext = ".rds")
  saveRDS(valid_dataset, good_file)

  loaded <- NULL
  expect_no_error(loaded <- nmr_dataset_load(good_file))
  expect_s3_class(loaded, "nmr_dataset_family")

  # Also exercise loading the original fixture path directly.
  loaded_fixture <- NULL
  expect_no_error(loaded_fixture <- nmr_dataset_load(fixture))
  expect_s3_class(loaded_fixture, "nmr_dataset_family")
})
