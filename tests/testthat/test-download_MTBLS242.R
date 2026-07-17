test_that("download_MTBLS242() rejects a zip archive with a path-traversal entry", {
    skip_if_not_installed("zip")
    skip_if_not_installed("fs")

    # download_MTBLS242() downloads a zip archive and extracts it locally. An
    # archive entry name is not trustworthy regardless of how the archive was
    # obtained: a maliciously (or accidentally) crafted entry can encode
    # "../" segments that, once combined with the extraction root, resolve
    # outside of it ("zip-slip"). Archive contents must never be allowed to
    # write outside the intended extraction directory, so download_MTBLS242()
    # checks each entry's resolved path against the extraction root before
    # calling zip::unzip() and aborts if any entry would escape it.
    #
    # This test builds such a malicious archive and confirms the guard
    # rejects it before any extraction happens.
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
    # The metadata file and the (here: malicious) per-sample zip archive are
    # served; requests for MetaboLights' canonical SHA-256 manifest
    # (HASHES/*.json) are deliberately left unhandled, so
    # download_MTBLS242() falls back to its local-only integrity check,
    # which is what this test is about.
    mock_curl_download_retry <- function(url, destfile, ...) {
        if (grepl("s_MTBLS242\\.txt$", url, ignore.case = TRUE)) {
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

    # download_MTBLS242() (with its inline zip-slip check) must refuse to
    # extract the malicious archive.
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

test_that("download_MTBLS242() extracts a benign zip archive normally", {
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

    # As above, HASHES/*.json requests are left unhandled on purpose so this
    # test exercises the local-only fallback rather than canonical
    # verification (covered separately below).
    mock_curl_download_retry <- function(url, destfile, ...) {
        if (grepl("s_MTBLS242\\.txt$", url, ignore.case = TRUE)) {
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

test_that("download_MTBLS242() verifies fresh downloads against MetaboLights' canonical SHA-256 manifest", {
    skip_if_not_installed("zip")
    skip_if_not_installed("fs")
    skip_if_not_installed("digest")
    skip_if_not_installed("jsonlite")

    # MetaboLights publishes canonical SHA-256 checksums for MTBLS242's
    # metadata and per-sample data files at
    # <url>/HASHES/{metadata,data}_sha256.json, served over HTTPS. When that
    # manifest is reachable, download_MTBLS242() verifies every freshly
    # downloaded file against it directly -- not just against a value it
    # pinned itself on a previous run -- so it also protects the very first
    # download of a file against tampering or corruption in transit.
    work_root <- withr::local_tempdir()
    dest_dir <- file.path(work_root, "mtbls_test")
    dst_rootdir <- file.path(dest_dir, "samples")
    dir.create(dst_rootdir, recursive = TRUE)

    meta_content <- c("Sample Name\tFactor Value[time point]", "0-0001-1\tpreop")
    meta_tmp <- withr::local_tempfile()
    writeLines(meta_content, meta_tmp)
    meta_hash <- digest::digest(meta_tmp, algo = "sha256", file = TRUE)

    benign_src <- file.path(work_root, "benign_src")
    dir.create(file.path(benign_src, "Obs0_0001s", "3"), recursive = TRUE)
    writeLines("spectrum-data", file.path(benign_src, "Obs0_0001s", "3", "1r"))
    benign_zip <- file.path(work_root, "benign.zip")
    withr::with_dir(benign_src, {
        zip::zip(zipfile = benign_zip, files = file.path("Obs0_0001s", "3", "1r"))
    })
    zip_hash <- digest::digest(benign_zip, algo = "sha256", file = TRUE)

    mock_curl_download_retry <- function(url, destfile, ...) {
        if (grepl("metadata_sha256\\.json$", url)) {
            jsonlite::write_json(list("s_MTBLS242.txt" = meta_hash), destfile, auto_unbox = TRUE)
        } else if (grepl("data_sha256\\.json$", url)) {
            jsonlite::write_json(list("FILES/Obs0_0001s.zip" = zip_hash), destfile, auto_unbox = TRUE)
        } else if (grepl("s_MTBLS242\\.txt$", url, ignore.case = TRUE)) {
            writeLines(meta_content, destfile)
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

test_that("download_MTBLS242() aborts a fresh download that doesn't match MetaboLights' canonical SHA-256", {
    skip_if_not_installed("zip")
    skip_if_not_installed("fs")
    skip_if_not_installed("digest")
    skip_if_not_installed("jsonlite")

    work_root <- withr::local_tempdir()
    dest_dir <- file.path(work_root, "mtbls_test")
    dst_rootdir <- file.path(dest_dir, "samples")
    dir.create(dst_rootdir, recursive = TRUE)

    wrong_hash <- strrep("0", 64) # deliberately does not match the served content

    mock_curl_download_retry <- function(url, destfile, ...) {
        if (grepl("metadata_sha256\\.json$", url)) {
            jsonlite::write_json(list("s_MTBLS242.txt" = wrong_hash), destfile, auto_unbox = TRUE)
        } else if (grepl("data_sha256\\.json$", url)) {
            jsonlite::write_json(list(), destfile)
        } else if (grepl("s_MTBLS242\\.txt$", url, ignore.case = TRUE)) {
            writeLines(
                c("Sample Name\tFactor Value[time point]", "0-0001-1\tpreop"),
                destfile
            )
        } else {
            stop("Unexpected curl_download_retry() call in test mock: ", url)
        }
        invisible(destfile)
    }
    testthat::local_mocked_bindings(
        curl_download_retry = mock_curl_download_retry,
        .package = "AlpsNMR"
    )

    expect_error(
        download_MTBLS242(
            dest_dir = dest_dir,
            force = TRUE,
            keep_only_CPMG_1r = TRUE,
            keep_only_preop_and_3months = TRUE,
            keep_only_complete_time_points = TRUE
        ),
        regexp = "published\\s+by MetaboLights"
    )
})

test_that("download_MTBLS242() falls back to pinning local SHA-256 checksums when the canonical manifest is unavailable", {
    skip_if_not_installed("zip")
    skip_if_not_installed("fs")
    skip_if_not_installed("digest")

    # If MetaboLights' canonical SHA-256 manifest cannot be fetched (network
    # issue, simulated here by simply not mocking the HASHES/*.json
    # requests), download_MTBLS242() falls back to a local-only safety net:
    # the SHA-256 of every downloaded file is pinned to
    # `<dest_dir>/SHA256SUMS` the first time it is saved, and re-verified on
    # every later call that reuses the cached file, so local
    # corruption/tampering between calls is still detected.
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
        if (grepl("s_MTBLS242\\.txt$", url, ignore.case = TRUE)) {
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

    # First download: the canonical manifest fetch fails (not mocked) and
    # checksums get pinned locally to SHA256SUMS instead.
    download_MTBLS242(
        dest_dir = dest_dir,
        force = TRUE,
        keep_only_CPMG_1r = TRUE,
        keep_only_preop_and_3months = TRUE,
        keep_only_complete_time_points = TRUE
    )
    manifest_file <- file.path(dest_dir, "SHA256SUMS")
    expect_true(file.exists(manifest_file))
    sample_zip <- file.path(dst_rootdir, "Obs0_0001s.zip")
    expected_hash <- digest::digest(sample_zip, algo = "sha256", file = TRUE)
    manifest <- readLines(manifest_file)
    expect_true(any(grepl(paste0("^", expected_hash, "  samples/Obs0_0001s\\.zip$"), manifest)))

    # A subsequent call that reuses the cached file (force = FALSE) must
    # verify it against the pinned checksum without error or re-download.
    expect_no_error(
        download_MTBLS242(
            dest_dir = dest_dir,
            force = FALSE,
            keep_only_CPMG_1r = TRUE,
            keep_only_preop_and_3months = TRUE,
            keep_only_complete_time_points = TRUE
        )
    )

    # Tampering with (or corrupting) the cached sample zip after it was
    # pinned must be caught and refused, not silently used.
    writeBin(as.raw(c(0, 1, 2, 3)), sample_zip)
    expect_error(
        download_MTBLS242(
            dest_dir = dest_dir,
            force = FALSE,
            keep_only_CPMG_1r = TRUE,
            keep_only_preop_and_3months = TRUE,
            keep_only_complete_time_points = TRUE
        ),
        regexp = "Checksum mismatch"
    )
})
