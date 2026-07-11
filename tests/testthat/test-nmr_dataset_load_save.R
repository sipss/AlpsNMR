test_that("nmr_dataset_load() rejects an .rds file that is not a valid AlpsNMR dataset", {
    # nmr_dataset_load() reads an arbitrary .rds file from disk and hands it
    # back to the caller as a dataset object. A corrupted or otherwise
    # unrelated .rds file must not be silently treated as a valid dataset --
    # it must be rejected with a clear error instead of returning garbage
    # that downstream code would then operate on.
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
