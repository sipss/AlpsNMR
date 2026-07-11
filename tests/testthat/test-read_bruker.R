test_that("nmr_read_samples_dir works", {
    dir_to_demo_dataset <- system.file("dataset-demo", package = "AlpsNMR")
    dataset <- nmr_read_samples_dir(dir_to_demo_dataset)
    expect_equal(dataset$num_samples, 3)
})

test_that("nmr_read_samples returns unique NMR experiments", {
    dir_to_demo_dataset <- system.file("dataset-demo", package = "AlpsNMR")
    dataset <- nmr_read_samples(
        c(
            file.path(dir_to_demo_dataset, "10.zip"),
            file.path(dir_to_demo_dataset, "10.zip")
        )
    )
    expect_equal(dataset$num_samples, 2)
    expect_false(any(duplicated(names(dataset))))
})

test_that("create_sample_names returns good unique guesses", {
    sample_names <- c("a", "b")
    expect_equal(create_sample_names(sample_names), sample_names)
    sample_names <- c("a.zip", "b.zip")
    expect_equal(create_sample_names(sample_names), c("a", "b"))
    sample_names <- c("bar/a.zip", "bar/b.zip")
    expect_equal(create_sample_names(sample_names), c("a", "b"))
    sample_names <- c("bar/a.zip", "foo/b.zip")
    expect_equal(create_sample_names(sample_names), c("a", "b"))
    sample_names <- c("bar/a.zip", "foo/a.zip")
    expect_equal(create_sample_names(sample_names), c("bar/a", "foo/a"))
})

test_that("read_orig_file() leaves an empty value for a key with no value", {
    dir <- withr::local_tempdir()
    writeLines(
        c("KEY1 value1", "SINGLETOKEN", "KEY2 value2 extra"),
        file.path(dir, "orig")
    )
    result <- read_orig_file(dir)

    # Regular key-value pairs are parsed normally
    expect_equal(result$KEY1, "value1")
    expect_equal(result$KEY2, "value2 extra")

    # A line with a single token (no value) parses to an empty string
    expect_equal(result$SINGLETOKEN, "")
    expect_false(grepl("NA", result$SINGLETOKEN, fixed = TRUE))
})

test_that("read_bin_data() reports the file-open error when the file does not exist", {
    nonexistent_file <- tempfile(pattern = "does-not-exist-")
    expect_false(file.exists(nonexistent_file))

    err <- tryCatch(
        read_bin_data(nonexistent_file, endian = "little"),
        error = function(e) e
    )
    expect_s3_class(err, "error")
    expect_false(grepl("object 'con' not found", conditionMessage(err), fixed = TRUE))
    expect_true(grepl("cannot open", conditionMessage(err), fixed = TRUE))
})

test_that("nmr_read_bruker_fid() reads acqus-driven byte order/dtype/TD/SW_h and returns a data frame", {
    dir <- withr::local_tempdir()
    writeLines(
        c(
            "##$BYTORDA= 0",
            "##$DTYPA= 0",
            "##$TD= 8",
            "##$SW_h= 4000"
        ),
        file.path(dir, "acqus")
    )
    writeBin(
        as.integer(c(10, 1, 20, 2, 30, 3, 40, 4)),
        file.path(dir, "fid"),
        size = 4,
        endian = "little"
    )

    result <- nmr_read_bruker_fid(dir)

    expect_s3_class(result, "data.frame")
    expect_equal(nrow(result), 4)
    expect_equal(Re(result$fid_complex), c(10, 20, 30, 40))
    expect_equal(Im(result$fid_complex), c(1, 2, 3, 4))
    expect_equal(result$time_s[1], 0)
    expect_equal(diff(result$time_s), rep(1 / 4000, 3))
})

test_that("nmr_read_bruker_fid() returns NULL when the sample has no fid file", {
    dir <- withr::local_tempdir()
    # No "fid" file created in this directory.
    result <- nmr_read_bruker_fid(dir)
    expect_null(result)
})
