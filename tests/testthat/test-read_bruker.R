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

test_that("nmr_read_samples_dir names samples by parent directory when the same EXPNO repeats across unrelated study folders", {
    # Regression test for https://github.com/sipss/AlpsNMR/issues/62: several
    # unrelated top-level directories (not sharing any common naming scheme)
    # each containing a single Bruker experiment folder named identically
    # (e.g. the default EXPNO "10"). The leaf name alone ("10") is
    # duplicated across all of them, so create_sample_names() must fall
    # back to prepending the (differing) parent directory name.
    dir_to_demo_dataset <- system.file("dataset-demo", package = "AlpsNMR")
    zip_files <- fs::dir_ls(dir_to_demo_dataset, glob = "*.zip")

    study_root <- withr::local_tempdir()
    study_dirs <- character(length(zip_files))
    for (i in seq_along(zip_files)) {
        study_dir <- file.path(study_root, paste0("Study", i))
        dir.create(study_dir)
        utils::unzip(zip_files[i], exdir = study_dir)
        expno_dir <- list.dirs(study_dir, recursive = FALSE)
        file.rename(expno_dir, file.path(study_dir, "10")) # force the same leaf name "10"
        study_dirs[i] <- study_dir
    }

    dataset <- nmr_read_samples_dir(study_dirs)

    expect_equal(dataset$num_samples, length(zip_files))
    expect_equal(sort(names(dataset)), sort(paste0("Study", seq_along(zip_files), "/10")))
    # None of the vctrs::vec_as_names "...N" disambiguation suffixes should
    # have been needed here, since the parent directory names already make
    # every sample unique:
    expect_false(any(grepl("\\.\\.\\.", names(dataset))))
})

test_that("nmr_read_samples_dir disambiguates samples whose collision goes deeper than one directory level", {
    # Regression test for https://github.com/sipss/AlpsNMR/issues/62: two
    # samples can share not just the same leaf EXPNO name ("10") but also
    # the same *immediate parent* name ("subject1"), while still being
    # distinguishable a level further up ("groupA" vs "groupB"). A naming
    # scheme that only prepends one parent directory level would still
    # collide on "subject1/10" for both and fall back to the unreadable
    # vctrs::vec_as_names() suffixes; create_sample_names() must walk up as
    # many levels as needed instead.
    dir_to_demo_dataset <- system.file("dataset-demo", package = "AlpsNMR")
    zip_files <- fs::dir_ls(dir_to_demo_dataset, glob = "*.zip")[1:2]

    study_root <- withr::local_tempdir()
    group_names <- c("groupA", "groupB")
    sample_dirs <- character(length(zip_files))
    for (i in seq_along(zip_files)) {
        subject_dir <- file.path(study_root, "study1", group_names[i], "subject1")
        dir.create(subject_dir, recursive = TRUE)
        utils::unzip(zip_files[i], exdir = subject_dir)
        expno_dir <- list.dirs(subject_dir, recursive = FALSE)
        file.rename(expno_dir, file.path(subject_dir, "10")) # force the same leaf name "10"
        sample_dirs[i] <- subject_dir
    }

    dataset <- nmr_read_samples_dir(sample_dirs)

    expect_equal(dataset$num_samples, length(zip_files))
    expect_equal(sort(names(dataset)), sort(paste0(group_names, "/subject1/10")))
    expect_false(any(grepl("\\.\\.\\.", names(dataset))))
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
    # A collision that only resolves two directory levels up: both the leaf
    # ("10") and the immediate parent ("subject1") repeat, but a
    # grandparent ("groupA"/"groupB") differs.
    sample_names <- c(
        "root/study1/groupA/subject1/10.zip",
        "root/study1/groupB/subject1/10.zip"
    )
    expect_equal(
        create_sample_names(sample_names),
        c("groupA/subject1/10", "groupB/subject1/10")
    )
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

    # file(file_name, "rb") lazily creates the connection object without
    # erroring; the actual open() attempt (and its "cannot open file"
    # warning) only happens on the first read below, and that warning isn't
    # caught by read_bin_data()'s own tryCatch (no warning handler there) or
    # by the error-only tryCatch here, so it would otherwise leak out of the
    # test and get flagged as a failure by CI (rcmdcheck's error_on =
    # "warning").
    err <- suppressWarnings(tryCatch(
        read_bin_data(nonexistent_file, endian = "little"),
        error = function(e) e
    ))
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
