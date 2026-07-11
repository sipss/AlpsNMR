# Regression tests for silent-data-corruption bugs found and fixed in a prior
# critical code review of the Bruker (R/bruker.R) and JDX (R/import_jdx.R)
# import code. Each test below targets one specific fix.

test_that("read_orig_file() does not inject a spurious NA for single-token lines", {
    dir <- withr::local_tempdir()
    writeLines(
        c("KEY1 value1", "SINGLETOKEN", "KEY2 value2 extra"),
        file.path(dir, "orig")
    )
    result <- read_orig_file(dir)

    # Regular key-value pairs are parsed normally
    expect_equal(result$KEY1, "value1")
    expect_equal(result$KEY2, "value2 extra")

    # The single-token line ("SINGLETOKEN", no value) used to produce
    # line[2:1] == c(NA, "SINGLETOKEN"), i.e. a value of "NA SINGLETOKEN".
    # It must now produce an empty string instead.
    expect_equal(result$SINGLETOKEN, "")
    expect_false(grepl("NA", result$SINGLETOKEN, fixed = TRUE))
})

test_that("read_bin_data() reports the real file-open failure, not a con-not-found error", {
    nonexistent_file <- tempfile(pattern = "does-not-exist-")
    expect_false(file.exists(nonexistent_file))

    # Previously, when file() failed to open the connection, `con` was never
    # assigned, and the `finally` block's `close(con)` raised
    # "object 'con' not found", masking the real error. Now `con` is
    # initialized to NULL and guarded, so the real file-open error surfaces.
    err <- tryCatch(
        read_bin_data(nonexistent_file, endian = "little"),
        error = function(e) e
    )
    expect_s3_class(err, "error")
    expect_false(grepl("object 'con' not found", conditionMessage(err), fixed = TRUE))
    expect_true(grepl("cannot open", conditionMessage(err), fixed = TRUE))
})

test_that("process_block() (X++(Y..Y)) parser does not inject NA for single-token data lines", {
    # A data block where one interior line has only one token (e.g. a stray
    # value with no accompanying Y value). Previously x[2:length(x)] with
    # length(x) == 1 evaluated to x[2:1] == c(NA, x[1]), injecting a spurious
    # NA (and a bogus extra value) into the Y data. It should now contribute
    # zero values instead.
    lines <- c(
        "##TITLE=Test3",
        "##FIRSTX=0",
        "##LASTX=10",
        "##NPOINTS=2",
        "##DATA TABLE= (X++(Y..Y))",
        "0 5",
        "999",
        "10 7",
        "##END="
    )
    result <- process_block(lines)
    data <- result$block[["DATA TABLE"]]

    expect_equal(nrow(data), 2)
    expect_false(anyNA(data$y))
    expect_equal(data$y, c(5, 7))
    expect_equal(data$x, c(0, 10))
})

test_that("process_block() raises a clear error when XUNITS is 'Hz' but .OBSERVE FREQUENCY is missing", {
    # Previously this either hard-errored with an unclear message, or (if the
    # NULL check was absent) silently divided the x-axis by NULL, wiping out
    # the data. It must now raise a clear, informative error instead.
    lines <- c(
        "##TITLE=Test4",
        "##XUNITS=HZ",
        "##FIRSTX=0",
        "##LASTX=10",
        "##NPOINTS=2",
        "##DATA TABLE= (X++(Y..Y))",
        "0 1 2",
        "##END="
    )
    expect_error(
        process_block(lines),
        regexp = "\\.OBSERVE FREQUENCY.*missing",
    )
})

test_that("process_block() works fine when XUNITS is 'Hz' and .OBSERVE FREQUENCY is present", {
    lines <- c(
        "##TITLE=Test5",
        "##XUNITS=HZ",
        "##.OBSERVE FREQUENCY=500",
        "##FIRSTX=0",
        "##LASTX=1000",
        "##NPOINTS=2",
        "##DATA TABLE= (X++(Y..Y))",
        "0 1 2",
        "##END="
    )
    result <- process_block(lines)
    data <- result$block[["DATA TABLE"]]
    # x-axis (in Hz, 0 to 1000) is divided by the observe frequency (500) to
    # convert to chemical shift / ppm.
    expect_equal(data$x, c(0, 2))
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
