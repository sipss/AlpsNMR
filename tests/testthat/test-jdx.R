test_that("comments are stripped", {
    expect_equal(
        strip_comments(c(
            "##TITLE=Hello $$ world",
            "$$ Hello world",
            "12343 $$ hello world"
        )),
        c("##TITLE=Hello ", "", "12343 ")
    )
})

test_that("comments are stripped", {
    expect_equal(
        strip_comments(c(
            "##TITLE=Hello $$ world",
            "$$ Hello world",
            "12343 $$ hello world"
        )),
        c("##TITLE=Hello ", "", "12343 ")
    )
})

test_that("process_block() (X++(Y..Y)) parser ignores single-token data lines", {
    # A data line with only one token (a stray value with no accompanying Y
    # value) contributes zero values to the parsed data instead of being
    # padded with NA.
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
