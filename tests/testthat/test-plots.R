make_large_nmr_dataset_1D <- function(n_samples) {
    ppm_axis <- seq(0, 1, length.out = 10)
    data_1r <- matrix(
        stats::runif(n_samples * length(ppm_axis)),
        nrow = n_samples,
        ncol = length(ppm_axis)
    )
    metadata <- list(
        external = data.frame(
            NMRExperiment = as.character(seq_len(n_samples)),
            stringsAsFactors = FALSE
        )
    )
    new_nmr_dataset_1D(
        ppm_axis = ppm_axis,
        data_1r = data_1r,
        metadata = metadata
    )
}

test_that("plot() notifies the user when more than 20 samples are subsampled to a random 10", {
    dataset <- make_large_nmr_dataset_1D(25)
    expect_true(dataset[["num_samples"]] > 20)

    expect_message(
        plot(dataset),
        regexp = "more than 20 samples.*random subset of 10"
    )
})

test_that("plot() does not emit a subsampling notice when there are 20 samples or fewer", {
    dataset <- make_large_nmr_dataset_1D(10)
    expect_true(dataset[["num_samples"]] <= 20)

    expect_no_message(
        plot(dataset),
        message = "more than 20 samples"
    )
})

test_that("plot_with_aes_string() also emits the subsampling notice", {
    dataset <- make_large_nmr_dataset_1D(25)

    # Passing a character (string) aes argument routes through the
    # deprecated aes_string code path inside plot.nmr_dataset_1D(), which
    # also warns about the deprecation itself.
    expect_message(
        suppressWarnings(plot(dataset, colour = "NMRExperiment")),
        regexp = "more than 20 samples.*random subset of 10"
    )
})

test_that("plot() rejects a chemshift_range that isn't length 2 or 3", {
    dataset <- make_large_nmr_dataset_1D(3)
    expect_error(
        plot(dataset, chemshift_range = c(1, 2, 3, 4)),
        "length 2 or 3"
    )
})

## is_using_aes_string ------------------------------------------------------

test_that("is_using_aes_string detects a character (deprecated aes_string) argument", {
    expect_false(is_using_aes_string())
    expect_true(is_using_aes_string(color = "Group"))
    expect_false(is_using_aes_string(color = quote(Group)))
})

## get_vars_from_aes_string / get_vars_from_aes -----------------------------

test_that("get_vars_from_aes_string extracts the variable names referenced in aes_string expressions", {
    expect_equal(get_vars_from_aes_string(character(0)), character(0))
    expect_equal(
        get_vars_from_aes_string(c("Group", "as.factor(Batch)")),
        c("Group", "Batch")
    )
})

test_that("get_vars_from_aes extracts the variable names referenced in aes expressions", {
    f <- function(...) get_vars_from_aes(...)
    expect_equal(f(), character(0))
    expect_equal(f(color = Group, y = log(Intensity)), c("Group", "Intensity"))
})

## prepare_aes ----------------------------------------------------------------

test_that("prepare_aes always sets x/y/group and defaults colour to NMRExperiment", {
    f <- function(...) prepare_aes(...)

    default_aes <- f()
    expect_setequal(names(default_aes), c("colour", "x", "y", "group"))

    with_color <- f(color = Group)
    expect_setequal(names(with_color), c("color", "x", "y", "group"))

    with_colour <- f(colour = Group)
    expect_setequal(names(with_colour), c("colour", "x", "y", "group"))
})

## decimate_axis --------------------------------------------------------------

test_that("decimate_axis keeps every point when there is no resolution target", {
    xaxis <- seq(0, 10, by = 0.1)
    expect_equal(sum(decimate_axis(xaxis, NULL)), length(xaxis))
})

test_that("decimate_axis restricts to a 2-value range without decimating", {
    xaxis <- seq(0, 10, by = 0.1)
    in_range <- decimate_axis(xaxis, c(2, 5))
    expect_equal(range(xaxis[in_range]), c(2, 5))
    expect_equal(sum(in_range), sum(xaxis >= 2 & xaxis <= 5))
})

test_that("decimate_axis subsamples to approximately the requested resolution when given a 3-value range", {
    xaxis <- seq(0, 10, by = 0.1) # 101 points in [0, 10]
    in_range <- decimate_axis(xaxis, c(0, 10, 1)) # want ~1 point per unit
    expect_lt(sum(in_range), 101)
    expect_lte(sum(in_range), 11)
})

test_that("decimate_axis fills in a missing range bound from the axis min/max", {
    xaxis <- seq(0, 10, by = 0.1)
    expect_equal(range(xaxis[decimate_axis(xaxis, c(NA, 5))]), c(0, 5))
    expect_equal(range(xaxis[decimate_axis(xaxis, c(5, NA))]), c(5, 10))
})

## tidy.nmr_dataset_1D --------------------------------------------------------

test_that("tidy.nmr_dataset_1D returns a long data frame with NMRExperiment/chemshift/intensity", {
    dataset <- make_large_nmr_dataset_1D(3)
    df <- tidy(dataset, chemshift_range = c(0.2, 0.6))

    expect_equal(colnames(df), c("NMRExperiment", "chemshift", "intensity"))
    expect_setequal(unique(df$NMRExperiment), names(dataset))
    expect_true(all(df$chemshift >= 0.2 & df$chemshift <= 0.6))
})

test_that("tidy.nmr_dataset_1D warns and excludes unknown NMRExperiment values when some are valid", {
    # Regression test for https://github.com/sipss/AlpsNMR/issues/69: an
    # unknown NMRExperiment used to silently produce NA-filled rows labelled
    # with the wrong (nonexistent) name instead of warning and dropping them.
    dataset <- make_large_nmr_dataset_1D(3)

    expect_warning(
        df <- tidy(dataset, NMRExperiment = c(names(dataset)[1], "does-not-exist")),
        "does-not-exist"
    )

    expect_equal(unique(df$NMRExperiment), names(dataset)[1])
    expect_false(anyNA(df$intensity))
})

test_that("tidy.nmr_dataset_1D errors when every given NMRExperiment value is unknown", {
    dataset <- make_large_nmr_dataset_1D(3)

    expect_error(
        tidy(dataset, NMRExperiment = c("does-not-exist-1", "does-not-exist-2")),
        "None of the given NMRExperiment values"
    )
})

test_that("tidy.nmr_dataset_1D's error/warning suggest some valid NMRExperiment values", {
    dataset <- make_large_nmr_dataset_1D(3)
    valid_examples <- paste(utils::head(names(dataset), 2), collapse = ", ")

    expect_error(
        tidy(dataset, NMRExperiment = "does-not-exist"),
        valid_examples,
        fixed = TRUE
    )
})

## plot_interactive / plot_webgl ----------------------------------------------

test_that("plot_interactive writes an html file and a lib/ folder, and aborts on a pre-existing lib/ unless told to overwrite", {
    skip_if_not_installed("plotly")
    skip_if_not_installed("htmltools")
    dataset <- make_large_nmr_dataset_1D(3)
    plt <- suppressMessages(plot(dataset))

    tmpdir <- withr::local_tempdir()
    html_file <- file.path(tmpdir, "plot.html")

    result <- suppressWarnings(plot_interactive(plt, html_file))
    expect_equal(result, html_file)
    expect_true(file.exists(html_file))
    expect_true(dir.exists(file.path(tmpdir, "lib")))

    # Calling again non-interactively with overwrite = NULL (the default)
    # must abort rather than silently clobber the existing lib/ folder:
    expect_error(
        suppressWarnings(plot_interactive(plt, html_file)),
        "already exists"
    )

    # overwrite = TRUE proceeds anyway:
    expect_true(file.exists(suppressWarnings(plot_interactive(plt, html_file, overwrite = TRUE))))
})

test_that("plot_webgl plots the dataset and writes it out via plot_interactive", {
    skip_if_not_installed("plotly")
    skip_if_not_installed("htmltools")
    dataset <- make_large_nmr_dataset_1D(3)

    tmpdir <- withr::local_tempdir()
    html_file <- file.path(tmpdir, "plot.html")

    result <- suppressWarnings(suppressMessages(plot_webgl(dataset, html_file)))

    expect_equal(result, html_file)
    expect_true(file.exists(html_file))
})
