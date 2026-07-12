make_dataset_1D <- function(with_row_names = TRUE) {
    data_1r <- matrix(
        1:9,
        nrow = 3,
        dimnames = if (with_row_names) list(c("10", "20", "30"), NULL) else NULL
    )
    new_nmr_dataset_1D(
        ppm_axis = c(1, 2, 3),
        data_1r = data_1r,
        metadata = list(external = data.frame(
            NMRExperiment = c("10", "20", "30"),
            Group = c("A", "B", "A")
        ))
    )
}

## new_nmr_dataset_1D / validate_nmr_dataset_1D --------------------------------

test_that("new_nmr_dataset_1D builds a valid object with an empty excluded_regions by default", {
    dataset_1D <- make_dataset_1D()

    expect_s3_class(dataset_1D, c("nmr_dataset_1D", "nmr_dataset_family"))
    expect_equal(dataset_1D$num_samples, 3)
    expect_equal(dataset_1D$axis, c(1, 2, 3))
    expect_equal(dataset_1D$excluded_regions, list())
})

test_that("validate_nmr_dataset_1D rejects an object of the wrong class", {
    dataset_1D <- unclass(make_dataset_1D())
    expect_error(validate_nmr_dataset_1D(dataset_1D), "Not an nmr_dataset_family object")
})

test_that("validate_nmr_dataset_1D requires axis and data_1r elements", {
    dataset_1D <- make_dataset_1D()

    no_axis <- dataset_1D
    no_axis[["axis"]] <- NULL
    expect_error(validate_nmr_dataset_1D(no_axis), "must have a ppm axis")

    no_data <- dataset_1D
    no_data[["data_1r"]] <- NULL
    expect_error(validate_nmr_dataset_1D(no_data), "must have a data_1r matrix")
})

test_that("validate_nmr_dataset_1D requires axis length to match ncol(data_1r)", {
    dataset_1D <- make_dataset_1D()
    dataset_1D$axis <- c(1, 2)

    expect_error(
        validate_nmr_dataset_1D(dataset_1D),
        "ppm axis does not have a length equal to ncol\\(data_1r\\)"
    )
})

test_that("validate_nmr_dataset_1D requires data_1r to be a numeric matrix", {
    dataset_1D <- make_dataset_1D()
    dataset_1D$data_1r <- matrix(letters[1:9], nrow = 3)

    expect_error(validate_nmr_dataset_1D(dataset_1D), "must be a numeric matrix")
})

test_that("validate_nmr_dataset_1D requires num_samples to match nrow(data_1r)", {
    dataset_1D <- make_dataset_1D()
    # Keep num_samples/metadata consistent with each other (so the
    # nmr_dataset_family-level check passes) but make data_1r itself have a
    # different number of rows, to reach the 1D-specific check.
    dataset_1D$data_1r <- dataset_1D$data_1r[1:2, , drop = FALSE]

    expect_error(
        validate_nmr_dataset_1D(dataset_1D),
        "num_samples value does not match nrow\\(data_1r\\)"
    )
})

test_that("validate_nmr_dataset_1D warns and backfills an empty excluded_regions when missing", {
    dataset_1D <- make_dataset_1D()
    dataset_1D$excluded_regions <- NULL

    expect_warning(
        result <- validate_nmr_dataset_1D(dataset_1D),
        "excluded_regions"
    )
    expect_equal(result$excluded_regions, list())
})

## is.nmr_dataset_1D -----------------------------------------------------------

test_that("is.nmr_dataset_1D distinguishes 1D datasets from other objects", {
    expect_true(is.nmr_dataset_1D(make_dataset_1D()))
    expect_false(is.nmr_dataset_1D(list()))
    expect_false(is.nmr_dataset_1D(1:3))
})

## print / format ---------------------------------------------------------------

test_that("format.nmr_dataset_1D summarizes the number of samples", {
    dataset_1D <- make_dataset_1D()
    expect_equal(format(dataset_1D), "An nmr_dataset_1D (3 samples)")
})

test_that("print.nmr_dataset_1D prints the format() text and returns x invisibly", {
    dataset_1D <- make_dataset_1D()
    expect_output(print(dataset_1D), "An nmr_dataset_1D \\(3 samples\\)")
    expect_identical(withr::with_output_sink(nullfile(), print(dataset_1D)), dataset_1D)
})

## [.nmr_dataset_1D ----------------------------------------------------------

test_that("[.nmr_dataset_1D subsets metadata, data_1r and num_samples together", {
    dataset_1D <- make_dataset_1D()

    subset_1D <- dataset_1D[1:2]

    expect_equal(subset_1D$num_samples, 2)
    expect_equal(dim(subset_1D$data_1r), c(2, 3))
    expect_equal(subset_1D$metadata$external$NMRExperiment, c("10", "20"))
})

test_that("[.nmr_dataset_1D also subsets data_1r_baseline when present", {
    dataset_1D <- make_dataset_1D()
    dataset_1D$data_1r_baseline <- matrix(100:108, nrow = 3)

    subset_1D <- dataset_1D[1:2]

    expect_equal(dim(subset_1D$data_1r_baseline), c(2, 3))
    expect_equal(subset_1D$data_1r_baseline, dataset_1D$data_1r_baseline[1:2, , drop = FALSE])
})

test_that("[.nmr_dataset_1D returns a still-valid object", {
    dataset_1D <- make_dataset_1D()
    expect_no_error(validate_nmr_dataset_1D(dataset_1D[c(1, 3)]))
})

## nmr_export_data_1r ------------------------------------------------------------

test_that("nmr_export_data_1r writes the spectra matrix as csv and returns the dataset unchanged", {
    dataset_1D <- make_dataset_1D()
    tmp <- withr::local_tempfile(fileext = ".csv")

    result <- nmr_export_data_1r(dataset_1D, tmp)

    expect_identical(result, dataset_1D)
    written <- utils::read.csv(tmp)
    expect_equal(unname(as.matrix(written)), unname(dataset_1D$data_1r))
})

test_that("nmr_export_data_1r rejects non-nmr_dataset_1D input", {
    expect_error(
        nmr_export_data_1r(list(), tempfile()),
        "An nmr_dataset_1D should be given"
    )
})

## SummarizedExperiment round trip ------------------------------------------------

test_that("nmr_data_1r_to_SummarizedExperiment rejects non-nmr_dataset_1D input", {
    skip_if_not_installed("SummarizedExperiment")
    expect_error(
        nmr_data_1r_to_SummarizedExperiment(list()),
        "An nmr_dataset_1D should be given"
    )
})

test_that("nmr_data_1r_to_SummarizedExperiment builds an SE with data_1r as its assay", {
    skip_if_not_installed("SummarizedExperiment")
    dataset_1D <- make_dataset_1D()

    se <- nmr_data_1r_to_SummarizedExperiment(dataset_1D)

    expect_s4_class(se, "SummarizedExperiment")
    expect_equal(
        unname(SummarizedExperiment::assay(se, "data_1r")),
        unname(dataset_1D$data_1r)
    )
})

test_that("SummarizedExperiment_to_nmr_data_1r round-trips the ppm axis, spectra and external metadata", {
    skip_if_not_installed("SummarizedExperiment")
    skip_if_not_installed("S4Vectors")
    dataset_1D <- make_dataset_1D()

    se <- nmr_data_1r_to_SummarizedExperiment(dataset_1D)
    dataset_1D_2 <- SummarizedExperiment_to_nmr_data_1r(se)

    expect_true(is.nmr_dataset_1D(dataset_1D_2))
    expect_equal(dataset_1D_2$axis, dataset_1D$axis)
    expect_equal(
        unname(dataset_1D_2$data_1r[c("10", "20", "30"), ]),
        unname(dataset_1D$data_1r)
    )
    expect_equal(dataset_1D_2$metadata$external$NMRExperiment, c("10", "20", "30"))
    expect_equal(dataset_1D_2$metadata$external$Group, c("A", "B", "A"))
})
