make_peak_table <- function() {
    peak_table <- matrix(
        c(1, 2, 3, 4, 5, 6),
        nrow = 3,
        dimnames = list(c("10", "20", "30"), c("ppm_1.4", "ppm_1.6"))
    )
    metadata <- list(external = data.frame(
        NMRExperiment = c("10", "20", "30"),
        Group = c("A", "B", "A")
    ))
    new_nmr_dataset_peak_table(peak_table, metadata)
}

## new_nmr_dataset_peak_table / validate_nmr_dataset_peak_table -------------

test_that("new_nmr_dataset_peak_table builds a valid object", {
    pt <- make_peak_table()

    expect_s3_class(pt, c("nmr_dataset_peak_table", "nmr_dataset_family"))
    expect_equal(pt$num_samples, 3)
    expect_equal(dim(pt$peak_table), c(3, 2))
    expect_equal(names(pt), c("10", "20", "30"))
})

test_that("validate_nmr_dataset_peak_table rejects an object of the wrong class", {
    pt <- unclass(make_peak_table())

    expect_error(
        validate_nmr_dataset_peak_table(pt),
        "Not an nmr_dataset_family object"
    )
})

test_that("validate_nmr_dataset_peak_table requires a peak_table element", {
    pt <- make_peak_table()
    pt[["peak_table"]] <- NULL

    expect_error(
        validate_nmr_dataset_peak_table(pt),
        "must have a peak_table matrix"
    )
})

test_that("validate_nmr_dataset_peak_table requires peak_table to be a numeric matrix", {
    pt <- make_peak_table()
    pt[["peak_table"]] <- matrix(letters[1:6], nrow = 3)

    expect_error(
        validate_nmr_dataset_peak_table(pt),
        "must be a numeric matrix"
    )
})

test_that("validate_nmr_dataset_peak_table requires num_samples to match nrow(peak_table)", {
    pt <- make_peak_table()
    # Keep num_samples/metadata consistent with each other (so the
    # nmr_dataset_family-level check passes) but make peak_table itself
    # have a different number of rows, to reach the peak-table-specific check.
    pt[["peak_table"]] <- pt[["peak_table"]][1:2, , drop = FALSE]

    expect_error(
        validate_nmr_dataset_peak_table(pt),
        "num_samples value does not match nrow\\(peak_table\\)"
    )
})

## is.nmr_dataset_peak_table --------------------------------------------------

test_that("is.nmr_dataset_peak_table distinguishes peak tables from other objects", {
    expect_true(is.nmr_dataset_peak_table(make_peak_table()))
    expect_false(is.nmr_dataset_peak_table(list()))
    expect_false(is.nmr_dataset_peak_table(1:3))
})

## print / format --------------------------------------------------------------

test_that("format.nmr_dataset_peak_table summarizes samples and peaks", {
    pt <- make_peak_table()
    expect_equal(format(pt), "An nmr_dataset_peak_table (3 samples, and 2 peaks)")
})

test_that("print.nmr_dataset_peak_table prints the format() text and returns x invisibly", {
    pt <- make_peak_table()
    expect_output(print(pt), "An nmr_dataset_peak_table \\(3 samples, and 2 peaks\\)")
    expect_identical(withr::with_output_sink(nullfile(), print(pt)), pt)
})

## [.nmr_dataset_peak_table ----------------------------------------------------

test_that("[.nmr_dataset_peak_table subsets samples in metadata, peak_table and num_samples together", {
    pt <- make_peak_table()

    subset_pt <- pt[1:2]

    expect_equal(subset_pt$num_samples, 2)
    expect_equal(rownames(subset_pt$peak_table), c("10", "20"))
    expect_equal(subset_pt$metadata$external$NMRExperiment, c("10", "20"))
    # The dropped sample is really gone, not just hidden:
    expect_false("30" %in% rownames(subset_pt$peak_table))
})

test_that("[.nmr_dataset_peak_table returns a still-valid object", {
    pt <- make_peak_table()
    subset_pt <- pt[c(1, 3)]
    expect_no_error(validate_nmr_dataset_peak_table(subset_pt))
})

## as.data.frame.nmr_dataset_peak_table ----------------------------------------

test_that("as.data.frame.nmr_dataset_peak_table combines external metadata and the peak table", {
    pt <- make_peak_table()

    df <- as.data.frame(pt)

    expect_equal(rownames(df), c("10", "20", "30"))
    expect_equal(df$NMRExperiment, c("10", "20", "30"))
    expect_equal(df$Group, c("A", "B", "A"))
    expect_equal(df$ppm_1.4, c(1, 2, 3))
    expect_equal(df$ppm_1.6, c(4, 5, 6))
})

## SummarizedExperiment round trip ---------------------------------------------

test_that("nmr_dataset_peak_table_to_SummarizedExperiment rejects non-peak-table input", {
    expect_error(
        nmr_dataset_peak_table_to_SummarizedExperiment(list()),
        "Not an nmr_dataset_peak_table"
    )
})

test_that("nmr_dataset_peak_table_to_SummarizedExperiment builds an SE with the peak table as its assay", {
    skip_if_not_installed("SummarizedExperiment")
    pt <- make_peak_table()

    se <- nmr_dataset_peak_table_to_SummarizedExperiment(pt)

    expect_s4_class(se, "SummarizedExperiment")
    expect_equal(dim(se), c(3, 2)) # 3 samples (rows) x 2 peaks (cols) in the nmr_dataset_peak_table sense
    expect_equal(
        unname(SummarizedExperiment::assay(se, "peak_table")),
        unname(pt$peak_table)
    )
})

test_that("SummarizedExperiment_to_nmr_dataset_peak_table round-trips peak values and external metadata", {
    skip_if_not_installed("SummarizedExperiment")
    skip_if_not_installed("S4Vectors")
    pt <- make_peak_table()

    se <- nmr_dataset_peak_table_to_SummarizedExperiment(pt)
    pt2 <- SummarizedExperiment_to_nmr_dataset_peak_table(se)

    expect_true(is.nmr_dataset_peak_table(pt2))
    expect_equal(pt2$num_samples, 3)
    expect_setequal(colnames(pt2$peak_table), colnames(pt$peak_table))
    expect_equal(
        pt2$peak_table[c("10", "20", "30"), colnames(pt$peak_table)],
        pt$peak_table[c("10", "20", "30"), ],
        ignore_attr = TRUE
    )
    expect_equal(pt2$metadata$external$NMRExperiment, c("10", "20", "30"))
    expect_equal(pt2$metadata$external$Group, c("A", "B", "A"))
})
