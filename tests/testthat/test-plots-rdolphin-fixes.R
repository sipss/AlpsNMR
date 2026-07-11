# Regression tests for fixes made to R/plots.R and R/to_rDolphin_blood.R
#
# 1. plots.R: plot_with_aes()/plot_with_aes_string() now emit a visible
#    cli::cli_inform() message when a dataset with more than 20 samples is
#    plotted without an explicit NMRExperiment (previously this silent
#    subsampling to 10 samples happened with no notice at all).
# 2. to_rDolphin_blood.R: files_to_rDolphin() now emits a message stating
#    the label -> numeric code mapping used to encode the Group/type column
#    (previously as.numeric(as.factor(...)) silently picked codes based on
#    alphabetical order with no indication to the user).

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

test_that("plot() notifies the user when subsampling to 10 samples (>20 samples)", {
    dataset <- make_large_nmr_dataset_1D(25)
    expect_true(dataset[["num_samples"]] > 20)

    expect_message(
        plot(dataset),
        regexp = "more than 20 samples.*random subset of 10"
    )
})

test_that("plot() does NOT emit the subsampling notice when <= 20 samples", {
    dataset <- make_large_nmr_dataset_1D(10)
    expect_true(dataset[["num_samples"]] <= 20)

    expect_no_message(
        plot(dataset),
        message = "more than 20 samples"
    )
})

test_that("plot_with_aes_string() (deprecated path) also emits the subsampling notice", {
    dataset <- make_large_nmr_dataset_1D(25)

    # Passing a character (string) aes argument routes through the
    # deprecated aes_string code path inside plot.nmr_dataset_1D(), which
    # also warns about the deprecation itself.
    expect_message(
        suppressWarnings(plot(dataset, colour = "NMRExperiment")),
        regexp = "more than 20 samples.*random subset of 10"
    )
})

test_that("files_to_rDolphin() reports the Group label -> numeric code mapping", {
    dataset <- new_nmr_dataset_1D(
        ppm_axis = 1:10,
        data_1r = matrix(
            stats::runif(9 * 10),
            nrow = 9
        ),
        metadata = list(
            external = data.frame(
                NMRExperiment = c("10", "20", "30", "40", "50", "60", "70", "80", "90"),
                stringsAsFactors = FALSE
            )
        )
    )
    dataset[["metadata"]][["external"]][["SubjectID"]] <-
        c("10", "20", "30", "40", "50", "60", "70", "80", "90")
    dataset[["metadata"]][["external"]][["Group"]] <-
        c("Control", "Control", "Control", "Control", "Control", "Patient", "Patient", "Patient", "Patient")

    messages <- testthat::capture_messages(
        result <- files_to_rDolphin(dataset, "blood")
    )
    mapping_message <- paste(messages, collapse = "\n")

    # The message must mention both original labels together with their
    # assigned numeric codes, so the user can see how the encoding was done.
    expect_match(mapping_message, "Control\\s*=\\s*1")
    expect_match(mapping_message, "Patient\\s*=\\s*2")

    # Sanity check that the returned data actually uses that same mapping.
    expect_equal(
        result[["meta_rDolphin"]][["type"]],
        c(1, 1, 1, 1, 1, 2, 2, 2, 2)
    )
})
