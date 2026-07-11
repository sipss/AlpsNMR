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
