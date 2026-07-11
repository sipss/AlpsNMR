#' Build a tiny, fast PLS-DA model for use as test fixture input to the
#' plot_plsda_* functions.
build_test_plsda_model <- function(seed = 1, n = 20, p = 4, ncomp = 1) {
    set.seed(seed)
    x <- matrix(rnorm(n * p), nrow = n, ncol = p)
    colnames(x) <- paste0("V", seq_len(p))
    y <- factor(rep(c("A", "B"), each = n / 2))
    x[y == "A", 1] <- x[y == "A", 1] + 5

    model <- plsda_build(x = x, y = y, identity = NULL, ncomp = ncomp)
    # plot_plsda_samples()/plot_plsda_multimodel() expect these attached,
    # matching what nmr_data_analysis()'s outer CV loop attaches (see
    # callback_plsda_auroc_vip() in R/plsda.R):
    model$X_test <- x
    model$Y_test <- y
    model
}

test_that("plot_plsda_samples closes its graphics device even when plotting fails", {
    # plot_plsda_samples() opens a hidden pdf(file = tempfile()) device
    # before calling mixOmics::plotIndiv(). The graphics device opened for
    # plotting must always be closed, even if plotting fails, or the device
    # (and its temp file) would be leaked.
    #
    # We force a failure *inside* mixOmics::plotIndiv() (deep enough that it
    # happens after the pdf device has already been opened) by mocking the
    # plotIndiv binding inside the mixOmics namespace (this affects all
    # mixOmics::plotIndiv() calls for the duration of the test only), then
    # confirm no extra open graphics device is left behind.
    skip_if_not_installed("mixOmics")
    model <- build_test_plsda_model(seed = 1, ncomp = 1)

    devices_before <- dev.list()

    testthat::local_mocked_bindings(
        plotIndiv = function(...) stop("forced failure for on.exit test"),
        .package = "mixOmics"
    )
    expect_error(plot_plsda_samples(model), "forced failure for on.exit test")

    expect_equal(dev.list(), devices_before)
})

test_that("plot_plsda_multimodel closes its graphics device even when plotting fails", {
    skip_if_not_installed("mixOmics")
    model <- list(
        outer_cv_results = list(
            list(model = build_test_plsda_model(seed = 1, ncomp = 1)),
            list(model = build_test_plsda_model(seed = 2, ncomp = 1))
        )
    )

    devices_before <- dev.list()

    testthat::local_mocked_bindings(
        plotIndiv = function(...) stop("forced failure for on.exit test (multimodel)"),
        .package = "mixOmics"
    )
    expect_error(
        plot_plsda_multimodel(model),
        "forced failure for on.exit test \\(multimodel\\)"
    )

    expect_equal(dev.list(), devices_before)
})

test_that("plot_plsda_samples leaves no open graphics device on the successful path", {
    skip_if_not_installed("mixOmics")
    model <- build_test_plsda_model(seed = 3, ncomp = 1)

    devices_before <- dev.list()
    plsda_plot <- plot_plsda_samples(model, plot = FALSE)
    devices_after <- dev.list()

    expect_s3_class(plsda_plot, "ggplot")
    expect_equal(devices_after, devices_before)
})
