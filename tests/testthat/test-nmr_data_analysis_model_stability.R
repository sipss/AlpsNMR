## models_stability_plot_plsda ----------------------------------------------

fake_plsda_model <- function(ncomp, loadings) {
    list(ncomp = ncomp, loadings = list(X = loadings))
}

test_that("models_stability_plot_plsda compares loadings pairwise across models and components", {
    skip_if_not_installed("ggplot2")
    set.seed(1)
    m1 <- fake_plsda_model(2, matrix(rnorm(4 * 2), nrow = 4, ncol = 2, dimnames = list(paste0("F", 1:4), NULL)))
    m2 <- fake_plsda_model(2, matrix(rnorm(4 * 2), nrow = 4, ncol = 2, dimnames = list(paste0("F", 1:4), NULL)))
    model <- list(outer_cv_results = list(list(model = m1), list(model = m2)))

    p <- expect_no_warning(models_stability_plot_plsda(model))

    expect_s3_class(p, "ggplot")
    # 2 models x 2 latent variables = 4 rows/cols, fully populated (no NAs to drop):
    expect_equal(nrow(p$data), 16)
})

test_that("models_stability_plot_plsda drops the cells for latent variables a model doesn't have", {
    skip_if_not_installed("ggplot2")
    set.seed(1)
    # model 1 has 2 components, model 2 only has 1: the "model 2 - LV2" row/col
    # is entirely NA and must be dropped rather than erroring.
    m1 <- fake_plsda_model(2, matrix(rnorm(4 * 2), nrow = 4, ncol = 2, dimnames = list(paste0("F", 1:4), NULL)))
    m2 <- fake_plsda_model(1, matrix(rnorm(4 * 1), nrow = 4, ncol = 1, dimnames = list(paste0("F", 1:4), NULL)))
    model <- list(outer_cv_results = list(list(model = m1), list(model = m2)))

    p <- expect_no_error(models_stability_plot_plsda(model))

    expect_s3_class(p, "ggplot")
    # 3 valid (model, LV) combinations remain -> 3x3 = 9 rows after dropping NAs:
    expect_equal(nrow(p$data), 9)
})

## models_stability_plot_bootstrap --------------------------------------------

test_that("models_stability_plot_bootstrap compares general_model loadings pairwise across folds", {
    skip_if_not_installed("ggplot2")
    set.seed(1)
    m1 <- fake_plsda_model(1, matrix(rnorm(4), nrow = 4, ncol = 1, dimnames = list(paste0("F", 1:4), NULL)))
    m2 <- fake_plsda_model(1, matrix(rnorm(4), nrow = 4, ncol = 1, dimnames = list(paste0("F", 1:4), NULL)))
    bp_results <- list(kfold_results = list(list(general_model = m1), list(general_model = m2)))

    p <- expect_no_warning(models_stability_plot_bootstrap(bp_results))

    expect_s3_class(p, "ggplot")
    # 2 folds x 1 latent variable = 2 rows/cols -> 4 rows after melt:
    expect_equal(nrow(p$data), 4)
})

## plot_bootstrap_multimodel --------------------------------------------------

build_multimodel_test_setup <- function(ncomp) {
    set.seed(1)
    n <- 12
    p <- 3
    y <- factor(rep(c("A", "B"), each = n / 2))
    x <- matrix(rnorm(n * p), nrow = n, ncol = p)
    colnames(x) <- paste0("V", seq_len(p))
    metadata <- data.frame(NMRExperiment = as.character(seq_len(n)), Condition = y)
    dataset <- new_nmr_dataset_peak_table(peak_table = x, metadata = list(external = metadata))

    idx1 <- 1:8
    idx2 <- 5:12
    m1 <- plsda_build(x[idx1, ], y[idx1], identity = NULL, ncomp = ncomp)
    m2 <- plsda_build(x[idx2, ], y[idx2], identity = NULL, ncomp = ncomp)

    list(
        dataset = dataset,
        bp_results = list(
            kfold_results = list(list(general_model = m1), list(general_model = m2)),
            kfold_index = list(idx1, idx2)
        )
    )
}

test_that("plot_bootstrap_multimodel builds a histogram of test-set scores for a single-component model", {
    skip_if_not_installed("mixOmics")
    setup <- build_multimodel_test_setup(ncomp = 1)

    p <- expect_no_warning(
        plot_bootstrap_multimodel(setup$bp_results, setup$dataset, "Condition", plot = FALSE)
    )

    expect_s3_class(p, "ggplot")
    expect_true(inherits(p$layers[[1]]$geom, "GeomBar")) # geom_histogram is built on GeomBar
})

test_that("plot_bootstrap_multimodel builds a scores scatterplot for a multi-component model", {
    # Regression test: geom_hline()/geom_vline() used a `size` argument, which
    # is deprecated in favor of `linewidth` since ggplot2 3.4.0 and emitted a
    # warning on every call.
    skip_if_not_installed("mixOmics")
    setup <- build_multimodel_test_setup(ncomp = 2)

    p <- expect_no_warning(
        plot_bootstrap_multimodel(setup$bp_results, setup$dataset, "Condition", plot = FALSE)
    )

    expect_s3_class(p, "ggplot")
    expect_true(inherits(p$layers[[1]]$geom, "GeomHline"))
})

test_that("plot_bootstrap_multimodel's plot argument doesn't change the returned object", {
    skip_if_not_installed("mixOmics")
    setup <- build_multimodel_test_setup(ncomp = 1)

    p_true <- plot_bootstrap_multimodel(setup$bp_results, setup$dataset, "Condition", plot = TRUE)
    p_false <- plot_bootstrap_multimodel(setup$bp_results, setup$dataset, "Condition", plot = FALSE)

    expect_equal(p_true$data, p_false$data)
})
