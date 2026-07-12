## permutation_test_model --------------------------------------------------

test_that("permutation_test_model returns an nPerm x 1 matrix of AUCs on label-permuted data", {
    skip_if_not_installed("mixOmics")
    skip_if_not_installed("BiocParallel")

    old_bpparam <- BiocParallel::bpparam()
    BiocParallel::register(BiocParallel::SerialParam())
    on.exit(BiocParallel::register(old_bpparam), add = TRUE)

    set.seed(1)
    n <- 16
    p <- 3
    y <- factor(rep(c("A", "B"), each = n / 2))
    x <- matrix(rnorm(n * p), nrow = n, ncol = p)
    colnames(x) <- paste0("V", seq_len(p))
    metadata <- data.frame(NMRExperiment = as.character(seq_len(n)), Condition = y)
    dataset <- new_nmr_dataset_peak_table(peak_table = x, metadata = list(external = metadata))

    methodology <- plsda_auroc_vip_method(ncomp = 1)
    nPerm <- 3

    permMatrix <- permutation_test_model(
        dataset,
        y_column = "Condition",
        identity_column = NULL,
        external_val = list(iterations = 1, test_size = 0.25),
        internal_val = list(iterations = 1, test_size = 0.25),
        data_analysis_method = methodology,
        nPerm = nPerm
    )

    expect_true(is.matrix(permMatrix))
    expect_equal(dim(permMatrix), c(nPerm, 1))
    expect_true(all(permMatrix >= 0 & permMatrix <= 1))
})

test_that("permutation_test_model permutes the y_column and leaves the rest of the dataset untouched", {
    # permutation_test_model() must reassign only the y_column values (a
    # shuffle of the original ones), not the row count/other metadata.
    skip_if_not_installed("BiocParallel")

    old_bpparam <- BiocParallel::bpparam()
    BiocParallel::register(BiocParallel::SerialParam())
    on.exit(BiocParallel::register(old_bpparam), add = TRUE)

    n <- 8
    y <- factor(rep(c("A", "B"), each = n / 2))
    metadata <- data.frame(
        NMRExperiment = as.character(seq_len(n)),
        Condition = y,
        Extra = seq_len(n)
    )
    dataset <- new_nmr_dataset_peak_table(
        peak_table = matrix(seq_len(n * 2), nrow = n, ncol = 2, dimnames = list(NULL, c("P1", "P2"))),
        metadata = list(external = metadata)
    )

    captured_y <- NULL
    fake_analysis <- function(dataset_perm, y_column, identity_column, external_val, internal_val,
        data_analysis_method, .enable_parallel) {
        captured_y[[length(captured_y) + 1L]] <<- nmr_meta_get_column(dataset_perm, column = y_column)
        list(outer_cv_results_digested = list(auroc = data.frame(auc = 0.5)))
    }
    testthat::local_mocked_bindings(nmr_data_analysis = fake_analysis)

    permutation_test_model(
        dataset,
        y_column = "Condition",
        identity_column = NULL,
        external_val = list(iterations = 1, test_size = 0.25),
        internal_val = list(iterations = 1, test_size = 0.25),
        data_analysis_method = NULL,
        nPerm = 2
    )

    expect_length(captured_y, 2)
    for (y_perm in captured_y) {
        expect_setequal(as.character(y_perm), as.character(y))
        expect_equal(length(y_perm), n)
    }
})

## permutation_test_plot ----------------------------------------------------

test_that("permutation_test_plot defaults xlim to c(0, 1)", {
    captured_hist_args <- NULL
    testthat::local_mocked_bindings(
        hist = function(...) {
            captured_hist_args <<- list(...)
            list(density = c(0.1, 0.2))
        },
        text = function(...) invisible(NULL),
        axis = function(...) invisible(NULL),
        lines = function(...) invisible(NULL)
    )

    permMatrix <- matrix(c(0.4, 0.5, 0.6), ncol = 1)
    model <- list(outer_cv_results_digested = list(auroc = data.frame(auc = 0.5)))

    permutation_test_plot(model, permMatrix)

    expect_equal(captured_hist_args$xlim, c(0, 1))
})

test_that("permutation_test_plot labels an extreme observed AUC as 'p<' the smallest achievable p-value", {
    captured_text <- NULL
    testthat::local_mocked_bindings(
        hist = function(...) list(density = c(0.1, 0.2, 0.3)),
        text = function(...) {
            captured_text <<- list(...)
            invisible(NULL)
        },
        axis = function(...) invisible(NULL),
        lines = function(...) invisible(NULL)
    )

    # The observed AUC (0.91) is higher than every permuted AUC:
    permMatrix <- matrix(c(0.48, 0.49, 0.5, 0.51, 0.52), ncol = 1)
    model <- list(outer_cv_results_digested = list(auroc = data.frame(auc = c(0.9, 0.92))))

    permutation_test_plot(model, permMatrix)

    expect_equal(captured_text$labels, paste0("p<", signif(1 / nrow(permMatrix), 4)))
})

test_that("permutation_test_plot labels a typical observed AUC with its two-sided empirical p-value", {
    captured_text <- NULL
    testthat::local_mocked_bindings(
        hist = function(...) list(density = c(0.1, 0.2, 0.3)),
        text = function(...) {
            captured_text <<- list(...)
            invisible(NULL)
        },
        axis = function(...) invisible(NULL),
        lines = function(...) invisible(NULL)
    )

    permMatrix <- matrix(seq(0.1, 0.95, length.out = 10), ncol = 1)
    # model_auc = 0.8: 8 of the 10 permuted values are <= 0.8, so the
    # empirical (right-tail) p-value is 1 - 0.8 = 0.2:
    model <- list(outer_cv_results_digested = list(auroc = data.frame(auc = 0.8)))

    permutation_test_plot(model, permMatrix)

    expect_equal(captured_text$labels, "p=0.2")
})

test_that("permutation_test_plot draws the observed AUC as a vertical line at its x-position", {
    captured_lines_args <- NULL
    testthat::local_mocked_bindings(
        hist = function(...) list(density = c(0.1, 0.2, 0.3)),
        text = function(...) invisible(NULL),
        axis = function(...) invisible(NULL),
        lines = function(...) {
            captured_lines_args <<- list(...)
            invisible(NULL)
        }
    )

    permMatrix <- matrix(c(0.4, 0.5, 0.6), ncol = 1)
    model <- list(outer_cv_results_digested = list(auroc = data.frame(auc = c(0.7, 0.75))))

    permutation_test_plot(model, permMatrix)

    expect_equal(captured_lines_args[[1]], rep(0.725, 2))
})

test_that("permutation_test_plot runs against the real base graphics device without leaving one open", {
    skip_on_cran()
    t <- tempfile(fileext = ".pdf")
    grDevices::pdf(t)
    devices_before <- grDevices::dev.list()
    on.exit(
        {
            grDevices::dev.off()
            unlink(t)
        },
        add = TRUE
    )

    permMatrix <- matrix(c(0.4, 0.45, 0.5, 0.55, 0.6), ncol = 1)
    model <- list(outer_cv_results_digested = list(auroc = data.frame(auc = 0.5)))

    expect_no_error(permutation_test_plot(model, permMatrix))
    expect_equal(grDevices::dev.list(), devices_before)
})
