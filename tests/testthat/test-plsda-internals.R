skip_if_not_installed("mixOmics")

build_plsda_xy <- function(seed = 1, n = 40, p = 4) {
    set.seed(seed)
    x <- matrix(rnorm(n * p), nrow = n, ncol = p)
    colnames(x) <- paste0("V", seq_len(p))
    y <- factor(rep(c("A", "B"), each = n / 2))
    x[y == "A", 1] <- x[y == "A", 1] + 5
    list(x = x, y = y)
}

test_that("plsda_build rejects a number of components larger than the number of features", {
    xy <- build_plsda_xy()
    expect_error(
        plsda_build(xy$x, xy$y, identity = NULL, ncomp = ncol(xy$x) + 1),
        "can't"
    )
})

test_that("plsda_build returns a mixOmics plsda model", {
    xy <- build_plsda_xy()
    model <- plsda_build(xy$x, xy$y, identity = NULL, ncomp = 2)
    expect_s3_class(model, "mixo_plsda")
    expect_equal(model$ncomp, 2)
})

test_that("plsda_vip returns the variable importance in projection matrix", {
    xy <- build_plsda_xy()
    model <- plsda_build(xy$x, xy$y, identity = NULL, ncomp = 2)
    vip <- plsda_vip(model)
    expect_equal(dim(vip), c(ncol(xy$x), 2))
    expect_equal(rownames(vip), colnames(xy$x))
})

test_that("plsda_vip returns NULL and informs the user when the underlying call fails", {
    testthat::local_mocked_bindings(
        vip = function(...) stop("forced vip failure"),
        .package = "mixOmics"
    )
    expect_message(result <- plsda_vip(list()), "VIP calculation failed")
    expect_null(result)
})

test_that("choose_best_nlv_impl returns integer(0) for an empty input", {
    empty <- data.frame(ncomp = integer(0), auc = numeric(0))
    expect_equal(choose_best_nlv_impl(empty, auc_threshold = 0.05), integer(0))
})

test_that("choose_best_nlv_impl returns the only ncomp when there is a single row", {
    single <- data.frame(ncomp = 3, auc = 0.7)
    expect_equal(choose_best_nlv_impl(single, auc_threshold = 0.05), 3)
})

test_that("choose_best_nlv_impl returns the first ncomp when all aucs are NA", {
    all_na <- data.frame(ncomp = c(1, 2, 3), auc = c(NA, NA, NA))
    expect_equal(choose_best_nlv_impl(all_na, auc_threshold = 0.05), 1)
})

test_that("choose_best_nlv_impl stops increasing complexity once the auc increment falls below the threshold", {
    ncomp_auc <- data.frame(ncomp = c(1, 2, 3, 4), auc = c(0.5, 0.6, 0.65, 0.66))
    # 1->2: +0.10 (keep going), 2->3: +0.05 (keep going, not < threshold), 3->4: +0.01 (< threshold, stop)
    expect_equal(choose_best_nlv_impl(ncomp_auc, auc_threshold = 0.05), 3)
})

test_that("choose_best_nlv_impl returns the highest ncomp when every increment clears the threshold", {
    ncomp_auc <- data.frame(ncomp = c(1, 2, 3), auc = c(0.5, 0.6, 0.8))
    expect_equal(choose_best_nlv_impl(ncomp_auc, auc_threshold = 0.05), 3)
})

build_inner_cv_results <- function() {
    list(
        "1_1" = list(auroc = data.frame(ncomp = c(1, 2, 3), auc = c(0.50, 0.60, 0.65))),
        "1_2" = list(auroc = data.frame(ncomp = c(1, 2, 3), auc = c(0.55, 0.62, 0.63))),
        "2_1" = list(auroc = data.frame(ncomp = c(1, 2, 3), auc = c(0.50, 0.70, 0.72))),
        "2_2" = list(auroc = data.frame(ncomp = c(1, 2, 3), auc = c(0.50, 0.65, 0.90)))
    )
}

test_that("choose_best_nlv picks a number of latent variables per outer iteration and returns diagnostics", {
    inner_cv_results <- build_inner_cv_results()
    result <- choose_best_nlv(inner_cv_results, auc_threshold = 0.05)

    expect_equal(result$num_latent_var$cv_outer_iteration, c(1, 2))
    expect_equal(result$train_evaluate_model_args, list(ncomp = result$num_latent_var$ncomp))
    expect_s3_class(result$diagnostic_plot, "ggplot")
    expect_s3_class(result$diagnostic_box_plot, "ggplot")
    expect_s3_class(result$model_performances, "data.frame")
    expect_setequal(
        colnames(result$model_performances),
        c("cv_outer_iteration", "cv_inner_iteration", "ncomp", "auc")
    )
    expect_equal(nrow(result$model_performances), 4 * 3)
})

test_that("fun_choose_best_ncomp_auc_threshold builds a closure that fixes the auc_threshold", {
    inner_cv_results <- build_inner_cv_results()
    choose_fn <- fun_choose_best_ncomp_auc_threshold(auc_threshold = 0.05)
    expect_equal(
        choose_fn(inner_cv_results)$num_latent_var,
        choose_best_nlv(inner_cv_results, auc_threshold = 0.05)$num_latent_var
    )
})

test_that("callback_outer_cv_auroc_vip keeps the auroc at the highest ncomp and computes VIP rank products", {
    outer_cv_results <- list(
        "1" = list(
            auroc = data.frame(ncomp = c(1, 2), auc = c(0.6, 0.7)),
            vip = matrix(c(1, 2, 3, 4, 5, 6), nrow = 3, ncol = 2, dimnames = list(c("f1", "f2", "f3"), NULL))
        ),
        "2" = list(
            auroc = data.frame(ncomp = c(1, 2), auc = c(0.55, 0.75)),
            vip = matrix(c(2, 1, 3, 5, 4, 6), nrow = 3, ncol = 2, dimnames = list(c("f1", "f2", "f3"), NULL))
        )
    )
    result <- callback_outer_cv_auroc_vip(outer_cv_results)

    expect_equal(result$auroc$ncomp, c(2, 2))
    expect_equal(result$auroc$auc, c(0.7, 0.75))
    expect_equal(result$vip_vectors[["1"]], c(f1 = 4, f2 = 5, f3 = 6))
    expect_equal(result$vip_vectors[["2"]], c(f1 = 5, f2 = 4, f3 = 6))
    # f3 is ranked worst (highest VIP -> best rank 1) in both folds -> geometric mean rank of 1
    expect_equal(unname(result$vip_rankproducts["f3"]), 1)
})

test_that("plsda_auroc_vip_method builds an nmr_data_analysis_method with the requested parameters", {
    method <- plsda_auroc_vip_method(ncomp = 3, auc_increment_threshold = 0.1)

    expect_s3_class(method, "nmr_data_analysis_method")
    expect_equal(method$train_evaluate_model_params_inner$ncomp, 3)
    expect_false(method$train_evaluate_model_params_inner$return_vip)
    expect_true(method$train_evaluate_model_params_outer$return_vip)
    expect_true(method$train_evaluate_model_params_outer$return_model)
    expect_true(is.function(method$choose_best_inner))
    expect_identical(method$train_evaluate_model_digest_outer, callback_outer_cv_auroc_vip)
})

test_that("plsda_auroc_vip_compare requires named arguments", {
    expect_error(plsda_auroc_vip_compare(list(auc = 1)), "should be named")
})

test_that("plsda_auroc_vip_compare builds a boxplot comparing the auc of each named model", {
    model1 <- list(outer_cv_results_digested = list(auroc = data.frame(auc = c(0.6, 0.7, 0.65))))
    model2 <- list(outer_cv_results_digested = list(auroc = data.frame(auc = c(0.5, 0.55, 0.52))))

    p <- plsda_auroc_vip_compare(model1 = model1, model2 = model2)

    expect_s3_class(p, "ggplot")
    expect_equal(as.character(unique(p$data$Group)), c("model1", "model2"))
})
