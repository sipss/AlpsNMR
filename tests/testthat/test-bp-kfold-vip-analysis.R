## bp_kfold_VIP_analysis ----------------------------------------------------
##
## These tests register BiocParallel::SerialParam() for the duration of the
## call, for the same reason as the bp_VIP_analysis() tests in
## test-nmr-data-analysis.R: bp_kfold_VIP_analysis() calls the
## (parallel-capable) bp_VIP_analysis() once per fold, and its bplapply()-based
## sample() calls need a deterministic, non-forking backend to be reproducible.

build_kfold_test_dataset <- function(seed = 1, n = 30, p = 4) {
    set.seed(seed)
    y <- factor(rep(c("A", "B"), times = n / 2))
    x <- matrix(rnorm(n * p), nrow = n, ncol = p)
    colnames(x) <- paste0("V", seq_len(p))
    # V1 is made strongly predictive of the class; the rest stay pure noise
    x[y == "A", 1] <- x[y == "A", 1] + 15

    metadata <- data.frame(NMRExperiment = as.character(seq_len(n)), Condition = y)
    new_nmr_dataset_peak_table(peak_table = x, metadata = list(external = metadata))
}

test_that("bp_kfold_VIP_analysis combines results across folds and ranks the predictive feature first", {
    skip_if_not_installed("mixOmics")
    skip_if_not_installed("BiocParallel")

    old_bpparam <- BiocParallel::bpparam()
    BiocParallel::register(BiocParallel::SerialParam())
    on.exit(BiocParallel::register(old_bpparam), add = TRUE)

    dataset <- build_kfold_test_dataset()

    set.seed(2)
    k <- 2
    result <- suppressWarnings(bp_kfold_VIP_analysis(
        dataset,
        y_column = "Condition",
        k = k,
        ncomp = 1,
        nbootstrap = 5
    ))

    expect_setequal(
        names(result),
        c("important_vips", "relevant_vips", "wilcoxon_vips", "vip_means", "vip_score_plot", "kfold_results", "kfold_index")
    )
    expect_length(result$kfold_results, k)
    expect_length(result$kfold_index, k)
    expect_s3_class(result$vip_score_plot, "ggplot")
    expect_setequal(names(result$vip_means), c("V1", "V2", "V3", "V4"))
    # The predictive feature has the highest mean VIP-difference score:
    expect_equal(names(result$vip_means)[1], "V1")
    expect_true(result$vip_means[["V1"]] > max(result$vip_means[setdiff(names(result$vip_means), "V1")]))
})

test_that("bp_kfold_VIP_analysis requires k > 1", {
    dataset <- build_kfold_test_dataset()
    expect_error(
        bp_kfold_VIP_analysis(dataset, y_column = "Condition", k = 1),
        "K must be integer greater than 1"
    )
})

test_that("bp_kfold_VIP_analysis requires at least two classes", {
    dataset <- build_kfold_test_dataset()
    dataset$metadata$external$Condition <- factor(rep("A", nrow(dataset$metadata$external)))
    expect_error(
        bp_kfold_VIP_analysis(dataset, y_column = "Condition", k = 2),
        "Only one class in data set"
    )
})
