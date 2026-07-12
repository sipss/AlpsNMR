## Prepare demo dataset
prepare_dataset <- function() {
    # 12 artificial samples created based on the 3 demo samples
    MeOH_plasma_extraction_dir <- system.file("dataset-demo", package = "AlpsNMR")
    MeOH_plasma_extraction_xlsx <- file.path(MeOH_plasma_extraction_dir, "dummy_metadata.xlsx")
    exp_subj_id <- readxl::read_excel(MeOH_plasma_extraction_xlsx, sheet = 1)

    zip_files <- fs::dir_ls(MeOH_plasma_extraction_dir, glob = "*.zip")

    dataset <- nmr_read_samples(sample_names = zip_files)
    dataset <- nmr_meta_add(dataset, metadata = exp_subj_id, by = "NMRExperiment")
    dataset <- nmr_interpolate_1D(dataset, axis = c(min = 3.7, max = 4.5, by = 2.3E-4))
    dataset <- nmr_baseline_removal(dataset, lambda = 6, p = 0.01)
    dataset <- nmr_normalize(dataset, method = "area")

    metadata <- nmr_meta_get(dataset, groups = "external")
    metadata$Group <- c("A", "B", "B")
    # Artificially create a larger dataset
    larger_metadata <- rbind(metadata, metadata, metadata, metadata, metadata)

    larger_metadata$NMRExperiment <- as.character(
        seq(from = 10, by = 10, length.out = nrow(larger_metadata))
    )
    data_matrix <- nmr_data(dataset)
    dataset <- new_nmr_dataset_1D(
        ppm_axis = dataset$axis,
        data_1r = rbind(data_matrix, data_matrix, data_matrix, data_matrix, data_matrix),
        metadata = list(external = larger_metadata)
    )
    dataset
}

## Dataset can be used

test_that("nmr_data_analysis works", {
    dataset <- prepare_dataset()
    methodology <- plsda_auroc_vip_method(ncomp = 2)
    set.seed(123L)
    out <- nmr_data_analysis(
        dataset,
        y_column = "Group",
        identity_column = NULL,
        external_val = list(iterations = 1, test_size = 0.25),
        internal_val = list(iterations = 2, test_size = 0.25),
        data_analysis_method = methodology
    )
    expect_false(is.null(out))
})

test_that("random subsampling works", {
    subject_id <- rep(c("Alice", "Bob", "Charlie", "Diana"), times = 2)
    replicate <- rep(c(1, 2), each = 4)
    set.seed(2563432L)
    sample_idx <- 1:8
    num_iterations <- 2L
    out <- random_subsampling(sample_idx,
        iterations = num_iterations, test_size = 0.25,
        keep_together = subject_id
    )
    expect_equal(length(out), num_iterations)
    expect_equal(length(out[[1]][["training"]]), 6L)
    expect_equal(length(out[[1]][["test"]]), 2L)
    # Subjects kept together in the split, no subject in train is present in test:
    expect_equal(
        length(
            intersect(
                subject_id[out[[1]][["test"]]],
                subject_id[out[[1]][["training"]]]
            )
        ),
        0L
    )
})

test_that("split_double_cv works", {
    nsamples <- 16L
    subject_id <- rep(c("Alice", "Bob", "Charlie", "Diana"), times = 4)
    replicate <- rep(c(1, 2), each = 8)
    metadata <- data.frame(
        NMRExperiment = as.character(seq(from = 10, by = 10, length.out = nsamples)),
        SubjectID = subject_id,
        Replicate = replicate
    )
    dataset <- new_nmr_dataset_1D(
        ppm_axis = 1:10,
        data_1r = matrix(sample(1:200, 10 * nsamples), ncol = 10, nrow = nsamples),
        metadata = list(external = metadata)
    )

    external_val_niter <- 2L
    internal_val_niter <- 4L
    external_test_size <- 0.25
    internal_test_size <- 0.34
    out <- split_double_cv(
        dataset = dataset,
        keep_together = "SubjectID",
        external_val = list(iterations = external_val_niter, test_size = external_test_size),
        internal_val = list(iterations = internal_val_niter, test_size = internal_test_size)
    )

    expect_equal(names(out), c("outer", "inner"))
    expect_equal(length(out[["outer"]]), external_val_niter)
    expect_equal(length(out[["inner"]]), external_val_niter * internal_val_niter)
    expected_samples_in_external_test <- floor(nsamples * external_test_size)
    expected_samples_in_train <- nsamples - expected_samples_in_external_test
    expected_samples_in_train_internal_test <- floor(expected_samples_in_train * internal_test_size)
    expected_samples_in_train_internal_train <- expected_samples_in_train - expected_samples_in_train_internal_test

    expect_equal(
        length(out$inner$`1_1`$inner_train_idx),
        expected_samples_in_train_internal_train
    )
})

## bp_VIP_analysis --------------------------------------------------------
##
## These tests register BiocParallel::SerialParam() for the duration of the
## call to bp_VIP_analysis() so that the sequence of sample() calls made
## inside its (parallel-capable) bplapply() loop is fully deterministic given
## a set.seed() call, and reproducible across machines/CI (the default
## MulticoreParam backend forks workers and its RNG behaviour is not a
## reliable basis for a reproducible unit test).

test_that("bp_VIP_analysis identifies a strongly predictive feature as relevant", {
    # bp_VIP_analysis() permutes each feature j *in place* to build a
    # null-importance baseline for that feature. We call bp_VIP_analysis()
    # directly (not a reimplementation of the permutation step) on a tiny
    # synthetic dataset with exactly one strongly predictive feature (V1)
    # and several pure-noise features. If permutation correctly shuffles
    # each feature's own column, V1 should reliably surface as the single
    # "relevant" VIP with a markedly higher importance score than the noise
    # features.
    skip_if_not_installed("mixOmics")
    skip_if_not_installed("BiocParallel")

    old_bpparam <- BiocParallel::bpparam()
    BiocParallel::register(BiocParallel::SerialParam())
    on.exit(BiocParallel::register(old_bpparam), add = TRUE)

    set.seed(1)
    n <- 30
    p <- 4
    y <- factor(rep(c("A", "B"), times = n / 2))
    x <- matrix(rnorm(n * p, sd = 1), nrow = n, ncol = p)
    colnames(x) <- paste0("V", seq_len(p))
    # V1 is made strongly predictive of the class; V2-V4 stay pure noise
    x[y == "A", 1] <- x[y == "A", 1] + 15

    metadata <- data.frame(
        NMRExperiment = as.character(seq_len(n)),
        Condition = y
    )
    dataset <- new_nmr_dataset_peak_table(
        peak_table = x,
        metadata = list(external = metadata)
    )

    train_index <- 1:20 # both classes present in train (1:20) and test (21:30)

    # bp_VIP_analysis() emits an informational cli_warn() when few VIPs
    # clear the "important" (stricter) threshold with so few bootstraps;
    # that is expected/benign here, so it is suppressed. nbootstrap = 30 (as
    # opposed to a handful) keeps the per-feature bootstrap standard
    # deviation used to standardize the VIP-score difference (see the paper's
    # Eq. (13)/step 6) from being dominated by noise on this tiny dataset.
    result <- suppressWarnings(bp_VIP_analysis(
        dataset,
        train_index,
        y_column = "Condition",
        ncomp = 1,
        nbootstrap = 30 # kept small for test speed
    ))

    # Expected shape: num_features x nbootstrap matrices
    expect_equal(dim(result$pls_vip), c(p, 30))
    expect_equal(dim(result$pls_vip_perm), c(p, 30))
    expect_setequal(rownames(result$pls_vip_means), colnames(x))

    # The informative feature is (the only feature) flagged as relevant:
    expect_equal(result$relevant_vips, "V1")

    # And its bootstrapped VIP-difference mean is clearly the largest,
    # well above every noise feature:
    means <- result$pls_vip_means[, 1]
    expect_gt(means["V1"], max(means[setdiff(names(means), "V1")]))
})

test_that("bp_VIP_analysis uses the ncomp-th (cumulative) VIP column, not a re-aggregation across components", {
    # mixOmics::vip() returns, in its column h, the VIP already computed
    # cumulatively over components 1..h (Eq. 9 of Afanador, Tran & Buydens,
    # 2013). bp_VIP_analysis() must therefore take column `ncomp` as-is; it
    # must not re-aggregate the per-component columns (e.g. via
    # sqrt(rowSums(x^2) / ncomp)), which would apply Eq. 9's normalization a
    # second time to already-normalized quantities.
    #
    # plsda_vip() is mocked to return an (almost) fixed, easily distinguished
    # matrix regardless of the fitted model (a small jitter is added so the
    # across-bootstrap variance used later isn't degenerately zero), so the
    # two candidate formulas produce numerically distinct, checkable
    # results: taking column `ncomp` (2) gives approximately c(10, 20, 30,
    # 40), while sqrt(rowSums(x^2) / 2) would instead give approximately
    # sqrt(c(101, 404, 909, 1616) / 2) ~= c(7.1, 14.2, 21.3, 28.4).
    skip_if_not_installed("mixOmics")
    skip_if_not_installed("BiocParallel")

    old_bpparam <- BiocParallel::bpparam()
    BiocParallel::register(BiocParallel::SerialParam())
    on.exit(BiocParallel::register(old_bpparam), add = TRUE)

    p <- 4
    base_vip <- matrix(
        c(1, 2, 3, 4, 10, 20, 30, 40),
        nrow = p, ncol = 2,
        dimnames = list(paste0("V", seq_len(p)), NULL)
    )
    # A little call-to-call jitter keeps each bootstrap replicate's
    # VIP-difference distribution non-degenerate (bp_VIP_analysis() divides
    # by its across-replicate standard deviation), without disturbing which
    # of the two candidate formulas the assertion below distinguishes.
    testthat::local_mocked_bindings(
        plsda_vip = function(plsda_model) base_vip + stats::rnorm(length(base_vip), sd = 0.01)
    )

    set.seed(1)
    n <- 20
    y <- factor(rep(c("A", "B"), each = n / 2))
    x <- matrix(rnorm(n * p), nrow = n, ncol = p)
    colnames(x) <- paste0("V", seq_len(p))

    metadata <- data.frame(
        NMRExperiment = as.character(seq_len(n)),
        Condition = y
    )
    dataset <- new_nmr_dataset_peak_table(
        peak_table = x,
        metadata = list(external = metadata)
    )
    train_index <- c(1:7, 11:17) # both classes in train (1:7, 11:17) and test (8:10, 18:20)

    result <- suppressWarnings(bp_VIP_analysis(
        dataset,
        train_index,
        y_column = "Condition",
        ncomp = 2,
        nbootstrap = 3
    ))

    expected <- base_vip[, 2]
    for (feature in names(expected)) {
        expect_equal(
            unname(result$pls_vip[feature, ]),
            rep(expected[[feature]], 3),
            tolerance = 0.1
        )
    }
})

test_that("bp_VIP_analysis's permutation score for feature j is feature j's own VIP in the j-permuted model", {
    # Per the paper's steps (4)-(5) (Section 2.4), the permutation baseline
    # used for feature j's importance must be feature j's own VIP value in
    # the model fit with (only) feature j permuted -- the diagonal entry
    # pls_vip_perm[j, j] -- not an average of feature j's VIP across the p
    # separate permuted-feature models.
    #
    # plsda_vip() is mocked so its return value depends only on the position
    # of the call within a bootstrap iteration: the first call (the
    # un-permuted model) returns a fixed vector, and the k-th call inside the
    # permutation loop (the model with feature k permuted) returns
    # 10*k + m for feature m, so the diagonal (m == k) is numerically
    # distinguishable from a column mean. A small jitter keeps the
    # across-bootstrap variance non-degenerate.
    skip_if_not_installed("mixOmics")
    skip_if_not_installed("BiocParallel")

    old_bpparam <- BiocParallel::bpparam()
    BiocParallel::register(BiocParallel::SerialParam())
    on.exit(BiocParallel::register(old_bpparam), add = TRUE)

    p <- 4
    call_count <- 0L
    testthat::local_mocked_bindings(
        plsda_vip = function(plsda_model) {
            call_count <<- call_count + 1L
            pos <- (call_count - 1L) %% (p + 1L)
            base_vip <- if (pos == 0L) {
                c(100, 200, 300, 400) # the un-permuted bootstrap model
            } else {
                10 * pos + seq_len(p) # model with feature `pos` permuted
            }
            vip <- base_vip + stats::rnorm(p, sd = 0.01)
            matrix(vip, nrow = p, ncol = 1, dimnames = list(paste0("V", seq_len(p)), NULL))
        }
    )

    set.seed(1)
    n <- 20
    y <- factor(rep(c("A", "B"), each = n / 2))
    x <- matrix(rnorm(n * p), nrow = n, ncol = p)
    colnames(x) <- paste0("V", seq_len(p))
    metadata <- data.frame(NMRExperiment = as.character(seq_len(n)), Condition = y)
    dataset <- new_nmr_dataset_peak_table(peak_table = x, metadata = list(external = metadata))
    train_index <- c(1:7, 11:17) # both classes in train (1:7, 11:17) and test (8:10, 18:20)

    result <- suppressWarnings(bp_VIP_analysis(
        dataset, train_index,
        y_column = "Condition", ncomp = 1, nbootstrap = 3
    ))

    expected_perm_score <- c(V1 = 11, V2 = 22, V3 = 33, V4 = 44)
    for (feature in names(expected_perm_score)) {
        expect_equal(
            unname(result$pls_vip_perm[feature, ]),
            rep(expected_perm_score[[feature]], 3),
            tolerance = 0.1
        )
    }
})

test_that("bp_VIP_analysis's significance threshold uses df = n_samples - 1, not df = nbootstrap - 1", {
    # Afanador, Tran & Buydens (2013) calculate the 95% confidence interval
    # for the (per-feature) standardized VIP-score difference by multiplying
    # its bootstrap standard deviation by "the appropriate quantile,
    # t_{1-alpha/2,n-1}" (text following Eqs. (12)-(13); the same quantile is
    # reused verbatim for the "Important" cut-off in the guidelines list of
    # Section 2.3). Throughout the paper, n denotes the number of *training
    # samples* the model was fit on (e.g. "VACCINE data set with n = 50" in
    # Section 3.6) -- not B, the number of bootstrap datasets (nbootstrap in
    # this code; B = 300 for the paper's own experiments). This is confirmed
    # numerically by the paper's own worked example: for the VACCINE
    # dataset (n = 50 training samples, B = 300 bootstraps), the text reports
    # the cut-off as "2.01 (t_{1-alpha/2,n-1})" -- which is qt(0.975, df = 49)
    # = 2.0096, not qt(0.975, df = 299) = 1.9679.
    #
    # Because `element <- pls_vip_score_diff[k, ] / sd(pls_vip_score_diff[k, ])`
    # standardizes each feature's difference vector to its own sample
    # standard deviation, boots_vip_sd[k] (recomputed from `element` with the
    # same n-1 estimator) is analytically always 1, so
    # result$error[k] == qt(0.975, df = <whatever df is used>) exactly,
    # independent of the (random) bootstrap data. This lets the two
    # candidate degrees of freedom be told apart without any mocking.
    skip_if_not_installed("mixOmics")
    skip_if_not_installed("BiocParallel")

    old_bpparam <- BiocParallel::bpparam()
    BiocParallel::register(BiocParallel::SerialParam())
    on.exit(BiocParallel::register(old_bpparam), add = TRUE)

    set.seed(1)
    n <- 20
    p <- 4
    y <- factor(rep(c("A", "B"), each = n / 2))
    x <- matrix(rnorm(n * p), nrow = n, ncol = p)
    colnames(x) <- paste0("V", seq_len(p))
    x[y == "A", 1] <- x[y == "A", 1] + 8

    metadata <- data.frame(
        NMRExperiment = as.character(seq_len(n)),
        Condition = y
    )
    dataset <- new_nmr_dataset_peak_table(
        peak_table = x,
        metadata = list(external = metadata)
    )
    train_index <- c(1:5, 11:15) # n_samples = 10, both classes in train and test
    nbootstrap <- 20 # deliberately != n_samples, so the two df choices disagree

    result <- suppressWarnings(bp_VIP_analysis(
        dataset,
        train_index,
        y_column = "Condition",
        ncomp = 1,
        nbootstrap = nbootstrap
    ))

    n_samples <- length(train_index)
    expected_error <- qt(0.975, df = n_samples - 1)
    wrong_error <- qt(0.975, df = nbootstrap - 1)
    expect_false(isTRUE(all.equal(expected_error, wrong_error)))
    expect_equal(unname(result$error[, 1]), rep(expected_error, p))
})

test_that("bp_VIP_analysis shuffles each feature's own column when building its permutation baseline", {
    # Each feature's own values must be shuffled independently, not swapped
    # for another feature's values, or the permutation-importance baseline
    # would not reflect that feature's own association with the outcome.
    body_txt <- paste(deparse(body(bp_VIP_analysis)), collapse = " ")

    expect_true(
        grepl(
            "x_train_boots_perm\\[, j\\] *<- *sample\\(x_train_boots\\[, *j\\]\\)",
            body_txt
        ),
        info = "Expected column j to be shuffled via sample(x_train_boots[, j])"
    )
    expect_false(
        grepl("random_pos", body_txt),
        info = "Column j must not be replaced by another (random) column's values"
    )
})

test_that("bp_VIP_analysis recovers from a degenerate single-class bootstrap resample", {
    # Inside the bootstrap loop of bp_VIP_analysis(), if a bootstrap
    # resample happens to contain only one class, the code redraws it (a
    # fresh sample-with-replacement) until it contains more than one class,
    # rather than patching a single element of the offending resample. To
    # exercise this code path with the real function (not a reimplementation
    # of it), we build a deliberately imbalanced train set (1 sample of
    # class A, 4 of class B). Bootstrap resampling with replacement from 5
    # elements, only one of which is class A, has a per-iteration
    # probability of ~(4/5)^5 = 32.8% of missing the lone A sample entirely
    # (a degenerate, single-class resample). With nbootstrap = 10 draws, the
    # probability of the degenerate branch firing at least once is
    # 1 - 0.672^10 ~= 98%. We pin a seed (verified across 15 candidate seeds
    # to all succeed) and assert the call completes without error.
    skip_if_not_installed("mixOmics")
    skip_if_not_installed("BiocParallel")

    old_bpparam <- BiocParallel::bpparam()
    BiocParallel::register(BiocParallel::SerialParam())
    on.exit(BiocParallel::register(old_bpparam), add = TRUE)

    set.seed(1)
    n <- 20
    p <- 4
    y <- factor(rep(c("A", "B"), each = n / 2))
    x <- matrix(rnorm(n * p), nrow = n, ncol = p)
    colnames(x) <- paste0("V", seq_len(p))
    x[y == "A", 1] <- x[y == "A", 1] + 8

    metadata <- data.frame(
        NMRExperiment = as.character(seq_len(n)),
        Condition = y
    )
    dataset <- new_nmr_dataset_peak_table(
        peak_table = x,
        metadata = list(external = metadata)
    )

    # 1 sample of class A (row 1) + 4 samples of class B (rows 11-14):
    train_index <- c(1L, 11L, 12L, 13L, 14L)

    result <- NULL
    expect_no_error(
        result <- suppressWarnings(bp_VIP_analysis(
            dataset,
            train_index,
            y_column = "Condition",
            ncomp = 1,
            nbootstrap = 10
        ))
    )
    expect_true(is.list(result))
    expect_true(all(c("pls_vip", "relevant_vips") %in% names(result)))
})

test_that("bp_VIP_analysis redraws a degenerate bootstrap resample instead of patching a single element", {
    # Patching a single element of a degenerate resample (e.g. always
    # overwriting position 1 with the first alternate-class sample found in
    # y_train) would bias that position towards a fixed, non-random value
    # instead of leaving the resample an unbiased draw with replacement.
    # Guard against that pattern reappearing.
    body_txt <- paste(deparse(body(bp_VIP_analysis)), collapse = " ")

    expect_false(
        grepl("x_train_boots\\[1, *\\] *<-", body_txt),
        info = "A degenerate resample must be redrawn, not fixed up by patching a single row"
    )
})

test_that("plot_vip_scores draws the importance threshold at df = n_samples - 1", {
    # Afanador, Tran & Buydens (2013) place the "Important" cut-off at the
    # quantile t_{1-alpha/2,n-1}, where n is the number of training samples the
    # models were fit on -- not the number of bootstrap datasets. The threshold
    # line must therefore use df = n_samples - 1.
    skip_if_not_installed("ggplot2")

    vip_means <- c(V1 = 3, V2 = 1, V3 = 0.2)
    error <- 2
    n_samples <- 12L

    p <- plot_vip_scores(vip_means, error, n_samples = n_samples, plot = FALSE)

    # Extract the y-intercept of the horizontal threshold line from the ggplot
    # (geom_hline stores its yintercept in that layer's data frame):
    hline_y <- NULL
    for (ly in p$layers) {
        d <- ly$data
        if (is.data.frame(d) && "yintercept" %in% names(d)) {
            hline_y <- d$yintercept
        }
    }
    expect_false(is.null(hline_y))
    expect_equal(hline_y, qt(0.975, df = n_samples - 1))
    # The old (wrong) API took an `nbootstrap` argument used for this df:
    expect_false("nbootstrap" %in% names(formals(plot_vip_scores)))
    expect_true("n_samples" %in% names(formals(plot_vip_scores)))
})

test_that("bp_kfold_VIP_analysis assigns folds at random, not by a fixed modulo split", {
    # The previous implementation shuffled a local copy of the data but then
    # built the folds with split(x, x %% k) over the *unshuffled* index vector,
    # so the partition was a deterministic modulo split and the shuffle was dead
    # code. Folds must instead be a genuine random partition of the samples.
    body_txt <- paste(deparse(body(bp_kfold_VIP_analysis)), collapse = " ")

    expect_false(
        grepl("%%", body_txt, fixed = TRUE),
        info = "Folds must be a random partition, not a deterministic x %% k modulo split"
    )
    expect_true(
        grepl("rep_len", body_txt, fixed = TRUE),
        info = "Random fold assignment is expected to be built from sample(rep_len(...))"
    )
})
