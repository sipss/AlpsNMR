## random_subsampling: balance_in_train branch ---------------------------

test_that("random_subsampling with balance_in_train keeps the training set balanced", {
    set.seed(42)
    sample_idx <- 1:12
    group <- rep(c("A", "B"), each = 6)

    out <- random_subsampling(sample_idx, iterations = 3, test_size = 0.25, balance_in_train = group)

    expect_length(out, 3)
    for (iter in out) {
        expect_equal(as.integer(table(group[iter$training])), c(4L, 4L))
        expect_equal(length(intersect(iter$training, iter$test)), 0L)
        expect_equal(sort(c(iter$training, iter$test)), sample_idx)
    }
})

test_that("random_subsampling errors when keep_together groups span more than one balance_in_train level", {
    expect_error(
        random_subsampling(
            1:4,
            iterations = 1, test_size = 0.25,
            keep_together = c("g1", "g1", "g2", "g2"),
            balance_in_train = c("A", "B", "A", "A")
        ),
        "more than one balancing group"
    )
})

test_that("random_subsampling with balance_in_train aborts when a split leaves an empty train or test set", {
    expect_error(
        random_subsampling(1:4, iterations = 1, test_size = 0, balance_in_train = c("A", "A", "B", "B")),
        "Too few samples"
    )
})

test_that("random_subsampling without balance_in_train aborts when a split leaves an empty train or test set", {
    expect_error(
        random_subsampling(1:4, iterations = 1, test_size = 0),
        "Too few samples"
    )
    expect_error(
        random_subsampling(1:4, iterations = 1, test_size = 1),
        "Too few samples"
    )
})

## split_build_perform ----------------------------------------------------

build_split_test_dataset <- function() {
    new_nmr_dataset_1D(
        ppm_axis = 1:5,
        data_1r = matrix(1:20, nrow = 4, ncol = 5),
        metadata = list(external = data.frame(
            NMRExperiment = as.character(1:4),
            Group = c("A", "A", "B", "B"),
            Subj = c("s1", "s1", "s2", "s2")
        ))
    )
}

test_that("split_build_perform splits x/y/identity into train/test and forwards extra args", {
    dataset <- build_split_test_dataset()
    captured <- NULL
    fake_train_evaluate <- function(x_train, y_train, identity_train, x_test, y_test, identity_test, ...) {
        captured <<- list(
            x_train = x_train, y_train = y_train, identity_train = identity_train,
            x_test = x_test, y_test = y_test, identity_test = identity_test,
            dots = list(...)
        )
        "result"
    }

    result <- split_build_perform(
        train_test_subset = list(c(1, 2), c(3, 4)),
        dataset = dataset,
        y_column = "Group",
        identity_column = "Subj",
        train_evaluate_model = fake_train_evaluate,
        extra_arg = 99
    )

    expect_equal(result, "result")
    expect_equal(nrow(captured$x_train), 2)
    expect_equal(captured$y_train, c("A", "A"))
    expect_equal(captured$identity_train, c("s1", "s1"))
    expect_equal(nrow(captured$x_test), 2)
    expect_equal(captured$y_test, c("B", "B"))
    expect_equal(captured$identity_test, c("s2", "s2"))
    expect_equal(captured$dots, list(extra_arg = 99))
})

test_that("split_build_perform leaves identity_train/identity_test NULL when identity_column is NULL", {
    dataset <- build_split_test_dataset()
    captured <- NULL
    fake_train_evaluate <- function(x_train, y_train, identity_train, x_test, y_test, identity_test, ...) {
        captured <<- list(identity_train = identity_train, identity_test = identity_test)
        NULL
    }

    split_build_perform(
        train_test_subset = list(c(1, 2), c(3, 4)),
        dataset = dataset,
        y_column = "Group",
        identity_column = NULL,
        train_evaluate_model = fake_train_evaluate
    )

    expect_null(captured$identity_train)
    expect_null(captured$identity_test)
})

test_that("split_build_perform errors when y_column does not exist in the dataset", {
    dataset <- build_split_test_dataset()
    expect_error(
        split_build_perform(
            train_test_subset = list(c(1, 2), c(3, 4)),
            dataset = dataset,
            y_column = "NoSuchColumn",
            identity_column = NULL,
            train_evaluate_model = function(...) NULL
        )
    )
})

## do_cv --------------------------------------------------------------------

test_that("do_cv runs train_evaluate_model once per train_test_subsets entry, iterating extra args and preserving names", {
    dataset <- new_nmr_dataset_1D(
        ppm_axis = 1:5,
        data_1r = matrix(1:30, nrow = 6, ncol = 5),
        metadata = list(external = data.frame(NMRExperiment = as.character(1:6), Group = rep(c("A", "B"), 3)))
    )
    train_test_subsets <- list(
        iter1 = list(c(1, 2, 3), c(4, 5, 6)),
        iter2 = list(c(4, 5, 6), c(1, 2, 3))
    )
    fake_train_evaluate <- function(x_train, y_train, identity_train, x_test, y_test, identity_test, ncomp) {
        list(ncomp = ncomp, n_train = nrow(x_train), n_test = nrow(x_test))
    }

    out <- do_cv(
        dataset = dataset,
        y_column = "Group",
        identity_column = NULL,
        train_evaluate_model = fake_train_evaluate,
        train_test_subsets = train_test_subsets,
        train_evaluate_model_args_iter = list(ncomp = c(1, 2)),
        .enable_parallel = FALSE
    )

    expect_equal(names(out), c("iter1", "iter2"))
    expect_equal(out$iter1$ncomp, 1)
    expect_equal(out$iter2$ncomp, 2)
    expect_equal(out$iter1$n_train, 3)
    expect_equal(out$iter1$n_test, 3)
})
