#' Random subsampling
#'
#' @param sample_idx Typically a numeric vector with sample index to be separated.
#'     A character vector with sample IDs could also be used
#' @param iterations An integer, the number of iterations in the random subsampling
#' @param test_size A number between 0 and 1. The samples to be included in the
#'    test set on each interation.
#' @param keep_together Either `NULL` or a factor with the same length as `sample_idx`.
#'     `keep_together` can be used to ensure that groups of samples are kept
#'     in together in all iterations (either on training or on test, but never split).
#'     A typical use case for this is when you have sample replicates and you want
#'     to keep all replicates together to prevent overoptimistic results (having
#'     one sample on the train subset and its replicate on the test subset would
#'     make the prediction easier to guess).
#'     Another use case for this is when you have a longitudinal study and you
#'     want to keep some subjects in the same train or test group, because you
#'     want to use some information in a longitudinal way (e.g. a multilevel plsda model).
#' @param balance_in_train Either `NULL` or a factor with the same length as `sample_idx`.
#'     `balance_in_train` can be used to force that on each iteration, the train
#'     partition contains the same number of samples of the given factor levels.
#'     For instance, if we have a dataset with 40 samples of class "A" and 20 samples
#'     of class "B", using a `test_size = 0.25`, we can force to always have 16
#'     samples of class "A" and 16 samples of class "B" in the training subset.
#'     This is beneficial to those algorithms that require that the training groups
#'     are balanced.
#'
#' @return A list of length equal to `iterations`. Each element of the list is
#' a list with two entries (`training` and `test`) containing the `sample_idx`
#' values that will belong to each subset.
#'
#' @examples
#' random_subsampling(1:100, iterations = 4, test_size = 0.25)
#'
#' subject_id <- c("Alice", "Bob", "Charlie", "Eve")
#' random_subsampling(1:4, iterations = 2, test_size = 0.25, keep_together = subject_id)
#'
#' @export
random_subsampling <- function(sample_idx,
    iterations = 10L,
    test_size = 0.25,
    keep_together = NULL,
    balance_in_train = NULL) {
    resample <- function(x, ...) x[sample.int(length(x), ...)]
    num_samples <- length(sample_idx)
    output <- vector("list", iterations)
    if (is.null(keep_together)) {
        keep_together <- sample_idx
    }
    keep_together <- as.factor(keep_together)
    stopifnot(length(keep_together) == num_samples)
    if (!is.null(balance_in_train)) {
        balance_in_train <- as.factor(balance_in_train)
        # Check keep_together & balance consistency
        stopifnot(length(balance_in_train) == num_samples)
        for (lev in levels(keep_together)) {
            samples_in_lev <- which(keep_together == lev)
            if (length(unique(balance_in_train[samples_in_lev])) > 1) {
                stop(glue::glue(
                    "Can't balance samples and keep samples together",
                    " for samples of keep_together == {lev}, because they",
                    " belong to more than one balancing groups: {paste0(unique(balance_in_train[samples_in_lev]))}"
                ))
            }
        }
        keep_balance <- data.frame(
            keep_together = keep_together,
            balance_in_train = balance_in_train,
            stringsAsFactors = TRUE
        )
        keep_balance <- unique(keep_balance)

        smaller_balancing_group <- min(table(keep_balance$balance_in_train))
        num_groups_train <- max(floor(smaller_balancing_group * (1 - test_size)), 1)

        for (i in seq_len(iterations)) {
            groups_train <- keep_balance %>%
                dplyr::group_by(balance_in_train) %>%
                dplyr::sample_n(num_groups_train) %>%
                dplyr::ungroup() %>%
                dplyr::pull(keep_together)
            train_samples <- sample_idx[keep_together %in% groups_train]
            test_samples <- base::setdiff(sample_idx, train_samples)
            output[[i]] <- list(
                training = train_samples,
                test = test_samples
            )
            if (length(train_samples) == 0 || length(test_samples) == 0) {
                cli::cli_abort(
                    message = c(
                        "Too few samples in a random subsampling split",
                        "i" = glue::glue("Train samples: {length(train_samples)}"),
                        "i" = glue::glue("Test samples: {length(test_samples)}"),
                        "i" = sprintf("Please %s the test_subset", ifelse(length(test_samples) == 0, "increase", "decrease"))
                    )
                )
            }
        }
    } else {
        groups_keep_together <- unique(keep_together)
        num_groups_test <- floor(length(groups_keep_together) * test_size)
        for (i in seq_len(iterations)) {
            groups_test <- resample(x = groups_keep_together, size = num_groups_test)
            test_samples <- sample_idx[keep_together %in% groups_test]
            train_samples <- base::setdiff(sample_idx, test_samples)
            output[[i]] <- list(
                training = train_samples,
                test = test_samples
            )
            if (length(train_samples) == 0 || length(test_samples) == 0) {
                cli::cli_abort(
                    message = c(
                        "Too few samples in a random subsampling split",
                        "i" = glue::glue("Train samples: {length(train_samples)}"),
                        "i" = glue::glue("Test samples: {length(test_samples)}"),
                        "i" = sprintf("Please %s the test_subset", ifelse(length(test_samples) == 0, "increase", "decrease"))
                    )
                )
            }
        }
    }
    return(output)
}

#' Split samples for double cross-validation
#'
#' @param dataset An [nmr_dataset_family] like object
#' @inheritParams random_subsampling
#' @inheritParams nmr_data_analysis
#'
#' @return A list with two elements: `outer` and `inner`. Each of those elements
#' is a list as long as the number of iterations given in `external_val` and `internal_val`
#' and it includes the indices of the samples that will belong to the train and test
#' subsets.
#' @noRd
split_double_cv <- function(dataset, keep_together = NULL,
    external_val = list(iterations = 16L, test_size = 0.25),
    internal_val = list(iterations = 10L, test_size = 0.25)) {
    if (is.null(keep_together)) {
        keep_together_data <- seq_len(dataset$num_samples)
    } else {
        keep_together_data <- as.factor(nmr_meta_get_column(dataset, keep_together))
    }
    sample_idx <- seq_len(dataset$num_samples)

    outer_cv_loop <- random_subsampling(
        sample_idx = sample_idx,
        iterations = external_val$iterations,
        test_size = external_val$test_size,
        keep_together = keep_together_data
    )
    outer_cv_iterations <- purrr::imap(outer_cv_loop, function(outer, outer_idx) {
        list(
            outer_iter = outer_idx,
            outer_train = outer$training,
            outer_test = outer$test
        )
    })
    names(outer_cv_iterations) <- as.character(purrr::map_int(outer_cv_iterations, "outer_iter"))

    inner_cv_iterations <- purrr::flatten(
        purrr::imap(outer_cv_loop, function(outer, outer_idx) {
            calib_idx <- outer$training
            blind_idx <- outer$test
            inner_cv <- random_subsampling(
                sample_idx = calib_idx,
                iterations = internal_val$iterations,
                test_size = internal_val$test_size,
                keep_together = keep_together_data[calib_idx]
            )
            purrr::imap(inner_cv, function(inner_iter, inner_iter_idx) {
                list(
                    outer_iter = outer_idx,
                    inner_iter = inner_iter_idx,
                    inner_train_idx = inner_iter$training,
                    inner_test_idx = inner_iter$test
                )
            })
        })
    )

    names(inner_cv_iterations) <-
        paste0(
            as.character(purrr::map_int(inner_cv_iterations, "outer_iter")),
            "_",
            as.character(purrr::map_int(inner_cv_iterations, "inner_iter"))
        )

    list(
        outer = outer_cv_iterations,
        inner = inner_cv_iterations
    )
}




#' Split dataset in matrices, and call a function
#'
#' @param train_test_subset A list with two elements. The first element has the
#' training index and the second the test indices
#' @inheritParams nmr_data_analysis
#' @inheritParams nmr_data_analysis_method
#'
#' @param ... Additional arguments passed to `train_evaluate_model`.
#'
#'
#' @return Whatever `train_evaluate_model` returns
#' @noRd
split_build_perform <- function(train_test_subset,
    dataset,
    y_column,
    identity_column = NULL,
    train_evaluate_model,
    ...) {
    x_all <- nmr_data(dataset)
    y_all <- nmr_meta_get_column(dataset, column = y_column)
    if (is.null(y_all)) {
        cli::cli_abort(
            c(
                "y_column not found",
                "x" = sprintf('Column "%s" does not exist in the dataset', y_column)
            )
        )
    }
    if (!is.null(identity_column)) {
        identity_all <- nmr_meta_get_column(dataset, column = identity_column)
    } else {
        identity_all <- NULL
    }

    # Split
    train_idx <- train_test_subset[[1]]
    test_idx <- train_test_subset[[2]]

    x_train <- x_all[train_idx, , drop = FALSE]
    y_train <- y_all[train_idx]
    if (!is.null(identity_all)) {
        identity_train <- identity_all[train_idx]
    } else {
        identity_train <- NULL
    }

    x_test <- x_all[test_idx, , drop = FALSE]
    y_test <- y_all[test_idx]
    if (!is.null(identity_all)) {
        identity_test <- identity_all[test_idx]
    } else {
        identity_test <- NULL
    }

    # Build & performance:
    train_evaluate_model(
        x_train = x_train, y_train = y_train, identity_train = identity_train,
        x_test = x_test, y_test = y_test, identity_test = identity_test, ...
    )
}



#' Do cross-validation
#' @param train_test_subsets A list of length the number of cross-validation iterations.
#'     Each list item should have two elements, first the train indices and second the test indices
#'
#' @inheritParams nmr_data_analysis
#' @inheritParams new_nmr_data_analysis_method
#' @param train_evaluate_model_args_iter A named list. The names are valid `train_evaluate_model` arguments. Each list element is a vector.
#'
#' @return A list of length equal to `train_test_subsets` with outputs from their corresponding `train_evaluate_model` outputs
#' @noRd
do_cv <- function(dataset, y_column, identity_column, train_evaluate_model,
    train_test_subsets, train_evaluate_model_args_iter = NULL, ..., .enable_parallel = TRUE) {
    if (.enable_parallel) {
        fun <- BiocParallel::bpmapply
    } else {
        fun <- mapply
    }
    output <- do.call(
        what = fun,
        args = c(
            list(
                FUN = split_build_perform,
                train_test_subset = train_test_subsets
            ),
            train_evaluate_model_args_iter,
            list(
                MoreArgs = list(
                    dataset = dataset,
                    y_column = y_column,
                    identity_column = identity_column,
                    train_evaluate_model = train_evaluate_model,
                    ...
                ),
                SIMPLIFY = FALSE
            )
        )
    )
    names(output) <- names(train_test_subsets)
    output
}
