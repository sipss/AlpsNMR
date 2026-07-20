#' Data analysis
#'
#' Data analysis on AlpsNMR can be performed on both [nmr_dataset_1D] full spectra
#' as well as [nmr_dataset_peak_table] peak tables.
#'
#' The workflow consists of a double cross validation strategy using random
#' subsampling for splitting into train and test sets. The classification model
#' and the metric to choose the best model can be customized (see
#' [new_nmr_data_analysis_method()]), but for now only a PLSDA classification
#' model with a best area under ROC curve metric is implemented (see
#' the examples here and [plsda_auroc_vip_method])
#'
#' @param dataset An [nmr_dataset_family] object
#' @param y_column A string with the name of the y column (present in the
#'    metadata of the dataset)
#' @param identity_column `NULL` or a string with the name of the identity column (present in the
#'    metadata of the dataset).
#' @param external_val,internal_val A list with two elements: `iterations` and `test_size`.
#'    See [random_subsampling] for further details
#' @param data_analysis_method An [nmr_data_analysis_method] object
#' @param .enable_parallel Set to `FALSE` to disable parallellization.
#' @return A list with the following elements:
#'
#' - `train_test_partitions`: A list with the indices used in train and test on each of the cross-validation iterations
#' - `inner_cv_results`: The output returned by `train_evaluate_model` on each inner cross-validation
#' - `inner_cv_results_digested`: The output returned by `choose_best_inner`.
#' - `outer_cv_results`: The output returned by `train_evaluate_model` on each outer cross-validation
#' - `outer_cv_results_digested`: The output returned by `train_evaluate_model_digest_outer`.
#' @examples
#' # Data analysis for a table of integrated peaks
#'
#' ## Generate an artificial nmr_dataset_peak_table:
#' ### Generate artificial metadata:
#' num_samples <- 32 # use an even number in this example
#' num_peaks <- 20
#' metadata <- data.frame(
#'     NMRExperiment = as.character(1:num_samples),
#'     Condition = rep(c("A", "B"), times = num_samples / 2)
#' )
#'
#' ### The matrix with peaks
#' peak_means <- runif(n = num_peaks, min = 300, max = 600)
#' peak_sd <- runif(n = num_peaks, min = 30, max = 60)
#' peak_matrix <- mapply(function(mu, sd) rnorm(num_samples, mu, sd),
#'     mu = peak_means, sd = peak_sd
#' )
#' colnames(peak_matrix) <- paste0("Peak", 1:num_peaks)
#'
#' ## Artificial differences depending on the condition:
#' peak_matrix[metadata$Condition == "A", "Peak2"] <-
#'     peak_matrix[metadata$Condition == "A", "Peak2"] + 70
#'
#' peak_matrix[metadata$Condition == "A", "Peak6"] <-
#'     peak_matrix[metadata$Condition == "A", "Peak6"] - 60
#'
#' ### The nmr_dataset_peak_table
#' peak_table <- new_nmr_dataset_peak_table(
#'     peak_table = peak_matrix,
#'     metadata = list(external = metadata)
#' )
#'
#' ## We will use a double cross validation, splitting the samples with random
#' ## subsampling both in the external and internal validation.
#' ## The classification model will be a PLSDA, exploring at maximum 3 latent
#' ## variables.
#' ## The best model will be selected based on the area under the ROC curve
#' methodology <- plsda_auroc_vip_method(ncomp = 3)
#' model <- nmr_data_analysis(
#'     peak_table,
#'     y_column = "Condition",
#'     identity_column = NULL,
#'     external_val = list(iterations = 3, test_size = 0.25),
#'     internal_val = list(iterations = 3, test_size = 0.25),
#'     data_analysis_method = methodology
#' )
#' ## Area under ROC for each outer cross-validation iteration:
#' model$outer_cv_results_digested$auroc
#' ## Rank Product of the Variable Importance in the Projection
#' ## (Lower means more important)
#' sort(model$outer_cv_results_digested$vip_rankproducts)
#'
#' @export
nmr_data_analysis <- function(dataset,
    y_column,
    identity_column,
    external_val,
    internal_val,
    data_analysis_method,
    .enable_parallel = TRUE) {
    train_evaluate_model <- data_analysis_method[["train_evaluate_model"]]
    train_evaluate_model_params_inner <- data_analysis_method[["train_evaluate_model_params_inner"]]
    choose_best_inner <- data_analysis_method[["choose_best_inner"]]
    train_evaluate_model_params_outer <- data_analysis_method[["train_evaluate_model_params_outer"]]
    train_evaluate_model_digest_outer <- data_analysis_method[["train_evaluate_model_digest_outer"]]

    # Prepare double cross-validation splits:
    train_test_blind_subsets <- split_double_cv(dataset,
        keep_together = identity_column,
        external_val = external_val,
        internal_val = internal_val
    )

    # These are ALL the inner cross-validation iterations:
    train_test_subsets_inner <- purrr::map(
        train_test_blind_subsets$inner,
        ~ .[c("inner_train_idx", "inner_test_idx")]
    )
    # We run the train_evaluate_model function for each of the inner CV.
    inner_cv_results <- do.call(
        what = do_cv,
        args = c(
            list(
                dataset = dataset,
                y_column = y_column,
                identity_column = identity_column,
                train_evaluate_model = train_evaluate_model,
                train_test_subsets = train_test_subsets_inner,
                .enable_parallel = .enable_parallel
            ),
            train_evaluate_model_params_inner
        )
    )

    # We choose the best hyper-parameters for each inner cross-validation:
    inner_cv_results_digested <- choose_best_inner(inner_cv_results)

    # Prepare the indices for the final outer cv models:
    train_test_subsets_outer <- purrr::map(
        train_test_blind_subsets$outer,
        ~ .[c("outer_train", "outer_test")]
    )

    # Compute the outer cv models
    outer_cv_results <- do.call(
        what = do_cv,
        args = c(
            list(
                dataset = dataset,
                y_column = y_column,
                identity_column = identity_column,
                train_evaluate_model = train_evaluate_model,
                train_test_subsets = train_test_subsets_outer,
                train_evaluate_model_args_iter = inner_cv_results_digested$train_evaluate_model_args,
                .enable_parallel = .enable_parallel
            ),
            train_evaluate_model_params_outer
        )
    )

    # Digest the results:
    outer_cv_results_dig <- train_evaluate_model_digest_outer(outer_cv_results)

    # Give output:
    list(
        train_test_partitions = train_test_blind_subsets,
        inner_cv_results = inner_cv_results,
        inner_cv_results_digested = inner_cv_results_digested,
        outer_cv_results = outer_cv_results,
        outer_cv_results_digested = outer_cv_results_dig
    )
}


get_test_accuracy <- function(model, x_test, y_test) {
    pred <- stats::predict(model, newdata = x_test, dist = "max.dist")
    y_test_pred <- pred$class$max.dist[,model$ncomp]
    conf_mat <- table(REAL = y_test, PRED = y_test_pred)
    accuracy <- diag(conf_mat)/sum(conf_mat)
    accuracy
}


train_models_with_only_vip_features <- function(x_train, y_train, x_test, y_test, ncomp, important_vips, relevant_vips, identity_train = NULL) {
    # important_vips is a more stringent subset of relevant_vips
    if (length(important_vips) == 0) {
        cli::cli_warn(
            c(
                "No VIPs are ranked as important",
                "i" = "Try increasing the number of bootstrap iterations"
            )
        )
        if (length(relevant_vips) == 0) {
            cli::cli_warn(
                c(
                    "bp_VIP_analysis: no relevant_vips found",
                    "i" = "You may try increasing the number of bootstraps"
                )
            )
        }
        return(list(vips_model = NULL, vips_CR = 0))
    }
    
    if (length(important_vips) < ncomp) {
        cli::cli_warn(
            c(
                "Number of important vips ({length(important_vips)}) smaller than requested ncomp ({ncomp})",
                "i" = "Attempting to use relevant vips (less stringent criteria)",
                "i" = "You may want to consider reducing ncomp"
            )
        )
        if (length(relevant_vips) < ncomp) {
            cli::cli_warn(
                c(
                    "Number of relevant vips ({length(relevant_vips)}) smaller than requested ncomp ({ncomp})",
                    "i" = "Can't compute sub-model using only VIPs",
                    "i" = "You may want to consider reducing ncomp"
                )
            )
            return(list(vips_model = NULL, vips_CR = 0))
        }
        x_train_reduced <- as.matrix(x_train[, relevant_vips, drop = FALSE])
        x_test_reduced <- as.matrix(x_test[, relevant_vips, drop = FALSE])
    } else {
        x_train_reduced <- as.matrix(x_train[, important_vips, drop = FALSE])
        x_test_reduced <- as.matrix(x_test[, important_vips, drop = FALSE])
    }
    
    # Fit PLS model
    vips_model <- plsda_build(
        x = x_train_reduced,
        y = y_train,
        identity = identity_train,
        ncomp = ncomp
    )
    
    # Measure the classification rate (CR) of the fold
    vips_CR <- get_test_accuracy(vips_model, x_test_reduced, y_test)
    list(vips_model = vips_model, vips_CR = vips_CR)
}


#' Create method for NMR data analysis
#'
#' @param train_evaluate_model A function. The `train_evaluate_model` must have the following signature:
#'
#'         function(x_train, y_train, identity_train, x_test, y_test, identity_test, ...)
#'
#' The `x_train` and `y_train` (and their test counterparts) are self-explanatory.
#'
#' The `identity_` arguments are expected to be factors. They can be used for
#' instance with a callback that uses [mixOmics::plsda] in a `multilevel` approach
#' for longitudinal studies. In those studies the `identity` would be an
#' identifier of the subject.
#'
#' The `...` arguments are free to be defined for each `train_evaluate_model`.
#'
#' @param train_evaluate_model_params_inner,train_evaluate_model_params_outer A list with additional
#'     arguments to pass to `train_evaluate_model` either in the inner cv loop or in the outer cv loop.
#'
#' @param choose_best_inner A function with a single argument:
#'
#'     function(inner_cv_results)
#'
#'    The argument is a list of `train_evaluate_model` outputs.
#'    The return value of must be a list with at least an element named `train_evaluate_model_args`.
#'    `train_evaluate_model_args` must be a named list.
#'
#'    - Each element must be named as one of the `train_evaluate_model` arguments.
#'    - Each element must be a vector as long as the number of outer cross-validations.
#'    - The values of each vector must be the values that the `train_evaluate_model`
#'    argument must take on each outer cross-validation iteration
#'    Additional list elements can be returned and will be given back to the user
#'
#' @param train_evaluate_model_digest_outer A function with a single argument:
#'
#'    function(outer_cv_results)
#'
#'    The argument is a list of `train_evaluate_model` outputs in outer cross-validation.
#'    The return value is returned by `nmr_data_analysis`
#'
#' @return An object encapsulating the method dependent functions that can be used with [nmr_data_analysis]
#' @name nmr_data_analysis_method
#' @export
#'
new_nmr_data_analysis_method <- function(train_evaluate_model,
    train_evaluate_model_params_inner,
    choose_best_inner,
    train_evaluate_model_params_outer,
    train_evaluate_model_digest_outer) {
    out <- list(
        train_evaluate_model = train_evaluate_model,
        train_evaluate_model_params_inner = train_evaluate_model_params_inner,
        choose_best_inner = choose_best_inner,
        train_evaluate_model_params_outer = train_evaluate_model_params_outer,
        train_evaluate_model_digest_outer = train_evaluate_model_digest_outer
    )
    class(out) <- "nmr_data_analysis_method"
    out
}
