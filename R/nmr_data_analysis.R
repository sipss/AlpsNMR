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
        keep_balance <- data.frame(keep_together = keep_together,
                                   balance_in_train = balance_in_train,
                                   stringsAsFactors = TRUE)
        keep_balance <- unique(keep_balance)
        
        smaller_balancing_group <- min(table(keep_balance$balance_in_train))
        num_groups_train <- max(floor(smaller_balancing_group*(1 - test_size)), 1)
        
        for (i in seq_len(iterations)) {
            groups_train <- keep_balance %>%
                dplyr::group_by(balance_in_train) %>% 
                dplyr::sample_n(num_groups_train) %>% 
                dplyr::ungroup() %>%
                dplyr::pull(keep_together)
            train_samples <- sample_idx[keep_together %in% groups_train]
            test_samples <- base::setdiff(sample_idx, train_samples)
            output[[i]] <- list(training = train_samples,
                                test = test_samples)
            if (length(train_samples) == 0 || length(test_samples) == 0) {
              stop("Error in random_subsampling, set of cero length: increase number of samples of the lowest set")
            }
        }
    } else {
        groups_keep_together <- unique(keep_together)
        num_groups_test <- floor(length(groups_keep_together)*test_size)
        for (i in seq_len(iterations)) {
            groups_test <- resample(x = groups_keep_together, size = num_groups_test)
            test_samples <- sample_idx[keep_together %in% groups_test]
            train_samples <- base::setdiff(sample_idx, test_samples)
            output[[i]] <- list(training = train_samples,
                                test = test_samples)
            if (length(train_samples) == 0 || length(test_samples) == 0) {
              stop("Error in random_subsampling, set of cero length: increase number of samples of the lowest set")
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
    
    outer_cv_loop <- random_subsampling(sample_idx = sample_idx,
                                        iterations = external_val$iterations,
                                        test_size = external_val$test_size,
                                        keep_together = keep_together_data)
    outer_cv_iterations <- purrr::imap(outer_cv_loop, function(outer, outer_idx) {
        list(outer_iter = outer_idx,
                 outer_train = outer$training,
                 outer_test = outer$test)
    })
    names(outer_cv_iterations) <- as.character(purrr::map_int(outer_cv_iterations, "outer_iter"))
    
    inner_cv_iterations <- purrr::flatten(
        purrr::imap(outer_cv_loop, function(outer, outer_idx) {
            calib_idx <- outer$training
            blind_idx <- outer$test
            inner_cv <- random_subsampling(sample_idx = calib_idx,
                                           iterations = internal_val$iterations,
                                           test_size = internal_val$test_size,
                                           keep_together = keep_together_data[calib_idx])
            purrr::imap(inner_cv, function(inner_iter, inner_iter_idx) {
                list(outer_iter = outer_idx,
                         inner_iter = inner_iter_idx,
                         inner_train_idx = inner_iter$training,
                         inner_test_idx = inner_iter$test)
            })
        }))
    
    names(inner_cv_iterations) <- 
        paste0(as.character(purrr::map_int(inner_cv_iterations, "outer_iter")),
                     "_", 
                     as.character(purrr::map_int(inner_cv_iterations, "inner_iter")))
    
    list(outer = outer_cv_iterations,
         inner = inner_cv_iterations)
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
      stop(y_column, "column not found in the dataset")
    }
    if (!is.null(identity_column)) {
        identity_all <- nmr_meta_get_column(dataset, column = identity_column)
    } else {
        identity_all <- NULL
    }
    
    # Split
    train_idx <- train_test_subset[[1]]
    test_idx <- train_test_subset[[2]]
    
    x_train <- x_all[train_idx,, drop = FALSE]
    y_train <- y_all[train_idx]
    if (!is.null(identity_all)) {
        identity_train <- identity_all[train_idx]
    } else {
        identity_train <- NULL
    }
    
    x_test <- x_all[test_idx,, drop = FALSE]
    y_test <- y_all[test_idx]
    if (!is.null(identity_all)) {
        identity_test <- identity_all[test_idx]
    } else {
        identity_test <- NULL
    }
    
    # Build & performance:
    train_evaluate_model(x_train = x_train, y_train = y_train, identity_train = identity_train,
                         x_test = x_test, y_test = y_test, identity_test = identity_test, ...)
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
                  train_test_subsets, train_evaluate_model_args_iter = NULL, ...) {
    if (show_progress_bar(length(train_test_subsets) > 5)) {
        prgrs <- TRUE
    } else {
        prgrs <- FALSE
    }

    output <- furrr::future_pmap(
        c(list(train_test_subset = train_test_subsets),
            train_evaluate_model_args_iter),
        split_build_perform,
        dataset = dataset,
        y_column = y_column,
        identity_column = identity_column,
        train_evaluate_model = train_evaluate_model,
        ...,
        .progress = prgrs,
        .options = furrr::future_options(globals = character(0),
                                         packages = character(0)))
    names(output) <- names(train_test_subsets)
    output
}


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
#'     Condition = rep(c("A", "B"), times = num_samples/2),
#'     stringsAsFactors = FALSE
#' )
#' 
#' ### The matrix with peaks
#' peak_means <- runif(n = num_peaks, min = 300, max = 600)
#' peak_sd <- runif(n = num_peaks, min = 30, max = 60)
#' peak_matrix <- mapply(function(mu, sd) rnorm(num_samples, mu, sd),
#'                                             mu = peak_means, sd = peak_sd)
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
#' print(model$outer_cv_results_digested$auroc)
#' ## Rank Product of the Variable Importance in the Projection
#' ## (Lower means more important)
#' print(sort(model$outer_cv_results_digested$vip_rankproducts))
#' 
#' @export
nmr_data_analysis <- function(dataset,
                              y_column,
                              identity_column, 
                              external_val,
                              internal_val,
                              data_analysis_method) {
    
    train_evaluate_model <- data_analysis_method[["train_evaluate_model"]]
    train_evaluate_model_params_inner <- data_analysis_method[["train_evaluate_model_params_inner"]]
    choose_best_inner <- data_analysis_method[["choose_best_inner"]]
    train_evaluate_model_params_outer <- data_analysis_method[["train_evaluate_model_params_outer"]]
    train_evaluate_model_digest_outer <- data_analysis_method[["train_evaluate_model_digest_outer"]]
    
    # Prepare double cross-validation splits:
    train_test_blind_subsets <- split_double_cv(dataset,
                                                keep_together = identity_column,
                                                external_val = external_val,
                                                internal_val = internal_val)
    
    # These are ALL the inner cross-validation iterations:
    train_test_subsets_inner <- purrr::map(train_test_blind_subsets$inner,
                                           ~ .[c("inner_train_idx", "inner_test_idx")])
    
    # We run the train_evaluate_model function for each of the inner CV.
    inner_cv_results <- do.call(
        what = do_cv,
        args = c(list(dataset = dataset,
                      y_column = y_column,
                      identity_column = identity_column,
                      train_evaluate_model = train_evaluate_model,
                      train_test_subsets = train_test_subsets_inner),
                 train_evaluate_model_params_inner))
    
    # We choose the best hyper-parameters for each inner cross-validation:
    inner_cv_results_digested <- choose_best_inner(inner_cv_results)
    
    # Prepare the indices for the final outer cv models:
    train_test_subsets_outer <- purrr::map(train_test_blind_subsets$outer,
                                           ~ .[c("outer_train", "outer_test")])
    
    # Compute the outer cv models
    outer_cv_results <- do.call(
        what = do_cv,
        args = c(list(dataset = dataset,
                      y_column = y_column,
                      identity_column = identity_column,
                      train_evaluate_model = train_evaluate_model,
                      train_test_subsets = train_test_subsets_outer,
                      train_evaluate_model_args_iter = inner_cv_results_digested$train_evaluate_model_args
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


#' Bootstrap and permutation over PLS-VIP
#' 
#' Bootstrap and permutation over PLS-VIP on AlpsNMR can be performed on both 
#' [nmr_dataset_1D] full spectra as well as [nmr_dataset_peak_table] peak tables.
#' 
#' Use of the bootstrap and permutation methods for a more robust
#' variable importance in the projection metric for partial least
#' squares regression
#' 
#' @param dataset An [nmr_dataset_family] object
#' @param train_index set of index used to generate the bootstrap datasets
#' @param y_column A string with the name of the y column (present in the
#'    metadata of the dataset)
#' @param ncomp number of components used in the plsda models
#' @param nbootstrap number of bootstrap dataset
#' @return A list with the following elements:
#' 
#' - `important_vips`: A list with the important vips selected
#' - `relevant_vips`: List of vips with some relevance
#' - `pls_vip`: Pls-VIPs of every bootstrap
#' - `pls_vip_perm`: Pls-VIPs of every bootstrap with permuted variables
#' - `pls_mean`: Pls-VIPs normaliced differences means
#' - `pls_vip_score_diff`: Differences of `pls_vip` and `pls_vip_perm`
#' - `error`: error spected in a t distribution
#' - `lower_bound`: lower bound of the confidence interval
#' - `upper_bound`: upper bound of the confidence interval
#' - `aucroc`: Auroc measures for validation of the method
#' 
#' @importFrom stats qt
#' @examples 
#' # Data analysis for a table of integrated peaks
#' 
#' ## Generate an artificial nmr_dataset_peak_table:
#' ### Generate artificial metadata:
#' num_samples <- 32 # use an even number in this example
#' num_peaks <- 20
#' metadata <- data.frame(
#'     NMRExperiment = as.character(1:num_samples),
#'     Condition = rep(c("A", "B"), times = num_samples/2),
#'     stringsAsFactors = FALSE
#' )
#' 
#' ### The matrix with peaks
#' peak_means <- runif(n = num_peaks, min = 300, max = 600)
#' peak_sd <- runif(n = num_peaks, min = 30, max = 60)
#' peak_matrix <- mapply(function(mu, sd) rnorm(num_samples, mu, sd),
#'                                             mu = peak_means, sd = peak_sd)
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
#'     external_val = list(iterations = 1, test_size = 0.25),
#'     internal_val = list(iterations = 3, test_size = 0.25),
#'     data_analysis_method = methodology
#' )
#' ## Area under ROC for each outer cross-validation iteration:
#' print(model$outer_cv_results_digested$auroc)
#' ## Rank Product of the Variable Importance in the Projection
#' ## (Lower means more important)
#' print(sort(model$outer_cv_results_digested$vip_rankproducts))
#' 
#' ## The number of components for the bootstrap models is selected 
#' ncomps <- model$outer_cv_results$`1`$model$ncomp
#' train_index <- model$train_test_partitions$outer$`1`$outer_train
#' 
#' # Bootstrap and permutation for VIP selection
#' bp_VIPS <- bp_VIP_analysis(peak_table, # Data to be analized
#'                            train_index,
#'                            y_column = "Condition", # Label
#'                            ncomp = ncomps,
#'                            nbootstrap = 100)
#'
#' aucs <- max(bp_VIPS$aucroc$aucs$auc)
#' message("AUC of the Bootstrap and permutation ", aucs)
#' @export
bp_VIP_analysis <- function(dataset,
                            train_index,
                            y_column,
                            ncomp,
                            nbootstrap = 300) {

    # Extract data and split for train and test
    x_all <- dataset$peak_table
    y_all <- nmr_meta_get_column(dataset, column = y_column)
    x_train <- x_all[train_index,, drop = FALSE]
    y_train <- y_all[train_index]
    x_test <- x_all[-train_index,, drop = FALSE]
    y_test <- y_all[-train_index]
    
    n <- dim(x_all)[2]
    names <- colnames(x_all)
    if (length(names) == 0) {
        stop("Error in bp_VIP_analysis, the dataset peak_table don't have colnames.")
    }    

    # # Bootstrap with replacement nbootstraps datasets
    # # Begin parallel processing
    # plan(multiprocess)
    # results <- furrr::future_map(seq_len(nbootstrap), function(i) {
    #     index <- sample(1:nrow(x_train),nrow(x_train), rep = TRUE)
    #     x_train_boots <- x_train[index,]
    #     y_train_boots <- y_train[index]
    # 
    #     # Fit PLS model
    #     model <-
    #         plsda_build(
    #             x = x_train_boots,
    #             y = y_train_boots,
    #             identity = NULL,
    #             ncomp = ncomp
    #         )
    #     # VIPs per component extraction
    #     pls_vip_comps <- plsda_vip(model)
    #     # Sum contributions of VIPs to each component
    #     pls_vip <- sqrt(rowSums(pls_vip_comps^2)/ncomp)
    # 
    #     # Permutation of variables
    #     for (j in seq_len(n)) {
    #         random_pos <- sample(seq_len(n), 1)
    #         x_train_boots_perm <- x_train_boots
    #         x_train_boots_perm[,j] <- x_train_boots[, random_pos]
    # 
    #         # Refit model with permuted variables
    #         model_perm <-
    #             plsda_build(
    #                 x = x_train_boots_perm,
    #                 y = y_train_boots,
    #                 identity = NULL,
    #                 ncomp = ncomp
    #             )
    #         # VIPs per component extraction
    #         pls_vip_comps_perm <- plsda_vip(model_perm)
    #         # Sum contributions of VIPs to each component
    #         pls_vip_perm[,j] <- sqrt(rowSums(pls_vip_comps_perm^2)/ncomp)
    #     }
    # 
    #     # bootsrapped and randomly permuted PLS-VIPs
    #     pls_vip_perm_score <- colSums(pls_vip_perm)/n
    #     # bootsrapped and randomly permuted standard desviation
    #     #pls_vip_perm_sd <- sqrt(sum((pls_vip_perm - pls_vip_perm_score) ^ 2) / (nbootstrap - 1))
    #     # bootsrapped and randomly permuted difference
    #     pls_vip_score_diff <- pls_vip - pls_vip_perm_score
    #     list(pls_vip = pls_vip,
    #          pls_vip_perm = pls_vip_perm,
    #          pls_vip_perm_score = pls_vip_perm_score,
    #          pls_vip_score_diff = pls_vip_score_diff)
    # })
    # 
    # # Transform list results to data frame
    # results_df <- as.data.frame(do.call(rbind, results), stringsAsFactors=FALSE)
    # pls_vip <- as.data.frame(do.call(rbind, results_df$pls_vip), stringsAsFactors=FALSE)
    # pls_vip_perm <- as.data.frame(do.call(rbind, results_df$pls_vip_perm), stringsAsFactors=FALSE)
    # pls_vip_perm_score <- as.data.frame(do.call(rbind, results_df$pls_vip_perm_score), stringsAsFactors=FALSE)
    # pls_vip_score_diff <- as.data.frame(do.call(rbind, results_df$pls_vip_score_diff), stringsAsFactors=FALSE)
    # 
    # # Normalization of the difference vector for each variable to
    # # its corresponding standard deviation and construct
    # # 95% confidence intervals around the differences
    # boots_vip <- matrix(nrow = n, dimnames = list(names))
    # boots_vip_sd <- matrix(nrow = n, dimnames = list(names))
    # error <- matrix(nrow = n, dimnames = list(names))
    # lower_bound <- matrix(nrow = n, dimnames = list(names))
    # upper_bound <- matrix(nrow = n, dimnames = list(names))
    # for (k in seq_len(n)){
    #     element <- pls_vip_score_diff[,k] / sd(pls_vip_score_diff[,k])
    #     boots_vip[k] <- sum(element)/nbootstrap
    #     boots_vip_sd[k] <- sqrt(sum((element - boots_vip[k])^2)/(nbootstrap-1))
    #     error[k] <- qt(0.975, df = nbootstrap - 1) * boots_vip_sd[k]
    #     lower_bound[k] <- boots_vip[k] - error[k]
    #     upper_bound[k] <- boots_vip[k] + error[k]
    # }
    # 
    # important_vips <- names[lower_bound > qt(0.975, df = nbootstrap - 1)]
    # relevant_vips <- names[lower_bound <= qt(0.975, df = nbootstrap - 1) & lower_bound > 0]
    # if (length(important_vips) == 0) {
    #     if (length(relevant_vips) == 0) {
    #         stop("Error in bp_VIP_analysis, none of the variables seems relevant:\n
    #          try increasing the number of bootstraps")
    #     }
    #     warning("No VIPs are ranked as important, use the relevant_vips or try again with more bootstraps")
    # }
    # 
    # # Building a model with just the important vips to check performance
    # # Spliting test
    # index <- sample(1:nrow(x_test),nrow(x_test)*0.75, rep = FALSE)
    # x_train_selected <- x_test[index, important_vips]
    # y_train_selected <- y_test[index]
    # x_test_selected <- x_test[-index, important_vips]
    # y_test_selected <- y_test[-index]
    # 
    # model_test <- plsda_build(
    #     x = x_train_selected,
    #     y = y_train_selected,
    #     identity = NULL,
    #     ncomp = ncomp
    # )
    # aucroc <- plsda_auroc(model_test, x_test_selected, y_test_selected, NULL)
    # 
    # # To return it ordered by mean of the normalized vectors
    # orden <- order(boots_vip, decreasing = TRUE)
    # #Return important vips and auc performance
    # list(important_vips = important_vips,
    #      relevant_vips = relevant_vips,
    #      pls_vip = pls_vip[orden,,drop=FALSE],
    #      pls_vip_perm = pls_vip_perm_score[orden,,drop=FALSE],
    #      pls_mean = boots_vip[orden,,drop=FALSE],
    #      pls_vip_score_diff = pls_vip_score_diff[orden,,drop=FALSE],
    #      error = error[orden,,drop=FALSE],
    #      lower_bound = lower_bound[orden,,drop=FALSE],
    #      upper_bound = upper_bound[orden,,drop=FALSE],
    #      aucroc = aucroc)
    
    
    pls_vip <- matrix(nrow = n, ncol = nbootstrap, dimnames = list(names, NULL))
    pls_vip_perm <- matrix(nrow = n, ncol = n, dimnames = list(names, NULL))
    pls_vip_perm_score <- matrix(nrow = n, ncol = nbootstrap, dimnames = list(names, NULL))
    pls_vip_perm_sd <- matrix(nrow = n, ncol = nbootstrap, dimnames = list(names, NULL))
    pls_vip_score_diff <- matrix(nrow = n, ncol = nbootstrap, dimnames = list(names, NULL)) # Bootstrap with replacement nbootstraps datasets
    
    # Bootstrap with replacement nbootstraps datasets
    for (i in seq_len(nbootstrap)) {
        index <- sample(1:nrow(x_train),nrow(x_train), rep = TRUE)
        x_train_boots <- x_train[index,]
        y_train_boots <- y_train[index]

        # Fit PLS model
        model <-
            plsda_build(
                x = x_train_boots,
                y = y_train_boots,
                identity = NULL,
                ncomp = ncomp
            )
        # VIPs per component extraction
        pls_vip_comps <- plsda_vip(model)
        # Sum contributions of VIPs to each component
        pls_vip[,i] <- sqrt(rowSums(pls_vip_comps^2)/ncomp)

        # Permutation of variables
        for (j in seq_len(n)) {
            random_pos <- sample(seq_len(n), 1)
            x_train_boots_perm <- x_train_boots
            x_train_boots_perm[,j] <- x_train_boots[, random_pos]

            # Refit model with permuted variables
            model_perm <-
                plsda_build(
                    x = x_train_boots_perm,
                    y = y_train_boots,
                    identity = NULL,
                    ncomp = ncomp
                )
            # VIPs per component extraction
            pls_vip_comps_perm <- plsda_vip(model_perm)
            # Sum contributions of VIPs to each component
            pls_vip_perm[,j] <- sqrt(rowSums(pls_vip_comps_perm^2)/ncomp)
        }

        # bootsrapped and randomly permuted PLS-VIPs
        pls_vip_perm_score[,i] <- colSums(pls_vip_perm)/n
        # bootsrapped and randomly permuted standard desviation
        #pls_vip_perm_sd[,i] <- sqrt(sum((pls_vip_perm - pls_vip_perm_score[,i]) ^ 2) / (nbootstrap - 1))
        # bootsrapped and randomly permuted difference
        pls_vip_score_diff[,i] <- pls_vip[,i] - pls_vip_perm_score[,i]
    }

    # Normalization of the difference vector for each variable to
    # its corresponding standard deviation and construct
    # 95% confidence intervals around the differences
    boots_vip <- matrix(nrow = n, dimnames = list(names))
    boots_vip_sd <- matrix(nrow = n, dimnames = list(names))
    error <- matrix(nrow = n, dimnames = list(names))
    lower_bound <- matrix(nrow = n, dimnames = list(names))
    upper_bound <- matrix(nrow = n, dimnames = list(names))
    for (k in seq_len(n)){
        element <- pls_vip_score_diff[k,] / sd(pls_vip_score_diff[k,])
        boots_vip[k] <- sum(element)/nbootstrap
        boots_vip_sd[k] <- sqrt(sum((element - boots_vip[k])^2)/(nbootstrap-1))
        error[k] <- qt(0.975, df = nbootstrap - 1) * boots_vip_sd[k]
        lower_bound[k] <- boots_vip[k] - error[k]
        upper_bound[k] <- boots_vip[k] + error[k]
    }
    
    important_vips <- names[lower_bound > qt(0.975, df = nbootstrap - 1)]
    relevant_vips <- names[lower_bound <= qt(0.975, df = nbootstrap - 1) & lower_bound > 0]
    if (length(important_vips) == 0) {
        if (length(relevant_vips) == 0) {
            stop("Error in bp_VIP_analysis, none of the variables seems relevant:\n
             try increasing the number of bootstraps")
        }
        warning("No VIPs are ranked as important, use the relevant_vips or try again with more bootstraps")
    }
    
    # Building a model with just the important vips to check performance
    # Spliting test
    index <- sample(1:nrow(x_test),nrow(x_test)*0.75, rep = FALSE)
    x_train_selected <- x_test[index, important_vips]
    y_train_selected <- y_test[index]
    x_test_selected <- x_test[-index, important_vips]
    y_test_selected <- y_test[-index]
    
    model_test <- plsda_build(
        x = x_train_selected,
        y = y_train_selected,
        identity = NULL,
        ncomp = ncomp
    )
    aucroc <- plsda_auroc(model_test, x_test_selected, y_test_selected, NULL)
    
    # To return it ordered by mean of the normalized vectors
    orden <- order(boots_vip, decreasing = TRUE)
    #Return important vips and auc performance
    list(important_vips = important_vips,
         relevant_vips = relevant_vips,
         pls_vip = pls_vip[orden,,drop=FALSE],
         pls_vip_perm = pls_vip_perm_score[orden,,drop=FALSE],
         pls_mean = boots_vip[orden,,drop=FALSE],
         pls_vip_score_diff = pls_vip_score_diff[orden,,drop=FALSE],
         error = error[orden,,drop=FALSE],
         lower_bound = lower_bound[orden,,drop=FALSE],
         upper_bound = upper_bound[orden,,drop=FALSE],
         aucroc = aucroc)
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
#' @examples 
#' help(new_nmr_data_analysis_method)
#' 
new_nmr_data_analysis_method <- function(train_evaluate_model,
                                         train_evaluate_model_params_inner,
                                         choose_best_inner,
                                         train_evaluate_model_params_outer,
                                         train_evaluate_model_digest_outer) {
    out <- list(train_evaluate_model = train_evaluate_model,
                            train_evaluate_model_params_inner = train_evaluate_model_params_inner,
                            choose_best_inner = choose_best_inner,
                            train_evaluate_model_params_outer = train_evaluate_model_params_outer,
                            train_evaluate_model_digest_outer = train_evaluate_model_digest_outer)
    class(out) <- "nmr_data_analysis_method"
    out
}

#' Permutation test
#'
#' Make permutations with data and default settings from an nmr_data_analysis_method
#' 
#' @param nPerm number of permutations
#' 
#' @inheritParams nmr_data_analysis
#' 
#' @return A permutation matrix with permuted values
#' @name permutation_test_model
#' @export
#' @examples
#' # Data analysis for a table of integrated peaks
#' 
#' ## Generate an artificial nmr_dataset_peak_table:
#' ### Generate artificial metadata:
#' num_samples <- 32 # use an even number in this example
#' num_peaks <- 20
#' metadata <- data.frame(
#'     NMRExperiment = as.character(1:num_samples),
#'     Condition = rep(c("A", "B"), times = num_samples/2),
#'     stringsAsFactors = FALSE
#' )
#' 
#' ### The matrix with peaks
#' peak_means <- runif(n = num_peaks, min = 300, max = 600)
#' peak_sd <- runif(n = num_peaks, min = 30, max = 60)
#' peak_matrix <- mapply(function(mu, sd) rnorm(num_samples, mu, sd),
#'                                             mu = peak_means, sd = peak_sd)
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
#' methodology <- plsda_auroc_vip_method(ncomp = 3)
#' model <- nmr_data_analysis(
#'     peak_table,
#'     y_column = "Condition",
#'     identity_column = NULL,
#'     external_val = list(iterations = 3, test_size = 0.25),
#'     internal_val = list(iterations = 3, test_size = 0.25),
#'     data_analysis_method = methodology
#' )
#' 
#' p = permutation_test_model(peak_table,
#'                            y_column = "Condition",
#'                            identity_column = NULL,
#'                            external_val = list(iterations = 3, test_size = 0.25),
#'                            internal_val = list(iterations = 3, test_size = 0.25),
#'                            data_analysis_method = methodology, 
#'                            nPerm = 10)
#'                            
permutation_test_model = function (dataset,
                                   y_column,
                                   identity_column, 
                                   external_val,
                                   internal_val,
                                   data_analysis_method,
                                   nPerm = 50)
{
    startTime=proc.time()[3]
    permMatrix=matrix(ncol=1,nrow=nPerm)
    #colnames(permMatrix)=c('Min','Mid','Max')
    dataset_perm <- dataset
    y_all <- nmr_meta_get_column(dataset, column = y_column)
    for (p in seq_len(nPerm)) {
        cat('\n permutation ',p,' of ',nPerm,'\n',sep = '')
        
        #Permutar columna y_colum del dataset
        YPerm=sample(y_all)
        dataset_perm[["metadata"]][["external"]][[y_column]] <- YPerm
        # print(sum(y_test!=y_all))
        permMod=nmr_data_analysis(dataset_perm,
                                  y_column = y_column,
                                  identity_column = identity_column,
                                  external_val = external_val,
                                  internal_val = internal_val,
                                  data_analysis_method = data_analysis_method)
        
        # I will use the mean of the auc of all the outer_cv for the test static
        test_stat = mean(permMod$outer_cv_results_digested$auroc$auc)

        permMatrix[p,1] = test_stat
        nowTime=proc.time()[3]
        timePerRep=(nowTime-startTime)/p
        timeLeft=(timePerRep*(nPerm-p))/60
        cat('\nEstimated time left:',timeLeft,'mins\n\n')
    }

    return (permMatrix)
}

#' Permutation test plot
#'
#' Plot permutation test using actual model and permutated models
#'
#' @param nmr_data_analysis_model A nmr_data_analysis_model
#' @param permMatrix A permutation fitness outcome from permutation_test_model
#' @param xlab optional xlabel
#' @param xlim optional x-range
#' @param ylim otional y-range
#' @param breaks optional custom histogram breaks (defaults to 'sturges')
#' @param main optional plot title (or TRUE for autoname)
#' 
#' @importFrom graphics axis hist lines text
#' @importFrom stats median pt sd ecdf na.omit
#' @return A plot with the comparison between the actual model versus the permuted models
#' @name permutation_test_plot
#' @export
#' @examples
#'# Data analysis for a table of integrated peaks
#' 
#' ## Generate an artificial nmr_dataset_peak_table:
#' ### Generate artificial metadata:
#' num_samples <- 32 # use an even number in this example
#' num_peaks <- 20
#' metadata <- data.frame(
#'     NMRExperiment = as.character(1:num_samples),
#'     Condition = rep(c("A", "B"), times = num_samples/2),
#'     stringsAsFactors = FALSE
#' )
#' 
#' ### The matrix with peaks
#' peak_means <- runif(n = num_peaks, min = 300, max = 600)
#' peak_sd <- runif(n = num_peaks, min = 30, max = 60)
#' peak_matrix <- mapply(function(mu, sd) rnorm(num_samples, mu, sd),
#'                                             mu = peak_means, sd = peak_sd)
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
#' methodology <- plsda_auroc_vip_method(ncomp = 3)
#' model <- nmr_data_analysis(
#'     peak_table,
#'     y_column = "Condition",
#'     identity_column = NULL,
#'     external_val = list(iterations = 3, test_size = 0.25),
#'     internal_val = list(iterations = 3, test_size = 0.25),
#'     data_analysis_method = methodology
#' )
#' 
#' p = permutation_test_model(peak_table,
#'                            y_column = "Condition",
#'                            identity_column = NULL,
#'                            external_val = list(iterations = 3, test_size = 0.25),
#'                            internal_val = list(iterations = 3, test_size = 0.25),
#'                            data_analysis_method = methodology, 
#'                            nPerm = 10)
#'                            
#' permutation_test_plot(model, p)
#' 
permutation_test_plot = function (nmr_data_analysis_model, 
                                  permMatrix,
                                  xlab = "AUCs",
                                  xlim,
                                  ylim = NULL,
                                  breaks = "Sturges",
                                  main = "Permutation test")
{
    
    h0=permMatrix[,1]
    if(missing(xlim)) {
        xlim=c(0,1)
    }
    h=hist(permMatrix,breaks,xlim=xlim,ylim=ylim,axes=FALSE,xlab=xlab,freq=FALSE,main=main)
    h2=max(h$density)*.75
    axis(1,pos=0)
    axis(2,pos=0,las=1)
    
    model_auc = mean(nmr_data_analysis_model$outer_cv_results_digested$auroc$auc)
    lines(rep(model_auc,2),c(0,h2))
    
    p=ecdf(h0)(model_auc) # Empirical
    #p1Stud=pt((model_auc-(mean(h0)/2))/sd(h0/2),(length(h0)-1)) # Students
    
    # warning: sometimes the auc of the permutation is NaN
    h0_median <- apply(permMatrix,2,function(x) median(na.omit(x)))
    
    pP=ifelse(model_auc<h0_median, p, 1-p)
    if(pP<1/length(h0)){
      text(h2,pos=2,labels=paste('p<',signif(1/length(h0),4),sep=''))
    } else {
      text(h2,pos=2,labels=paste('p=',signif(pP,4),sep=''))
    }
}
