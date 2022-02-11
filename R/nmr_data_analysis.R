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
      rlang::abort(
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

    warn_future_to_biocparallel()
    output <- do.call(
        what = BiocParallel::bpmapply,
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
        args = c(
            list(
                dataset = dataset,
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
#' @name bp_VIP_analysis
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
#' - `pls_vip_means`: Pls-VIPs normaliced differences means
#' - `pls_vip_score_diff`: Differences of `pls_vip` and `pls_vip_perm`
#' - `pls_models`: pls models of the diferent bootstraps
#' - `pls_perm_models`: pls permuted models of the diferent bootstraps
#' - `classif_rate`: classification rate of the bootstrap models
#' - `general_model`: pls model trained with all train data
#' - `general_CR`: classification rate of the `general_model`
#' - `vips_model`: pls model trained with vips selection over all train data
#' - `vips_CR`: classification rate of the `vips_model`
#' - `error`: error spected in a t distribution
#' - `lower_bound`: lower bound of the confidence interval
#' - `upper_bound`: upper bound of the confidence interval
#' 
#' @importFrom stats qt
#' @importFrom stats wilcox.test
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
#' model$outer_cv_results_digested$auroc
#' 
#' ## The number of components for the bootstrap models is selected 
#' ncomps <- model$outer_cv_results$`1`$model$ncomp
#' train_index <- model$train_test_partitions$outer$`1`$outer_train
#' 
#' # Bootstrap and permutation for VIP selection
#' bp_VIPS <- bp_VIP_analysis(peak_table, # Data to be analized
#'                            train_index,
#'                            y_column = "Condition",
#'                            ncomp = ncomps,
#'                            nbootstrap = 10)
#'
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
    # For check performance
    x_test <- x_all[-train_index,, drop = FALSE]
    y_test <- y_all[-train_index]
    
    if (length(unique(y_train)) == 1) {
        stop("Only one class in train set, increase number of samples")
    }
    if (length(unique(y_test)) == 1) {
        stop("Only one class in test set, increase number of samples")
    }
    
    n <- dim(x_all)[2]
    names <- colnames(x_all)
    #some checks
    if (length(names) == 0) {
        stop("Error in bp_VIP_analysis, the dataset peak_table don't have colnames.")
    }
    
    pls_vip <- matrix(nrow = n, ncol = nbootstrap, dimnames = list(names, NULL))
    pls_vip_perm <- matrix(nrow = n, ncol = n, dimnames = list(names, NULL))
    pls_vip_perm_score <- matrix(nrow = n, ncol = nbootstrap, dimnames = list(names, NULL))
    pls_vip_perm_sd <- matrix(nrow = n, ncol = nbootstrap, dimnames = list(names, NULL))
    pls_vip_score_diff <- matrix(nrow = n, ncol = nbootstrap, dimnames = list(names, NULL)) # Bootstrap with replacement nbootstraps datasets
    pls_models <- list()
    pls_perm_models <- list()
    CR <- list()
    
    # Bootstrap with replacement nbootstraps datasets
    for (i in seq_len(nbootstrap)) {
        index <- sample(seq_len(nrow(x_train)),nrow(x_train), replace = TRUE)
        x_train_boots <- x_train[index,]
        y_train_boots <- y_train[index]
        
        #if y_train_boots only have one class, plsda models give error
        if (length(unique(y_train_boots)) == 1) {
            #Replace the first element for the first element of another class
            for(i in seq_len(y_train)){
                if (y_train[i] != y_train_boots[1]){
                    y_train_boots[1] = y_train[i]
                    x_train_boots[1] = x_train[i]
                    break
                }
            }
        }
        
        # Rename rownames, because if they are repeated, plsda fails
        rownames(x_train_boots) <- paste0("Sample", seq_len(dim(x_train_boots)[1]))
        
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
        # Measure the classification rate (CR) of the bootstrap model
        perf <- mixOmics::perf(model, newdata = x_test)
        CR[[i]] <- 1 - perf$error.rate$overall[1]
                 
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
        # bootsrapped and randomly permuted difference
        pls_vip_score_diff[,i] <- pls_vip[,i] - pls_vip_perm_score[,i]
        pls_models[[i]] <- model
        pls_perm_models[[i]] <- model_perm
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
    relevant_vips <- names[lower_bound > 0]
    
    # Chequing performance
    # Fit PLS model
    general_model <-
        plsda_build(
            x = x_train,
            y = y_train,
            identity = NULL,
            ncomp = ncomp
        )
    # Measure the classification rate (CR) of the fold
    perf <- mixOmics::perf(general_model, newdata = x_test)
    general_CR <- 1 - perf$error.rate$overall[1]
    
    
    if (length(important_vips) == 0) {
        if (length(relevant_vips) == 0) {
            warning(
                "Error in bp_VIP_analysis, none of the variables seems relevant:\n
             try increasing the number of bootstraps"
            )
        }
        warning(
            "No VIPs are ranked as important, use the relevant_vips or try again with more bootstraps"
        )
        vips_model <- NULL
        vips_CR <- 0
    } else {
        # Chequing performance of selected vips
        if (length(important_vips) == 1) {
            if (length(relevant_vips) == 1) {
                warning(
                    "Only one VIP ranked as important and relevant, you can try again with more bootstraps"
                )
                vips_model <- NULL
                vips_CR <- 0
            } else {
                # if the are only one importan vip, we use relevants instead
                x_train_reduced <-
                    as.matrix(x_all[train_index, relevant_vips, drop = FALSE])
                x_test_reduced <-
                    as.matrix(x_all[-train_index, relevant_vips, drop = FALSE])
                
                # Fit PLS model
                vips_model <-
                    plsda_build(
                        x = x_train_reduced,
                        y = y_train,
                        identity = NULL,
                        ncomp = ncomp
                    )
                
                # Measure the classification rate (CR) of the fold
                perf <-
                    mixOmics::perf(vips_model, newdata = x_test_reduced)
                vips_CR  <- 1 - perf$error.rate$overall[1]
            }
        } else {
            x_train_reduced <-
                as.matrix(x_all[train_index, important_vips, drop = FALSE])
            x_test_reduced <-
                as.matrix(x_all[-train_index, important_vips, drop = FALSE])
            
            # Fit PLS model
            vips_model <-
                plsda_build(
                    x = x_train_reduced,
                    y = y_train,
                    identity = NULL,
                    ncomp = ncomp
                )
            
            # Measure the classification rate (CR) of the fold
            perf <- mixOmics::perf(vips_model, newdata = x_test_reduced)
            vips_CR  <- 1 - perf$error.rate$overall[1]
        }
    }
    
    # To return it ordered by mean of the normalized vectors
    orden <- order(boots_vip, decreasing = TRUE)
    #Return important vips and auc performance
    list(important_vips = important_vips,
         relevant_vips = relevant_vips,
         pls_vip = pls_vip[orden,,drop=FALSE],
         pls_vip_perm = pls_vip_perm_score[orden,,drop=FALSE],
         pls_vip_means = boots_vip[orden,,drop=FALSE],
         pls_vip_score_diff = pls_vip_score_diff[orden,,drop=FALSE],
         pls_models = pls_models,
         pls_perm_models = pls_perm_models,
         classif_rate = CR,
         general_model = general_model,
         general_CR = general_CR,
         vips_model = vips_model,
         vips_CR = vips_CR,
         error = error[orden,,drop=FALSE],
         lower_bound = lower_bound[orden,,drop=FALSE],
         upper_bound = upper_bound[orden,,drop=FALSE])
}


#' K-fold bootstrap and permutation over PLS-VIP
#' 
#' Bootstrap and permutation over PLS-VIP on AlpsNMR can be performed on both 
#' [nmr_dataset_1D] full spectra as well as [nmr_dataset_peak_table] peak tables.
#' 
#' Use of the bootstrap and permutation methods for a more robust
#' variable importance in the projection metric for partial least
#' squares regression, in a k-fold cross validation
#' 
#' @name bp_kfold_VIP_analysis
#' @param dataset An [nmr_dataset_family] object
#' @param y_column A string with the name of the y column (present in the
#'    metadata of the dataset)
#' @param k Number of folds, recomended between 4 to 10
#' @param ncomp number of components for the bootstrap models
#' @param nbootstrap number of bootstrap dataset
#' @return A list with the following elements:
#' 
#' - `important_vips`: A list with the important vips selected
#' - `relevant_vips`: List of vips with some relevance
#' - `wilcoxon_vips`: List of vips that pass a wilcoxon test
#' - `vip_means`: Means of the vips scores
#' - `vip_score_plot`: plot of the vips scores
#' - `kfold_resuls`: results of the k [bp_VIP_analysis]
#' - `kfold_index`: list of index of partitions of the folds
#' 
#' @export
#' @examples 
#' # Data analysis for a table of integrated peaks
#' 
#' ## Generate an artificial nmr_dataset_peak_table:
#' ### Generate artificial metadata:
#' num_samples <- 64 # use an even number in this example
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
#' rownames(peak_matrix) <- paste0("Sample", 1:num_samples)
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
#' ## We will use bootstrap and permutation method for VIPs selection 
#' ## in a a k-fold cross validation 
#' bp_results <- bp_kfold_VIP_analysis(peak_table, # Data to be analized
#'                            y_column = "Condition", # Label
#'                            k = 3,
#'                            nbootstrap = 10)
#'
#' message("Selected VIPs are: ", bp_results$importarn_vips)
#' 
bp_kfold_VIP_analysis <- function(dataset,
                            y_column,
                            k = 4,
                            ncomp = 3,
                            nbootstrap = 300) {

    if (k <= 1) {
        stop("K must be integer greater than 1")
    }
    
    # Extract data and split for train and test
    x_all <- dataset$peak_table
    y_all <- nmr_meta_get_column(dataset, column = y_column)
    if (length(unique(y_all)) == 1) {
        stop("Only one class in data set, at least two needed")
    }
    
    # Random and spliting
    x <- seq_len(length(y_all))
    index <- sample(x, replace = FALSE)
    x_all <- x_all[index,]
    y_all <- y_all[index]
    
    # Split data for k-fold
    k_fold_split <- split(x, x%%k)
    k_fold_index <- list()
    for(i in seq_len(k)){
        k_fold_index[[i]] <- seq_len(length(y_all))[-k_fold_split[[i]]]
    }
    
    results <- BiocParallel::bplapply(
        k_fold_index, function(index, dataset = dataset, y_column = y_column,
                               ncomp = ncomp, nbootstrap = nbootstrap) {
        bp_VIP_analysis(
            dataset,
            index,
            y_column = y_column,
            ncomp = ncomp,
            nbootstrap = nbootstrap
        )}, dataset = dataset, y_column = y_column, 
        ncomp = ncomp, nbootstrap = nbootstrap)

    # Mean of the vips of the different folds for the plot
    means <- purrr::map(results, "pls_vip_means")
    #means <- sapply(results, "[", "pls_vip_means")
    names_order <- sort(rownames(means[[1]]))
    ordered_means <- matrix(nrow = k, ncol = length(names_order), dimnames = list(NULL, names_order))
    for(i in seq_len(k)){
        ordered_means[i,] <- means[[i]][order(rownames(means[[i]]))]
    }
    ordered_means <- colSums(ordered_means)/k
    vip_means <- ordered_means[order(ordered_means, decreasing = TRUE)]
    error <- results[[1]]$error
    
    # Selection based on the means (deprecated, now ussing intersection of the vips)
    # important_vips <- vip_means[vip_means-2*error > 0]
    # relevant_vips <- vip_means[vip_means-error > 0]
    
    ## Wilcoxon test
    num_var <- dim(results[[1]]$pls_vip)[1]
    wt <- matrix(nrow = k, ncol = num_var)
    wt_vips <- list()
    for(i in seq_len(k)) {
        for (j in seq_len(num_var)) {
            x <- results[[i]]$pls_vip[j, ]
            y <- results[[i]]$pls_vip_perm[j, ]
            #wt_object <- wilcox.test(x, y, paired = TRUE, alternative = "two.sided")
            wt_object <- wilcox.test(x, y, paired = TRUE, alternative = "greater")
            wt[i,j] <- wt_object$p.value
        }
        wt_vips[[i]] <- rownames(results[[i]]$pls_vip)[wt[i,]<0.05]
    }
    
    # Plot of the scores
    x <- seq_len(length(vip_means))
    p <- ggplot2::ggplot() +
        ggplot2::geom_segment(
            mapping = ggplot2::aes(
                x = x,
                y = vip_means - error,
                xend = x,
                yend = vip_means + error
            ),
            arrow = NULL
        ) +
        ggplot2::geom_point(mapping = ggplot2::aes(x = x, y = vip_means),
                            shape = 21,
                            fill = "white") +
        ggplot2::geom_hline(yintercept = qt(0.975, df = nbootstrap - 1)) +
        ggplot2::ggtitle("BP-VIP") +
        ggplot2::labs(x = "Variables", y = "Scores") +
        ggplot2::theme_bw()
    
    # Reporting the intersection of the vips of the different folds
    list(important_vips = Reduce(intersect, purrr::map(results, "important_vips")),
         relevant_vips = Reduce(intersect, purrr::map(results, "relevant_vips")),
         wilcoxon_vips = unique(unlist(wt_vips)),
         vip_means = vip_means,
         vip_score_plot = p,
         kfold_results = results,
         kfold_index = k_fold_index)    
}


#' Plot vip scores of bootstrap
#'
#' @param vip_means vips means values of bootstraps
#' @param error error tolerated, calculated in the bootstrap
#' @param nbootstrap number of bootstraps realiced
#' @param plot A boolean that indicate if results are plotted or not
#'
#' @return A plot of the results or a ggplot object
#' @export
#' @examples
#' # Data analysis for a table of integrated peaks
#' 
#' ## Generate an artificial nmr_dataset_peak_table:
#' ### Generate artificial metadata:
#' num_samples <- 64 # use an even number in this example
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
#' ## We will use bootstrap and permutation method for VIPs selection 
#' ## in a a k-fold cross validation 
#' #bp_results <- bp_kfold_VIP_analysis(peak_table, # Data to be analized
#' #                           y_column = "Condition", # Label
#' #                           k = 3,
#' #                           nbootstrap = 10)
#'
#' #message("Selected VIPs are: ", bp_results$importarn_vips)
#' 
#' #plot_vip_scores(bp_results$kfold_results[[1]]$vip_means, 
#' #                bp_results$kfold_results[[1]]$error[1],
#' #                nbootstrap = 10)
#' 
plot_vip_scores <- function(vip_means, error, nbootstrap, plot = TRUE) {
    
    # Plot of the scores
    x <- seq_len(length(vip_means))
    vip_score_plot <- ggplot2::ggplot() +
        ggplot2::geom_segment(
            mapping = ggplot2::aes(
                x = x,
                y = vip_means - error,
                xend = x,
                yend = vip_means + error
            ),
            arrow = NULL
        ) +
        ggplot2::geom_point(mapping = ggplot2::aes(x = x, y = vip_means),
                            shape = 21,
                            fill = "white") +
        ggplot2::geom_hline(yintercept = qt(0.975, df = nbootstrap - 1)) +
        ggplot2::ggtitle("BP-VIP") +
        ggplot2::labs(x = "Variables", y = "Scores") +
        ggplot2::theme_bw()
    if(plot){
        vip_score_plot
    } else {
        return(vip_score_plot)
    }
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

#' Models stability plot
#'
#' Plot stability among models of the external cross validation
#'
#' @param model A nmr_data_analysis_model
#' 
#' @return A plot of models stability
#' @name models_stability_plot_plsda
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
#' #models_stability_plot_plsda(model)
#' 
models_stability_plot_plsda = function (model)
{
    # Loadings of the models
    ex_n = length(model$outer_cv_results)
    max_ncomp = 0
    for (n in seq_len(ex_n)) {
        if (model$outer_cv_results[[n]]$model$ncomp > max_ncomp) {
            max_ncomp = model$outer_cv_results[[n]]$model$ncomp
        }
    }
    
    loadings_sp <-
        matrix(nrow = ex_n * max_ncomp, ncol = ex_n * max_ncomp)
    labelsX <- list()
    for (n in seq_len(max_ncomp)) {
        for (m in seq_len(max_ncomp)) {
            for (i in seq_len(ex_n)) {
                for (j in seq_len(ex_n)) {
                    tryCatch({
                        loadings_sp[(n - 1) * ex_n + i, (m - 1) * ex_n + j] <-
                            (
                                model$outer_cv_results[[i]]$model$loadings$X[, n] %*% model$outer_cv_results[[j]]$model$loadings$X[, m]
                            )[1, 1]
                        labelsX[[(n - 1) * ex_n + i]] <-
                            paste("M", i, "- LV", n)
                    }, error = function(e) {
                    })
                }
            }
        }
    }
    #deleting na cows and cols
    loadings_sp <-
        loadings_sp[, colSums(is.na(loadings_sp)) != nrow(loadings_sp)]
    loadings_sp <-
        loadings_sp[rowSums(is.na(loadings_sp)) != ncol(loadings_sp), ]
    #rotate results
    loadings_sp <- t(loadings_sp[nrow(loadings_sp):1, , drop = FALSE])
    melted_loadings_sp <- reshape2::melt(loadings_sp)
    
    ggplot2::ggplot(melted_loadings_sp, ggplot2::aes(x = .data[["Var1"]], y = .data[["Var2"]], fill =
                                                         .data[["value"]])) +
        ggplot2::geom_tile() + ggplot2::coord_equal() + ggplot2::theme_bw() +
        ggplot2::scale_fill_distiller(palette = "Blues",
                                      direction = 1,
                                      na.value = "white") +
        ggplot2::guides(fill = FALSE) + # removing legend for `fill`
        ggplot2::labs(title = "Stability among models") + # using a title instead
        ggplot2::geom_text(
            ggplot2::aes(
                label = round(.data[["value"]], digits = 2),
                color = ifelse(.data[["value"]] > 0.5, 1, 0)
            ),
            size = 2.5,
            show.legend = FALSE
        ) +
        ggplot2::scale_x_discrete(limits = unlist(labelsX)) +
        ggplot2::scale_y_discrete(limits = rev(unlist(labelsX))) +
        ggplot2::labs(x = "Model - Latent variable", 
                      y = "Model - Latent variable") +
        ggplot2::theme(axis.text.x = ggplot2::element_text(
            angle = 45,
            vjust = 1,
            hjust = 1
        ),
        text = ggplot2::element_text(size = 10))
}

#' Models stability plot
#'
#' Plot stability among models of the external cross validation
#'
#' @param bp_results bp_kfold_VIP_analysis results
#' 
#' @return A plot of models stability
#' @name models_stability_plot_bootstrap
#' @export
#' @examples
#'# Data analysis for a table of integrated peaks
#' 
#' ## Generate an artificial nmr_dataset_peak_table:
#' ### Generate artificial metadata:
#' num_samples <- 64 # use an even number in this example
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
#' ## We will use bootstrap and permutation method for VIPs selection 
#' ## in a a k-fold cross validation 
#' #bp_results <- bp_kfold_VIP_analysis(peak_table, # Data to be analized
#' #                           y_column = "Condition", # Label
#' #                           k = 3,
#' #                           nbootstrap = 10)
#'
#' #message("Selected VIPs are: ", bp_results$importarn_vips)
#' 
#' #models_stability_plot_bootstrap(bp_results)
#' 
models_stability_plot_bootstrap = function (bp_results)
{
    # Loadings of the models
    n_models = length(bp_results$kfold_results)
    max_ncomp = bp_results$kfold_results[[1]]$general_model$ncomp
    
    loadings_sp <-
        matrix(nrow = n_models * max_ncomp, ncol = n_models * max_ncomp)
    labelsX <- list()
    for (n in seq_len(max_ncomp)) {
        for (m in seq_len(max_ncomp)) {
            for (i in seq_len(n_models)) {
                for (j in seq_len(n_models)) {
                    tryCatch({
                        loadings_sp[(n - 1) * n_models + i, (m - 1) * n_models + j] <-
                            (
                                bp_results$kfold_results[[i]]$general_model$loadings$X[, n] %*% bp_results$kfold_results[[j]]$general_model$loadings$X[, m]
                            )[1, 1]
                        labelsX[[(n - 1) * n_models + i]] <-
                            paste("M", i, "- LV", n)
                    }, error = function(e) {
                    })
                }
            }
        }
    }
    #deleting na cows and cols
    loadings_sp <-
        loadings_sp[, colSums(is.na(loadings_sp)) != nrow(loadings_sp)]
    loadings_sp <-
        loadings_sp[rowSums(is.na(loadings_sp)) != ncol(loadings_sp), ]
    #rotate results
    loadings_sp <- t(loadings_sp[nrow(loadings_sp):1, , drop = FALSE])
    melted_loadings_sp <- reshape2::melt(loadings_sp)
    
    ggplot2::ggplot(melted_loadings_sp, ggplot2::aes(x = .data[["Var1"]], y = .data[["Var2"]], fill =
                                                         .data[["value"]])) +
        ggplot2::geom_tile() + ggplot2::coord_equal() + ggplot2::theme_bw() +
        ggplot2::scale_fill_distiller(palette = "Blues",
                                      direction = 1,
                                      na.value = "white") +
        ggplot2::guides(fill = FALSE) + # removing legend for `fill`
        ggplot2::labs(title = "Stability among models") + # using a title instead
        ggplot2::geom_text(
            ggplot2::aes(
                label = round(.data[["value"]], digits = 2),
                color = ifelse(.data[["value"]] > 0.5, 1, 0)
            ),
            size = 2.5,
            show.legend = FALSE
        ) +
        ggplot2::scale_x_discrete(limits = unlist(labelsX)) +
        ggplot2::scale_y_discrete(limits = rev(unlist(labelsX))) +
        ggplot2::labs(x = "Model - Latent variable", 
                      y = "Model - Latent variable") +
        ggplot2::theme(axis.text.x = ggplot2::element_text(
            angle = 45,
            vjust = 1,
            hjust = 1
        ),
        text = ggplot2::element_text(size = 10))
}

#' Bootstrap plot predictions
#'
#' @param bp_results bp_kfold_VIP_analysis results
#' @param dataset An [nmr_dataset_family] object
#' @param y_column A string with the name of the y column (present in the
#' metadata of the dataset)
#' @param plot A boolean that indicate if results are plotted or not
#'
#' @name plot_bootstrap_multimodel
#' @return A plot of the results or a ggplot object
#' @importFrom stats predict
#' @importFrom mixOmics mixOmics
#' @export
#' @examples
#' # Data analysis for a table of integrated peaks
#' 
#' ## Generate an artificial nmr_dataset_peak_table:
#' ### Generate artificial metadata:
#' num_samples <- 64 # use an even number in this example
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
#' ## We will use bootstrap and permutation method for VIPs selection 
#' ## in a a k-fold cross validation 
#' #bp_results <- bp_kfold_VIP_analysis(peak_table, # Data to be analized
#' #                           y_column = "Condition", # Label
#' #                           k = 3,
#' #                           nbootstrap = 10)
#'
#' #message("Selected VIPs are: ", bp_results$importarn_vips)
#' 
#' #plot_bootstrap_multimodel(bp_results, peak_table, "Condition")
#' 
plot_bootstrap_multimodel <- function(bp_results, dataset, y_column, plot = TRUE) {
    
    n_models = length(bp_results$kfold_results)
    ncomp = bp_results$kfold_results[[1]]$general_model$ncomp
    # Extract data and split for train and test
    x_all <- dataset$peak_table
    y_all <- nmr_meta_get_column(dataset, column = y_column)
    te_data <- data.frame()
    for(i in seq_len(n_models)){
        # Predictions of test set
        predictions <- predict(bp_results$kfold_results[[i]]$general_model, newdata = x_all[-bp_results$kfold_index[[i]],, drop = FALSE])
        # Individuals plot
        if(ncomp == 1){
            te_data <- rbind(te_data, data.frame(x = predictions$variates[,1],
                                                 label = paste("test ", y_all[-bp_results$kfold_index[[i]]])))
        } else {
            te_y <- predictions$variates[, 2]
            te_data <- rbind(te_data, data.frame(x = predictions$variates[,1],
                                                 y = te_y,
                                                 label= y_all[-bp_results$kfold_index[[i]]],
                                                 group = "test "))
        }
    }
    
    # Individuals plot
    if(ncomp == 1){
        # This is needed if the model only have one component
        plsda_plot <- ggplot2::ggplot(data = te_data, ggplot2::aes(.data[["x"]], fill = .data[["label"]])) +
            ggplot2::geom_histogram(alpha = .5, bins = 10,
                                    position="identity") +
            ggplot2::ggtitle("PLS-DA") +
            ggplot2::labs(x = "Latent variable 1") +
            ggplot2::theme_bw()
    } else {
        plsda_plot <- ggplot2::ggplot(data = te_data,
                                      ggplot2::aes(shape = .data[["group"]],
                                                   col = .data[["label"]]
                                      )) +
            ggplot2::geom_hline(yintercept=0, linetype="dashed", 
                                color = "black", size=0.5) +
            ggplot2::geom_vline(xintercept=0, linetype="dashed", 
                                color = "black", size=0.5) +
            ggplot2::geom_point(ggplot2::aes(.data[["x"]], .data[["y"]]), size = 1.5) +
            ggplot2::ggtitle("PLS-DA") +
            ggplot2::labs(y = "Latent variable 2",
                          x = "Latent variable 1") +
            ggplot2::theme_bw()  
    }
    if(plot){
        plsda_plot
    } else {
        return(plsda_plot)
    }
}
