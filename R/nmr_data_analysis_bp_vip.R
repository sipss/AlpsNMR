#' Summarize bootstrap-permutation VIP differences into a per-feature CI and importance call
#'
#' Implements Afanador, Tran & Buydens (2013)'s normalization and importance
#' guidelines (Section 2.4, step 6): each feature's difference vector
#' (bootstrapped VIP minus permuted-bootstrapped VIP, across replicates) is
#' standardized to its own standard deviation, and the mean and 95%
#' confidence interval of those standardized differences are used to classify
#' the feature (lower-bound > `t_{1-alpha/2,n-1}` = important; > 0 = merely
#' relevant/marginally important; otherwise not important).
#'
#' `pls_vip_score_diff`'s replicates can be a single fold's bootstrap
#' iterations, or several folds' bootstrap iterations pooled together
#' (`bp_kfold_VIP_analysis()` does the latter, combining every fold's raw
#' differences before this single thresholding step, instead of thresholding
#' within each fold and then intersecting the resulting feature sets).
#'
#' @param pls_vip_score_diff A features x replicates matrix (`dimnames` gives
#'   feature names in rows) of `pls_vip - pls_vip_perm` differences.
#' @param n_samples Number of training samples backing the difference vector,
#'   used for the `t_{1-alpha/2,n-1}` quantile (n = training samples, not the
#'   number of replicates, per the paper).
#' @return A list with `boots_vip`, `boots_vip_sd`, `error`, `lower_bound`,
#'   `upper_bound` (each a one-column matrix over features, `dimnames` giving
#'   feature names), and `important_vips`/`relevant_vips` (character vectors
#'   of feature names).
#' @noRd
summarize_vip_differences <- function(pls_vip_score_diff, n_samples) {
    names <- rownames(pls_vip_score_diff)
    num_features <- nrow(pls_vip_score_diff)
    n_replicates <- ncol(pls_vip_score_diff)

    boots_vip <- matrix(nrow = num_features, dimnames = list(names))
    boots_vip_sd <- matrix(nrow = num_features, dimnames = list(names))
    error <- matrix(nrow = num_features, dimnames = list(names))
    lower_bound <- matrix(nrow = num_features, dimnames = list(names))
    upper_bound <- matrix(nrow = num_features, dimnames = list(names))
    for (k in seq_len(num_features)) {
        element <- pls_vip_score_diff[k, ] / sd(pls_vip_score_diff[k, ])
        boots_vip[k] <- sum(element) / n_replicates
        boots_vip_sd[k] <- sqrt(sum((element - boots_vip[k])^2) / (n_replicates - 1))
        error[k] <- qt(0.975, df = n_samples - 1) * boots_vip_sd[k]
        lower_bound[k] <- boots_vip[k] - error[k]
        upper_bound[k] <- boots_vip[k] + error[k]
    }

    # The "important" threshold is the same t_{1-alpha/2, n-1} quantile used
    # to build the confidence interval above (n = number of training
    # samples, not the number of replicates): a variable is "important" when
    # its entire two-sided (1-alpha) confidence interval lies above that
    # quantile, and merely "marginally important" (relevant_vips) when its
    # lower bound clears zero.
    important_vips <- names[lower_bound > qt(0.975, df = n_samples - 1)]
    relevant_vips <- names[lower_bound > 0]

    list(
        boots_vip = boots_vip,
        boots_vip_sd = boots_vip_sd,
        error = error,
        lower_bound = lower_bound,
        upper_bound = upper_bound,
        important_vips = important_vips,
        relevant_vips = relevant_vips
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
#' @param identity_column `NULL` or a string with the name of the identity column
#'    (present in the metadata of the dataset). When given, bootstrap resamples
#'    are drawn by resampling whole `identity_column` groups (e.g. subjects)
#'    with replacement, rather than individual rows, so that every repeated
#'    measurement of a resampled subject is kept together, and a multilevel
#'    (repeated-measures) `plsda` model is fitted (see [mixOmics::plsda]'s
#'    `multilevel` argument), matching what [nmr_data_analysis()] does for the
#'    same `identity_column`.
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
#' bp_VIPS <- bp_VIP_analysis(peak_table, # Data to be analyzed
#'     train_index,
#'     y_column = "Condition",
#'     ncomp = ncomps,
#'     nbootstrap = 10
#' )
#'
#' @export
bp_VIP_analysis <- function(dataset,
    train_index,
    y_column,
    identity_column = NULL,
    ncomp,
    nbootstrap = 300) {

    # Extract data and split for train and test
    x_all <- dataset$peak_table
    y_all <- nmr_meta_get_column(dataset, column = y_column)
    x_train <- x_all[train_index, , drop = FALSE]
    y_train <- y_all[train_index]
    n_samples <- nrow(x_train)
    # For check performance
    x_test <- x_all[-train_index, , drop = FALSE]
    y_test <- y_all[-train_index]

    if (!is.null(identity_column)) {
        identity_all <- nmr_meta_get_column(dataset, column = identity_column)
        identity_train <- identity_all[train_index]
    } else {
        identity_train <- NULL
    }

    if (length(unique(y_train)) == 1) {
        stop("Only one class in train set, please increase number of samples")
    }
    if (length(unique(y_test)) == 1) {
        stop("Only one class in test set, please increase number of samples")
    }

    num_features <- ncol(x_all)
    
    if (ncomp > num_features) {
        cli::cli_abort("ncomp ({ncomp}) can't be larger than num_features ({num_features})")
    }
    names <- colnames(x_all)
    # some checks
    if (length(names) == 0) {
        stop("Error in bp_VIP_analysis, the dataset peak_table doesn't have colnames.")
    }

    pls_vip <- matrix(nrow = num_features, ncol = nbootstrap, dimnames = list(names, NULL))
    pls_vip_perm_score <- matrix(nrow = num_features, ncol = nbootstrap, dimnames = list(names, NULL))
    pls_vip_score_diff <- matrix(nrow = num_features, ncol = nbootstrap, dimnames = list(names, NULL)) # Bootstrap with replacement nbootstraps datasets
    pls_models <- list()
    pls_perm_models <- list()
    CR <- list()

    # Bootstrap with replacement nbootstraps datasets
    res <- BiocParallel::bplapply(
        seq_len(nbootstrap),
        function(i, x_train, y_train, identity_train, ncomp) {
            num_features <- ncol(x_train)
            # A bootstrap resample can happen to draw only one class (plsda
            # models require at least two); when that happens, redraw rather
            # than patching a single element of the offending resample, which
            # would bias every degenerate resample towards the same fixed
            # (first-encountered) replacement sample instead of leaving the
            # resample an unbiased draw with replacement.
            if (is.null(identity_train)) {
                repeat {
                    index <- sample(seq_len(nrow(x_train)), nrow(x_train), replace = TRUE)
                    x_train_boots <- x_train[index, ]
                    y_train_boots <- y_train[index]
                    if (length(unique(y_train_boots)) > 1) {
                        break
                    }
                }
                identity_train_boots <- NULL
            } else {
                # Resample whole identity_train groups (e.g. subjects) with
                # replacement, instead of individual rows, so that every
                # repeated measurement of a resampled subject stays together
                # (required for a multilevel plsda model).
                groups <- unique(identity_train)
                repeat {
                    drawn_groups <- sample(groups, length(groups), replace = TRUE)
                    index <- unlist(
                        lapply(drawn_groups, function(g) which(identity_train == g)),
                        use.names = FALSE
                    )
                    x_train_boots <- x_train[index, ]
                    y_train_boots <- y_train[index]
                    if (length(unique(y_train_boots)) > 1) {
                        break
                    }
                }
                identity_train_boots <- identity_train[index]
            }

            # Rename rownames, because if they are repeated, plsda fails
            rownames(x_train_boots) <- paste0("Sample", seq_len(nrow(x_train_boots)))

            # Fit PLS model
            model <-
                plsda_build(
                    x = x_train_boots,
                    y = y_train_boots,
                    identity = identity_train_boots,
                    ncomp = ncomp
                )
            # VIPs per component extraction. mixOmics::vip() already returns,
            # in column h, the cumulative VIP (Eq. 9 of Afanador et al. 2013)
            # computed over components 1..h, so the last column (h = ncomp)
            # is the VIP of the fitted model. Indexed by the model's actual
            # fitted ncomp (model$ncomp), not the requested ncomp: mixOmics can
            # silently fit fewer components than requested on degenerate data
            # (e.g. a resample with a near-constant feature), in which case
            # plsda_vip()'s matrix has fewer than `ncomp` columns.
            pls_vip_comps <- plsda_vip(model)
            pls_vip <- pls_vip_comps[, model$ncomp]
            # Measure the classification rate (CR) of the bootstrap model
            CR <- get_test_accuracy(model, x_test, y_test)

            # bootsrapped and randomly permuted PLS-VIPs: for each feature j,
            # take its own VIP from the model fit with feature j (and only
            # feature j) permuted, rather than averaging feature j's VIP
            # across all num_features permuted-feature models.
            pls_vip_perm_score <- stats::setNames(numeric(num_features), names)
            for (j in seq_len(num_features)) {
                x_train_boots_perm <- x_train_boots
                # Shuffle column j's own values (breaks its row-order
                # association with the outcome) to obtain the permutation
                # null distribution for feature j's importance.
                x_train_boots_perm[, j] <- sample(x_train_boots[, j])

                # Refit model with permuted variables
                model_perm <-
                    plsda_build(
                        x = x_train_boots_perm,
                        y = y_train_boots,
                        identity = identity_train_boots,
                        ncomp = ncomp
                    )
                # VIPs per component extraction (see note above: take the
                # cumulative VIP through component ncomp, not a re-aggregation
                # across components, and index by the model's actual fitted
                # ncomp).
                pls_vip_comps_perm <- plsda_vip(model_perm)
                pls_vip_perm_score[j] <- pls_vip_comps_perm[j, model_perm$ncomp]
            }

            # bootsrapped and randomly permuted difference
            pls_vip_score_diff <- pls_vip - pls_vip_perm_score

            list(
                pls_vip = pls_vip,
                pls_vip_perm_score = pls_vip_perm_score,
                pls_vip_score_diff = pls_vip_score_diff,
                model = model,
                model_perm = model_perm,
                CR = CR
            )
        },
        x_train = x_train,
        y_train = y_train,
        identity_train = identity_train,
        ncomp = ncomp
    )
    pls_vip <- do.call(cbind, purrr::map(res, "pls_vip"))
    rownames(pls_vip) <- names
    pls_vip_perm_score <- do.call(cbind, purrr::map(res, "pls_vip_perm_score"))
    rownames(pls_vip_perm_score) <- names
    pls_vip_score_diff <- do.call(cbind, purrr::map(res, "pls_vip_score_diff"))
    rownames(pls_vip_score_diff) <- names
    pls_models <- purrr::map(res, "model")
    pls_perm_models <- purrr::map(res, "model_perm")
    CR <- purrr::map(res, "CR")

    # Normalization of the difference vector for each variable to its
    # corresponding standard deviation, construction of 95% confidence
    # intervals around the differences, and the importance guidelines
    # (Afanador et al. 2013, Section 2.4 step 6). See summarize_vip_differences().
    vip_summary <- summarize_vip_differences(pls_vip_score_diff, n_samples = n_samples)
    boots_vip <- vip_summary$boots_vip
    error <- vip_summary$error
    lower_bound <- vip_summary$lower_bound
    upper_bound <- vip_summary$upper_bound
    important_vips <- vip_summary$important_vips
    relevant_vips <- vip_summary$relevant_vips

    # Checking performance
    # Fit PLS model
    general_model <-
        plsda_build(
            x = x_train,
            y = y_train,
            identity = identity_train,
            ncomp = ncomp
        )
    # Measure the classification rate (CR) of the fold
    general_CR <- get_test_accuracy(general_model, x_test, y_test)

    vips_results <- train_models_with_only_vip_features(
        x_train, y_train, x_test, y_test,
        ncomp, important_vips, relevant_vips,
        identity_train = identity_train
    )
    vips_model <- vips_results$vips_model
    vips_CR <- vips_results$vips_CR
    
    # To return it ordered by mean of the normalized vectors
    orden <- order(boots_vip, decreasing = TRUE)
    # Return important vips and auc performance
    list(
        important_vips = important_vips,
        relevant_vips = relevant_vips,
        pls_vip = pls_vip[orden, , drop = FALSE],
        pls_vip_perm = pls_vip_perm_score[orden, , drop = FALSE],
        pls_vip_means = boots_vip[orden, , drop = FALSE],
        pls_vip_score_diff = pls_vip_score_diff[orden, , drop = FALSE],
        pls_models = pls_models,
        pls_perm_models = pls_perm_models,
        classif_rate = CR,
        general_model = general_model,
        general_CR = general_CR,
        vips_model = vips_model,
        vips_CR = vips_CR,
        error = error[orden, , drop = FALSE],
        lower_bound = lower_bound[orden, , drop = FALSE],
        upper_bound = upper_bound[orden, , drop = FALSE]
    )
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
#' `important_vips` and `relevant_vips` are derived by pooling every fold's
#' raw bootstrap-permutation VIP differences into a single collection of
#' replicates and applying Afanador, Tran & Buydens (2013)'s own mean/SD/CI
#' procedure to that pooled collection exactly once (see
#' [bp_VIP_analysis()]), rather than by classifying each fold on its own and
#' intersecting the resulting feature sets: a single fold's estimate falling
#' just short of the "important" cut-off no longer veto's a feature that
#' every other fold flagged as important, since the folds' evidence is
#' combined before any thresholding happens.
#'
#' @name bp_kfold_VIP_analysis
#' @param dataset An [nmr_dataset_family] object
#' @param y_column A string with the name of the y column (present in the
#'    metadata of the dataset)
#' @param identity_column `NULL` or a string with the name of the identity column
#'    (present in the metadata of the dataset). When given, whole
#'    `identity_column` groups (e.g. subjects) are assigned to the same fold
#'    together, instead of assigning individual samples to folds, and the
#'    bootstrap models fitted within each fold use a multilevel `plsda` (see
#'    [bp_VIP_analysis()]), matching what [nmr_data_analysis()] does for the
#'    same `identity_column`.
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
#' set.seed(42)
#' ## Generate an artificial nmr_dataset_peak_table:
#' ### Generate artificial metadata:
#' num_samples <- 64 # use an even number in this example
#' num_peaks <- 10
#' metadata <- data.frame(
#'     NMRExperiment = as.character(1:num_samples),
#'     Condition = sample(rep(c("A", "B"), times = num_samples / 2), num_samples)
#' )
#'
#' ### The matrix with peaks
#' peak_means <- runif(n = num_peaks, min = 300, max = 600)
#' peak_sd <- runif(n = num_peaks, min = 30, max = 60)
#' peak_matrix <- mapply(function(mu, sd) rnorm(num_samples, mu, sd),
#'     mu = peak_means, sd = peak_sd
#' )
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
#' bp_results <- bp_kfold_VIP_analysis(peak_table, # Data to be analyzed
#'     y_column = "Condition", # Label
#'     k = 2,
#'     ncomp = 1,
#'     nbootstrap = 5
#' )
#'
#' message("Selected VIPs are: ", bp_results$important_vips)
#'
bp_kfold_VIP_analysis <- function(dataset,
    y_column,
    identity_column = NULL,
    k = 4,
    ncomp = 3,
    nbootstrap = 300) {
    if (k <= 1) {
        stop("K must be integer greater than 1")
    }

    # Extract data and split for train and test
    y_all <- nmr_meta_get_column(dataset, column = y_column)
    if (length(unique(y_all)) == 1) {
        stop("Only one class in data set, at least two needed")
    }

    # Randomly assign each sample to one of the k folds. `rep_len()` spreads the
    # fold labels as evenly as possible and `sample()` shuffles them, so the
    # partition is random and (nearly) balanced. The indices refer to the
    # original dataset ordering, because bp_VIP_analysis() below reads the
    # samples straight from `dataset` using these indices. Each fold in turn is
    # held out as the test set, so k_fold_index[[i]] is its training set (every
    # sample not in fold i).
    n_all <- length(y_all)
    if (is.null(identity_column)) {
        fold_of_sample <- sample(rep_len(seq_len(k), n_all))
    } else {
        # Assign whole identity_column groups (e.g. subjects) to folds
        # together, so that no group is ever split across a fold's train and
        # test sides.
        identity_all <- as.factor(nmr_meta_get_column(dataset, column = identity_column))
        groups <- levels(identity_all)
        fold_of_group <- sample(rep_len(seq_len(k), length(groups)))
        names(fold_of_group) <- groups
        fold_of_sample <- unname(fold_of_group[as.character(identity_all)])
    }
    k_fold_index <- list()
    for (i in seq_len(k)) {
        k_fold_index[[i]] <- which(fold_of_sample != i)
    }

    # bp_VIP_analysis is already parallellized.
    results <- lapply(
        k_fold_index, function(index, dataset = dataset, y_column = y_column,
    identity_column = identity_column, ncomp = ncomp, nbootstrap = nbootstrap) {
            bp_VIP_analysis(
                dataset,
                index,
                y_column = y_column,
                identity_column = identity_column,
                ncomp = ncomp,
                nbootstrap = nbootstrap
            )
        },
        dataset = dataset, y_column = y_column, identity_column = identity_column,
        ncomp = ncomp, nbootstrap = nbootstrap
    )

    # Pool every fold's raw bootstrap-permutation VIP differences into a
    # single features x (k * nbootstrap) matrix, then apply Afanador et al.
    # (2013)'s own mean/SD/CI procedure (summarize_vip_differences()) exactly
    # once, on the pooled replicates, instead of thresholding within each
    # fold first and then intersecting the resulting important/relevant
    # feature sets. With an intersection, a single fold falling just short of
    # the "important" bar veto'd a feature entirely, even if every other
    # fold cleared it; pooling combines the evidence before any thresholding
    # happens, matching how the paper itself combines all of its bootstrap
    # replicates before applying a single threshold.
    #
    # bp_VIP_analysis() returns each fold's pls_vip_score_diff reordered by
    # that fold's own importance ranking, so rows must be realigned by
    # feature name (not position) before pooling -- row i of one fold's
    # matrix is not necessarily the same feature as row i of another's.
    names_order <- rownames(results[[1]]$pls_vip_score_diff)
    pooled_diff <- do.call(cbind, lapply(results, function(r) {
        r$pls_vip_score_diff[names_order, , drop = FALSE]
    }))

    # `n_samples` for the paper's t_{1-alpha/2,n-1} quantile is the training
    # sample size, not the number of (pooled) bootstrap replicates. Folds
    # have (nearly) equal training sizes by construction (`rep_len()` above
    # spreads samples/groups across folds as evenly as possible), so the
    # first fold's size is used as a representative value.
    n_samples_fold1 <- length(k_fold_index[[1]])
    pooled_summary <- summarize_vip_differences(pooled_diff, n_samples = n_samples_fold1)
    boots_vip_vec <- stats::setNames(as.numeric(pooled_summary$boots_vip), rownames(pooled_summary$boots_vip))
    error_vec <- stats::setNames(as.numeric(pooled_summary$error), rownames(pooled_summary$error))
    vip_order <- order(boots_vip_vec, decreasing = TRUE)
    vip_means <- boots_vip_vec[vip_order]
    error <- error_vec[vip_order]
    important_vips <- pooled_summary$important_vips
    relevant_vips <- pooled_summary$relevant_vips

    ## Wilcoxon test
    num_var <- dim(results[[1]]$pls_vip)[1]
    wt <- matrix(nrow = k, ncol = num_var)
    wt_vips <- list()
    for (i in seq_len(k)) {
        for (j in seq_len(num_var)) {
            x <- results[[i]]$pls_vip[j, ]
            y <- results[[i]]$pls_vip_perm[j, ]
            # wt_object <- wilcox.test(x, y, paired = TRUE, alternative = "two.sided")
            wt_object <- wilcox.test(x, y, paired = TRUE, alternative = "greater")
            wt[i, j] <- wt_object$p.value
        }
        wt_vips[[i]] <- rownames(results[[i]]$pls_vip)[wt[i, ] < 0.05]
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
        ggplot2::geom_point(
            mapping = ggplot2::aes(x = x, y = vip_means),
            shape = 21,
            fill = "white"
        ) +
        ggplot2::geom_hline(yintercept = qt(0.975, df = n_samples_fold1 - 1)) +
        ggplot2::ggtitle("BP-VIP") +
        ggplot2::labs(x = "Variables", y = "Scores") +
        ggplot2::theme_bw()

    list(
        important_vips = important_vips,
        relevant_vips = relevant_vips,
        wilcoxon_vips = unique(unlist(wt_vips)),
        vip_means = vip_means,
        vip_score_plot = p,
        kfold_results = results,
        kfold_index = k_fold_index
    )
}


#' Plot vip scores of bootstrap
#'
#' @param vip_means vips means values of bootstraps
#' @param error error tolerated, calculated in the bootstrap
#' @param n_samples number of training samples the bootstrap models were fit on.
#'    It sets the importance threshold line at `qt(0.975, df = n_samples - 1)`,
#'    following Afanador, Tran & Buydens (2013), where the cut-off quantile
#'    `t_{1-alpha/2,n-1}` uses n = number of samples (not the number of bootstraps).
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
#' ## We will use bootstrap and permutation method for VIPs selection
#' ## in a a k-fold cross validation
#' # bp_results <- bp_kfold_VIP_analysis(peak_table, # Data to be analized
#' #                           y_column = "Condition", # Label
#' #                           k = 3,
#' #                           ncomp = 1,
#' #                           nbootstrap = 10)
#'
#' # message("Selected VIPs are: ", bp_results$importarn_vips)
#'
#' # plot_vip_scores(bp_results$kfold_results[[1]]$vip_means,
#' #                bp_results$kfold_results[[1]]$error[1],
#' #                n_samples = length(bp_results$kfold_index[[1]]))
#'
plot_vip_scores <- function(vip_means, error, n_samples, plot = TRUE) {

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
        ggplot2::geom_point(
            mapping = ggplot2::aes(x = x, y = vip_means),
            shape = 21,
            fill = "white"
        ) +
        ggplot2::geom_hline(yintercept = qt(0.975, df = n_samples - 1)) +
        ggplot2::ggtitle("BP-VIP") +
        ggplot2::labs(x = "Variables", y = "Scores") +
        ggplot2::theme_bw()
    if (plot) {
        vip_score_plot
    } else {
        return(vip_score_plot)
    }
}
