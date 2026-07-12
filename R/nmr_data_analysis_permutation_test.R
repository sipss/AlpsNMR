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
#' p <- permutation_test_model(peak_table,
#'     y_column = "Condition",
#'     identity_column = NULL,
#'     external_val = list(iterations = 3, test_size = 0.25),
#'     internal_val = list(iterations = 3, test_size = 0.25),
#'     data_analysis_method = methodology,
#'     nPerm = 10
#' )
#'
permutation_test_model <- function(dataset, y_column, identity_column, external_val, internal_val,
    data_analysis_method, nPerm = 50) {
    permMatrix <- BiocParallel::bplapply(
        seq_len(nPerm),
        function(p, dataset, y_column, identity_column, external_val, internal_val, data_analysis_method) {
            dataset_perm <- dataset
            y_all <- nmr_meta_get_column(dataset, column = y_column)
            YPerm <- sample(y_all)
            dataset_perm[["metadata"]][["external"]][[y_column]] <- YPerm
            # print(sum(y_test!=y_all))
            permMod <- nmr_data_analysis(
                dataset_perm,
                y_column = y_column,
                identity_column = identity_column,
                external_val = external_val,
                internal_val = internal_val,
                data_analysis_method = data_analysis_method,
                .enable_parallel = FALSE
            )

            # I will use the mean of the auc of all the outer_cv for the test static
            test_stat <- mean(permMod$outer_cv_results_digested$auroc$auc)
            test_stat
        },
        dataset = dataset,
        y_column = y_column,
        identity_column = identity_column,
        external_val = external_val,
        internal_val = internal_val,
        data_analysis_method = data_analysis_method
    )
    permMatrix <- matrix(unlist(permMatrix), ncol = 1)
    return(permMatrix)
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
#' p <- permutation_test_model(peak_table,
#'     y_column = "Condition",
#'     identity_column = NULL,
#'     external_val = list(iterations = 3, test_size = 0.25),
#'     internal_val = list(iterations = 3, test_size = 0.25),
#'     data_analysis_method = methodology,
#'     nPerm = 10
#' )
#'
#' permutation_test_plot(model, p)
#'
permutation_test_plot <- function(nmr_data_analysis_model,
    permMatrix,
    xlab = "AUCs",
    xlim,
    ylim = NULL,
    breaks = "Sturges",
    main = "Permutation test") {
    h0 <- permMatrix[, 1]
    if (missing(xlim)) {
        xlim <- c(0, 1)
    }
    h <- hist(permMatrix, breaks, xlim = xlim, ylim = ylim, axes = FALSE, xlab = xlab, freq = FALSE, main = main)
    h2 <- max(h$density) * .75
    axis(1, pos = 0)
    axis(2, pos = 0, las = 1)

    model_auc <- mean(nmr_data_analysis_model$outer_cv_results_digested$auroc$auc)
    lines(rep(model_auc, 2), c(0, h2))

    p <- ecdf(h0)(model_auc) # Empirical
    # p1Stud=pt((model_auc-(mean(h0)/2))/sd(h0/2),(length(h0)-1)) # Students

    # warning: sometimes the auc of the permutation is NaN
    h0_median <- apply(permMatrix, 2, function(x) median(na.omit(x)))

    pP <- ifelse(model_auc < h0_median, p, 1 - p)
    if (pP < 1 / length(h0)) {
        text(h2, pos = 2, labels = paste("p<", signif(1 / length(h0), 4), sep = ""))
    } else {
        text(h2, pos = 2, labels = paste("p=", signif(pP, 4), sep = ""))
    }
}
