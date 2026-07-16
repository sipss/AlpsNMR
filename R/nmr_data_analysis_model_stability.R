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
#' # models_stability_plot_plsda(model)
#'
models_stability_plot_plsda <- function(model) {
    # Loadings of the models
    ex_n <- length(model$outer_cv_results)
    max_ncomp <- 0
    for (n in seq_len(ex_n)) {
        if (model$outer_cv_results[[n]]$model$ncomp > max_ncomp) {
            max_ncomp <- model$outer_cv_results[[n]]$model$ncomp
        }
    }

    loadings_sp <-
        matrix(nrow = ex_n * max_ncomp, ncol = ex_n * max_ncomp)
    labelsX <- list()
    for (n in seq_len(max_ncomp)) {
        for (m in seq_len(max_ncomp)) {
            for (i in seq_len(ex_n)) {
                for (j in seq_len(ex_n)) {
                    tryCatch(
                        {
                            loadings_sp[(n - 1) * ex_n + i, (m - 1) * ex_n + j] <-
                                (
                                    model$outer_cv_results[[i]]$model$loadings$X[, n] %*% model$outer_cv_results[[j]]$model$loadings$X[, m]
                                )[1, 1]
                            labelsX[[(n - 1) * ex_n + i]] <-
                                paste("M", i, "- LV", n)
                        },
                        error = function(e) {
                        }
                    )
                }
            }
        }
    }
    # deleting na cows and cols
    loadings_sp <-
        loadings_sp[, colSums(is.na(loadings_sp)) != nrow(loadings_sp)]
    loadings_sp <-
        loadings_sp[rowSums(is.na(loadings_sp)) != ncol(loadings_sp), ]
    # rotate results
    loadings_sp <- t(loadings_sp[nrow(loadings_sp):1, , drop = FALSE])
    melted_loadings_sp <- reshape2::melt(loadings_sp)

    ggplot2::ggplot(melted_loadings_sp, ggplot2::aes(
        x = .data[["Var1"]], y = .data[["Var2"]], fill =
            .data[["value"]]
    )) +
        ggplot2::geom_tile() +
        ggplot2::coord_equal() +
        ggplot2::theme_bw() +
        ggplot2::scale_fill_distiller(
            palette = "Blues",
            direction = 1,
            na.value = "white"
        ) +
        ggplot2::guides(fill = "none") + # removing legend for `fill`
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
        ggplot2::labs(
            x = "Model - Latent variable",
            y = "Model - Latent variable"
        ) +
        ggplot2::theme(
            axis.text.x = ggplot2::element_text(
                angle = 45,
                vjust = 1,
                hjust = 1
            ),
            text = ggplot2::element_text(size = 10)
        )
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
#' #                           nbootstrap = 10)
#'
#' # message("Selected VIPs are: ", bp_results$importarn_vips)
#'
#' # models_stability_plot_bootstrap(bp_results)
#'
models_stability_plot_bootstrap <- function(bp_results) {
    # Loadings of the models
    n_models <- length(bp_results$kfold_results)
    max_ncomp <- bp_results$kfold_results[[1]]$general_model$ncomp

    loadings_sp <-
        matrix(nrow = n_models * max_ncomp, ncol = n_models * max_ncomp)
    labelsX <- list()
    for (n in seq_len(max_ncomp)) {
        for (m in seq_len(max_ncomp)) {
            for (i in seq_len(n_models)) {
                for (j in seq_len(n_models)) {
                    tryCatch(
                        {
                            loadings_sp[(n - 1) * n_models + i, (m - 1) * n_models + j] <-
                                (
                                    bp_results$kfold_results[[i]]$general_model$loadings$X[, n] %*% bp_results$kfold_results[[j]]$general_model$loadings$X[, m]
                                )[1, 1]
                            labelsX[[(n - 1) * n_models + i]] <-
                                paste("M", i, "- LV", n)
                        },
                        error = function(e) {
                        }
                    )
                }
            }
        }
    }
    # deleting na cows and cols
    loadings_sp <-
        loadings_sp[, colSums(is.na(loadings_sp)) != nrow(loadings_sp)]
    loadings_sp <-
        loadings_sp[rowSums(is.na(loadings_sp)) != ncol(loadings_sp), ]
    # rotate results
    loadings_sp <- t(loadings_sp[nrow(loadings_sp):1, , drop = FALSE])
    melted_loadings_sp <- reshape2::melt(loadings_sp)

    ggplot2::ggplot(melted_loadings_sp, ggplot2::aes(
        x = .data[["Var1"]], y = .data[["Var2"]], fill =
            .data[["value"]]
    )) +
        ggplot2::geom_tile() +
        ggplot2::coord_equal() +
        ggplot2::theme_bw() +
        ggplot2::scale_fill_distiller(
            palette = "Blues",
            direction = 1,
            na.value = "white"
        ) +
        ggplot2::guides(fill = "none") + # removing legend for `fill`
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
        ggplot2::labs(
            x = "Model - Latent variable",
            y = "Model - Latent variable"
        ) +
        ggplot2::theme(
            axis.text.x = ggplot2::element_text(
                angle = 45,
                vjust = 1,
                hjust = 1
            ),
            text = ggplot2::element_text(size = 10)
        )
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
#' #                           nbootstrap = 10)
#'
#' # message("Selected VIPs are: ", bp_results$importarn_vips)
#'
#' # plot_bootstrap_multimodel(bp_results, peak_table, "Condition")
#'
plot_bootstrap_multimodel <- function(bp_results, dataset, y_column, plot = TRUE) {
    n_models <- length(bp_results$kfold_results)
    ncomp <- bp_results$kfold_results[[1]]$general_model$ncomp
    # Extract data and split for train and test
    x_all <- dataset$peak_table
    y_all <- nmr_meta_get_column(dataset, column = y_column)
    te_data <- data.frame()
    for (i in seq_len(n_models)) {
        # Predictions of test set
        predictions <- predict(bp_results$kfold_results[[i]]$general_model, newdata = x_all[-bp_results$kfold_index[[i]], , drop = FALSE])
        # Individuals plot
        if (ncomp == 1) {
            te_data <- rbind(te_data, data.frame(
                x = predictions$variates[, 1],
                label = paste("test ", y_all[-bp_results$kfold_index[[i]]])
            ))
        } else {
            te_y <- predictions$variates[, 2]
            te_data <- rbind(te_data, data.frame(
                x = predictions$variates[, 1],
                y = te_y,
                label = y_all[-bp_results$kfold_index[[i]]],
                group = "test "
            ))
        }
    }

    # Individuals plot
    if (ncomp == 1) {
        # This is needed if the model only have one component
        plsda_plot <- ggplot2::ggplot(data = te_data, ggplot2::aes(.data[["x"]], fill = .data[["label"]])) +
            ggplot2::geom_histogram(
                alpha = .5, bins = 10,
                position = "identity"
            ) +
            ggplot2::ggtitle("PLS-DA") +
            ggplot2::labs(x = "Latent variable 1") +
            ggplot2::theme_bw()
    } else {
        plsda_plot <- ggplot2::ggplot(
            data = te_data,
            ggplot2::aes(
                shape = .data[["group"]],
                col = .data[["label"]]
            )
        ) +
            ggplot2::geom_hline(
                yintercept = 0, linetype = "dashed",
                color = "black", linewidth = 0.5
            ) +
            ggplot2::geom_vline(
                xintercept = 0, linetype = "dashed",
                color = "black", linewidth = 0.5
            ) +
            ggplot2::geom_point(ggplot2::aes(.data[["x"]], .data[["y"]]), size = 1.5) +
            ggplot2::ggtitle("PLS-DA") +
            ggplot2::labs(
                y = "Latent variable 2",
                x = "Latent variable 1"
            ) +
            ggplot2::theme_bw()
    }
    if (plot) {
        plsda_plot
    } else {
        return(plsda_plot)
    }
}
