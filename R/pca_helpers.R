#' Build a PCA on for an nmr_dataset
#'
#' This function builds a PCA model with all the NMR spectra. Regions with
#' zero values (excluded regions) or near-zero variance regions are automatically
#' excluded from the analysis.
#'
#' @param nmr_dataset a [nmr_dataset_1D] object
#' @inheritParams mixOmics::pca
#' @param ... Additional arguments passed on to [mixOmics::pca]
#' @family PCA related functions
#' @return
#' A PCA model as given by [mixOmics::pca] with two additional attributes:
#'    - `nmr_data_axis` containing the full ppm axis
#'    - `nmr_included` with the data points included in the model
#' These attributes are used internally by AlpsNMR to create loading plots
#'
#' @export
#' @examples
#' dir_to_demo_dataset <- system.file("dataset-demo", package = "AlpsNMR")
#' dataset <- nmr_read_samples_dir(dir_to_demo_dataset)
#' dataset_1D <- nmr_interpolate_1D(dataset, axis = c(min = -0.5, max = 10, by = 2.3E-4))
#' model <- nmr_pca_build_model(dataset_1D)
#'
nmr_pca_build_model <- function(nmr_dataset,
                                ncomp = NULL,
                                center = TRUE,
                                scale = FALSE,
                                ...) {
    UseMethod("nmr_pca_build_model")
}


#' @rdname nmr_pca_build_model
#' @family nmr_dataset_1D functions
#' @export
#' @examples
#' dir_to_demo_dataset <- system.file("dataset-demo", package = "AlpsNMR")
#' dataset <- nmr_read_samples_dir(dir_to_demo_dataset)
#' dataset_1D <- nmr_interpolate_1D(dataset, axis = c(min = -0.5, max = 10, by = 2.3E-4))
#' model <- nmr_pca_build_model(dataset_1D)
#'
nmr_pca_build_model.nmr_dataset_1D <- function(nmr_dataset,
                                               ncomp = NULL,
                                               center = TRUE,
                                               scale = FALSE,
                                               ...) {
    data_1r <- nmr_dataset$data_1r
    rownames(data_1r) <-
        nmr_meta_get_column(nmr_dataset, column = "NMRExperiment")
    pca_model <-
        mixOmics::pca(
            X = data_1r,
            ncomp = ncomp,
            center = center,
            scale = scale,
            ...
        )
    # These attributes are used by nmr_pca_loadingplot:
    attr(pca_model, "nmr_data_axis") <- nmr_dataset$axis
    attr(pca_model, "nmr_included") <-
        seq_along(nmr_dataset$axis)
    pca_model
}

#' Plotting functions for PCA
#'
#' @param nmr_dataset an [nmr_dataset_1D] object
#' @param pca_model A PCA model trained with [nmr_pca_build_model]
#' @param comp Components to represent
#' @param ... Additional aesthetics passed on to [ggplot2::aes] (use bare unquoted names)
#'
#' @family PCA related functions
#' @name nmr_pca_plots
#' @return Plot of PCA
#' @examples
#' dir_to_demo_dataset <- system.file("dataset-demo", package = "AlpsNMR")
#' dataset <- nmr_read_samples_dir(dir_to_demo_dataset)
#' dataset_1D <- nmr_interpolate_1D(dataset, axis = c(min = -0.5, max = 10, by = 2.3E-4))
#' model <- nmr_pca_build_model(dataset_1D)
#' nmr_pca_plot_variance(model)
#' nmr_pca_scoreplot(dataset_1D, model)
#' nmr_pca_loadingplot(model, 1)
#'
NULL

#' @rdname nmr_pca_plots
#' @export
nmr_pca_plot_variance <- function(pca_model) {
    cum_var_percent <-
        100 * cumsum(pca_model$sdev ^ 2 / pca_model$var.tot)
    ggplot2::qplot(x = seq_along(cum_var_percent),
                   y = cum_var_percent,
                   geom = "line") +
        ggplot2::xlab("Number of Principal Components") +
        ggplot2::ylab("Cummulated Explained Variance (%)")
}

#' @rdname nmr_pca_plots
#' @export
nmr_pca_scoreplot <- function(nmr_dataset,
                              pca_model,
                              comp = seq_len(2), ...) {
    nmr_metadata <- nmr_meta_get(nmr_dataset)
    scores_df <- as.data.frame(pca_model$X)
    colnames(scores_df) <- paste0("PC", seq_len(ncol(scores_df)))
    scores <- dplyr::left_join(
        cbind(NMRExperiment = rownames(pca_model$X), scores_df),
        nmr_metadata,
        by = "NMRExperiment"
    )
    var_percent <- 100 * pca_model$sdev ^ 2 / pca_model$var.tot
    axis_labels <-
        paste0("PC",
               seq_along(var_percent),
               " (",
               round(var_percent, 2),
               "%)")
    names(axis_labels) <- paste0("PC", seq_along(var_percent))
    if (length(comp) == 2) {
        xy <- paste0("PC", comp)
        gplt <- ggplot2::ggplot(
            data = scores,
            mapping = ggplot2::aes(
                x = !!rlang::sym(xy[1]),
                y = !!rlang::sym(xy[2]),
                ...
            )
        ) +
            ggplot2::geom_point() +
            ggplot2::xlab(axis_labels[xy[1]]) +
            ggplot2::ylab(axis_labels[xy[2]])
    } else {
        if (!requireNamespace("GGally", quietly = TRUE)) {
            rlang::abort("Please install the GGally package to plot more than two components at once")
        }
        gplt <- GGally::ggpairs(
            data = scores,
            mapping = ggplot2::aes(...),
            columns = comp + 1,
            # +1 because the NMRExperiment becomes the first column
            progress = FALSE,
            labeller = ggplot2::as_labeller(axis_labels)
        )
    }
    gplt
}

#' @rdname nmr_pca_plots
#' @export
nmr_pca_loadingplot <- function(pca_model, comp) {
    ppm_axis <- attr(pca_model, "nmr_data_axis")
    loadings <-
        matrix(0, nrow = length(ppm_axis), ncol = length(comp))
    loadings[attr(pca_model, "nmr_included"),] <-
        pca_model$loadings$X[, comp, drop = FALSE]
    # loadings[,1] # first loading
    loadings <- as.data.frame(loadings)
    loadings$ppm <- ppm_axis
    loadings_long <-
        reshape2::melt(
            loadings,
            id.vars = "ppm",
            variable.name = "component",
            value.name = "loading"
        )
    ggplot2::ggplot(loadings_long,
                    ggplot2::aes_string(x = "ppm", y = "loading", group = "component")) +
        ggplot2::geom_line() +
        ggplot2::scale_x_reverse()
}

#' Compute PCA residuals and score distance for each sample
#'
#' @param nmr_dataset An [nmr_dataset_1D] object
#' @param pca_model A pca model returned by [nmr_pca_build_model]
#' @param ncomp Number of components to use. Use `NULL` for 90% of the variance
#' @param quantile_critical critical quantile
#'
#' @family PCA related functions
#' @family outlier detection functions
#' @family nmr_dataset_1D functions
#' @return
#'
#' A list with:
#'
#'    - outlier_info: A data frame with the NMRExperiment, the Q residuals and T scores
#'    - ncomp: Number of components used to compute Q and T
#'    - Tscore_critical, QResidual_critical: Critical values, given a quantile, for both Q and T.
#'
#' @export
#' @examples
#' dir_to_demo_dataset <- system.file("dataset-demo", package = "AlpsNMR")
#' dataset <- nmr_read_samples_dir(dir_to_demo_dataset)
#' dataset_1D <- nmr_interpolate_1D(dataset, axis = c(min = -0.5, max = 10, by = 2.3E-4))
#' model <- nmr_pca_build_model(dataset_1D)
#' outliers_info <- nmr_pca_outliers(dataset_1D, model)
#'
nmr_pca_outliers <- function(nmr_dataset,
                             pca_model,
                             ncomp = NULL,
                             quantile_critical = 0.975) {
    nmr_dataset <- validate_nmr_dataset_1D(nmr_dataset)
    
    if (is.null(ncomp)) {
        cum_var_percent <-
            100 * cumsum(pca_model$sdev ^ 2 / pca_model$var.tot)
        ncomp <- which(cum_var_percent > 90)[1]
    }
    
    # T scores:
    scores <-
        pca_model$variates$X[, seq_len(ncomp), drop = FALSE]
    variances <- utils::head(pca_model$sdev ^ 2, ncomp)
    loadings <-
        pca_model$loadings$X[, seq_len(ncomp), drop = FALSE]
    Tscore <-
        sqrt(apply(scores ^ 2 / rep(variances, each = nrow(scores)), 1, sum))
    
    # Q residuals
    Xs <-
        scale(nmr_dataset$data_1r,
              center = pca_model$center,
              scale = pca_model$scale)
    residuals <- Xs - scores %*% t(loadings)
    Qres <- sqrt(apply(residuals ^ 2, 1, sum))
    
    # compute critical values
    Tscore_critical <-
        sqrt(stats::qchisq(quantile_critical, ncomp))
    QResidual_critical <-
        (stats::median(Qres ^ (2 / 3)) + stats::mad(Qres ^ (2 / 3)) * stats::qnorm(quantile_critical)) ^
        (3 / 2)
    
    outlier_info <- nmr_dataset %>%
        nmr_meta_get(columns = "NMRExperiment") %>%
        dplyr::mutate(Tscores = Tscore,
                      QResiduals = Qres)
    list(
        outlier_info = outlier_info,
        ncomp = ncomp,
        Tscore_critical = Tscore_critical,
        QResidual_critical = QResidual_critical
    )
}


#' Outlier detection through robust PCA
#'
#' @param nmr_dataset An nmr_dataset_1D object
#' @param ncomp Number of rPCA components to use
#'
#' We have observed that the statistical test used as a threshold for
#' outlier detection usually flags as outliers too many samples, due possibly
#' to a lack of gaussianity
#'
#' As a workaround, a heuristic method has been implemented: We know that in the
#' Q residuals vs T scores plot from [nmr_pca_outliers_plot()] outliers are
#' on the right or on the top of the plot, and quite separated from non-outlier
#' samples.
#'
#' To determine the critical value, both for Q and T, we find the biggest gap
#' between samples in the plot and use as critical value the center of the gap.
#'
#' This approach seems to work well when there are outliers, but it fails when there
#' isn't any outlier. For that case, the gap would be placed anywhere and that is
#' not desirable as many samples would be incorrectly flagged. The
#' second assumption that we use is that no more than 10\% of
#' the samples may pass our critical value. If more than 10\% of the samples
#' pass the critical value, then we assume that our heuristics are not reasonable
#' and we don't set any critical limit.
#'
#'
#' @return A list similar to [nmr_pca_outliers]
#' @export
#' @family PCA related functions
#' @family outlier detection functions
#' @family nmr_dataset_1D functions
#' @examples
#' dir_to_demo_dataset <- system.file("dataset-demo", package = "AlpsNMR")
#' dataset <- nmr_read_samples_dir(dir_to_demo_dataset)
#' dataset_1D <- nmr_interpolate_1D(dataset, axis = c(min = -0.5, max = 10, by = 2.3E-4))
#' outliers_info <- nmr_pca_outliers_robust(dataset_1D)
#'
nmr_pca_outliers_robust <- function(nmr_dataset, ncomp = 5) {
    nmr_dataset <- validate_nmr_dataset_1D(nmr_dataset)
    
    Xprep_rob <- scale(
        nmr_dataset$data_1r,
        center = apply(nmr_dataset$data_1r, 2, function(x_col)
            stats::median(x_col, na.rm = TRUE)),
        scale = apply(nmr_dataset$data_1r, 2, function(x_col)
            stats::mad(x_col, na.rm = TRUE))
    )
    
    if (missing(ncomp) && nmr_dataset$num_samples <= ncomp) {
        ncomp <- nmr_dataset$num_samples - 1
    }
    
    pca_model <- pcaPP::PCAgrid(Xprep_rob, k = ncomp, scale = NULL)
    
    
    # T scores:
    scores <- pca_model$scores[, seq_len(ncomp), drop = FALSE]
    variances <- utils::head(pca_model$sdev ^ 2, ncomp)
    loadings <- pca_model$loadings[, seq_len(ncomp), drop = FALSE]
    Tscore <-
        sqrt(apply(scores ^ 2 / rep(variances, each = nrow(scores)), 1, sum))
    
    # Q residuals
    residuals <- Xprep_rob - scores %*% t(loadings)
    Qres <- sqrt(apply(residuals ^ 2, 1, sum))
    
    # compute critical values
    find_crit_thres <- function(x) {
        if (length(x) == 0) {
            return(numeric(0))
        }
        xsorted <- sort(x)
        crit_thresh_ind <- which.max(diff(xsorted))
        proposed_thresh <-
            (xsorted[crit_thresh_ind] + xsorted[crit_thresh_ind + 1]) / 2
        outliers_detected <- sum(x > proposed_thresh) / length(x)
        if (outliers_detected > 0.1) {
            return(Inf)
        } else {
            return(proposed_thresh)
        }
    }
    
    Tscore_critical <-    find_crit_thres(Tscore)
    QResidual_critical <- find_crit_thres(Qres)
    
    outlier_info <- nmr_dataset %>%
        nmr_meta_get(columns = "NMRExperiment") %>%
        dplyr::mutate(Tscores = Tscore,
                      QResiduals = Qres)
    
    list(
        outlier_info = outlier_info,
        ncomp = ncomp,
        Tscore_critical = Tscore_critical,
        QResidual_critical = QResidual_critical
    )
}


#' Plot for outlier detection diagnostic
#'
#' @param nmr_dataset An [nmr_dataset_1D] object
#' @param pca_outliers The output from [nmr_pca_outliers()]
#' @param ... Additional parameters passed on to [ggplot2::aes_string()]
#'
#' @return A plot for the outlier detection
#' @export
#'
#' @family PCA related functions
#' @family outlier detection functions
#' @family nmr_dataset_1D functions
#' @importFrom rlang .data
#' @examples
#' #dir_to_demo_dataset <- system.file("dataset-demo", package = "AlpsNMR")
#' #dataset <- nmr_read_samples_dir(dir_to_demo_dataset)
#' #dataset_1D <- nmr_interpolate_1D(dataset, axis = c(min = -0.5, max = 10, by = 2.3E-4))
#' #model <- nmr_pca_build_model(dataset_1D)
#' #outliers_info <- nmr_pca_outliers(dataset_1D, model)
#' #nmr_pca_outliers_plot(dataset_1D, outliers_info)
#' 
nmr_pca_outliers_plot <- function(nmr_dataset, pca_outliers, ...) {
    has_ggrepel <- requireNamespace("ggrepel", quietly = TRUE)
    if (!has_ggrepel) {
        rlang::warn(
            message = "Please install ggrepel to avoid overlap of text labels in the plot",
            .frequency = "once", 
            .frequency_id = "install_ggrepel"
        )
    }
    outlier_info <- pca_outliers[["outlier_info"]]
    tscore_crit <- pca_outliers[["Tscore_critical"]]
    qres_crit <- pca_outliers[["QResidual_critical"]]
    ncomp <- pca_outliers[["ncomp"]]
    
    pca_outliers_with_meta <-
        dplyr::left_join(nmr_meta_get(nmr_dataset, groups = "external"),
                         outlier_info,
                         by = "NMRExperiment")
    
    pca_outliers_with_meta_only_out <- dplyr::filter(
        pca_outliers_with_meta,
        .data$Tscores > tscore_crit | .data$QResiduals > qres_crit
    )
    
    gplt <- ggplot2::ggplot(
        pca_outliers_with_meta,
        ggplot2::aes_string(x = "Tscores", y = "QResiduals", label = "NMRExperiment")
    ) +
        ggplot2::geom_point(ggplot2::aes_string(...))
    if (has_ggrepel) {
        gplt <- gplt +
            ggrepel::geom_text_repel(data = pca_outliers_with_meta_only_out)
    } else {
        gplt <- gplt +
            ggplot2::geom_text(data = pca_outliers_with_meta_only_out)
    }
    gplt <- gplt +
        ggplot2::geom_vline(xintercept = tscore_crit,
                            colour = "red",
                            linetype = "dashed") +
        ggplot2::geom_hline(yintercept = qres_crit,
                            colour = "red",
                            linetype = "dashed") +
        ggplot2::ggtitle(glue::glue_data(
            list(ncomp = ncomp),
            "PCA Residuals and Score distance ({ncomp} components)"
        ))
    gplt
}

#' Exclude outliers
#'
#' @inheritParams nmr_pca_outliers_plot
#'
#' @return An [nmr_dataset_1D] without the detected outliers
#' @export
#' @family PCA related functions
#' @family outlier detection functions
#' @family nmr_dataset_1D functions
#' @family subsetting functions
#' @examples
#' dir_to_demo_dataset <- system.file("dataset-demo", package = "AlpsNMR")
#' dataset <- nmr_read_samples_dir(dir_to_demo_dataset)
#' dataset_1D <- nmr_interpolate_1D(dataset, axis = c(min = -0.5, max = 10, by = 2.3E-4))
#' model <- nmr_pca_build_model(dataset_1D)
#' outliers_info <- nmr_pca_outliers(dataset_1D, model)
#' dataset_whitout_outliers <- nmr_pca_outliers_filter(dataset_1D, outliers_info)
#' 
nmr_pca_outliers_filter <- function(nmr_dataset, pca_outliers) {
    outlier_info <- pca_outliers[["outlier_info"]]
    tscore_crit <- pca_outliers[["Tscore_critical"]]
    qres_crit <- pca_outliers[["QResidual_critical"]]
    
    nmrexp_to_keep <- outlier_info %>%
        dplyr::filter(.data$Tscores < tscore_crit &
                          .data$QResiduals < qres_crit) %>%
        dplyr::pull("NMRExperiment")
    
    dplyr::filter(nmr_dataset, .data$NMRExperiment %in% nmrexp_to_keep)
}
