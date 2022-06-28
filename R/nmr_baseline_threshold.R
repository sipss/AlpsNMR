#' Threshold estimation for peak detection
#'
#' Estimates a robust baseline threshold value for peak detection on an [nmr_dataset_1D] object.
#'    This is performed computing the median and the mean absolute deviation of each spectrum
#'    in the 9.5-10 ppm range. The threshold is computed as the `median+3*mad`.
#'
#' @family peak detection functions
#' @param nmr_dataset An [nmr_dataset_1D].
#' @param range_without_peaks A vector with two doubles describing a range without peaks suitable for baseline detection
#' @return Numerical. A threshold value in intensity below that no peak is detected.
#' @export
#' @examples 
#' dir_to_demo_dataset <- system.file("dataset-demo", package = "AlpsNMR")
#' dataset <- nmr_read_samples_dir(dir_to_demo_dataset)
#' dataset_1D <- nmr_interpolate_1D(dataset, axis = c(min = -0.5, max = 10, by = 2.3E-4))
#' bl_threshold <- nmr_baseline_threshold(dataset_1D)
#' 
nmr_baseline_threshold <- function(nmr_dataset, range_without_peaks = c(9.5, 10)) {
    # FIXME: Maybe a whole baseline would be better, so we can cope with slowly changing baselines better
    if (length(range_without_peaks) != 2) {
        rlang::abort("range_without_peaks must have length 2")
    }
    r_start <- min(range_without_peaks)
    r_end <- max(range_without_peaks)
    threshold_ind <- nmr_dataset$axis >= r_start & nmr_dataset$axis < r_end
    out <- rep(NA_real_, nmr_dataset$num_samples)
    for (i in seq_len(nmr_dataset$num_samples)) {
        if ("data_1r_baseline" %in% names(unclass(nmr_dataset))) {
            spec_region <- nmr_dataset$data_1r[i, threshold_ind] - nmr_dataset$data_1r_baseline[i, threshold_ind]
        } else {
            spec_region <- nmr_dataset$data_1r[i, threshold_ind]
        }
        out[i] <- stats::median(spec_region) + 3*stats::mad(spec_region)
    }
    names(out) <- names(nmr_dataset)
    out
}


#' Plot the baseline thresholds
#'
#' If you have a lot of samples you can make the plot window bigger (or
#' use "` ```{r fig.height=10, fig.width=10}`" in notebooks), or choose some NMRExperiments.
#' 
#' @inheritParams plot.nmr_dataset_1D
#' @param nmr_dataset An [nmr_dataset_1D] object
#' @param thresholds A named vector. The values are baseline thresholds. The names are NMRExperiments.
#' @param NMRExperiment The NMRExperiments to plot (Use `"all"` to plot all of them)
#' @param chemshift_range The range to plot, as a first check use the `range_without_peaks` from [nmr_baseline_threshold]
#'
#' @return A plot.
#' @export
#' @examples 
#' dir_to_demo_dataset <- system.file("dataset-demo", package = "AlpsNMR")
#' dataset <- nmr_read_samples_dir(dir_to_demo_dataset)
#' dataset_1D <- nmr_interpolate_1D(dataset, axis = c(min = -0.5, max = 10, by = 2.3E-4))
#' baselineThresh <- nmr_baseline_threshold(dataset_1D)
#' nmr_baseline_threshold_plot(dataset_1D, baselineThresh)
nmr_baseline_threshold_plot <- function(nmr_dataset, thresholds, NMRExperiment = "all", chemshift_range = c(9.5, 10), ...) {
    if (is.null(NMRExperiment)) {
        if (nmr_dataset$num_samples > 20) {
            NMRExperiment <- sample(names(nmr_dataset), size = 10)
        } else {
            NMRExperiment <- names(nmr_dataset)
        }
    } else if (identical(NMRExperiment, "all")) {
        NMRExperiment <- names(nmr_dataset)
    }
    if (!identical(NMRExperiment, "all")) {
        thresholds <- thresholds[NMRExperiment]
    }
    
    aes_str <- as.character(list(...))
    columns_to_request <- c("NMRExperiment", get_vars_from_aes_string(aes_str))
    to_plot <- tidy(
        nmr_dataset,
        chemshift_range = chemshift_range,
        NMRExperiment = NMRExperiment,
        columns = columns_to_request,
        matrix_name = "data_1r"
    )
    
    if ("data_1r_baseline" %in% names(unclass(nmr_dataset))) {
        to_plot_threshold <- tidy(
            nmr_dataset,
            chemshift_range = chemshift_range,
            NMRExperiment = NMRExperiment,
            columns = columns_to_request,
            matrix_name = "data_1r_baseline"
        )
        to_plot_threshold <- dplyr::left_join(
            to_plot_threshold,
            tibble::enframe(
                thresholds,
                name = "NMRExperiment",
                value = "threshold"
            ),
            by = "NMRExperiment"
        )
        to_plot_threshold$intensity <- to_plot_threshold$intensity + to_plot_threshold$threshold
    } else {
        to_plot_threshold <-  dplyr::left_join(
            to_plot,
            tibble::enframe(
                thresholds,
                name = "NMRExperiment",
                value = "threshold"
            ),
            by = "NMRExperiment"
        )
        to_plot_threshold$intensity <- to_plot_threshold$threshold
    }
    
    dotdotdot_aes <- list(...)
    fixed_aes <- list(
        x = "chemshift",
        y = "intensity",
        group = "NMRExperiment"
    )
    all_aes <- c(fixed_aes, dotdotdot_aes)
    if (!"color" %in% names(all_aes) && !"colour" %in% names(all_aes)) {
        all_aes <- c(all_aes, list(color = "NMRExperiment"))
    }
    ymax <- 1.5*max(to_plot_threshold$intensity)

    gplt <- ggplot2::ggplot() +
        ggplot2::geom_line(mapping = do.call(ggplot2::aes_string, all_aes), data = to_plot) +
        ggplot2::geom_line(mapping = do.call(ggplot2::aes_string, all_aes), data = to_plot_threshold, linetype = "dashed", color = "black") +
        ggplot2::labs(x = "Chemical Shift (ppm)", y = "Intensity (a.u.)") +
        ggplot2::scale_x_reverse(limits = rev(chemshift_range[seq_len(2)])) +
        ggplot2::scale_y_continuous(labels = scales::label_number(scale_cut = scales::cut_si("")), limits = c(0, ymax)) +
        ggplot2::facet_wrap(~factor(NMRExperiment, levels = unique(NMRExperiment))) +
        ggplot2::theme(legend.position = "none")
    gplt
}