#' Threshold estimation for peak detection
#'
#' Estimates the threshold value for peak detection on an [nmr_dataset_1D] object by examining
#' a range without peaks, by default the 9.5 - 10 ppm range.
#'
#' Two methods can be used:
#'
#' - "mean3sd": The mean3sd method computes the mean and the standard deviation of each spectrum
#' in the 9.5 - 10 ppm range. The mean spectrum and the mean standard deviation are both vectors
#' of length equal to the number of points in the given range. The mean of the mean spectrum
#  and the mean of the standard deviations are used to summarize the center and dispersion of
#' the noise. The threshold is defined as `center + 3 dispersion`, and it is one single threshold
#' for the whole dataset. This is the default for backwards compatibility.
#'
#' - "median3mad": First we take the data matrix. If we have estimated a baseline already,
#'   subtract it. In the defined region without peaks, estimate the median of each sample and
#'   its median absolute deviation. Return a vector of length equal to the number of samples
#'   with the `median+3mad` for each sample. This is a new more robust method.
#'
#' @family peak detection functions
#' @param nmr_dataset An [nmr_dataset_1D].
#' @param method Either "mean3sd" or the more robust "median3mad". See the details.
#' @param range_without_peaks A vector with two doubles describing a range without peaks suitable for baseline detection
#' @return Numerical. A threshold value in intensity below that no peak is detected.
#' @export
#' @examples
#' dir_to_demo_dataset <- system.file("dataset-demo", package = "AlpsNMR")
#' dataset <- nmr_read_samples_dir(dir_to_demo_dataset)
#' dataset_1D <- nmr_interpolate_1D(dataset, axis = c(min = -0.5, max = 10, by = 2.3E-4))
#' bl_threshold <- nmr_baseline_threshold(dataset_1D)
#'
nmr_baseline_threshold <- function(nmr_dataset, range_without_peaks = c(9.5, 10), method = c("mean3sd", "median3mad")) {
    # FIXME: Maybe a whole baseline would be better, so we can cope with slowly changing baselines better
    method <- match.arg(method)
    if (length(range_without_peaks) != 2) {
        rlang::abort("range_without_peaks must have length 2")
    }
    r_start <- min(range_without_peaks)
    r_end <- max(range_without_peaks)
    threshold_ind <- nmr_dataset$axis >= r_start & nmr_dataset$axis < r_end
    if (method == "mean3sd") {
        cent <- mean(apply(nmr_dataset$data_1r[, threshold_ind], 2, mean))
        disp <- mean(apply(nmr_dataset$data_1r[, threshold_ind], 2, stats::sd))
        return(cent + 3 * disp)
    } else if (method == "median3mad") {
        out <- rep(NA_real_, nmr_dataset$num_samples)
        for (i in seq_len(nmr_dataset$num_samples)) {
            if ("data_1r_baseline" %in% names(unclass(nmr_dataset))) {
                spec_region <- nmr_dataset$data_1r[i, threshold_ind] - nmr_dataset$data_1r_baseline[i, threshold_ind]
            } else {
                spec_region <- nmr_dataset$data_1r[i, threshold_ind]
            }
            out[i] <- stats::median(spec_region) + 3 * stats::mad(spec_region)
        }
        names(out) <- names(nmr_dataset)
        return(out)
    } else {
        stop("Unexpected method", method)
    }
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
    tidy_data <- tidy_spectra_baseline_and_threshold(
        dataset = nmr_dataset,
        thresholds = thresholds,
        chemshift_range = chemshift_range,
        NMRExperiment = NMRExperiment,
        columns = columns_to_request
    )
    to_plot <- tidy_data$spectra
    to_plot_baseline <- tidy_data$baselines
    to_plot_threshold <- tidy_data$thresholds


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
    ymax <- 1.5 * max(to_plot_threshold$intensity)

    gplt <- ggplot2::ggplot() +
        # The spectra:
        ggplot2::geom_line(mapping = do.call(ggplot2::aes_string, all_aes), data = to_plot)

    if (!is.null(to_plot_baseline)) {
        # The baseline:
        gplt <- gplt + ggplot2::geom_line(mapping = do.call(ggplot2::aes_string, all_aes), data = to_plot_baseline, linetype = "dashed")
    }
    gplt <- gplt +
        # The threshold
        ggplot2::geom_line(mapping = do.call(ggplot2::aes_string, all_aes), data = to_plot_threshold, linetype = "dashed", color = "black") +
        # Other plotting options
        ggplot2::labs(x = "Chemical Shift (ppm)", y = "Intensity (a.u.)") +
        ggplot2::scale_x_reverse(limits = rev(chemshift_range[seq_len(2)])) +
        ggplot2::scale_y_continuous(labels = scales::label_number(scale_cut = scales::cut_si("")), limits = c(0, ymax)) +
        ggplot2::facet_wrap(~ factor(NMRExperiment, levels = unique(NMRExperiment))) +
        ggplot2::theme(legend.position = "none")
    gplt
}

tidy_spectra_baseline_and_threshold <- function(dataset, thresholds, chemshift_range, NMRExperiment, columns = character(0L)) {
    to_plot <- tidy(
        dataset,
        chemshift_range = chemshift_range,
        NMRExperiment = NMRExperiment,
        columns = columns,
        matrix_name = "data_1r"
    )
    if ("data_1r_baseline" %in% names(unclass(dataset))) {
        to_plot_baseline <- tidy(
            dataset,
            chemshift_range = chemshift_range,
            NMRExperiment = NMRExperiment,
            columns = columns,
            matrix_name = "data_1r_baseline"
        )
        if (is.null(thresholds)) {
            to_plot_threshold <- NULL
        } else {
            to_plot_threshold <- dplyr::left_join(
                to_plot_baseline,
                tibble::enframe(
                    thresholds,
                    name = "NMRExperiment",
                    value = "threshold"
                ),
                by = "NMRExperiment"
            )
            to_plot_threshold$intensity <- to_plot_threshold$intensity + to_plot_threshold$threshold
        }
    } else {
        to_plot_baseline <- NULL
        if (is.null(thresholds)) {
            to_plot_threshold <- NULL
        } else {
            to_plot_threshold <- dplyr::left_join(
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
    }
    list(
        spectra = to_plot,
        baselines = to_plot_baseline,
        thresholds = to_plot_threshold
    )
}
