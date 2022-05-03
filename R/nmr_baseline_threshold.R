#' Threshold estimation for peak detection
#'
#' Estimates a robust baseline threshold value for peak detection on an [nmr_dataset_1D] object.
#'    This is performed computing the median and the mean absolute deviation of each spectrum
#'    in the 9.5-10 ppm range. The threshold is computed as the `median+3*mad`.
#'
#' @family peak detection functions
#' @family nmr_dataset_1D functions
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
        spec_region <- nmr_dataset$data_1r[i, threshold_ind]
        out[i] <- stats::median(spec_region) + 3*stats::mad(spec_region)
    }
    names(out) <- nmr_meta_get_column(nmr_dataset, column = "NMRExperiment")
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
#'
nmr_baseline_threshold_plot <- function(nmr_dataset, thresholds, NMRExperiment = "all", chemshift_range = c(9.5, 10), ...) {
    if (NMRExperiment != "all") {
        thresholds <- thresholds[NMRExperiment]
    }
    plot(nmr_dataset, chemshift_range = chemshift_range, NMRExperiment = NMRExperiment, ...) +
        ggplot2::geom_hline(
            ggplot2::aes(yintercept = .data$threshold),
            data = tibble::enframe(thresholds, name = "NMRExperiment", value = "threshold")
        ) +
        ggplot2::facet_wrap(~NMRExperiment) +
        ggplot2::theme(legend.position = "none")
}