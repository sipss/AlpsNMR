#' Threshold estimation for peak detection
#'
#' Estimates the threshold value for peak detection on an [nmr_dataset_1D] object.
#'    This is performed computing the mean and the standard deviation of each spectrum
#'    beyond 9.5 ppm. The threshold is then averaged of means and adding 3 times the
#'    mean of the standard deviations
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
    if (length(range_without_peaks) != 2) {
        rlang::abort("range_without_peaks must have length 2")
    }
    r_start <- min(range_without_peaks)
    r_end <- max(range_without_peaks)
    threshold_ind <- nmr_dataset$axis >= r_start & nmr_dataset$axis < r_end
    cent <- mean(apply(nmr_dataset$data_1r[, threshold_ind], 2, mean))
    disp <- mean(apply(nmr_dataset$data_1r[, threshold_ind], 2, stats::sd))
    cent + 3*disp
}
