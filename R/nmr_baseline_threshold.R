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
#' @return Numerical. A threshold value in intensity below that no peak is detected.
#' @export
#'

nmr_baseline_threshold <- function(nmr_dataset) {
    range_noise_ppm = c(9.5, 10)
    threshold_ind = which((nmr_dataset$axis > range_noise_ppm[1]) |
                                                    (nmr_dataset$axis < range_noise_ppm[2]))
    cent <- mean(apply(nmr_dataset$data_1r[, threshold_ind], 2, mean))
    disp <-
        3 * mean(apply(nmr_dataset$data_1r[, threshold_ind], 2, stats::sd))
    baselineThresh <- cent + disp
    baselineThresh
}
