#' Baseline Removal NMR
#'
#' Removes the baseline on an [nmr_dataset_1D] object, using [baseline::baseline.als].
#'
#' @family baseline removal functions
#' @seealso [baseline::baseline.als]
#' @param nmr_dataset An [nmr_dataset_1D].
#' @inheritParams baseline::baseline.als
#' @return The same [nmr_dataset_1D] object after baseline removal.
#' @export
#'
#' @examples 
#' dir_to_demo_dataset <- system.file("dataset-demo", package = "AlpsNMR")
#' dataset <- nmr_read_samples_dir(dir_to_demo_dataset)
#' dataset_1D <- nmr_interpolate_1D(dataset, axis = c(min = -0.5, max = 10, by = 2.3E-4))
#' dataset_no_base_line <- nmr_baseline_removal(dataset_1D, lambda = 6, p = 0.01)
#' 
nmr_baseline_removal <- function(nmr_dataset,
                                 lambda = 6,
                                 p = 0.05,
                                 maxit = 20) {
    results <-
        baseline::baseline.als(
            nmr_dataset$data_1r,
            lambda = lambda,
            p = p,
            maxit = maxit
        )
    nmr_dataset$data_1r <- results$corrected
    nmr_dataset
}
