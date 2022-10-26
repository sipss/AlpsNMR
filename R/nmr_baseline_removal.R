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
#' dataset_1D <- nmr_dataset_load(system.file("extdata", "nmr_dataset.rds", package = "AlpsNMR"))
#' dataset_no_base_line <- nmr_baseline_removal(dataset_1D, lambda = 6, p = 0.01)
#'
nmr_baseline_removal <- function(nmr_dataset,
    lambda = 6,
    p = 0.05,
    maxit = 20) {
    results <- baseline::baseline(
        nmr_dataset$data_1r,
        method = "als",
        lambda = lambda,
        p = p,
        maxit = maxit
    )
    nmr_dataset$data_1r <- baseline::getCorrected(results)
    nmr_dataset
}


#' Estimate the baseline on an [nmr_dataset_1D] object, using [baseline::baseline.als].
#'
#' @family baseline removal functions
#' @seealso [baseline::baseline.als]
#' @param nmr_dataset An [nmr_dataset_1D].
#' @inheritParams baseline::baseline.als
#' @return The same [nmr_dataset_1D] object with the `data_1r_baseline` element.
#' @export
#'
#' @examples
#' dataset_1D <- nmr_dataset_load(system.file("extdata", "nmr_dataset.rds", package = "AlpsNMR"))
#' dataset_1D <- nmr_baseline_estimation(dataset_1D, lambda = 9, p = 0.01)
#'
nmr_baseline_estimation <- function(nmr_dataset,
    lambda = 9,
    p = 0.05,
    maxit = 20) {
    results <- baseline::baseline(
        nmr_dataset$data_1r,
        method = "als",
        lambda = lambda,
        p = p,
        maxit = maxit
    )
    nmr_dataset$data_1r_baseline <- baseline::getBaseline(results)
    nmr_dataset
}
