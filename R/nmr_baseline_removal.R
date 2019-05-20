#' Baseline Removal NMR
#'
#' Removes the baseline on an [nmr_dataset_1D] object, using [baseline::baseline.als].
#' 
#' @family baseline removal functions
#' @family nmr_dataset_1D functions
#' @seealso [baseline::baseline.als]
#' @param nmr_dataset An [nmr_dataset_1D].
#' @inheritParams baseline::baseline.als
#' @return The same [nmr_dataset_1D] object after baseline removal.
#' @export
#'
#'
nmr_baseline_removal <- function(nmr_dataset,
                             lambda = 6,
                             p = 0.05,
                             maxit = 20){
  
 results <- baseline::baseline.als(nmr_dataset$data_1r, lambda = lambda, p = p, maxit = maxit)
 nmr_dataset$data_1r <- results$corrected
 nmr_dataset
}
