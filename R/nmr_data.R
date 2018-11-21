
#' Return the full spectra matrix
#'
#' @param dataset 
#' @param ... 
#'
#' @return a matrix
#' @export
#'
nmr_data <- function(dataset, ...) {
  UseMethod("nmr_data")
}

#' @noRd
#' @export
nmr_data.nmr_dataset_1D <- function(nmr_dataset, ...) {
  data_1r <- nmr_dataset$data_1r
  rownames(data_1r) <- nmr_meta_get_column(nmr_dataset, "NMRExperiment")
  colnames(data_1r) <- nmr_dataset$axis
  data_1r
}