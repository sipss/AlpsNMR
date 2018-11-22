
#' Set/Return the full spectra matrix
#'
#' @param nmr_dataset The NMR object to get the raw data from
#' @param ... Unused and left for future compatibility
#'
#' @return a matrix
#' @export
#'
nmr_data <- function(nmr_dataset, ...) {
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


#' @rdname nmr_data
#' @param value A matrix
#'
#' @return The given nmr_dataset
#' @export
#'
"nmr_data<-" <- function(nmr_dataset, value) {
  UseMethod("nmr_data<-")
}

#' @rdname nmr_data
#' @export
"nmr_data<-.nmr_dataset_1D" <- function(nmr_dataset, value) {
  value <- as.matrix(value)
  stopifnot(nrow(value) == nmr_dataset$num_samples)
  rownames(value) <- colnames(value) <- NULL
  nmr_dataset$data_1r <- value
  nmr_dataset
}