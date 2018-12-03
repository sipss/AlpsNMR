#' Functions to load and save nmr_dataset objects
#' 
#' @name load_and_save_functions
#' @param file_name The file name to load or save to
#' @param nmr_dataset An object from the [nmr_dataset_family]
#' @param ... Additional arguments passed to [saveRDS].
#' @family nmr_dataset functions
#' @family nmr_dataset_1D functions
#' @family nmr_dataset_peak_table functions
#' @family import/export functions
NULL

#' @rdname load_and_save_functions
#' @export
nmr_dataset_load <- function(file_name) {
  return(readRDS(file_name))
}

#' @rdname load_and_save_functions
#' @export
nmr_dataset_save <- function(nmr_dataset, file_name, ...) {
  nmr_diagnose(nmr_dataset) <- NULL
  saveRDS(nmr_dataset, file_name)
  return(nmr_dataset)
}
