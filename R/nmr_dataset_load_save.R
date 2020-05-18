#' @return Functions to load and save nmr_dataset objects
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
#' @return load nmr dataset
#' @examples 
#' # dataset <- nmr_dataset_load("test")
#' nmr_dataset <- nmr_dataset_load(system.file("extdata", "nmr_dataset.rds", package = "AlpsNMR"))
#' 
nmr_dataset_load <- function(file_name) {
    return(readRDS(file_name))
}

#' @rdname load_and_save_functions
#' @export
#' @return save nmr dataset
#' @examples 
#' dir_to_demo_dataset <- system.file("dataset-demo", package = "AlpsNMR")
#' dataset <- nmr_read_samples_dir(dir_to_demo_dataset)
#' # nmr_dataset_save(dataset, "test")
#' 
nmr_dataset_save <- function(nmr_dataset, file_name, ...) {
    attr(nmr_dataset, "normalize_extra_info") <- NULL
    saveRDS(nmr_dataset, file_name)
    return(nmr_dataset)
}
