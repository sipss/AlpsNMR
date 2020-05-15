#' NMR file lister
#'
#' The function lists samples from the chosen folder required to import and
#' create a [nmr_dataset_1D] object. The function is based on the [fs::dir_ls()]
#' function.
#' @return lists of samples from the chosen folder
#' @param dataset_path_nmr A character vector of the path where samples are.
#' @param glob A wildcard or globbing pattern common for the samples to be read,
#'   for example ending with *0 (spectra acquired by a NOESY sequence often end
#'   by 0: 10, 20, 30...) or *s (for example, samples from the tutorial in this
#'   package) passed on to `grep()` to filter paths.
#'
#'
#' @family nmr_dataset_1D functions
#' @export
#'
#' @examples
#' dir_to_demo_dataset <- system.file("dataset-demo", package = "AlpsNMR")
#' lists_of_samples <- file_lister(dir_to_demo_dataset, "*0")
#' 
file_lister <- function(dataset_path_nmr, glob) {
    as.character(fs::dir_ls(dataset_path_nmr, glob = glob))
}
