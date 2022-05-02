#' AlpsNMR: Automated spectraL Processing System for NMR
#'
#' AlpsNMR allows you to import NMR spectra into R and provides automated and efficient signal processing for untargeted NMR metabolomics.
#'
#' The following functions can be combined with the pipe. They
#' create or modify the [nmr_dataset] object.
#'
#' - [nmr_read_samples_dir()] or [nmr_read_samples()]
#' - [nmr_interpolate_1D()]
#' - [nmr_exclude_region()]
#' - [nmr_normalize()]
#' - [plot()][plot.nmr_dataset_1D()]
#' 
#' There are also functions to extract the metadata and submit the samples to
#' irods, see the example below.
#'
#' The [nmr_dataset] object is essentially a list, so it is easy to access
#' its components for further analysis.
#'
#' @importFrom magrittr %>%
#' @examples
#' dir_to_demo_dataset <- system.file("dataset-demo", package = "AlpsNMR")
#' dataset <- nmr_read_samples_dir(dir_to_demo_dataset)
#' my_nmr_dataset <- dataset %>%
#'   nmr_interpolate_1D(axis = c(0.4, 10)) %>%
#'   nmr_exclude_region(exclude = list(water = c(4.6, 5))) %>%
#'   nmr_normalize(method = "pqn") %>%
#'   plot
"_PACKAGE"
