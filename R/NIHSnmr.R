#' NIHSnmr: Routines for importing and processing NMR data
#'
#' NIHSnmr allows you to import NMR samples into R and process them.
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
#' @examples
#' \dontrun{
#' library(NIHSnmr)
#' sample_names <- c("dataset/10/", "dataset/20/", "dataset/30/")
#' my_nmr_dataset <- nmr_read_samples(sample_names) %>%
#'   nmr_interpolate_1D(axis = c(0.4, 10)) %>%
#'   nmr_exclude_region(exclude = list(water = c(4.6, 5))) %>%
#'   nmr_normalize(method = "pqn") %>%
#'   plot
#' }
#' \dontrun{
#' library(NIHSnmr)
#' sample_names <- c("dataset/10/", "dataset/20/", "dataset/30/")
#' nmrdata <- nmr_read_samples(sample_names, metadata_only = TRUE)
#' irods_metadata <- nmr_get_irods_meta(nmrdata)
#' # explore and confirm everything is in irods_metadata
#' # Now, let's compress the samples to zip files for irods
#' irods_metadata <- nmr_prepare_zip_files(meta_irods = irods_metadata,
#'                                         workdir = 'my_dataset_zip_files')
#' # Finally push the data to the irods directory
#' nmr_push_to_irods(irods_metadata, "/NIHSData/DUND-XXXXXX/study/Metabolomics")
#' }
"_PACKAGE"
