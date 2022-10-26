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
#' @importFrom glue glue glue_collapse
#' @examples
#' dir_to_demo_dataset <- system.file("dataset-demo", package = "AlpsNMR")
#' dataset <- nmr_read_samples_dir(dir_to_demo_dataset)
#' my_nmr_dataset <- dataset %>%
#'     nmr_interpolate_1D(axis = c(0.4, 10)) %>%
#'     nmr_exclude_region(exclude = list(water = c(4.6, 5))) %>%
#'     nmr_normalize(method = "pqn") %>%
#'     plot()
"_PACKAGE"


#' @importFrom generics tidy
#' @export
generics::tidy


#' @importFrom dplyr rename
#' @export
dplyr::rename


#' @importFrom dplyr filter
#' @export
dplyr::filter

#' @importFrom magrittr %>%
#' @export
magrittr::`%>%`


#' @importFrom utils .DollarNames
#' @export
utils::.DollarNames



.refresh_nmr_dataset_rds <- function() {
    # Run this private function to refresh the extdata/nmr_dataset.rds file
    # if you change the demo data or the nmr_dataset object
    dir_to_demo_dataset <- system.file("dataset-demo", package = "AlpsNMR")
    dataset <- nmr_read_samples_dir(dir_to_demo_dataset)
    dataset_1D <- nmr_interpolate_1D(dataset, axis = c(min = 1, max = 5, by = 0.02))
    # annot_fn <- system.file("dataset-demo", "dummy_metadata.xlsx", package = "AlpsNMR")
    # annotations <- readxl::read_excel(annot_fn)
    annotations <- 
        structure(
            list(NMRExperiment = c("10", "20", "30"),
                 SubjectID = c("Ana", "Ana", "Elia"),
                 TimePoint = c("baseline", "3 months", "baseline")
            ),
            class = c("tbl_df", "tbl", "data.frame"),
            row.names = c(NA, -3L)
        )
    dataset_1D <- nmr_meta_add(dataset_1D, metadata = annotations, by = "NMRExperiment")
    nmr_dataset_save(nmr_dataset = dataset_1D, "inst/extdata/nmr_dataset.rds")
}