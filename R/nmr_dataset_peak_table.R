#' nmr_dataset_peak_table (S3 class)
#'
#' An `nmr_dataset_peak_table` represents a peak table with metadata.
#' It is defined as an S3 class, and it can be treated as a regular list.
#'
#' - `metadata`: A list of data frames. Each data frame contains metadata. Usually
#' the list only has one data frame named "external".
#'
#' - `peak_table`: A matrix with one sample on each row and the peaks in the
#' columns
#'
#' @name nmr_dataset_peak_table
NULL

#' Validate nmr_dataset_peak_table objects
#' @param nmr_dataset_peak_table An [nmr_dataset_peak_table] object
#' @return The [nmr_dataset_peak_table] unchanged
#'
#' This function is useful for its side-effects: Stopping in case of error
#'
#' @family nmr_dataset_peak_table functions
#' @family class helper functions
#' @export
#' @examples 
#' dir_to_demo_dataset <- system.file("dataset-demo", package = "AlpsNMR")
#' nmr_dataset <- nmr_read_samples_dir(dir_to_demo_dataset)
#' nmr_dataset <- nmr_interpolate_1D(dataset, axis = c(min = -0.5, max = 10, by = 2.3E-4))
#'
#' # 1.Peak detection in the dataset.
#' peak_data <- nmr_detect_peaks(nmr_dataset,
#'                               nDivRange_ppm = 0.1, # Size of detection segments
#'                               scales = seq(1, 16, 2),
#'                               baselineThresh = 0, # Minimum peak intensity
#'                               SNR.Th = 4) # Signal to noise ratio
#'
#' # 2.Find the reference spectrum to align with.
#' NMRExp_ref <- nmr_align_find_ref(nmr_dataset, peak_data)
#'
#' # 3.Spectra alignment using the ref spectrum and a maximum alignment shift
#' nmr_dataset <- nmr_align(nmr_dataset, # the dataset
#'                          peak_data, # detected peaks
#'                          NMRExp_ref = NMRExp_ref, # ref spectrum
#'                          maxShift_ppm = 0.0015, # max alignment shift
#'                          acceptLostPeak = FALSE) # lost peaks
#'
#' # 4.PEAK INTEGRATION (please, consider previous normalization step).
#' # First we take the peak table from the reference spectrum
#' peak_data_ref <- filter(peak_data, NMRExperiment == NMRExp_ref)
#'
#' # Then we integrate spectra considering the peaks from the ref spectrum
#' nmr_peak_table <- nmr_integrate_peak_positions(
#'                       samples = nmr_dataset,
#'                       peak_pos_ppm = peak_data_ref$ppm,
#'                       peak_width_ppm = NULL)
#' 
#' validate_nmr_dataset_peak_table(nmr_peak_table)
#'
validate_nmr_dataset_peak_table <-
    function(nmr_dataset_peak_table) {
        validate_nmr_dataset_family(nmr_dataset_peak_table)
        assert_that(inherits(nmr_dataset_peak_table, "nmr_dataset_peak_table"),
                                msg = "Not an nmr_dataset_peak_table")
        
        assert_that("peak_table" %in% names(nmr_dataset_peak_table),
                                msg = "nmr_dataset_peak_table must have a peak_table matrix")
        
        peak_table <- nmr_dataset_peak_table[["peak_table"]]
        assert_that(is.matrix(peak_table) && is.numeric(peak_table),
                                msg = "peak_table must be a numeric matrix")
        
        num_samples <- nrow(peak_table)
        
        assert_that(num_samples == nmr_dataset_peak_table[["num_samples"]],
                                msg = "The num_samples value does not match nrow(peak_table)")
        
        
        
        nmr_dataset_peak_table
    }

#' Creates a new nmr_dataset_peak_table object from scratch
#' 
#' @param peak_table A numeric matrix with one NMR spectrum on each row
#' @param metadata A list of data frames with at least the `NMRExperiment` column
#' @return Creates a new nmr_dataset_peak_table object from scratch
#' @name new_nmr_dataset_peak_table 
#' @importFrom assertthat assert_that
#' @importFrom glue glue
#' @family nmr_dataset_peak_table functions
#' @family class helper functions
#' @export
#' @examples
#' metadata_1D <- list(external = data.frame(NMRExperiment = c("10", "20")))
#' # Sample 10 and Sample 20 can have different lengths (due to different setups)
#' data_fields_1D <- list(data_1r = list(runif(16), runif(32)))
#' # Each sample has its own axis list, with one element (because this example is 1D)
#' axis_1D <- list(list(1:16), list(1:32))
#' nmr_dataset <- new_nmr_dataset(metadata_1D, data_fields_1D, axis_1D)
#'
#' # 1.Peak detection in the dataset.
#' peak_data <- nmr_detect_peaks(nmr_dataset,
#'                               nDivRange_ppm = 0.1, # Size of detection segments
#'                               scales = seq(1, 16, 2),
#'                               baselineThresh = 0, # Minimum peak intensity
#'                               SNR.Th = 4) # Signal to noise ratio
#'
#' # 2.Find the reference spectrum to align with.
#' NMRExp_ref <- nmr_align_find_ref(nmr_dataset, peak_data)
#'
#' # 3.Spectra alignment using the ref spectrum and a maximum alignment shift
#' nmr_dataset <- nmr_align(nmr_dataset, # the dataset
#'                          peak_data, # detected peaks
#'                          NMRExp_ref = NMRExp_ref, # ref spectrum
#'                          maxShift_ppm = 0.0015, # max alignment shift
#'                          acceptLostPeak = FALSE) # lost peaks
#'
#' # 4.PEAK INTEGRATION (please, consider previous normalization step).
#' # First we take the peak table from the reference spectrum
#' peak_data_ref <- filter(peak_data, NMRExperiment == NMRExp_ref)
#'
#' # Then we integrate spectra considering the peaks from the ref spectrum
#' nmr_peak_table <- nmr_integrate_peak_positions(
#'                       samples = nmr_dataset,
#'                       peak_pos_ppm = peak_data_ref$ppm,
#'                       peak_width_ppm = NULL)
#'
#' new_nmr_dataset_peak_table(nmr_peak_table, metadata_1D)
#'
new_nmr_dataset_peak_table <- function(peak_table, metadata) {
    samples <- list()
    samples[["metadata"]] <- metadata
    samples[["peak_table"]] <- as.matrix(peak_table)
    samples[["num_samples"]] <- nrow(peak_table)
    class(samples) <-
        c("nmr_dataset_peak_table", "nmr_dataset_family")
    validate_nmr_dataset_peak_table(samples)
    samples
}

#' Object is of [nmr_dataset_peak_table] class
#' @param x An object
#' @return `TRUE` if the object is an `nmr_dataset_peak_table`, `FALSE` otherwise
#' @export
#' @family nmr_dataset_peak_table functions
#' @family class helper functions
is.nmr_dataset_peak_table <-
    function(x)
        inherits(x, "nmr_dataset_peak_table")

#' @export
#' @family nmr_dataset_peak_table functions
#' @family class helper functions
print.nmr_dataset_peak_table <- function(x, ...) {
    cat(format(x, ...), "\n")
    invisible(x)
}

#' @export
#' @family nmr_dataset_peak_table functions
#' @family class helper functions
format.nmr_dataset_peak_table <- function(x, ...) {
    paste0(
        "An nmr_dataset_peak_table (",
        x$num_samples,
        " samples, and ",
        ncol(x$peak_table),
        " peaks)"
    )
}

#' Extract parts of an nmr_dataset_peak_table
#' @param x an [nmr_dataset_peak_table] object
#' @param i indices of the samples to keep
#' @return an nmr_dataset_peak_table with the extracted samples
#' @family subsetting functions
#' @family nmr_dataset_peak_table functions
#' @export
`[.nmr_dataset_peak_table` <- function(x, i) {
    output <- x
    output$metadata <- purrr::map(output$metadata, function(metad) {
        metad[i, , drop = FALSE]
    })
    output[["peak_table"]] <-
        output[["peak_table"]][i, , drop = FALSE]
    output$num_samples <- nrow(output$metadata[[1]])
    validate_nmr_dataset_peak_table(output)
    return(output)
}
