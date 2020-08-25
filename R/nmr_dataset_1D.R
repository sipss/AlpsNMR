#' nmr_dataset_1D (S3 class)
#'
#' An `nmr_dataset_1D` represents a set of 1D interpolated NMR samples.
#' It is defined as an S3 class, and it can be treated as a regular list.
#'
#' It currently has the following elements:
#'
#'
#' - `metadata`: A list of data frames. Each data frame contains metadata of
#' a given area (acquisition parameters, preprocessing parameters, general sample information...)
#'
#' - `axis`: A numeric vector with the chemical shift axis in ppm.
#'
#' - `data_1r`: A matrix with one sample on each row and the chemical
#' shifts in the columns.
#'
#'
#' @name nmr_dataset_1D
#' @family AlpsNMR dataset objects
NULL

#' Validate 1D nmr datasets
#' @name validate_nmr_dataset
#' @param nmr_dataset_1D An [nmr_dataset_1D] object
#' @return The [nmr_dataset_1D] unchanged
#' 
#' This function is useful for its side-effects. Stopping in case of error
#'
#' @family class helper functions
#' @family nmr_dataset_1D functions
#' @export
#' @examples 
#' dir_to_demo_dataset <- system.file("dataset-demo", package = "AlpsNMR")
#' dataset <- nmr_read_samples_dir(dir_to_demo_dataset)
#' dataset_1D <- nmr_interpolate_1D(dataset, axis = c(min = -0.5, max = 10, by = 2.3E-4))
#' dataset_1D_validated <- validate_nmr_dataset_1D(dataset_1D)
#' 
validate_nmr_dataset_1D <- function(nmr_dataset_1D) {
    validate_nmr_dataset_family(nmr_dataset_1D)
    assert_that(inherits(nmr_dataset_1D, "nmr_dataset_1D"),
                            msg = "Not an nmr_dataset_1D")
    
    assert_that("axis" %in% names(nmr_dataset_1D),
                            msg = "nmr_dataset_1D must have a ppm axis")
    assert_that("data_1r" %in% names(nmr_dataset_1D),
                            msg = "nmr_dataset_1D must have a data_1r matrix")
    
    ppm_axis <- nmr_dataset_1D[["axis"]]
    data_1r <- nmr_dataset_1D[["data_1r"]]
    assert_that(is.vector(ppm_axis) && is.numeric(ppm_axis),
                            msg = "axis must be a numeric vector")
    assert_that(is.matrix(data_1r) && is.numeric(data_1r),
                            msg = "data_1r must be a numeric matrix")
    
    assert_that(length(ppm_axis) == ncol(data_1r),
                            msg = "ppm axis does not have a length equal to ncol(data_1r)")
    num_samples <- nrow(data_1r)
    
    assert_that(num_samples == nmr_dataset_1D[["num_samples"]],
                            msg = "The num_samples value does not match nrow(data_1r)")
    
    nmr_dataset_1D
}

#' Creates a new 1D nmr_dataset object from scratch
#'
#' @name new_nmr_dataset_1D
#' @param ppm_axis A numeric vector with the ppm values for the columns of data_1r
#' @param data_1r A numeric matrix with one NMR spectrum on each row
#' @param metadata A list of data frames with at least the `NMRExperiment` column
#'
#' @importFrom assertthat assert_that
#' @importFrom glue glue
#' @family class helper functions
#' @family nmr_dataset_1D functions
#' @export
#' @return Creates a new 1D nmr_dataset object from scratch
#' 
#' @examples
#' # Create a random spectra matrix
#' nsamp <- 12
#' npoints <- 20
#' dummy_ppm_axis <- seq(from = 0.2, to = 10, length.out = npoints)
#' dummy_spectra_matrix <- matrix(runif(nsamp*npoints), nrow = nsamp, ncol = npoints)
#' metadata <- list(external = data.frame(NMRExperiment = paste0("Sample", 1:12),
#'                                        DummyClass = c("a", "b"),
#'                                        stringsAsFactors = FALSE))
#' dummy_nmr_dataset_1D <- new_nmr_dataset_1D(ppm_axis = dummy_ppm_axis,
#'                                            data_1r = dummy_spectra_matrix,
#'                                            metadata = metadata)
#'                                                  
#'                                                  
new_nmr_dataset_1D <- function(ppm_axis, data_1r, metadata) {
    samples <- list()
    samples[["metadata"]] <- metadata
    samples[["data_1r"]] <- data_1r
    samples[["axis"]] <- ppm_axis
    samples[["num_samples"]] <- nrow(data_1r)
    class(samples) <- c("nmr_dataset_1D", "nmr_dataset_family")
    validate_nmr_dataset_1D(samples)
    samples
}

#' Object is of [nmr_dataset_1D] class
#' @param x an [nmr_dataset_1D] object
#' @return `TRUE` if the object is an [nmr_dataset_1D], `FALSE` otherwise
#' @export
#' @family class helper functions
#' @family nmr_dataset_1D functions
#' @examples 
#' dir_to_demo_dataset <- system.file("dataset-demo", package = "AlpsNMR")
#' dataset <- nmr_read_samples_dir(dir_to_demo_dataset)
#' dataset_1D <- nmr_interpolate_1D(dataset, axis = c(min = -0.5, max = 10, by = 2.3E-4))
#' result <- is(dataset_1D)
is.nmr_dataset_1D <- function(x)
    inherits(x, "nmr_dataset_1D")

#' print for nmr_dataset_1D
#' @param x an [nmr_dataset_1D] object
#' @param ... for future use
#' @family class helper functions
#' @family nmr_dataset_1D functions
#' @return print for nmr_dataset_1D
#' @export
#' @examples 
#' dir_to_demo_dataset <- system.file("dataset-demo", package = "AlpsNMR")
#' dataset <- nmr_read_samples_dir(dir_to_demo_dataset)
#' dataset_1D <- nmr_interpolate_1D(dataset, axis = c(min = -0.5, max = 10, by = 2.3E-4))
#' print(dataset_1D)
print.nmr_dataset_1D <- function(x, ...) {
    cat(format(x, ...), "\n")
    invisible(x)
}

#' format for nmr_dataset_1D
#' @param x an [nmr_dataset_1D] object
#' @param ... for future use
#' @family class helper functions
#' @family nmr_dataset_1D functions
#' @return format for nmr_dataset_1D
#' @export
#' @examples 
#' dir_to_demo_dataset <- system.file("dataset-demo", package = "AlpsNMR")
#' dataset <- nmr_read_samples_dir(dir_to_demo_dataset)
#' dataset_1D <- nmr_interpolate_1D(dataset, axis = c(min = -0.5, max = 10, by = 2.3E-4))
#' format(dataset_1D)
format.nmr_dataset_1D <- function(x, ...) {
    paste0("An nmr_dataset_1D (", x$num_samples, " samples)")
}

#' Extract parts of an nmr_dataset_1D
#' @param x an [nmr_dataset_1D] object
#' @param i indices of the samples to keep
#' @return an nmr_dataset_1D with the extracted samples
#' @family subsetting functions
#' @family nmr_dataset_1D functions
#' @export
#' @examples 
#' dir_to_demo_dataset <- system.file("dataset-demo", package = "AlpsNMR")
#' dataset <- nmr_read_samples_dir(dir_to_demo_dataset)
#' dataset_1D <- nmr_interpolate_1D(dataset, axis = c(min = -0.5, max = 10, by = 2.3E-4))
#' dataset_1D[0]
`[.nmr_dataset_1D` <- function(x, i) {
    output <- x
    output$metadata <- purrr::map(output$metadata, function(metad) {
        metad[i, , drop = FALSE]
    })
    output[["data_1r"]] <- output[["data_1r"]][i, , drop = FALSE]
    output$num_samples <- nrow(output$metadata[[1]])
    validate_nmr_dataset_1D(output)
    return(output)
}

#' Export 1D NMR data to a CSV file
#'
#' @param nmr_dataset An [nmr_dataset_1D] object
#' @param filename The csv filename
#'
#' @return The nmr_dataset object (unmodified)
#' @export 
#' @examples 
#' dir_to_demo_dataset <- system.file("dataset-demo", package = "AlpsNMR")
#' dataset <- nmr_read_samples_dir(dir_to_demo_dataset)
#' dataset_1D <- nmr_interpolate_1D(dataset, axis = c(min = -0.5, max = 10, by = 2.3E-4))
#' #nmr_export_data_1r(dataset_1D, "exported_nmr_dataset")
nmr_export_data_1r <- function(nmr_dataset, filename) {
    # FIXME: remove me (nmr_data() covers for this)
    assert_that(is.nmr_dataset_1D(nmr_dataset), msg = "An nmr_dataset_1D should be given")
    data_1r <- nmr_data(nmr_dataset)
    utils::write.csv(data_1r, file = filename, row.names = FALSE)
    nmr_dataset
}
