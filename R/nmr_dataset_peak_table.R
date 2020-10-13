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
#' @rdname Peak_detection
validate_nmr_dataset_peak_table <- function(nmr_dataset_peak_table) {
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
#' dir_to_demo_dataset <- system.file("dataset-demo", package = "AlpsNMR")
#' dataset <- nmr_read_samples_dir(dir_to_demo_dataset)
#' dataset_1D <- nmr_interpolate_1D(dataset, axis = c(min = -0.5, max = 10, by = 2.3E-4))
#' meta <- file.path(dir_to_demo_dataset, "dummy_metadata.xlsx")
#' metadata <- readxl::read_excel(meta, sheet = 1)
#' dataset_1D <- nmr_meta_add(dataset_1D, metadata = metadata, by = "NMRExperiment")
#' metadata <- list(external = dataset_1D[["metadata"]][["external"]])
#' peak_table <- nmr_data(dataset_1D)
#' new <- new_nmr_dataset_peak_table(peak_table, metadata)
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
#' @param x an [nmr_dataset_peak_table] object
#' @return `TRUE` if the object is an `nmr_dataset_peak_table`, `FALSE` otherwise
#' @export
#' @family nmr_dataset_peak_table functions
#' @family class helper functions
#' @examples
#' dir_to_demo_dataset <- system.file("dataset-demo", package = "AlpsNMR")
#' dataset <- nmr_read_samples_dir(dir_to_demo_dataset)
#' dataset_1D <- nmr_interpolate_1D(dataset, axis = c(min = -0.5, max = 10, by = 2.3E-4))
#' meta <- file.path(dir_to_demo_dataset, "dummy_metadata.xlsx")
#' metadata <- readxl::read_excel(meta, sheet = 1)
#' dataset_1D <- nmr_meta_add(dataset_1D, metadata = metadata, by = "NMRExperiment")
#' metadata <- list(external = dataset_1D[["metadata"]][["external"]])
#' peak_table <- nmr_data(dataset_1D)
#' new <- new_nmr_dataset_peak_table(peak_table, metadata)
#' is(new)
#' 
is.nmr_dataset_peak_table <-
    function(x)
        inherits(x, "nmr_dataset_peak_table")

#' print for nmr_dataset_peak_table
#' @param x an [nmr_dataset_peak_table] object
#' @param ... for future use
#' @export
#' @return print for nmr_dataset_peak_table
#' @family nmr_dataset_peak_table functions
#' @family class helper functions
#' @examples
#' dir_to_demo_dataset <- system.file("dataset-demo", package = "AlpsNMR")
#' dataset <- nmr_read_samples_dir(dir_to_demo_dataset)
#' dataset_1D <- nmr_interpolate_1D(dataset, axis = c(min = -0.5, max = 10, by = 2.3E-4))
#' meta <- file.path(dir_to_demo_dataset, "dummy_metadata.xlsx")
#' metadata <- readxl::read_excel(meta, sheet = 1)
#' dataset_1D <- nmr_meta_add(dataset_1D, metadata = metadata, by = "NMRExperiment")
#' metadata <- list(external = dataset_1D[["metadata"]][["external"]])
#' peak_table <- nmr_data(dataset_1D)
#' new <- new_nmr_dataset_peak_table(peak_table, metadata)
#' new
print.nmr_dataset_peak_table <- function(x, ...) {
    cat(format(x, ...), "\n")
    invisible(x)
}

#' Format for nmr_dataset_peak_table
#' @param x an [nmr_dataset_peak_table] object
#' @param ... for future use
#' @export
#' @return Format for nmr_dataset_peak_table
#' @family nmr_dataset_peak_table functions
#' @family class helper functions
#' @examples
#' dir_to_demo_dataset <- system.file("dataset-demo", package = "AlpsNMR")
#' dataset <- nmr_read_samples_dir(dir_to_demo_dataset)
#' dataset_1D <- nmr_interpolate_1D(dataset, axis = c(min = -0.5, max = 10, by = 2.3E-4))
#' meta <- file.path(dir_to_demo_dataset, "dummy_metadata.xlsx")
#' metadata <- readxl::read_excel(meta, sheet = 1)
#' dataset_1D <- nmr_meta_add(dataset_1D, metadata = metadata, by = "NMRExperiment")
#' metadata <- list(external = dataset_1D[["metadata"]][["external"]])
#' peak_table <- nmr_data(dataset_1D)
#' new <- new_nmr_dataset_peak_table(peak_table, metadata)
#' format(new)
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
#' @examples
#' dir_to_demo_dataset <- system.file("dataset-demo", package = "AlpsNMR")
#' dataset <- nmr_read_samples_dir(dir_to_demo_dataset)
#' dataset_1D <- nmr_interpolate_1D(dataset, axis = c(min = -0.5, max = 10, by = 2.3E-4))
#' meta <- file.path(dir_to_demo_dataset, "dummy_metadata.xlsx")
#' metadata <- readxl::read_excel(meta, sheet = 1)
#' dataset_1D <- nmr_meta_add(dataset_1D, metadata = metadata, by = "NMRExperiment")
#' metadata <- list(external = dataset_1D[["metadata"]][["external"]])
#' peak_table <- nmr_data(dataset_1D)
#' new <- new_nmr_dataset_peak_table(peak_table, metadata)
#' new[0]
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
