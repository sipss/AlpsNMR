#' Set/Return the full spectra matrix
#'
#' @param nmr_dataset An object from the [nmr_dataset_family] to get the raw data from
#' @param ... Unused and left for future compatibility
#'
#' @return a matrix
#' @export
#' @family import/export functions
#' @examples
#' dir_to_demo_dataset <- system.file("dataset-demo", package = "AlpsNMR")
#' dataset <- nmr_read_samples_dir(dir_to_demo_dataset)
#' dataset_1D <- nmr_interpolate_1D(dataset, axis = c(min = -0.5, max = 10, by = 2.3E-4))
#' dataset_data <- nmr_data(dataset_1D)
nmr_data <- function(nmr_dataset, ...) {
    UseMethod("nmr_data")
}

#' @rdname nmr_data
#' @param what What data do we want to get (default: `data_1r`)
#' @export
#' @examples
#' dir_to_demo_dataset <- system.file("dataset-demo", package = "AlpsNMR")
#' dataset <- nmr_read_samples_dir(dir_to_demo_dataset)
#' dataset_1D <- nmr_interpolate_1D(dataset, axis = c(min = -0.5, max = 10, by = 2.3E-4))
#' dataset_1D_data <- nmr_data(dataset_1D)
nmr_data.nmr_dataset_1D <- function(nmr_dataset, what = "data_1r", ...) {
    mat <- nmr_dataset[[what]]
    rownames(mat) <- names(nmr_dataset)
    colnames(mat) <- nmr_dataset$axis
    mat
}

#' @noRd
#' @export
#' @examples
#' dir_to_demo_dataset <- system.file("dataset-demo", package = "AlpsNMR")
#' dataset <- nmr_read_samples_dir(dir_to_demo_dataset)
#' dataset_1D <- nmr_interpolate_1D(dataset, axis = c(min = -0.5, max = 10, by = 2.3E-4))
#' peak_table <- nmr_data(dataset_1D)
nmr_data.nmr_dataset_peak_table <- function(nmr_dataset, ...) {
    peak_table <- nmr_dataset$peak_table
    rownames(peak_table) <- names(nmr_dataset)
    peak_table
}


#' @rdname nmr_data
#' @param value A matrix
#' @param ... Passed on to methods for compatibility
#'
#' @return The given nmr_dataset
#' @export
"nmr_data<-" <- function(nmr_dataset, ..., value) {
    UseMethod("nmr_data<-")
}

#' @rdname nmr_data
#' @export
"nmr_data<-.nmr_dataset_1D" <- function(nmr_dataset, what = "data_1r", ..., value) {
    if (!is.null(value)) {
        value <- as.matrix(value)
        stopifnot(nrow(value) == nmr_dataset$num_samples)
        rownames(value) <- colnames(value) <- NULL
    }
    nmr_dataset[[what]] <- value
    nmr_dataset
}
