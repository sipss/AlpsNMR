#' Set/Return the full spectra matrix
#'
#' @param nmr_dataset An object from the [nmr_dataset_family] to get the raw data from
#' @param ... Unused and left for future compatibility
#'
#' @return a matrix
#' @export
#' @family import/export functions
#' @examples
#' dataset_rds <- system.file("extdata", "nmr_dataset.rds", package = "AlpsNMR")
#' dataset_1D <- nmr_dataset_load(dataset_rds)
#' dataset_data <- nmr_data(dataset_1D)
nmr_data <- function(nmr_dataset, ...) {
    UseMethod("nmr_data")
}

#' @rdname nmr_data
#' @param what What data do we want to get (default: `data_1r`)
#' @export
#' @examples
#' dataset_rds <- system.file("extdata", "nmr_dataset.rds", package = "AlpsNMR")
#' dataset_1D <- nmr_dataset_load(dataset_rds)
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
#' peak_table_exp <- matrix(1:6, nrow = 2)
#' metadata <- list(external = data.frame(NMRExperiment = letters[1:2]))
#' dataset_peak_table <- new_nmr_dataset_peak_table(peak_table, metadata)
#' peak_table <- nmr_data(dataset_peak_table)
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
