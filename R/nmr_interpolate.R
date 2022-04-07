verify_dimensionality <- function(samples, valid_dimensions) {
    # 1. Check that each of the loaded samples has the same dimensionality
    dimensions_per_sample <- vapply(samples$axis, length, numeric(1))
    if (any(dimensions_per_sample != dimensions_per_sample[1])) {
        stop("Samples were expected to have the same dimensionality")
    }
    if (!all(dimensions_per_sample %in% valid_dimensions)) {
        stop("Samples have unexpected dimensionality")
    }
    dimensions_per_sample
} 

verify_axisn <- function(axisn, samples) {
    if (is.null(axisn)) {
        axisn <- c(
          min = max(purrr::map_dbl(samples$axis, ~ min(.[[1]]))),
          max = min(purrr::map_dbl(samples$axis, ~ max(.[[1]]))),
          by = min(purrr::map_dbl(samples$axis, ~ stats::median(abs(diff(.[[1]])))))
        )
    } else {
        if (length(axisn) == 2) {
            axisn <- c(axisn, min(purrr::map_dbl(samples$axis, ~ stats::median(abs(diff(.[[1]]))))))
            names(axisn) <- c("min", "max", "by")
        }
        if (is.null(names(axisn))) {
            names(axisn) <- c("min", "max", "by")
        }
    }
    return(axisn)
}

#' Interpolate a set of 1D NMR Spectra
#' 
#' @param samples An NMR dataset
#' @param axis The ppm axis range and optionally the ppm step. Set it to `NULL` for autodetection
#' @return Interpolate a set of 1D NMR Spectra
#' @name nmr_interpolate_1D
#' @export
#' @examples 
#' dir_to_demo_dataset <- system.file("dataset-demo", package = "AlpsNMR")
#' dataset <- nmr_read_samples_dir(dir_to_demo_dataset)
#' dataset_1D <- nmr_interpolate_1D(dataset, axis = c(min = -0.5, max = 10, by = 2.3E-4))
#' 
nmr_interpolate_1D <- function(samples, axis = c(min = 0.2, max = 10, by = 0.0008)) {
    UseMethod("nmr_interpolate_1D")
}

#' @rdname nmr_interpolate_1D
#' @family nmr_dataset functions
#' @family nmr_dataset_1D functions
#' @export 
#' @examples 
#' dir_to_demo_dataset <- system.file("dataset-demo", package = "AlpsNMR")
#' dataset <- nmr_read_samples_dir(dir_to_demo_dataset)
#' dataset_1D <- nmr_interpolate_1D(dataset, axis = c(min = -0.5, max = 10, by = 2.3E-4))
#' 
nmr_interpolate_1D.nmr_dataset <- function(samples, axis = c(min = 0.2, max = 10, by = 0.0008)) {
    # Check if we can interpolate:
    verify_dimensionality(samples, valid_dimensions = 1)
    
    # 2. Check that we have the interpolation axis
    axis <- verify_axisn(axis, samples)
    
    axis_full <- seq(from = axis["min"], to = axis["max"], by = axis["by"])
    if (show_progress_bar(samples$num_samples > 5)) {
        message("Interpolating data_1r...")
    }
    data_1r <- interpolate_1d(list_of_ppms = purrr::map(samples[["axis"]], 1),
                                                        list_of_1r = samples[["data_1r"]],
                                                        ppm_axis = axis_full)

    new_nmr_dataset_1D(ppm_axis = axis_full,
                       data_1r = data_1r,
                        metadata = samples$metadata)
}


interpolate_1d <- function(list_of_ppms, list_of_1r, ppm_axis) {
    num_samples <- length(list_of_1r)
    data_matr <- matrix(NA, nrow = num_samples, ncol = length(ppm_axis))
    tryCatch({
        pb <- NULL
        if (show_progress_bar(num_samples > 5)) {
          pb <- progress_bar_new(
            name = "Interpolating samples",
            total = num_samples
          )
        }
        for (i in seq_len(num_samples)) {
            data_matr[i, ] <- signal::interp1(x = list_of_ppms[[i]],
                                              y = list_of_1r[[i]],
                                              xi = ppm_axis,
                                              method = "spline")
            progress_bar_update(pb)
        }
    }, finally = {
        progress_bar_end(pb)
    })
    data_matr
}
