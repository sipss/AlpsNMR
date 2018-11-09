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

verify_axisn <- function(axisn, one_sample_axis) {
  if (is.null(axisn)) {
    # axisn will be built on sample_axis:
    axisn <- c(min = min(one_sample_axis),
               max = max(one_sample_axis),
               by = abs(one_sample_axis[2] - one_sample_axis[1]))
  } else {
    if (length(axisn) == 2) {
      axisn <- c(axisn, abs(one_sample_axis[2] - one_sample_axis[1]))
      names(axisn) <- c("min", "max", "by")
    }
    if (is.null(names(axisn))) {
      names(axisn) <- c("min", "max", "by")
    }
  }
  return(axisn)
}

#' Interpolate a set of 1D NMR Spectra
#' @param samples An NMR dataset
#' @param axis1 The desired ppm axis range and optionally the ppm step
#' @export
nmr_interpolate_1D <- function(samples, axis1 = c(min = 0.2, max = 10, by = 0.0008)) {
  UseMethod("nmr_interpolate_1D")
}

#' @rdname nmr_interpolate_1D
#' @export
nmr_interpolate_1D.nmr_dataset <- function(samples, axis1 = c(min = 0.2, max = 10, by = 0.0008)) {
  # Check if we can interpolate:
  verify_dimensionality(samples, valid_dimensions = 1)
  
  # 2. Check that we have the interpolation axis
  axis1 <- verify_axisn(axis1, samples[["axis"]][[1]][[1]])
  
  axis1_full <- seq(from = axis1["min"], to = axis1["max"], by = axis1["by"])
  if (show_progress_bar(samples$num_samples > 5)) {
    message("Interpolating data_1r...")
  }
  data_1r <- interpolate_1d(list_of_ppms = purrr::map(samples[["axis"]], 1),
                            list_of_1r = samples[["data_1r"]],
                            ppm_axis = axis1_full)

  new_nmr_dataset_1D(ppm_axis = axis1_full,
                     data_1r = data_1r,
                     metadata = samples$metadata)
}


interpolate_1d <- function(list_of_ppms, list_of_1r, ppm_axis) {
  num_samples <- length(list_of_1r)
  data_matr <- matrix(NA, nrow = num_samples, ncol = length(ppm_axis))
  tryCatch({
    pb <- NULL
    if (show_progress_bar(num_samples > 5)) {
      pb <- utils::txtProgressBar(min = 0, max = num_samples, style = 3)
    }
    for (i in seq_len(num_samples)) {
      data_matr[i, ] <- signal::interp1(x = list_of_ppms[[i]],
                                        y = list_of_1r[[i]],
                                        xi = ppm_axis,
                                        method = "spline")
      if (!is.null(pb)) {
        utils::setTxtProgressBar(pb, i)
      }
    }
  }, finally = {
    if (!is.null(pb)) {
      close(pb)
    }
  })
  data_matr
}
