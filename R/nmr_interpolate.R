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

#' Interpolate samples
#'
#' Interpolate the samples so they share the same axis and can be placed
#' in the same matrix
#'
#' For 2-D samples, we use two 1-D linear interpolations.
#'
#'
#' @param samples A [nmr_dataset] object
#' @param axis1,axis2 Axis given as `c(minimum, maximum, step)`
#' @return The [nmr_dataset] object, with interpolated samples
#' @export
nmr_interpolate <- function(samples,
                            axis1=c(min = 0.2, max = 10, by = 0.0008),
                            axis2=NULL) {
  # Check if we can interpolate:
  dimensions_per_sample <- verify_dimensionality(samples, valid_dimensions = c(1, 2))
  
  if (dimensions_per_sample[1] == 1) {
    nmr_interpolate_1D_legacy(samples, axis1 = axis1)
  } else {
    nmr_interpolate_2D_legacy(samples, axis1 = axis1, axis2 = axis2)
  }
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

# _legacy, because it returns an nmr_dataset instead of an nmr_dataset_1D
# once the whole package works with the new structure, drop this
nmr_interpolate_1D_legacy <- function(samples, axis1 = c(min = 0.2, max = 10, by = 0.0008)) {
  # Check if we can interpolate:
  verify_dimensionality(samples, valid_dimensions = 1)
  
  # 2. Check that we have the interpolation axis
  axis1 <- verify_axisn(axis1, samples[["axis"]][[1]][[1]])
  
  data_fields <- names(samples)[grepl(pattern = "^data_.*", x = names(samples))]
  
  axis1_full <- seq(from = axis1["min"], to = axis1["max"], by = axis1["by"])
  for (data_field in data_fields) {
    if (show_progress_bar(samples$num_samples > 5)) {
      message("Interpolating ", data_field, "...")
    }
    data_matr <- interpolate_1d(list_of_ppms = purrr::map(samples[["axis"]], 1),
                                list_of_1r = samples[[data_field]],
                                ppm_axis = axis1_full)
    samples[[data_field]] <- data_matr
  }
  samples[["axis"]] <- list(axis1_full)
  samples[["processing"]][["interpolation"]] <- TRUE
  samples
}

nmr_interpolate_2D_legacy <- function(samples,
                                      axis1=c(min = 0.4, max = 10, by = 0.0008),
                                      axis2=NULL) {
  # Check if we can interpolate:
  verify_dimensionality(samples, valid_dimensions = 2)
  
  # 2. Check that we have the interpolation axis
  axis1 <- verify_axisn(axis1, samples[["axis"]][[1]][[1]])
  
  # Do interpolation:
  data_fields <- names(samples)[grepl(pattern = "^data_.*", x = names(samples))]
  # axis2 will be based on sample 1, dimension 2
  axis2 <- verify_axisn(axis2, samples[["axis"]][[1]][[2]])
  
  axis1_full <- seq(from = axis1["min"], to = axis1["max"], by = axis1["by"])
  axis2_full <- seq(from = axis2["min"], to = axis2["max"], by = axis2["by"])
  
  for (data_field in data_fields) {
    if (show_progress_bar()) {
      message("Interpolating ", data_field, "...")
    }
    data_matr <- interpolate_2d(list_of_f1 = purrr::map(samples[["axis"]], 1),
                                list_of_f2 = purrr::map(samples[["axis"]], 2),
                                list_of_data = samples[[data_field]],
                                f1_axis = axis1_full,
                                f2_axis = axis2_full)
    samples[[data_field]] <- data_matr
  }
  samples[["axis"]] <- list(axis1_full, axis2_full)
  samples[["processing"]][["interpolation"]] <- TRUE
  samples
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

interpolate_2d <- function(list_of_f1, list_of_f2, list_of_data,
                           f1_axis, f2_axis) {
  # For 2-D samples, we use two 1-D linear interpolations instead of a bilinear one. (FIXME)
  # Akima needs similar axis for interpolation. We will scale the axes.
  # rescale <- function(x0, y0, x1, y1) {
  #   slope <- (y1 - y0)/(x1 - x0)
  #   intercept <- y0 - slope * x0
  #   rescale_fun <- function(x) {
  #     return(slope*x + intercept)
  #   }
  #   return(rescale_fun)
  # }
  num_samples <- length(list_of_data)
  data_matr <- array(data = NA, dim = c(num_samples, length(f1_axis), length(f2_axis)))
  tryCatch({
    pb <- NULL
    if (show_progress_bar()) {
      pb <- utils::txtProgressBar(min = 0, max = 2*num_samples, style = 3)
    }
    for (i in seq_len(num_samples)) {
      data_interp_1ax <- apply(list_of_data[[i]], 2,
                               function(col) {
                                 signal::interp1(x = list_of_f1[[i]],
                                                 y = col,
                                                 xi = f1_axis,
                                                 method = "spline")
                               })
      if (!is.null(pb)) {
        utils::setTxtProgressBar(pb, 2*i - 1)
      }
      data_interp_1ax <- t(apply(data_interp_1ax, 1,
                                 function(row) {
                                   signal::interp1(x = list_of_f2[[i]],
                                                   y = row,
                                                   xi = f2_axis,
                                                   method = "spline")
                                 }))
      if (!is.null(pb)) {
        utils::setTxtProgressBar(pb, 2*i)
      }
      data_matr[i, , ] <- data_interp_1ax
      # I could not make this work with the time I had
      # orig_axis_i_1_scaled <- axis1_scalefun(orig_axis[[i]][[1]])
      # orig_axis_i_2_scaled <- axis2_scalefun(orig_axis[[i]][[2]])
      #
      # x <- rep(orig_axis_i_1_scaled, times = length(orig_axis_i_2_scaled))
      # y <- rep(orig_axis_i_2_scaled, each = length(orig_axis_i_1_scaled))
      # z <- samples[[data_field]][[i]]
      # z <- as.vector(z)
      # z0 <- akima::interp(x = x, y = y, z = z,
      #                     xo = axis1_scaled, yo = axis2_scaled,
      #                     linear = TRUE)
      # data_matr[i, , ] <- z0$z
    }
  }, finally = {
    if (!is.null(pb)) {
      close(pb)
    }
  })
  data_matr
}