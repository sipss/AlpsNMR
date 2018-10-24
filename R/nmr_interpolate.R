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
                            axis1=c(min = 0.4, max = 10, by = 0.0008),
                            axis2=NULL) {
  # Check if we can interpolate:
  
  # 1. Check that each of the loaded samples has the same dimensionality
  dimensions_per_sample <- vapply(samples$axis, length, numeric(1))
  if (any(dimensions_per_sample != dimensions_per_sample[1])) {
    # Have different dimensionality, can't merge
    stop("Samples have different dimensionality. Cant't interpolate.")
  }
  dimensions_per_sample <- dimensions_per_sample[1]
  if (!dimensions_per_sample %in% c(1, 2)) {
    stop("Interpolation not implemented for dimensions different than 1 or 2.")
  }
  
  # 2. Check that we have the interpolation axis
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
  if (dimensions_per_sample >= 1) {
    # axis1 will be based on sample 1, dimension 1
    axis1 <- verify_axisn(axis1, samples[["axis"]][[1]][[1]])
  }
  if (dimensions_per_sample >= 2) {
    # axis2 will be based on sample 1, dimension 2
    axis2 <- verify_axisn(axis2, samples[["axis"]][[1]][[2]])
  }
  if (dimensions_per_sample >= 3) {
    warning("Interpolation not implemented for dimensions > 2. Skipping interpolation.")
    return(samples)
  }
  
  # Do interpolation:
  orig_axis <- samples[["axis"]]
  data_fields <- names(samples)[grepl(pattern = "^data_.*", x = names(samples))]
  num_samples <- samples[["num_samples"]]
  if (dimensions_per_sample == 1) {
    axis1_full <- seq(from = axis1["min"], to = axis1["max"], by = axis1["by"])
    samples[["axis"]] <- list(axis1_full)
    for (data_field in data_fields) {
      if (show_progress_bar(samples$num_samples > 5)) {
        message("Interpolating ", data_field, "...")
      }
      data_matr <- matrix(NA, nrow = num_samples, ncol = length(axis1_full))
      tryCatch({
        pb <- NULL
        if (show_progress_bar(samples$num_samples > 5)) {
          pb <- utils::txtProgressBar(min = 0, max = samples$num_samples, style = 3)
        }
        for (i in seq_len(samples$num_samples)) {
          data_matr[i, ] <- signal::interp1(x = orig_axis[[i]][[1]],
                                            y = samples[[data_field]][[i]],
                                            xi = axis1_full, method = "spline")
          if (!is.null(pb)) {
            utils::setTxtProgressBar(pb, i)
          }
        }
      }, finally = {
        if (!is.null(pb)) {
          close(pb)
        }
      })
      samples[[data_field]] <- data_matr
    }
  } else if (dimensions_per_sample == 2) {
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
    axis1_full <- seq(from = axis1["min"], to = axis1["max"], by = axis1["by"])
    axis2_full <- seq(from = axis2["min"], to = axis2["max"], by = axis2["by"])
    # axis1_scalefun <- rescale(axis1["min"], 0, axis1["max"], 1)
    # axis1_scaled <- axis1_scalefun(axis1_full)
    # axis2_scalefun <- rescale(axis2["min"], 0, axis2["max"], axis1["by"]/axis2["by"])
    # axis2_scaled <- axis2_scalefun(axis2_full)
    
    for (data_field in data_fields) {
      if (show_progress_bar()) {
        message("Interpolating ", data_field, "...")
      }
      data_matr <- array(data = NA, dim = c(num_samples, length(axis1_full), length(axis2_full)))
      tryCatch({
        pb <- NULL
        if (show_progress_bar(samples$num_samples > 0)) {
          pb <- utils::txtProgressBar(min = 0, max = 2*samples$num_samples, style = 3)
        }
        for (i in seq_len(samples$num_samples)) {
          data_interp_1ax <- apply(samples[[data_field]][[i]], 2,
                                   function(col) {
                                     signal::interp1(x = orig_axis[[i]][[1]],
                                                     y = col,
                                                     xi = axis1_full,
                                                     method = "spline")
                                   })
          if (!is.null(pb)) {
            utils::setTxtProgressBar(pb, 2*i - 1)
          }
          data_interp_1ax <- t(apply(data_interp_1ax, 1,
                                     function(row) {
                                       signal::interp1(x = orig_axis[[i]][[2]],
                                                       y = row,
                                                       xi = axis2_full,
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
      samples[["axis"]] <- list(axis1_full, axis2_full)
      samples[[data_field]] <- data_matr
    }
  }
  samples[["processing"]][["interpolation"]] <- TRUE
  return(samples)
}
