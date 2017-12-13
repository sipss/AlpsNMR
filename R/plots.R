#' Plot an nmr_dataset
#' @param x a \code{\link{nmr_dataset}} object
#' @inheritParams nmr_get_long_df
#' @param sample_idx numeric vector with the sample indices to include. Use "all" to include all
#' @param quantile_plot If \code{TRUE} plot the 10\%, 50\%, 90\% percentiles of the spectra as reference.
#'                      If two numbers between 0 and 1 are given then a custom percentile can be plotted
#' @param ... arguments passed to ggplot2::aes.
#' @param interactive if \code{TRUE} return an interactive plotly plot, otherwise return a ggplot one.
#' @return The plot
#' @method plot nmr_dataset
#' @export
plot.nmr_dataset <- function(x, sample_idx = NULL,
                             chemshift_range = NULL,
                             quantile_plot = FALSE,
                             interactive = FALSE, ...) {

  dimension <- unique(x$metadata$info_dimension)
  if (length(dimension) > 1) {
    stop("Samples with different dimensionality")
  }

  if (dimension == 1) {
    output <- plot_nmr_dataset_1D(x = x, sample_idx = sample_idx,
                                  chemshift_range = chemshift_range,
                                  interactive = interactive,
                                  quantile_plot = quantile_plot, ...)
  } else if (dimension == 2) {
    output <- plot_nmr_dataset_2D(x = x, sample_idx = sample_idx,
                                  chemshift_range = chemshift_range,
                                  interactive = interactive, ...)
  } else {
    stop("Can't plot samples of dimension ", as.character(dimension))
  }
  return(output)
}

plot_nmr_dataset_1D <- function(x, sample_idx = NULL,
                                chemshift_range = NULL,
                                interactive = FALSE,
                                quantile_plot = NULL,
                                quantile_colors = NULL,
                                ...) {
  if (!isTRUE(x$processing$data_loaded)) {
    stop("Can't plot an nmr_dataset without data")
  }
  if (x$processing$interpolation == FALSE) {
    stop("The dataset needs to be interpolated. Use nmr_interpolate for that.")
  }

  if (is.null(chemshift_range)) {
    chemshift_range <- range(x$axis[[1]])
  }

  if (is.null(sample_idx)) {
    if (x$num_samples > 20) {
      sample_idx <- sample(1:x$num_samples, size = 10)
    } else {
      sample_idx <- 1:x$num_samples
    }
  } else if (any(sample_idx == "all")) {
    sample_idx <- 1:x$num_samples
  }
  longdf <- nmr_get_long_df(nmr_data = x, sample_idx = sample_idx,
                            chemshift_range = chemshift_range)
  fixed_aes <- list(x = "chemshift", y = "intensity", group = "injection_id")
  dotdotdot_aes <- list(...)
  all_aes <- c(fixed_aes, dotdotdot_aes)
  if (!"color" %in% names(all_aes)) {
    all_aes <- c(all_aes, list(color = "injection_id"))
  }

  linera <- NULL
  if (is.null(quantile_plot) || identical(quantile_plot, FALSE)) {
    quantile_plot <- NULL
  } else {
    if (isTRUE(quantile_plot)) {
      quantile_plot <- c(0.1, 0.5, 0.9)
    }
    if (is.null(quantile_colors)) {
      if (length(quantile_plot) %% 2 == 0) {
        # even number of quantiles.
        tmp <- grDevices::gray.colors(length(quantile_plot)/2, start = 0.5, end = 0.8)
        quantile_colors <- c(rev(tmp), tmp)
      } else {
        tmp <- grDevices::gray.colors(ceiling(length(quantile_plot)/2), start = 0.5, end = 0.8)
        quantile_colors <- c(rev(tmp[2:length(tmp)]), tmp)
      }
    }
    stopifnot(length(quantile_plot) == length(quantile_colors))
    decimate_qspectra <- decimate_axis(xaxis = x$axis[[1]],
                                       xrange = chemshift_range)
    q_spectra <- apply(x$data_1r[,decimate_qspectra], 2, function(x) stats::quantile(x, quantile_plot))
    linera <- tibble::tibble(injection_id = rep(paste0("Quantile ", 100*quantile_plot, "%"),
                                                each = sum(decimate_qspectra)),
                             color = rep(quantile_colors, each = sum(decimate_qspectra)),
                             chemshift = rep(x$axis[[1]][decimate_qspectra],
                                             times = length(quantile_plot)),
                             intensity = as.numeric(t(q_spectra)))
  }

  gplt <- ggplot2::ggplot(longdf)

  if (!is.null(quantile_plot)) {
    gplt <- gplt +
      ggplot2::geom_line(data = linera,
                         ggplot2::aes_string(x = "chemshift", y = "intensity", group = "injection_id"),
                         color = linera$color, # out of aes so it does not show up in the legend
                         size = 1, linetype = "dashed")
  }

  gplt <- gplt +
    ggplot2::geom_line(do.call(ggplot2::aes_string, all_aes)) +
    ggplot2::scale_x_reverse(name = "Chemical Shift (ppm)", limits = rev(chemshift_range[1:2])) +
    ggplot2::scale_y_continuous(name = "Intensity (a.u.)")

  if (interactive) {
    if (!requireNamespace("plotly", quietly = TRUE)) {
      stop("plotly needed for this plot to work. Please install it or use `interactive = FALSE`.",
           call. = FALSE)
    }
    output <- plotly::ggplotly(gplt)
  } else {
    output <- gplt
  }
  return(output)
}

plot_nmr_dataset_2D <- function(x, sample_idx = NULL,
                                chemshift_range = NULL,
                                interactive = FALSE, ...) {
  if (!isTRUE(x$processing$data_loaded)) {
    stop("Can't plot an nmr_dataset without data")
  }
  if (x$processing$interpolation == FALSE) {
    stop("The dataset needs to be interpolated. Use nmr_interpolate for that.")
  }
  if (is.null(sample_idx)) {
    sample_idx <- sample(1:x$num_samples, size = 1)
  }
  if (length(sample_idx) > 1) {
    warning("In 2D only one sample at a time can be plotted")
    sample_idx <- sample_idx[1]
  }
  if (sample_idx == "all") {
    warning("Plotting sample_idx = 1")
    sample_idx <- 1
  }

  if (!isTRUE(x$processing$data_loaded)) {
    stop("Can't plot an nmr_dataset without data")
  }
  if (x$processing$interpolation == FALSE) {
    stop("The dataset needs to be interpolated. Use nmr_interpolate for that.")
  }
  if (is.null(chemshift_range)) {
    chemshift_range <- list()
  }
  if (!"axis1" %in% names(chemshift_range)) {
    chemshift_range[["axis1"]] <- range(x$axis[[1]])
  }
  if (!"axis2" %in% names(chemshift_range)) {
    chemshift_range[["axis2"]] <- range(x$axis[[2]])
  }

  chemshift_in_range1 <- decimate_axis(xaxis = x$axis[[1]],
                                       xrange = chemshift_range[["axis1"]])

  chemshift_in_range2 <- decimate_axis(xaxis = x$axis[[2]],
                                       xrange = chemshift_range[["axis2"]])


  data <- x$data_2rr[sample_idx, chemshift_in_range1, chemshift_in_range2]
  xaxis <- x$axis[[1]][chemshift_in_range1]
  yaxis <- x$axis[[2]][chemshift_in_range2]

  sample_long <- reshape2::melt(data)
  sample_long$Var1 <- xaxis[sample_long$Var1]
  sample_long$Var2 <- yaxis[sample_long$Var2]
  colnames(sample_long) <- c("axis1", "axis2", "intensity")

  gplt <- ggplot2::ggplot(sample_long, ggplot2::aes_string(x = "axis1", y = "axis2")) +
    ggplot2::geom_raster(ggplot2::aes_string(fill = "intensity")) +
    ggplot2::scale_fill_gradient(trans = "log10_minus_min") +
    ggplot2::ggtitle(paste("Sample:", x$metadata$injection_id[sample_idx]))

  if (interactive) {
    if (!requireNamespace("plotly", quietly = TRUE)) {
      stop("plotly needed for this plot to work. Please install it or use `interactive = FALSE`.",
           call. = FALSE)
    }
    output <- plotly::ggplotly(gplt)
  } else {
    output <- gplt
  }
  return(output)
}


#' Custom scale transformation for 2D plots of nmr_dataset objects
#' @export
log10_minus_min_trans = function() {
  num_breaks <- 5 # number of breaks in the plot
  log_base <- 10
  # First we define our transformation:
  kept_minx <- NULL
  transform <-  function(x) {
    # (Here x is our data)
    # This will store the minimum so we can compute the inverse later:
    kept_minx <<- min(x)
    return(log(x - kept_minx + 1, base = log_base))
  }
  # And the inverse transformation
  inverse <- function(y) {
    # Here y is the transformed data
    return(10 ^ y - 1 + kept_minx)
  }

  # This function will define the breaks:
  breaks <- function(x) {
    # here x will be the limits of the data to be plotted
    minx <- min(x, na.rm = TRUE)
    trans_x_lim <- log(range(x, na.rm = TRUE) - minx + 1, base = log_base)
    min_tr_x <- floor(trans_x_lim[1])
    max_tr_x <- ceiling(trans_x_lim[2])
    if (min_tr_x == max_tr_x) {
      return(log_base ^ min)
    }
    by <- floor((max_tr_x - min_tr_x)/num_breaks) + 1
    return(log_base ^ seq(min_tr_x, max_tr_x, by = by))
  }
  # Finally we create the scale:
  scales::trans_new(name = "log10_minus_min",
                    transform = transform,
                    inverse = inverse,
                    breaks = breaks)
}



#' Get a long data frame from nmr_data object
#'
#' This dataframe is useful for plotting with ggplot, although it may be very
#' long and therefore eat a lot of RAM.
#'
#' @param nmr_data a \code{\link{nmr_dataset}} object
#' @param sample_idx numeric vector with the sample indices to include
#' @param chemshift_range range of the chemical shifts to be included. Can be of length 3
#'        to include the resolution in the third element (e.g. \code{c(0.2, 0.8, 0.05)})
#' @export
nmr_get_long_df <- function(nmr_data, sample_idx = NULL, chemshift_range = NULL) {
  if (is.null(sample_idx)) {
    sample_idx <- 1:nmr_data$num_samples
  }
  chemshift_in_range <- decimate_axis(xaxis = nmr_data$axis[[1]],
                                      xrange = chemshift_range)
  meta_df <- nmr_get_metadata(nmr_data)
  injection_ids <- meta_df$injection_id[sample_idx]
  chemshifts <- nmr_data$axis[[1]][chemshift_in_range]
  raw_data <- reshape2::melt(nmr_data$data_1r[sample_idx, chemshift_in_range, drop = FALSE])
  raw_data$Var1 <- injection_ids[raw_data$Var1]
  raw_data$Var2 <- chemshifts[raw_data$Var2]
  colnames(raw_data) <- c("injection_id", "chemshift", "intensity")
  result <- dplyr::left_join(raw_data, meta_df, by = "injection_id")
  return(result)
}

decimate_axis <- function(xaxis, xrange = NULL) {
  if (is.null(xrange)) {
    xrange <- range(xaxis)
  }
  # make sure we have min-max
  if (is.na(xrange[1])) {
    xrange[1] <- min(xaxis)
  }
  if (is.na(xrange[2])) {
    xrange[2] <- max(xaxis)
  }
  xrange_minmax <- range(xrange[1:2])
  x_in_range <- xaxis >= xrange_minmax[1] & xaxis <= xrange_minmax[2]
  points_in_range <- sum(x_in_range)
  if (length(xrange) == 3 && !is.na(xrange[3])) {
    max_points <- ceiling(diff(xrange_minmax)/xrange[3])
  } else {
    max_points <- points_in_range
  }
  if (points_in_range > max_points) {
    decimate_factor <- ceiling(points_in_range / max_points)
  } else {
    decimate_factor <- 1
  }
  decimate_vector <- rep(c(TRUE, rep(FALSE, decimate_factor - 1)),
                         length.out = length(xaxis))
  decimated_axis_bool <- x_in_range & decimate_vector
  return(decimated_axis_bool)
}

