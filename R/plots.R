#' Plot an nmr_dataset_1D
#' @family plotting functions
#' @param x a [nmr_dataset_1D] object
#' @param chemshift_range range of the chemical shifts to be included. Can be of length 3
#'                to include the resolution in the third element (e.g. `c(0.2, 0.8, 0.005)`)
#' @param NMRExperiment A character vector with the NMRExperiments to include. Use "all" to include all experiments.
#' @param quantile_plot If `TRUE` plot the 10\%, 50\%, 90\% percentiles of the spectra as reference.
#'                                            If two numbers between 0 and 1 are given then a custom percentile can be plotted
#' @param quantile_colors A vector with the colors for each of the quantiles
#' @param ... arguments passed to [ggplot2::aes_string].
#' @param interactive if `TRUE` return an interactive plotly plot, otherwise return a ggplot one.
#' @return The plot
#' @importFrom graphics plot
#' @export
#' @examples 
#' dir_to_demo_dataset <- system.file("dataset-demo", package = "AlpsNMR")
#' #dataset <- nmr_read_samples_dir(dir_to_demo_dataset)
#' #dataset_1D <- nmr_interpolate_1D(dataset, axis = c(min = -0.5, max = 10, by = 2.3E-4))
#' #plot(dataset_1D)
#' 
plot.nmr_dataset_1D <- function(x,
                                NMRExperiment = NULL,
                                chemshift_range = NULL,
                                interactive = FALSE,
                                quantile_plot = NULL,
                                quantile_colors = NULL,
                                ...) {
    if (interactive) {
        require_pkgs("plotly", msgs = c("i" = "Otherwise, you can set interactive=FALSE."))
    }
    if (is.null(chemshift_range)) {
        chemshift_range <- range(x$axis)
    } else if (length(chemshift_range) == 2) {
        chemshift_range <- range(chemshift_range)
    } else if (length(chemshift_range) == 3) {
        chemshift_range <-
            c(range(chemshift_range[seq_len(2)]), chemshift_range[3])
    } else {
        stop("chemshift_range should be a numeric vector of length 2 or 3.")
    }
    
    if (is.null(NMRExperiment)) {
        if (x$num_samples > 20) {
            NMRExperiment <-
                sample(nmr_meta_get_column(x, "NMRExperiment"), size = 10)
        } else {
            NMRExperiment <- nmr_meta_get_column(x, "NMRExperiment")
        }
    } else if (identical(NMRExperiment, "all")) {
        NMRExperiment <- nmr_meta_get_column(x, "NMRExperiment")
    }
    sample_idx <-
        which(nmr_meta_get_column(x, "NMRExperiment") %in% NMRExperiment)
    
    longdf <- tidy(
        x,
        sample_idx = sample_idx,
        chemshift_range = chemshift_range
    )
    fixed_aes <-
        list(x = "chemshift", y = "intensity", group = "NMRExperiment")
    dotdotdot_aes <- list(...)
    all_aes <- c(fixed_aes, dotdotdot_aes)
    if (!"color" %in% names(all_aes) && !"colour" %in% names(all_aes)) {
        all_aes <- c(all_aes, list(color = "NMRExperiment"))
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
                tmp <-
                    grDevices::gray.colors(length(quantile_plot) / 2,
                                           start = 0.5,
                                           end = 0.8)
                quantile_colors <- c(rev(tmp), tmp)
            } else {
                tmp <-
                    grDevices::gray.colors(ceiling(length(quantile_plot) / 2),
                                           start = 0.5,
                                           end = 0.8)
                quantile_colors <- c(rev(tmp[2:length(tmp)]), tmp)
            }
        }
        stopifnot(length(quantile_plot) == length(quantile_colors))
        decimate_qspectra <- decimate_axis(xaxis = x$axis,
                                           xrange = chemshift_range)
        q_spectra <-
            apply(x$data_1r[, decimate_qspectra], 2, function(x)
                stats::quantile(x, quantile_plot))
        linera <-
            tibble::tibble(
                NMRExperiment = rep(
                    paste0("Quantile ", 100 * quantile_plot, "%"),
                    each = sum(decimate_qspectra)
                ),
                color = rep(quantile_colors, each = sum(decimate_qspectra)),
                chemshift = rep(x$axis[decimate_qspectra],
                                times = length(quantile_plot)),
                intensity = as.numeric(t(q_spectra))
            )
    }
    
    gplt <- ggplot2::ggplot(longdf)
    
    if (!is.null(quantile_plot)) {
        gplt <- gplt +
            ggplot2::geom_line(
                data = linera,
                ggplot2::aes_string(
                    x = "chemshift",
                    y = "intensity",
                    group = "NMRExperiment"
                ),
                color = linera$color,
                # out of aes so it does not show up in the legend
                size = 1,
                linetype = "dashed"
            )
    }
    
    gplt <- gplt +
        ggplot2::geom_line(do.call(ggplot2::aes_string, all_aes)) +
        ggplot2::labs(x = "Chemical Shift (ppm)", y = "Intensity (a.u.)") +
        ggplot2::scale_x_reverse(limits = rev(chemshift_range[seq_len(2)])) +
        ggplot2::scale_y_continuous(labels = scales::label_number(scale_cut = scales::cut_si("")))

    
    if (interactive) {
        output <- plotly::ggplotly(gplt)
    } else {
        output <- gplt
    }
    return(output)
}

#' Get a tidy data frame from nmr_data object
#'
#' This dataframe is useful for plotting with ggplot, although it may be very
#' long and therefore use a lot of RAM.
#'
#' @param nmr_data a \code{\link{nmr_dataset}} object
#' @param sample_idx numeric vector with the sample indices to include
#' @param columns A character vector with the metadata columns to get, by default gets all of them
#' @param matrix_name A string with the matrix name, typically "data_1r"
#' @param axis_name A string with the axis name, for now "axis" is the only valid option
#' @inheritParams plot.nmr_dataset_1D
#' @importFrom generics tidy
#' @export
#' @examples 
#' dir_to_demo_dataset <- system.file("dataset-demo", package = "AlpsNMR")
#' dataset <- nmr_read_samples_dir(dir_to_demo_dataset)
#' dataset_1D <- nmr_interpolate_1D(dataset, axis = c(min = -1.0, max = 1.6, by = 2.3E-4))
#' dummy_metadata <- system.file("dataset-demo", "dummy_metadata.xlsx", package = "AlpsNMR")
#' NMRExp_SubjID <- readxl::read_excel(dummy_metadata, sheet = 1)
#' dataset_1D <- nmr_meta_add(dataset_1D, NMRExp_SubjID)
#' df_for_ggplot <- tidy(dataset_1D, chemshift_range = c(1.2, 1.4), columns = "SubjectID")
#' head(df_for_ggplot)
tidy.nmr_dataset_1D <-
    function(x,
             sample_idx = NULL,
             chemshift_range = NULL,
             columns = NULL,
             matrix_name = "data_1r",
             axis_name = "axis",
             ...
             ) {
        if (is.null(sample_idx)) {
            sample_idx <- seq_len(x$num_samples)
        }
        chemshift_in_range <- decimate_axis(xaxis = x[[axis_name]],
                                            xrange = chemshift_range)
        meta_df <- nmr_meta_get(x, columns = columns)
        NMRExperiments <- meta_df$NMRExperiment[sample_idx]
        chemshifts <- x[[axis_name]][chemshift_in_range]
        raw_data <-
            reshape2::melt(x[[matrix_name]][sample_idx, chemshift_in_range, drop = FALSE])
        raw_data$Var1 <- NMRExperiments[raw_data$Var1]
        raw_data$Var2 <- chemshifts[raw_data$Var2]
        colnames(raw_data) <-
            c("NMRExperiment", "chemshift", "intensity")
        result <-
            dplyr::left_join(raw_data, meta_df, by = "NMRExperiment")
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
    xrange_minmax <- range(xrange[seq_len(2)])
    x_in_range <-
        xaxis >= xrange_minmax[1] & xaxis <= xrange_minmax[2]
    points_in_range <- sum(x_in_range)
    if (length(xrange) == 3 && !is.na(xrange[3])) {
        max_points <- ceiling(diff(xrange_minmax) / xrange[3])
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

#' Plot a dataset into a HTML file
#'
#' Uses WebGL for performance
#'
#' @family plotting nmr datasets
#' @param nmr_dataset An [nmr_dataset_1D]
#' @param html_filename The output HTML filename to be created
#' @inheritDotParams plot.nmr_dataset_1D
#'
#' @return the html filename created
#' @export
#' @examples 
#' dir_to_demo_dataset <- system.file("dataset-demo", package = "AlpsNMR")
#' #dataset <- nmr_read_samples_dir(dir_to_demo_dataset)
#' #dataset_1D <- nmr_interpolate_1D(dataset, axis = c(min = -0.5, max = 10, by = 2.3E-4))
#' #html_plot <- plot_webgl(dataset_1D, "html_plot.html")
#' 
plot_webgl <- function(nmr_dataset, html_filename, ...) {
    plt <- plot(nmr_dataset, ...)
    plot_interactive(plt = plt, html_filename = html_filename)
    html_filename
}

#' Plots in WebGL
#'
#' @family plotting functions
#' @param plt A plot created with plotly or ggplot2
#' @param html_filename The file name where the plot will be saved
#' @param overwrite Overwrite the lib/ directory (use `NULL` to prompt the user)
#'
#' @return The html_filename
#' @export
#' @examples 
#' dir_to_demo_dataset <- system.file("dataset-demo", package = "AlpsNMR")
#' dataset <- nmr_read_samples_dir(dir_to_demo_dataset)
#' dataset_1D <- nmr_interpolate_1D(dataset, axis = c(min = -0.5, max = 10, by = 2.3E-4))
#' # plot <- plot(dataset_1D)
#' # html_plot_interactive <- plot_interactive(plot, "html_plot_interactive.html")
#' 
plot_interactive <- function(plt, html_filename, overwrite = NULL) {
    #Check if lib folder exists
    basedir <- dirname(html_filename)
    libdir <- file.path(basedir, "lib")
    libdir_exists <- dir.exists(libdir)
    if (is.null(overwrite) && libdir_exists) {
        if (interactive()) {
            # warning user before some contents of lib folder could be destroyed
            rlang::inform("{libdir} folder already exists, plot_interactive will replace it. Continue? [y/n]:")
            response <- scan("stdin", character(), n=1)
            overwrite <- response %in% c("y", "Y")
        } else {
            overwrite <- FALSE
        }
    }
    if (libdir_exists && isFALSE(overwrite)) {
        rlang::abort(message = c("plot_interactive aborted", "x" = "{libdir} folder already exists, use overwrite=TRUE"))
    }
    suppressMessages(utils::capture.output({
        htmltools::save_html(
            html = htmltools::as.tags(x = plotly::toWebGL(plt),
                                      standalone = TRUE),
            file = html_filename
        )
    }))
    html_filename
}
