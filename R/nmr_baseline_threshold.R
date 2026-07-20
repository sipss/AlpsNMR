#' Threshold estimation for peak detection
#'
#' Estimates the threshold value for peak detection on an [nmr_dataset_1D] object by examining
#' a range without peaks that you must provide (there is no ppm range guaranteed to be free of
#' peaks for every sample type).
#'
#' Two methods can be used:
#'
#' - "mean3sd": The mean3sd method computes the mean and the standard deviation of each spectrum
#' in the given range. The mean spectrum and the mean standard deviation are both vectors
#' of length equal to the number of points in the given range. The mean of the mean spectrum
#  and the mean of the standard deviations are used to summarize the center and dispersion of
#' the noise. The threshold is defined as `center + 3 dispersion`, and it is one single threshold
#' for the whole dataset. This is the default for backwards compatibility.
#'
#' - "median3mad": First we take the data matrix. If we have estimated a baseline already,
#'   subtract it. In the defined region without peaks, estimate the median of each sample and
#'   its median absolute deviation. Return a vector of length equal to the number of samples
#'   with the `median+3mad` for each sample. This is a new more robust method.
#'
#' @family peak detection functions
#' @param nmr_dataset An [nmr_dataset_1D].
#' @param method Either "mean3sd" or the more robust "median3mad". See the details.
#' @param range_without_peaks A vector with two doubles describing a range without peaks suitable for baseline detection.
#' There is no such a range that works for every sample type, so you must inspect your spectra and provide one.
#' @return Numerical. A threshold value in intensity below that no peak is detected.
#' @export
#' @examples
#' ppm_axis <- seq(from = 0, to = 10, length.out = 1000)
#' data_1r <- matrix(runif(1000, 0, 10), nrow = 1) + 100
#' dataset_1D <- new_nmr_dataset_1D(
#'     ppm_axis = ppm_axis,
#'     data_1r = data_1r,
#'     metadata = list(external=data.frame(NMRExperiment = "10"))
#' )
#' bl_threshold <- nmr_baseline_threshold(dataset_1D, range_without_peaks = c(9.5,10))
#'
nmr_baseline_threshold <- function(nmr_dataset, range_without_peaks = NULL, method = c("mean3sd", "median3mad")) {
    # FIXME: Maybe a whole baseline would be better, so we can cope with slowly changing baselines better
    method <- match.arg(method)
    if (is.null(range_without_peaks)) {
        cli::cli_abort(
            message = c(
                "range_without_peaks must be given",
                "i" = "There is no ppm range guaranteed to be free of peaks for every sample type.",
                "i" = "Inspect your spectra and pass a two-value ppm range (e.g. c(9.5, 10)) known to be free of peaks for your samples."
            )
        )
    }
    if (length(range_without_peaks) != 2) {
        cli::cli_abort("range_without_peaks must have length 2")
    }
    r_start <- min(range_without_peaks)
    r_end <- max(range_without_peaks)
    threshold_ind <- nmr_dataset$axis >= r_start & nmr_dataset$axis < r_end
    if (sum(threshold_ind) < 10) {
        cli::cli_abort(
            message = c(
                "Can't estimate a baseline threshold reliably",
                "i" = glue("The selected range_without_peaks [{r_start},{r_end}] ppm contains {sum(threshold_ind)} points."),
                "i" = glue("Either change the interpolation axis limits or choose a different range here")
            )
        )
    }
    if (method == "mean3sd") {
        if (nmr_dataset$num_samples > 1) {
            cent <- mean(apply(nmr_dataset$data_1r[, threshold_ind, drop = FALSE], 2, mean))
            disp <- mean(apply(nmr_dataset$data_1r[, threshold_ind, drop = FALSE], 2, stats::sd))
        } else {
            cent <- mean(as.numeric(nmr_dataset$data_1r[, threshold_ind]))
            disp <- stats::sd(as.numeric(nmr_dataset$data_1r[, threshold_ind]))
        }
        return(cent + 3 * disp)
    } else if (method == "median3mad") {
        out <- rep(NA_real_, nmr_dataset$num_samples)
        for (i in seq_len(nmr_dataset$num_samples)) {
            if ("data_1r_baseline" %in% names(unclass(nmr_dataset))) {
                spec_region <- nmr_dataset$data_1r[i, threshold_ind, drop = FALSE] - nmr_dataset$data_1r_baseline[i, threshold_ind]
            } else {
                spec_region <- nmr_dataset$data_1r[i, threshold_ind, drop = FALSE]
            }
            out[i] <- stats::median(spec_region) + 3 * stats::mad(spec_region)
        }
        names(out) <- names(nmr_dataset)
        return(out)
    } else {
        stop("Unexpected method", method)
    }
}


#' Plot the baseline thresholds
#'
#' With many samples, a single page can't legibly show one facet per sample;
#' `nrow`/`ncol`/`page` paginate the facets instead of cramming (or silently
#' subsampling) them all onto one page.
#'
#' @inheritParams plot.nmr_dataset_1D
#' @param nmr_dataset An [nmr_dataset_1D] object
#' @param thresholds A named vector. The values are baseline thresholds. The names are NMRExperiments.
#' @param NMRExperiment The NMRExperiments to plot. `NULL` (the default) plots every sample
#'    (paginated via `nrow`/`ncol`/`page`); `"all"` is a synonym for `NULL`; or pass a character
#'    vector to filter to a specific subset of samples.
#' @param chemshift_range The range to plot, as a first check use the `range_without_peaks` from [nmr_baseline_threshold]
#' @param nrow,ncol Number of rows/columns of facets per page. `NULL` (the default) picks a
#'    snug grid for the number of samples requested: 1x`n` for fewer than 4 samples, 2x2 for 4,
#'    2x3 for 5-6, and a fixed 3x3 (paginated via `page`) for 7 or more.
#' @param page Which page of facets to plot (1-indexed). Requesting a page beyond the number
#'    available is an error.
#'
#' @return A plot.
#' @export
#' @examples
#'
#' ppm_axis <- seq(from = 0, to = 10, length.out = 1000)
#' data_1r <- matrix(runif(1000, 0, 10), nrow = 1) + 100
#' dataset_1D <- new_nmr_dataset_1D(
#'     ppm_axis = ppm_axis,
#'     data_1r = data_1r,
#'     metadata = list(external=data.frame(NMRExperiment = "10"))
#' )
#' bl_threshold <- nmr_baseline_threshold(dataset_1D, range_without_peaks = c(9.5,10))
#' nmr_baseline_threshold_plot(dataset_1D, bl_threshold, chemshift_range = c(9.5, 10))
nmr_baseline_threshold_plot <- function(nmr_dataset, thresholds, NMRExperiment = NULL, chemshift_range = NULL,
    nrow = NULL, ncol = NULL, page = 1, ...) {
    if (is.null(chemshift_range)) {
        cli::cli_abort(
            message = c(
                "chemshift_range must be given",
                "i" = "There is no ppm range guaranteed to be free of peaks for every sample type.",
                "i" = "Pass the same two-value ppm range you used as range_without_peaks in nmr_baseline_threshold()."
            )
        )
    }
    if (is.null(NMRExperiment) || identical(NMRExperiment, "all")) {
        NMRExperiment <- names(nmr_dataset)
    }
    num_samples <- length(NMRExperiment)
    if (is.null(nrow) && is.null(ncol)) {
        if (num_samples >= 7) {
            nrow <- 3
            ncol <- 3
        } else if (num_samples %in% c(5, 6)) {
            nrow <- 2
            ncol <- 3
        } else if (num_samples == 4) {
            nrow <- 2
            ncol <- 2
        } else {
            nrow <- 1
            ncol <- max(num_samples, 1)
        }
    }
    samples_per_page <- nrow * ncol
    total_pages <- max(ceiling(num_samples / samples_per_page), 1)
    if (page < 1 || page > total_pages) {
        cli::cli_abort(
            message = c(
                "page ({page}) is out of bounds",
                "i" = "There {cli::qty(total_pages)} {?is/are} {total_pages} page{?s} for {num_samples} sample{?s} with nrow = {nrow}, ncol = {ncol} ({samples_per_page} per page)."
            )
        )
    }
    page_start <- (page - 1) * samples_per_page + 1
    page_end <- min(page * samples_per_page, num_samples)
    NMRExperiment <- NMRExperiment[page_start:page_end]
    if (length(thresholds) == 1L) {
        thresholds <- rep(thresholds, length = length(NMRExperiment))
        names(thresholds) <- NMRExperiment
    }
    thresholds <- thresholds[NMRExperiment]

    is_aes_string <- is_using_aes_string(...)
    
    if (is_aes_string) {
        cli::cli_warn(
            c(
                "!" = "Passing aes_string arguments to nmr_baseline_threshold_plot(nmr_dataset, ...) is deprecated.",
                "i" = "Please pass aes arguments instead"
            ),
            .frequency = "regularly",
            .frequency_id = "nmr_baseline_threshold_plot_plotting_with_aes_string",
        )
        aes_str <- as.character(list(...))
        columns_to_request <- c("NMRExperiment", get_vars_from_aes_string(aes_str))
    } else {
        columns_to_request <- c("NMRExperiment", get_vars_from_aes(...))
    }
    tidy_data <- tidy_spectra_baseline_and_threshold(
        dataset = nmr_dataset,
        thresholds = thresholds,
        chemshift_range = chemshift_range,
        NMRExperiment = NMRExperiment,
        columns = columns_to_request
    )
    to_plot <- tidy_data$spectra
    to_plot_baseline <- tidy_data$baselines
    to_plot_threshold <- tidy_data$thresholds

    ymax <- 1.5 * max(to_plot_threshold$intensity)

    if (is_aes_string) {
        return(
            nmr_baseline_threshold_plot_aes_string(
                to_plot,
                to_plot_baseline,
                to_plot_threshold,
                chemshift_range,
                ymax,
                NMRExperiment,
                nrow = nrow,
                ncol = ncol,
                ...
            )
        )
    }
    dots_aes_args <- prepare_aes(...)
    
    gplt <- ggplot2::ggplot() +
        # The spectra:
        ggplot2::geom_line(mapping = ggplot2::aes(!!!dots_aes_args), data = to_plot)
    
    if (!is.null(to_plot_baseline)) {
        # The baseline:
        gplt <- gplt +
            ggplot2::geom_line(
                mapping = ggplot2::aes(!!!dots_aes_args),
                data = to_plot_baseline,
                linetype = "dashed"
            )
    }

    gplt <- gplt +
        # The threshold
        ggplot2::geom_line(
            mapping = ggplot2::aes(!!!dots_aes_args),
            data = to_plot_threshold,
            linetype = "dashed", 
            color = "black"
        ) +
        # Other plotting options
        ggplot2::labs(x = "Chemical Shift (ppm)", y = "Intensity (a.u.)") +
        ggplot2::scale_x_reverse(limits = rev(chemshift_range[seq_len(2)])) +
        ggplot2::scale_y_continuous(labels = scales::label_number(scale_cut = scales::cut_si("")), limits = c(0, ymax)) +
        ggplot2::facet_wrap(~ factor(NMRExperiment, levels = unique(NMRExperiment)), nrow = nrow, ncol = ncol) +
        ggplot2::theme(legend.position = "none")
    gplt
}

# deprecated
nmr_baseline_threshold_plot_aes_string <- function(to_plot, to_plot_baseline, to_plot_threshold, chemshift_range, ymax, NMRExperiment, nrow = NULL, ncol = NULL, ...) {
    dotdotdot_aes <- list(...)
    fixed_aes <- list(
        x = "chemshift",
        y = "intensity",
        group = "NMRExperiment"
    )
    all_aes <- c(fixed_aes, dotdotdot_aes)
    if (!"color" %in% names(all_aes) && !"colour" %in% names(all_aes)) {
        all_aes <- c(all_aes, list(color = "NMRExperiment"))
    }

    gplt <- ggplot2::ggplot() +
        # The spectra:
        ggplot2::geom_line(mapping = do.call(ggplot2::aes_string, all_aes), data = to_plot)
    
    if (!is.null(to_plot_baseline)) {
        # The baseline:
        gplt <- gplt + ggplot2::geom_line(mapping = do.call(ggplot2::aes_string, all_aes), data = to_plot_baseline, linetype = "dashed")
    }
    gplt <- gplt +
        # The threshold
        ggplot2::geom_line(mapping = do.call(ggplot2::aes_string, all_aes), data = to_plot_threshold, linetype = "dashed", color = "black") +
        # Other plotting options
        ggplot2::labs(x = "Chemical Shift (ppm)", y = "Intensity (a.u.)") +
        ggplot2::scale_x_reverse(limits = rev(chemshift_range[seq_len(2)])) +
        ggplot2::scale_y_continuous(labels = scales::label_number(scale_cut = scales::cut_si("")), limits = c(0, ymax)) +
        ggplot2::facet_wrap(~ factor(NMRExperiment, levels = unique(NMRExperiment)), nrow = nrow, ncol = ncol) +
        ggplot2::theme(legend.position = "none")
    gplt
}

tidy_spectra_baseline_and_threshold <- function(dataset, thresholds, chemshift_range, NMRExperiment, columns = character(0L)) {
    to_plot <- tidy(
        dataset,
        chemshift_range = chemshift_range,
        NMRExperiment = NMRExperiment,
        columns = columns,
        matrix_name = "data_1r"
    )
    if ("data_1r_baseline" %in% names(unclass(dataset))) {
        to_plot_baseline <- tidy(
            dataset,
            chemshift_range = chemshift_range,
            NMRExperiment = NMRExperiment,
            columns = columns,
            matrix_name = "data_1r_baseline"
        )
        if (is.null(thresholds)) {
            to_plot_threshold <- NULL
        } else {
            to_plot_threshold <- dplyr::left_join(
                to_plot_baseline,
                tibble::enframe(
                    thresholds,
                    name = "NMRExperiment",
                    value = "threshold"
                ),
                by = "NMRExperiment"
            )
            to_plot_threshold$intensity <- to_plot_threshold$intensity + to_plot_threshold$threshold
        }
    } else {
        to_plot_baseline <- NULL
        if (is.null(thresholds)) {
            to_plot_threshold <- NULL
        } else {
            to_plot_threshold <- dplyr::left_join(
                to_plot,
                tibble::enframe(
                    thresholds,
                    name = "NMRExperiment",
                    value = "threshold"
                ),
                by = "NMRExperiment"
            )
            to_plot_threshold$intensity <- to_plot_threshold$threshold
        }
    }
    list(
        spectra = to_plot,
        baselines = to_plot_baseline,
        thresholds = to_plot_threshold
    )
}
