#' Peak detection for NMR
#'
#' @name Peak_detection
#' @examples
#' dir_to_demo_dataset <- system.file("dataset-demo", package = "AlpsNMR")
#' nmr_dataset <- nmr_read_samples_dir(dir_to_demo_dataset)
#' dataset_1D <- nmr_interpolate_1D(nmr_dataset, axis = c(min = -0.5, max = 10, by = 2.3E-4))
#'
#' sample_10 <- filter(dataset_1D, NMRExperiment == "10")
#'
#' # 1.Peak detection in the dataset.
#' peak_data <- nmr_detect_peaks(dataset_1D,
#'                               nDivRange_ppm = 0.1, # Size of detection segments
#'                               scales = seq(1, 16, 2),
#'                               baselineThresh = 0, # Minimum peak intensity
#'                               SNR.Th = 4) # Signal to noise ratio
#'
#' #nmr_detect_peaks_plot(sample_10, peak_data, "NMRExp_ref")
#' 
#' peaks_detected <- nmr_detect_peaks_tune_snr(
#'   sample_10, 
#'   SNR_thresholds = seq(from = 2, to = 3, by = 0.5),
#'   nDivRange_ppm = 0.03,
#'   scales = seq(1, 16, 2),
#'   baselineThresh = 0
#' )
#'
#' 
#' # 2.Find the reference spectrum to align with.
#' NMRExp_ref <- nmr_align_find_ref(dataset_1D, peak_data)
#'
#' # 3.Spectra alignment using the ref spectrum and a maximum alignment shift
#' nmr_dataset <- nmr_align(dataset_1D, # the dataset
#'                          peak_data, # detected peaks
#'                          NMRExp_ref = NMRExp_ref, # ref spectrum
#'                          maxShift_ppm = 0.0015, # max alignment shift
#'                          acceptLostPeak = FALSE) # lost peaks
#'                          
#' # 4.PEAK INTEGRATION (please, consider previous normalization step).
#' # First we take the peak table from the reference spectrum
#' peak_data_ref <- filter(peak_data, NMRExperiment == NMRExp_ref)
#'
#' # Then we integrate spectra considering the peaks from the ref spectrum
#' nmr_peak_table <- nmr_integrate_peak_positions(
#'                       samples = nmr_dataset,
#'                       peak_pos_ppm = peak_data_ref$ppm,
#'                       peak_width_ppm = NULL)
#'
#' validate_nmr_dataset_peak_table(nmr_peak_table)
#' 
#' #If you wanted the final peak table before machine learning you can run
#' nmr_peak_table_completed <- get_integration_with_metadata(nmr_peak_table)
#' 
NULL


#' Peak detection for NMR
#'
#' The function detects peaks on an [nmr_dataset_1D] object, using
#' [speaq::detectSpecPeaks]. `detectSpecPeaks` divides the whole spectra into
#' smaller segments and uses [MassSpecWavelet::peakDetectionCWT] for peak
#' detection.
#' 
#' Afterwards, the peak apex and the peak inflection points are used to efficiently
#' adjust a lorentzian to each peak, and compute the peak area and width, as well as
#' the error of the fit. These peak features can be used afterwards to reject false
#' detections.
#'
#' @family peak detection functions
#' @family nmr_dataset_1D functions
#' @seealso [nmr_align] for peak alignment with the detected peak table
#' @param nmr_dataset An [nmr_dataset_1D].
#' @param nDivRange_ppm Segment size, in ppms, to divide the spectra and search
#'     for peaks.
#' @param baselineThresh All peaks with intensities below the thresholds are excluded. Either:
#'   - A numeric vector of length the number of samples. Each number is a threshold for that sample
#'   - A single number. All samples use this number as baseline threshold.
#'   - A function, that takes the `nmr_dataset` and `range_without_peaks` and returns a numeric vector with the thresholds.
#'   - `NULL`. If that's the case, a default function is used ([nmr_baseline_threshold()])
#' @inheritParams speaq::detectSpecPeaks
#' @param range_without_peaks A numeric vector of length two with a region without peaks, only used when `baselineThresh = NULL`
#' @param verbose Logical (`TRUE` or `FALSE`). Show informational messages, such as the estimated baseline
#' @return A data frame with the NMRExperiment, the sample index, the position
#'     in ppm and index and the peak intensity
#' @seealso Peak_detection
#' @export
nmr_detect_peaks <- function(nmr_dataset,
                             nDivRange_ppm = 0.1,
                             scales = seq(1, 16, 2),
                             baselineThresh = NULL,
                             SNR.Th = 3,
                             range_without_peaks = c(9.5, 10),
                             verbose = FALSE) {
    nmr_dataset <- validate_nmr_dataset_1D(nmr_dataset)
    
    # Convert ppm to number of data points
    ppm_resolution <- stats::median(diff(nmr_dataset$axis))
    nDivRange <- round(nDivRange_ppm / ppm_resolution)
    
    # Computes the Baseline Threshold
    baselineThresh_fun <- NULL
    if (is.null(baselineThresh)) {
        baselineThresh_fun <- nmr_baseline_threshold
    } else if (rlang::is_function(baselineThresh)) {
        baselineThresh_fun <- baselineThresh
    }
    if (!is.null(baselineThresh_fun)) {
        baselineThresh <- rlang::exec(baselineThresh_fun, nmr_dataset, range_without_peaks = range_without_peaks)
        if (isTRUE(verbose)) {
            rlang::inform(
                message = c(
                    "i" = glue::glue("Using baselineThresh={baselineThresh}", baselineThresh=baselineThresh)
                )
            )
            rlang::inform(
                message = c(
                    "i" = glue::glue(
                        "You may plot(<your-dataset>, chemshift_range=c({rmin}, {rmax})) ",
                        "to assess that the baselineThresh was estimated correctly",
                        rmin = range_without_peaks[1], rmax = range_without_peaks[2]
                    )
                ),
                .frequency = "regularly",
                .frequency_id = "suggest-plot-for-baseline-threshold-estimation"
            )
        }
    }
    
    if (!is.numeric(baselineThresh)) {
        rlang::abort(glue::glue("The baseline threshold is not numeric: {baselineThresh}"))
    }
    if (length(baselineThresh) == 1) {
        baselineThresh <- rep(baselineThresh, times = nmr_dataset$num_samples)
    }

    data_matrix_to_list <-
        lapply(seq_len(nrow(nmr_dataset$data_1r)),
               function(i)
                   matrix(nmr_dataset$data_1r[i, ], nrow = 1))

    warn_future_to_biocparallel()
    peakList <- BiocParallel::bpmapply(
        FUN = function(spec, thresh, ...) {
            speaq::detectSpecPeaks(
                spec,
                baselineThresh = thresh,
                ...
            )[[1]]
        },
        data_matrix_to_list,
        baselineThresh,
        MoreArgs = list(
            nDivRange = nDivRange,
            scales = scales,
            SNR.Th = SNR.Th,
            verbose = FALSE
        ),
        SIMPLIFY = FALSE
    )
    peak_data <- peakList_to_dataframe(nmr_dataset, peakList)
    peaklist_fit_lorentzians(peak_data, nmr_dataset)
}

#' Overview of the peak detection results
#' 
#' This plot allows to explore the performance of the peak detection across all the samples, by summarizing
#' how many peaks are detected on each sample at each chemical shift range.
#' 
#' You can use this plot to find differences in the number of detected peaks across your dataset, and then use
#' [nmr_detect_peaks_plot()] to have a finer look at specific samples and chemical shifts, and assess graphically that the
#' peak detection results that you have are correct.
#' 
#'
#' @param peak_data The output of [nmr_detect_peaks()]
#' @param ppm_breaks A numeric vector with the breaks that will be used to count the number of the detected peaks.
#' @param accepted_only If `peak_data` contains a logical column named `accepted`, only those with `accepted=TRUE` will be counted.
#'
#' @return A scatter plot, with samples on one axis and chemical shift bins in the other axis. The size of each dot
#'   represents the number of peaks found on a sample within a chemical shift range.
#' @export
#' @seealso Peak_detection
#' @family peak detection functions
nmr_detect_peaks_plot_overview <- function(peak_data, ppm_breaks = pretty(range(peak_data$ppm), n = 20), accepted_only = TRUE) {
    to_plot <- peak_data
    if (accepted_only && "accepted" %in% colnames(to_plot)) {
        to_plot <- to_plot[to_plot$accepted, , drop = FALSE]
    }
    to_plot <- dplyr::mutate(to_plot, ppm_grp = cut(.data$ppm, breaks = !!ppm_breaks))
    to_plot <- to_plot[!is.na(to_plot$ppm_grp),, drop = FALSE]
    to_plot <- dplyr::group_by(to_plot, .data$NMRExperiment, .data$ppm_grp)
    to_plot <- dplyr::summarize(to_plot, num_peaks = dplyr::n(), .groups = "drop")
    to_plot$ppm_grp <- factor(to_plot$ppm_grp, levels = rev(levels(to_plot$ppm_grp)))
    to_plot$NMRExperiment <- factor(
        to_plot$NMRExperiment,
        levels = stringr::str_sort(unique(to_plot$NMRExperiment), numeric = TRUE)
    )
    
    gplt <- ggplot2::ggplot(to_plot) +
        ggplot2::geom_point(
            ggplot2::aes(
                x = .data$ppm_grp,
                y = .data$NMRExperiment,
                size = .data$num_peaks
            )
        ) +
        ggplot2::labs(x = "chemical shift bin", y = "NMRExperiment", size = "Number of peaks detected") +
        ggplot2::theme(legend.position = "bottom", axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
        ggplot2::coord_flip()
    gplt
}


#' Plot peak detection results
#'
#' @family peak detection functions
#' @inheritParams nmr_detect_peaks
#' @param peak_data The peak table returned by [nmr_detect_peaks]
#' @param NMRExperiment a single NMR experiment to plot
#' @param peak_id A character vector. If given, plot only that peak id.
#' @param ... Arguments passed to [plot.nmr_dataset_1D] (`chemshift_range`, `...`)
#' @param accepted_only If `peak_data` contains a logical column named `accepted`, only those with `accepted=TRUE` will be counted. By default, `accepted_only = TRUE`, unless a `peak_id` is given
#' @export
#' @return Plot peak detection results
#' 
#' @seealso Peak_detection nmr_detect_peaks
#' @family peak detection functions
#' @family nmr_dataset_1D functions
nmr_detect_peaks_plot <- function(nmr_dataset,
                                  peak_data,
                                  NMRExperiment = NULL,
                                  peak_id = NULL,
                                  accepted_only = NULL,
                                  ...) {
    
    if (!rlang::is_scalar_character(NMRExperiment) && is.null(peak_id)) {
        stop("NMRExperiment should be a string or peak_id should be a peak from one experiment")
    }
    peak_data_to_show <- peak_data
    if (is.null(NMRExperiment)) {
        peak_data_to_show <- peak_data_to_show[peak_data_to_show$peak_id %in% peak_id,,drop=FALSE]
        NMRExperiment <- unique(peak_data_to_show$NMRExperiment)
        if (!rlang::is_scalar_character(NMRExperiment)) {
            stop(glue::glue(
                "Peak ids `{peakids}` should belong to only one NMRExperiment ({exper})",
                peakids = paste0(peak_id, collapse = ", "),
                exper = paste0(NMRExperiment, collapse = ", ")
                ))
        }
    }
    if (is.null(accepted_only)) {
        if (!is.null(peak_id)) {
            accepted_only <- FALSE
        } else {
            accepted_only <- TRUE
        }
    }
    # If we plot only a subset of the data, we also plot only the required
    # vertical lines:
    dots <- list(...)
    if ("chemshift_range" %in% names(dots)) {
        chemshift_range <- dots[["chemshift_range"]]
        chemshift_range[seq_len(2)] <- range(chemshift_range[seq_len(2)])
    } else {
        chemshift_range <- range(c(peak_data_to_show$ppm, peak_data_to_show$ppm_infl_min-0.01, peak_data_to_show$ppm_infl_max+0.01))
    }
    dots[["chemshift_range"]] <- chemshift_range
    peak_data_to_show <- dplyr::filter(
        peak_data_to_show,
        .data$NMRExperiment == !!NMRExperiment,
        .data$ppm > chemshift_range[1] &
            .data$ppm < chemshift_range[2]
    )
    if (accepted_only && "accepted" %in% colnames(peak_data_to_show)) {
        peak_data_to_show <- peak_data_to_show[peak_data_to_show$accepted,,drop=FALSE]
    }
    if (!is.null(peak_id)) {
        peak_data_to_show <- peak_data_to_show[peak_data_to_show$peak_id %in% peak_id,,drop=FALSE]
    }
    # We can't make it interactive here
    interactive <- "interactive" %in% names(dots) && dots[["interactive"]]
    dots[["interactive"]] <- FALSE
    plt <- rlang::exec(
        plot,
        nmr_dataset,
        NMRExperiment = NMRExperiment,
        !!!dots
    )
    plt <- plt +
        ggplot2::geom_vline(
            data = peak_data_to_show,
            ggplot2::aes_string(xintercept = "ppm"),
            color = "black",
            linetype = "dashed"
        )

    if (interactive) {
        plotly::ggplotly(plt)
    } else {
        plt
    }
}

signif_transformer <- function(digits = 3) {
    force(digits)
    function(text, envir) {
        x <- glue::identity_transformer(text, envir)
        if (is.numeric(x)) {
            signif(x, digits = digits)
        } else {
            x
        }
    }
}


#' Plot multiple peaks from a peak list
#'
#' @param nmr_dataset The `nmr_dataset_1D` object with the spectra
#' @param peak_data A data frame, the peak list
#' @param peak_ids The peak ids to plot
#' @param caption The caption for each subplot
#'
#' @return A plot object
#' @export
#'
nmr_detect_peaks_plot_peaks <- function(
        nmr_dataset, 
        peak_data,
        peak_ids,
        caption = paste("{peak_id}", "(NMRExp. {NMRExperiment},", "\u03B3 = {gamma},",
                        "\narea = {area},", "nrmse = {norm_rmse})")) {
    
    require_pkgs(pkg = c("cowplot", "gridExtra"))
    force(nmr_dataset)
    force(peak_data)
    # Workaround https://github.com/r-lib/roxygen2/issues/1342
    caption <- gsub("NMRExp. ", "NMRExp.\u00A0", caption)
    caption <- gsub(" = ", "\u00a0=\u00a0", caption)
    plots <- purrr::map(peak_ids, function(peak_id) {
        peak_metadata <- peak_data[peak_data$peak_id == peak_id, , drop = FALSE]
        nmr_detect_peaks_plot(nmr_dataset, peak_data, peak_id = peak_id) +
            ggplot2::labs(caption = glue::glue_data(.x = peak_metadata, caption, .transformer = signif_transformer(3))) +
            ggplot2::theme(legend.position = "none", axis.title = ggplot2::element_blank())
    })
    all_plots <- cowplot::plot_grid(plotlist = plots)
    
    axis_title_x <- grid::textGrob("Chemical shift (ppm)", gp=grid::gpar(fontsize=14))
    axis_title_y <- grid::textGrob("Intensity (a.u.)", gp=grid::gpar(fontsize=14), rot=90)

    gridExtra::grid.arrange(gridExtra::arrangeGrob(all_plots, left = axis_title_y, bottom = axis_title_x))
}

#' Convert a speaq::detectSpecPeaks peak list to an interpretable data frame
#' @noRd
#' @param nmr_dataset The [nmr_dataset_1D] used to detect the peaks
#' @param peakList The peakList as returned by speaq::detectSpecPeaks
#' @return A data frame with NMRExperiment, ppm and intensity, among other columns
#' @keywords internal
peakList_to_dataframe <- function(nmr_dataset, peakList) {
    NMRExperiment <- nmr_meta_get_column(nmr_dataset, "NMRExperiment")
    out <- purrr::imap_dfr(peakList, function(peak_idx,
                                       sample_idx,
                                       nmr_dataset,
                                       NMRExperiment) {
        num_of_peaks_in_sample <- length(peak_idx)
        spec <- as.numeric(nmr_dataset$data_1r[sample_idx, ])
        data.frame(
            peak_id = "",
            NMRExperiment = rep(NMRExperiment[sample_idx], num_of_peaks_in_sample),
            sample_idx = rep(sample_idx, num_of_peaks_in_sample),
            ppm = nmr_dataset$axis[peak_idx],
            pos = peak_idx,
            intensity = spec[peak_idx],
            stringsAsFactors = FALSE
        )
    }, nmr_dataset = nmr_dataset, NMRExperiment = NMRExperiment)
    num_digits <- nchar(format(nrow(out)))
    fmt_str <- sprintf("Peak%%0%dd", num_digits)
    out$peak_id <- sprintf(fmt_str, seq_len(nrow(out)))
    out
}

#' Convert the data frame created by [peakList_to_dataframe] back to a peakList
#' This is required because speaq::dohCluster in nmr_align needs the peakList
#' @noRd
#' @param peak_data The peak list returned by [peakList_to_dataframe]
#' @return a peakList
#' @importFrom rlang .data
#' @keywords internal
peak_data_to_peakList <- function(nmr_dataset, peak_data) {
    peakList <- rep(list(numeric(0)), nmr_dataset$num_samples)
    sample_idx_peaks <- peak_data %>%
        dplyr::arrange(.data$sample_idx, .data$pos) %>%
        dplyr::group_by(.data$sample_idx) %>%
        dplyr::summarise(peak_pos = list(.data$pos)) %>%
        dplyr::ungroup()
    peakList[sample_idx_peaks$sample_idx] <-
        sample_idx_peaks$peak_pos
    peakList
}

#' Diagnose SNR threshold in peak detection
#'
#' @param ds An [nmr_dataset_1D] dataset
#' @param NMRExperiment A string with the single NMRExperiment used explore the SNR thresholds. If not given, use the first one.
#' @param SNR_thresholds A numeric vector with the SNR thresholds to explore
#' @inheritDotParams nmr_detect_peaks
#'
#' @return A list with the following elements:
#'    - `peaks_detected`: A data frame with the columns from the [nmr_detect_peaks] output and an additional column
#'                                        `SNR_threshold` with the threshold used on each row.
#'    - `num_peaks_per_region`: A summary of the `peaks_detected` table, with the number of peaks detected on
#'                                                        each chemical shift region
#'    - `plot_num_peaks_per_region`: A visual representation of `num_peaks_per_region`
#'    - `plot_spectrum_and_detections`: A visual representation of the spectrum and the peaks detected with each
#'         SNR threshold. Use [plotly::ggplotly] or [plot_interactive] on this to zoom and explore the results.
#' @family peak detection functions
#' @family nmr_dataset_1D functions
#' @export
#' @seealso nmr_detect_peaks
nmr_detect_peaks_tune_snr <- function(
    ds,
    NMRExperiment = NULL,
    SNR_thresholds = seq(from = 2, to = 6, by = 0.1),
    ...
) {
        if (is.null(NMRExperiment)) {
            NMRExperiment <-
                utils::head(nmr_meta_get_column(ds, column = "NMRExperiment"), n = 1)
        }
        ds1 <- filter(ds, NMRExperiment == !!NMRExperiment)
        names(SNR_thresholds) <- SNR_thresholds

        warn_future_to_biocparallel()
        peaks_detected_list <- BiocParallel::bplapply(
            X = SNR_thresholds,
            FUN = function(SNR.Th, nmr_dataset, ...) {
                nmr_detect_peaks(
                    nmr_dataset = ds1,
                    SNR.Th = SNR.Th,
                    ...
                )
            },
            ...
        )
        peaks_detected <- dplyr::bind_rows(peaks_detected_list, .id = "SNR_threshold")
        
        peaks_detected$SNR_threshold <-
            as.numeric(peaks_detected$SNR_threshold)
        peaks_per_region <- peaks_detected %>%
            dplyr::mutate(ppm_region = round(.data$ppm/0.5)*0.5) %>%
            dplyr::group_by(.data$SNR_threshold, .data$ppm_region) %>%
            dplyr::summarize(num_peaks = dplyr::n()) %>%
            dplyr::ungroup()
        
        
        gplt <- ggplot2::ggplot(peaks_per_region) +
            ggplot2::geom_tile(
                ggplot2::aes(
                    x = .data$ppm_region,
                    y = .data$SNR_threshold,
                    fill = .data$num_peaks
                )
            ) +
            ggplot2::scale_x_reverse(name = "ppm region", limits = rev(range(ds1$axis))) +
            ggplot2::scale_y_continuous(name = "Signal to Noise Ratio Threshold") +
            ggplot2::scale_fill_viridis_c(name = "Number of peaks detected") +
            ggplot2::ggtitle(
                "Number of peaks detected on each ppm region for several SNR thresholds",
                "Too low SNR thresholds may have false peaks, too high may miss actual peaks"
            )
        
        df1 <- data.frame(ppm = ds1$axis,
                          intensity = as.numeric(ds1$data_1r[1, ]))
        
        ord_thresh <- sort(unique(peaks_detected$SNR_threshold))
        num_thresholds <- length(ord_thresh)
        
        peak_tuning_plt <- peaks_detected %>%
            dplyr::mutate(intensity_scaled = .data$intensity * (match(.data$SNR_threshold, ord_thresh) - 1) /
                              num_thresholds)
        
        gplt2 <- ggplot2::ggplot() +
            ggplot2::geom_line(
                data = df1,
                mapping = ggplot2::aes(x = .data$ppm, y = .data$intensity),
                color = "red"
            ) +
            ggplot2::geom_segment(
                data = peak_tuning_plt,
                mapping = ggplot2::aes(
                    x = .data$ppm,
                    xend = .data$ppm,
                    y = .data$intensity_scaled,
                    yend = .data$intensity_scaled + .data$intensity / num_thresholds,
                    color = .data$SNR_threshold
                )
            ) +
            ggplot2::scale_x_reverse(name = "Chemical shift (ppm)") +
            ggplot2::scale_y_continuous("Intensity (a.u.)") +
            ggplot2::scale_color_continuous(name = "SNR thresholds") +
            ggplot2::ggtitle(
                "Peak detection for several SNR thresholds",
                "Zoom onto a peak to see the SNR thresholds were it is detected"
            )
        
        list(
            peaks_detected = peaks_detected,
            num_peaks_per_region = peaks_per_region,
            plot_num_peaks_per_region = gplt,
            plot_spectrum_and_detections = gplt2
        )
    }

#' Align NMR spectra
#'
#' This function is based on [speaq::dohCluster].
#'
#' @family alignment functions
#' @param nmr_dataset An [nmr_dataset_1D]
#' @param peak_data The detected peak data given by [nmr_detect_peaks].
#' @param maxShift_ppm The maximum shift allowed, in ppm
#' @param NMRExp_ref NMRExperiment of the reference to use for alignment
#' @inheritParams speaq::dohCluster
#' 
#' @return An [nmr_dataset_1D], with the spectra aligned
#' @export
#' @family peak alignment functions
#' @family nmr_dataset_1D functions
nmr_align <- function(nmr_dataset,
                      peak_data,
                      NMRExp_ref = NULL,
                      maxShift_ppm = 0.0015,
                      acceptLostPeak = FALSE) {
    nmr_dataset <- validate_nmr_dataset_1D(nmr_dataset)
    maxShift <-
        round(maxShift_ppm / nmr_ppm_resolution(nmr_dataset))
    if (is.null(NMRExp_ref)) {
        NMRExp_ref <- nmr_align_find_ref(nmr_dataset, peak_data)
    }
    NMRExp <- nmr_meta_get_column(nmr_dataset, "NMRExperiment")
    refInd <- which(NMRExp == NMRExp_ref)
    if (length(refInd) != 1) {
        stop("Wrong NMRExperiment as align_ref? Please check.")
    }
    peakList <- peak_data_to_peakList(nmr_dataset, peak_data)
    nmr_dataset$data_1r <- speaq::dohCluster(
        nmr_dataset$data_1r,
        peakList = peakList,
        refInd = refInd,
        maxShift    = maxShift,
        acceptLostPeak = acceptLostPeak,
        verbose = FALSE
    )
    nmr_dataset
}


#' Find alignment reference
#'
#' @inheritParams nmr_align
#' @family alignment functions
#'
#' @return The NMRExperiment of the reference sample
#' @export
#'
#' @family peak alignment functions
#' @family nmr_dataset_1D functions
nmr_align_find_ref <- function(nmr_dataset, peak_data) {
    peakList <- peak_data_to_peakList(nmr_dataset, peak_data)
    resFindRef <- speaq::findRef(peakList)
    NMRExperiment <-
        nmr_meta_get_column(nmr_dataset, "NMRExperiment")
    c(NMRExperiment = NMRExperiment[resFindRef$refInd])
}

#' Build list of regions for peak integration
#'
#' @param peak_pos_ppm The peak positions, in ppm
#' @param peak_width_ppm The peak widths (or a single peak width for all peaks)
#'
#' @return A list of regions suitable for [nmr_integrate_regions]
#'
#' @family peak detection functions
regions_from_peak_table <- function(peak_pos_ppm, peak_width_ppm) {
    if (length(peak_width_ppm) == 1) {
        peak_width_ppm <- rep(peak_width_ppm, length(peak_pos_ppm))
    }
    purrr::map2(peak_pos_ppm, peak_width_ppm,
                function(ppm, peak_width) {
                    c(ppm - peak_width / 2, ppm + peak_width / 2)
                })
}


#' PPM resolution of the spectra
#'
#' The function gets the ppm resolution of the dataset using the median of the
#' difference of data points.
#'
#' @param nmr_dataset An object containing NMR samples
#'
#' @return Numeric (the ppm resolution, measured in ppms)
#' @examples
#' nmr_dataset <- nmr_dataset_load(system.file("extdata", "nmr_dataset.rds", package = "AlpsNMR"))
#' nmr_ppm_resolution(nmr_dataset)
#' message("the ppm resolution of this dataset is ", nmr_ppm_resolution(nmr_dataset), " ppm")
#'
#' @export
nmr_ppm_resolution <- function(nmr_dataset) {
    UseMethod("nmr_ppm_resolution")
}

#' @rdname nmr_ppm_resolution
#' @family nmr_dataset functions
#' @export 
#' @examples
#' nmr_dataset <- nmr_dataset_load(system.file("extdata", "nmr_dataset.rds", package = "AlpsNMR"))
#' nmr_ppm_resolution(nmr_dataset)
#' message("the ppm resolution of this dataset is ", nmr_ppm_resolution(nmr_dataset), " ppm")
#' 
nmr_ppm_resolution.nmr_dataset <- function(nmr_dataset) {
    # For each sample:
    purrr::map(nmr_dataset$axis, function(axis_sample) {
        # For each dimension of each sample:
        purrr::map_dbl(axis_sample, function(axis_dimension) {
            stats::median(abs(diff(axis_dimension)))
        })
    })
}

#' @rdname nmr_ppm_resolution
#' @family nmr_dataset_1D functions
#' @export 
#' @examples 
#' nmr_dataset <- nmr_dataset_load(system.file("extdata", "nmr_dataset.rds", package = "AlpsNMR"))
#' nmr_ppm_resolution(nmr_dataset)
#' message("the ppm resolution of this dataset is ", nmr_ppm_resolution(nmr_dataset), " ppm")
nmr_ppm_resolution.nmr_dataset_1D <- function(nmr_dataset) {
    stats::median(abs(diff(nmr_dataset$axis)))
}

#' Unlisted PPM resolution
#'
#' A wrapper to unlist the output from the function
#' `nmr_ppm_resolution(nmr_dataset)` when no interpolation has been applied.
#'
#' @param nmr_dataset An object containing NMR samples
#'
#' @return A number (the ppm resolution, measured in ppms)
#' @examples
#' nmr_dataset <- nmr_dataset_load(system.file("extdata", "nmr_dataset.rds", package = "AlpsNMR"))
#' nmr_ppm_resolution(nmr_dataset)
#' @return Numeric (the ppm resolution, measured in ppms)
#' @export 
#'
ppm_resolution <- function (nmr_dataset) {
    unlist(nmr_ppm_resolution(nmr_dataset[1]))
}
