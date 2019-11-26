#' Peak detection for NMR
#'
#' The function detects peaks on an [nmr_dataset_1D] object, using
#' [speaq::detectSpecPeaks]. `detectSpecPeaks` divides the whole spectra into
#' smaller segments and uses [MassSpecWavelet::peakDetectionCWT] for peak
#' detection.
#' 
#' @family peak detection functions
#' @family nmr_dataset_1D functions
#' @seealso [nmr_align] for peak alignment with the detected peak table
#' @param nmr_dataset An [nmr_dataset_1D].
#' @param nDivRange_ppm Segment size, in ppms, to divide the spectra and search
#'   for peaks.
#' @param baselineThresh  It will remove all peaks under an intensity set by
#'   baselineThresh. If you set it to 'NULL', nmr_detect_peaks will
#'   automatically compute an aproximate value considering baseline between 9.5
#'   and 10.0 ppm (automatically calculation using `baselineThresh = NULL` will
#'   not work if spectra were not interpolated up to 10.0 ppm)
#' @inheritParams speaq::detectSpecPeaks
#' @return A data frame with the NMRExperiment, the sample index, the position
#'   in ppm and index and the peak intensity
#' @examples 
#' \dontrun{
#' 
#' # 0. Multiprocess (parallelization) to set the number of cores working in your PC
#' plan(multiprocess, workers = 12)
#' 
#' # 1.Peak detection in the dataset.
#' peak_data <- nmr_detect_peaks(nmr_dataset,
#'                               nDivRange_ppm = 0.1, # Size of detection segments
#'                               scales = seq(1, 16, 2),
#'                               baselineThresh = 0, # Minimum peak intensity
#'                               SNR.Th = 4) # Signal to noise ratio
#'                               
#' # 2.Find the reference spectrum to align with.
#' NMRExp_ref <- nmr_align_find_ref(nmr_dataset, peak_data)
#' 
#' # 3.Spectra alignment using the ref spectrum and a maximum alignment shift
#' nmr_dataset <- nmr_align(nmr_dataset, # the dataset
#'                          peak_data, # detected peaks
#'                          NMRExp_ref = NMRExp_ref, # ref spectrum
#'                          maxShift_ppm = 0.0015, # max alignment shift
#'                          acceptLostPeak = FALSE) # lost peaks
#'                          
#' # 4.Set sequential working to finish parallelization
#' plan(sequential)
#' 
#' # 5.Peak integration (please, consider previous normalization step).
#' # First we take the peak table from the reference spectrum
#' peak_data_ref <- filter(peak_data, NMRExperiment == NMRExp_ref)
#' 
#' # Then we integrate spectra considering the peaks from the ref spectrum
#' nmr_peak_table <- nmr_integrate_peak_positions(
#' samples = nmr_dataset,
#' peak_pos_ppm = peak_data_ref$ppm,
#' peak_width_ppm = NULL) 
#' 
#' #If you wanted the final peak table before machine learning you can run
#' nmr_peak_table_completed <- get_integration_with_metadata(nmr_peak_table)
#' 
#' }
#' @export
#'
nmr_detect_peaks <- function(nmr_dataset, nDivRange_ppm = 0.1,
                             scales = seq(1, 16, 2),
                             baselineThresh = NULL, SNR.Th = 3) {
  validate_nmr_dataset_1D(nmr_dataset)
  
  # Convert ppm to number of data points
  ppm_resolution <- stats::median(diff(nmr_dataset$axis))
  nDivRange <- round(nDivRange_ppm/ppm_resolution)
  
  # Computes the Baseline Threshold ###NEW###
  if(is.null(baselineThresh)){
    baselineThresh <- nmr_baseline_threshold(nmr_dataset)
  }
  
  
  # A dependency of the speaq package uses partial matching of arguments.
  #
  # ## What is partial argument matching?
  # partial matching is an R feature that allows to call functions without
  # fully specifying the arguments, for instance:
  # nmr_detect_peaks(n = dataset)
  # nmr_detect_peaks(nmr_dat = dataset)
  # nmr_detect_peaks(nmr_datase = dataset)
  # These can work thanks to partial matching of arguments.
  #
  # ## Why is partial argument matching a problem?
  # Partial Arg matching follows the idea that "what you write is *probably* what you want"
  # A better approach is to use the full argument name.
  #
  # ## Can I be warned if I used partial arguments?
  # Yes, you can use the warnPartialMatchArgs option
  # 
  # ## Why am I reading this?
  #
  # If you enable warnPartialMatchArgs, you will get several warnings when
  # calling speaq::detectSpecPeaks() because MassSpecWavelet::cwt() (a dependency)
  # uses partial argument matching by mistake on a seq() function call.
  # Wrong:
  #   MassSpecWavelet/R/cwt.R:   psi_xval <- seq(-8, 8, length = 1024)
  # Right:
  #   MassSpecWavelet/R/cwt.R:   psi_xval <- seq(-8, 8, length.out = 1024)
  #
  # As we don't control the MassSpecWavelet code, we disable warnPartialMatchArgs
  # in this function if you had it enabled, and we restore it when this function finishes
  # 
  warnPartialMatchArgs <- getOption("warnPartialMatchArgs")
  if (isTRUE(warnPartialMatchArgs)) {
    options(warnPartialMatchArgs = FALSE)
    on.exit({options(warnPartialMatchArgs = warnPartialMatchArgs)})
  }
  
  prgrs <- show_progress_bar(nrow(nmr_dataset$data_1r) > 5)
  
  data_matrix_to_list <- lapply(
    seq_len(nrow(nmr_dataset$data_1r)),
    function(i) matrix(nmr_dataset$data_1r[i,], nrow = 1)
  )
  
  peakList <- furrr::future_map(
    data_matrix_to_list,
    function(spec, ...) {
      speaq::detectSpecPeaks(spec, ...)[[1]]
    },
    nDivRange = nDivRange,
    scales = scales,
    baselineThresh = baselineThresh,
    SNR.Th = SNR.Th,
    verbose = FALSE,
    .progress = prgrs,
    .options = furrr::future_options(globals = character(0L),
                                     packages = character(0L))
  )
  
  # peakList <- speaq::detectSpecPeaks(
  #   nmr_dataset$data_1r,
  #   nDivRange = nDivRange,
  #   scales = scales,
  #   baselineThresh = baselineThresh,
  #   SNR.Th = SNR.Th,
  #   verbose = FALSE
  # )
  
  peakList_to_dataframe(nmr_dataset, peakList)
}

#' Plot peak detection results
#' 
#' @family peak detection functions
#' @inheritParams nmr_detect_peaks
#' @param peak_data The peak table returned by [nmr_detect_peaks]
#' @param NMRExperiment a single NMR experiment to plot
#' @param ... Arguments passed to [plot.nmr_dataset_1D] (`chemshift_range`, `...`)
#' @export
#' @examples 
#' \dontrun{
#' 
#' # 0. Multiprocess (parallelization) to set the number of cores working in your PC
#' plan(multiprocess, workers = 12)
#' 
#' # 1.Peak detection in the dataset.
#' peak_data <- nmr_detect_peaks(nmr_dataset,
#'                               nDivRange_ppm = 0.1, # Size of detection segments
#'                               scales = seq(1, 16, 2),
#'                               baselineThresh = 0, # Minimum peak intensity
#'                               SNR.Th = 4) # Signal to noise ratio
#'                               
#' # 2.Find the reference spectrum to align with.
#' NMRExp_ref <- nmr_align_find_ref(nmr_dataset, peak_data)
#' 
#' # 3.Spectra alignment using the ref spectrum and a maximum alignment shift
#' nmr_dataset <- nmr_align(nmr_dataset, # the dataset
#'                          peak_data, # detected peaks
#'                          NMRExp_ref = NMRExp_ref, # ref spectrum
#'                          maxShift_ppm = 0.0015, # max alignment shift
#'                          acceptLostPeak = FALSE) # lost peaks
#'                          
#' # 4.Set sequential working to finish parallelization
#' plan(sequential)
#' 
#' # 5.Peak integration (please, consider previous normalization step).
#' # First we take the peak table from the reference spectrum
#' peak_data_ref <- filter(peak_data, NMRExperiment == NMRExp_ref)
#' 
#' # Then we integrate spectra considering the peaks from the ref spectrum
#' nmr_peak_table <- nmr_integrate_peak_positions(
#' samples = nmr_dataset,
#' peak_pos_ppm = peak_data_ref$ppm,
#' peak_width_ppm = NULL) 
#' 
#' #If you wanted the final peak table before machine learning you can run
#' nmr_peak_table_completed <- get_integration_with_metadata(nmr_peak_table)
#' 
#' }
#' @family peak detection functions
#' @family nmr_dataset_1D functions
nmr_detect_peaks_plot <- function(nmr_dataset, peak_data, NMRExperiment, ...) {
  if (!rlang::is_scalar_character(NMRExperiment)) {
    stop("NMRExperiment should be a string")
  }
  # If we plot only a subset of the data, we also plot only the required 
  # vertical lines:
  dots <- list(...)
  if ("chemshift_range" %in% names(dots)) {
    chemshift_range <- dots[["chemshift_range"]]
    chemshift_range[1:2] <- range(chemshift_range[1:2])
  } else {
    chemshift_range <- range(peak_data$ppm)
  }
  peak_data_to_show <- dplyr::filter(
    peak_data,
    .data$NMRExperiment == !!NMRExperiment,
    .data$ppm > chemshift_range[1] & .data$ppm < chemshift_range[2])
  # Plot:
  plot(nmr_dataset, NMRExperiment = NMRExperiment, ..., interactive = FALSE) +
    ggplot2::geom_vline(data = peak_data_to_show,
                        ggplot2::aes_string(xintercept = "ppm"), color = "black", linetype = "dashed")
}

#' Convert a speaq::detectSpecPeaks peak list to an interpretable data frame
#' @noRd
#' @param nmr_dataset The [nmr_dataset_1D] used to detect the peaks
#' @param peakList The peakList as returned by speaq::detectSpecPeaks
#' @return A data frame with NMRExperiment, ppm and intensity, among other columns
#' @keywords internal
peakList_to_dataframe <- function(nmr_dataset, peakList) {
  NMRExperiment <- nmr_meta_get_column(nmr_dataset, "NMRExperiment")
  purrr::imap_dfr(peakList, function(peak_idx, sample_idx, nmr_dataset, NMRExperiment) {
    num_of_peaks_in_sample <- length(peak_idx)
    spec <- as.numeric(nmr_dataset$data_1r[sample_idx,])
    data.frame(NMRExperiment = rep(NMRExperiment[sample_idx], num_of_peaks_in_sample),
               sample_idx = rep(sample_idx, num_of_peaks_in_sample),
               ppm = nmr_dataset$axis[peak_idx],
               pos = peak_idx,
               intensity = spec[peak_idx],
               stringsAsFactors = FALSE)
  }, nmr_dataset = nmr_dataset, NMRExperiment = NMRExperiment)
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
  peakList[sample_idx_peaks$sample_idx] <- sample_idx_peaks$peak_pos
  peakList
}

#' Diagnose SNR threshold in peak detection
#' 
#' @section Parallelization:
#'
#' This function accepts parallellization with future strategies. 
#' You can use `plan(multiprocess)` or `plan(sequential)` before calling this function to determine
#' if it should be parallellized or not.
#' 
#' @param ds An [nmr_dataset_1D] dataset
#' @param NMRExperiment A string with the single NMRExperiment used explore the SNR thresholds. If not given, use the first one.
#' @param SNR_thresholds A numeric vector with the SNR thresholds to explore
#'
#' @return A list with the following elements:
#'  - `peaks_detected`: A data frame with the columns from the [nmr_detect_peaks] output and an additional column
#'                    `SNR_threshold` with the threshold used on each row.
#'  - `num_peaks_per_region`: A summary of the `peaks_detected` table, with the number of peaks detected on
#'                            each chemical shift region
#'  - `plot_num_peaks_per_region`: A visual representation of `num_peaks_per_region`
#'  - `plot_spectrum_and_detections`: A visual representation of the spectrum and the peaks detected with each
#'     SNR threshold. Use [plotly::ggplotly] or [plot_interactive] on this to zoom and explore the results.
#' @family peak detection functions
#' @family nmr_dataset_1D functions
#' @export
nmr_detect_peaks_tune_snr <- function(ds, NMRExperiment = NULL, SNR_thresholds = seq(from = 2, to = 6, by = 0.1)) {
  
  if (is.null(NMRExperiment)) {
    NMRExperiment <- utils::head(nmr_meta_get_column(ds, column = "NMRExperiment"), n = 1)
  }
  ds1 <- filter(ds, NMRExperiment == !!NMRExperiment)
  names(SNR_thresholds) <- SNR_thresholds
  
  peaks_detected <- furrr::future_map_dfr(
    SNR_thresholds,
    ~nmr_detect_peaks(ds1,
                      nDivRange_ppm = 0.03, 
                      scales = seq(1, 16, 2),
                      baselineThresh = NULL, SNR.Th = .), #baselineThresh = 0 antes
    .id = "SNR_threshold",
    .options = furrr::future_options(globals = character(), packages = character())
  )
  
  peaks_detected$SNR_threshold <- as.numeric(peaks_detected$SNR_threshold)
  peaks_per_region <- peaks_detected %>%
    dplyr::mutate(ppm_region = plyr::round_any(.data$ppm, 0.5)) %>%
    dplyr::group_by(.data$SNR_threshold, .data$ppm_region) %>%
    dplyr::summarize(num_peaks = dplyr::n()) %>%
    dplyr::ungroup()
  
  
  gplt <- ggplot2::ggplot(peaks_per_region) + 
    ggplot2::geom_tile(ggplot2::aes(x = .data$ppm_region, y = .data$SNR_threshold, fill = .data$num_peaks)) +
    ggplot2::scale_x_reverse(name = "ppm region", limits = rev(range(ds1$axis))) +
    ggplot2::scale_y_continuous(name = "Signal to Noise Ratio Threshold") +
    ggplot2::scale_fill_viridis_c(name = "Number of peaks detected") +
    ggplot2::ggtitle("Number of peaks detected on each ppm region for several SNR thresholds",
                     "Too low SNR thresholds may have false peaks, too high may miss actual peaks")
  
  df1 <- data.frame(
    ppm = ds1$axis,
    intensity = as.numeric(ds1$data_1r[1,])
  )
  
  ord_thresh <- sort(unique(peaks_detected$SNR_threshold))
  num_thresholds <- length(ord_thresh)
  
  peak_tuning_plt <- peaks_detected %>%
    dplyr::mutate(intensity_scaled = .data$intensity*(match(.data$SNR_threshold, ord_thresh) - 1)/num_thresholds)
  
  gplt2 <- ggplot2::ggplot() + 
    ggplot2::geom_line(
      data = df1,
      mapping = ggplot2::aes(x = .data$ppm, y = .data$intensity), color = "red"
    ) +
    ggplot2::geom_segment(
      data = peak_tuning_plt, 
      mapping = ggplot2::aes(
        x = .data$ppm, xend = .data$ppm,
        y = .data$intensity_scaled,
        yend = .data$intensity_scaled + .data$intensity/num_thresholds,
        color = .data$SNR_threshold)
    ) +
    ggplot2::scale_x_reverse(name = "Chemical shift (ppm)") +
    ggplot2::scale_y_continuous("Intensity (a.u.)") +
    ggplot2::scale_color_continuous(name = "SNR thresholds") +
    ggplot2::ggtitle(
      "Peak detection for several SNR thresholds",
      "Zoom onto a peak to see the SNR thresholds were it is detected"
    )
  
  list(peaks_detected = peaks_detected,
       num_peaks_per_region = peaks_per_region,
       plot_num_peaks_per_region = gplt,
       plot_spectrum_and_detections = gplt2)
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
#' @examples 
#' \dontrun{
#' 
#' # 0. Multiprocess (parallelization) to set the number of cores working in your PC
#' plan(multiprocess, workers = 12)
#' 
#' # 1.Peak detection in the dataset.
#' peak_data <- nmr_detect_peaks(nmr_dataset,
#'                               nDivRange_ppm = 0.1, # Size of detection segments
#'                               scales = seq(1, 16, 2),
#'                               baselineThresh = 0, # Minimum peak intensity
#'                               SNR.Th = 4) # Signal to noise ratio
#'                               
#' # 2.Find the reference spectrum to align with.
#' NMRExp_ref <- nmr_align_find_ref(nmr_dataset, peak_data)
#' 
#' # 3.Spectra alignment using the ref spectrum and a maximum alignment shift
#' nmr_dataset <- nmr_align(nmr_dataset, # the dataset
#'                          peak_data, # detected peaks
#'                          NMRExp_ref = NMRExp_ref, # ref spectrum
#'                          maxShift_ppm = 0.0015, # max alignment shift
#'                          acceptLostPeak = FALSE) # lost peaks
#'                          
#' # 4.Set sequential working to finish parallelization
#' plan(sequential)
#' 
#' # 5.Peak integration (please, consider previous normalization step).
#' # First we take the peak table from the reference spectrum
#' peak_data_ref <- filter(peak_data, NMRExperiment == NMRExp_ref)
#' 
#' # Then we integrate spectra considering the peaks from the ref spectrum
#' nmr_peak_table <- nmr_integrate_peak_positions(
#' samples = nmr_dataset,
#' peak_pos_ppm = peak_data_ref$ppm,
#' peak_width_ppm = NULL) 
#' 
#' #If you wanted the final peak table before machine learning you can run
#' nmr_peak_table_completed <- get_integration_with_metadata(nmr_peak_table)
#' 
#' }
#' @return An [nmr_dataset_1D], with the spectra aligned
#' @export
#' @family peak alignment functions
#' @family nmr_dataset_1D functions
nmr_align <- function(nmr_dataset, peak_data, NMRExp_ref = NULL,
                      maxShift_ppm = 0.0015, acceptLostPeak = FALSE) {
  validate_nmr_dataset_1D(nmr_dataset)
  maxShift <- round(maxShift_ppm/nmr_ppm_resolution(nmr_dataset))
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
    maxShift  = maxShift,
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
  NMRExperiment <- nmr_meta_get_column(nmr_dataset, "NMRExperiment")
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
               c(ppm - peak_width/2, ppm + peak_width/2)
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


