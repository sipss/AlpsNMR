#' Peak detetion for NMR
#'
#' @param nmr_dataset An [nmr_dataset_1D].
#' @inheritParams speaq::detectSpecPeaks
#' @return A data frame with the NMRExperiment, the sample index, the position in ppm and index and the peak intensity
#' @export
#'
nmr_detect_peaks <- function(nmr_dataset, nDivRange = 128, scales = seq(1, 16, 2),
                             baselineThresh = 0.01, SNR.Th = -1) {
  validate_nmr_dataset_1D(nmr_dataset)
  # The speaq package uses partial matching of arguments.
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
  # calling speaq::detectSpecPeaks because speaq v.2.3.3 uses partial argument
  # matching. :-(
  #
  # As we don't control the speaq code, we disable warnPartialMatchArgs
  # in this function if you had it enabled and we restore it when this function finishes
  # 
  warnPartialMatchArgs <- getOption("warnPartialMatchArgs")
  if (isTRUE(warnPartialMatchArgs)) {
    options(warnPartialMatchArgs = FALSE)
    on.exit({options(warnPartialMatchArgs = warnPartialMatchArgs)})
  }
  # Also speaq::detectSpecPeaks can print a lot of messages to stdout, even
  # when verbose = FALSE. This can lead to problems, especially with large
  # datasets in Rmd documents, so we wrap the detectSpecPeaks with capture.output
  # to prevent the output of the function to leak.
  utils::capture.output({
    peakList <- speaq::detectSpecPeaks(
      nmr_dataset$data_1r,
      nDivRange = nDivRange,
      scales = scales,
      baselineThresh = baselineThresh,
      SNR.Th = SNR.Th,
      verbose = FALSE
    )})
  
  peakList_to_dataframe(nmr_dataset, peakList)
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

#' Align NMR spectra
#'
#' @param nmr_dataset An [nmr_dataset_1D]
#' @param peak_data The detected peak data given by [nmr_detect_peaks].
#' @param NMRExp_ref NMRExperiment of the reference to use for alignment
#' @inheritParams speaq::dohCluster
#'
#' @return An [nmr_dataset_1D], with the spectra aligned
#' @export
#'
nmr_align <- function(nmr_dataset, peak_data, NMRExp_ref = NULL, maxShift = 3, acceptLostPeak = FALSE) {
  validate_nmr_dataset_1D(nmr_dataset)
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
#'
#' @return The NMRExperiment of the reference sample
#' @export
#'
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
#' @return A list of regions suitable for nmr_integrate_regions
#'
regions_from_peak_table <- function(peak_pos_ppm, peak_width_ppm) {
  if (length(peak_width_ppm) == 1) {
    peak_width_ppm <- rep(peak_width_ppm, length(peak_pos_ppm))
  }
  purrr::map2(peak_pos_ppm, peak_width_ppm,
             function(ppm, peak_width) {
               c(ppm - peak_width/2, ppm + peak_width/2)
             })
}
