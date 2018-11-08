#' Peak detetion for NMR
#'
#' @param nmr_dataset An [nmr_dataset_1D].
#' @inheritParams speaq::detectSpecPeaks
#' @return A data frame with the NMRExperiment, the sample index, the position in ppm and index and the peak intensity
#' @export
#'
nmr_detect_peaks <- function(nmr_dataset, nDivRange = 128, scales = seq(1, 16, 2),
                             baselineThresh = 0.01, SNR.Th = -1) {
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
  NMRExperiment <- nmr_get_metadata(nmr_dataset, columns = "NMRExperiment")$NMRExperiment
  purrr::imap_dfr(peakList, function(peak_idx, sample_idx) {
    data.frame(NMRExperiment = NMRExperiment[sample_idx],
               sample_idx = sample_idx,
               ppm = nmr_dataset$axis[peak_idx],
               pos = peak_idx,
               intensity = as.numeric(nmr_dataset$data_1r[sample_idx, peak_idx]),
               stringsAsFactors = FALSE)
  })
}

#' Convert the data frame created by [peakList_to_dataframe] back to a peakList
#' This is required because speaq::dohCluster in nmr_align needs the peakList
#' @noRd 
#' @param peak_data The peak list returned by [peakList_to_dataframe]
#' @return a peakList
#' @keywords internal
peak_data_to_peakList <- function(peak_data) {
  peak_data %>%
    dplyr::arrange(sample_idx, pos) %>%
    dplyr::group_by(sample_idx) %>%
    dplyr::summarize(peak_pos = list(pos)) %>%
    dplyr::ungroup() %>%
    dplyr::pull(peak_pos)
}

#' Align NMR spectra
#'
#' @param nmr_dataset An [nmr_dataset_1D]
#' @param peak_data The detected peak data given by [nmr_detect_peaks].
#' @inheritParams speaq::dohCluster
#'
#' @return An [nmr_dataset_1D], with the spectra aligned
#' @export
#'
nmr_align <- function(nmr_dataset, peak_data, maxShift = 3, acceptLostPeak = FALSE) {
  peakList <- peak_data_to_peakList(peak_data)
  nmr_dataset_align <- nmr_dataset
  resFindRef <- speaq::findRef(peakList)
  nmr_dataset_align$data_1r <- speaq::dohCluster(
    nmr_dataset$data_1r,
    peakList = peakList,
    refInd = resFindRef$refInd,
    maxShift  = 3,
    acceptLostPeak = FALSE, 
    verbose = FALSE
  )
  nmr_dataset_align
}
