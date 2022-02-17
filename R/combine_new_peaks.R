#' Combining peaks across the spectra
#' 
#' The function combines a peak list from the reference spectrum with detected
#' peaks in spectra other that the reference spectrum with a given window.This
#' ppm threshold should be the same value set in the [nmr_detect_peaks]
#' function.
#'
#' @param nmr_dataset An object containing NMR samples
#' @param peak_data The peak table returned by [nmr_detect_peaks]
#' @param NMRExp_ref NMRExperiment of the reference to use for alignment,
#'   returned by [nmr_align_find_ref]
#' @param maxShift_ppm The maximum shift allowed, in ppm, to consider a
#'   different peak
#'
#' @return A peak list including peaks others than the ones found in the ref
#'   spectrum.
#' @export
#'
#' @examples
#' \dontrun{
#' plan(multiprocess, workers = 12)
#' 
#' # 1.Peak detection in the dataset.
#' peak_data <- nmr_detect_peaks(nmr_dataset,
#'                               nDivRange_ppm = 0.1,
#'                               scales = seq(1, 16, 2),
#'                               baselineThresh = 0,
#'                               SNR.Th = 4) 
#'                               
#' # 2.Find the reference spectrum to align with.
#' NMRExp_ref <- nmr_align_find_ref(nmr_dataset, peak_data)
#' 
#' # 3.Spectra alignment using the ref spectrum and a maximum alignment shift
#' nmr_dataset <- nmr_align(nmr_dataset, 
#'                          peak_data, 
#'                          NMRExp_ref = NMRExp_ref, 
#'                          maxShift_ppm = 0.01, # max alignment shift
#'                          acceptLostPeak = FALSE) # lost peaks
#' 
#' # 4.Adding new peaks others than ref spectrum
#' new_peak_data <- combine_new_peaks(nmr_dataset,
#'                                    peak_data, 
#'                                    NMRExp_ref,
#'                                    maxShift_ppm = 0.01)
#' }
combine_new_peaks <- function(nmr_dataset, peak_data, NMRExp_ref, maxShift_ppm = 0.01) {
    refvec <- dplyr::filter(peak_data, NMRExperiment == !!NMRExp_ref)$ppm
    message("Number of peaks before adding new peaks: ", length(refvec))
    
    peak_data_to_peakList <- function(nmr_dataset, peak_data) {
        peakList <- rep(list(numeric(0)), nmr_dataset$num_samples)
        sample_idx_peaks <- peak_data %>%
            dplyr::arrange(.data$sample_idx, .data$ppm) %>%
            dplyr::group_by(.data$sample_idx) %>%
            dplyr::summarise(peak_ppm = list(.data$ppm)) %>%
            dplyr::ungroup()
        peakList[sample_idx_peaks$sample_idx] <- sample_idx_peaks$peak_ppm
        peakList
    }
    
    peakListX2 <- peak_data_to_peakList(nmr_dataset, peak_data)
    
    for (i in 1:nrow(nmr_data(nmr_dataset))){
        if (i!=NMRExp_ref) {
            peakvec = peakListX2[[i]]
            for (j in 1:length(peakvec)) {
                if (min(abs(peakvec[j]-refvec)) > (maxShift_ppm + 0.001)) refvec = c(refvec, peakvec[j])
            }
        }
    }
    refvec=sort(refvec) 
    
    refvec_df <- as.data.frame(refvec)
    colnames(refvec_df) <- "ppm"
    
    # And the original peak data 
    peak_data_original <- dplyr::filter(peak_data, 
                                        NMRExperiment == !!NMRExp_ref)
    
    message("Final number of peaks: ", dim(refvec_df)[1])
    
    return(refvec_df)
}