#' Peak width estimation for integration
#'
#' Estimates the peak width (ppm width) to perform peak integration using
#' `nmr_integrate_peak_positions`. for this purpose, the full width at half maximum
#' of a peak from alanine doublet is considered.
#'
#' @family peak integration functions
#' @family nmr_dataset_1D functions
#' @param nmr_dataset An [nmr_dataset_1D].
#' @return Numerical. A peak width (ppm) that may be set in `nmr_integrate_peak_positions`

library(AlpsNMR)
library(gridExtra)
library(ggplot2)
library(plotly)
library(GGally)


get_peak_bounds <- function(up, down, length_sgf, pos){
  up <- up[up - pos > 0]
  down <- down[pos - down > 0]
  if(length(up) == 0 || length(down) == 0){
    return(c("left" = 0, "right" = 0))
  }
  right <- up[1]
  left <- down[length(down)]
  numr <- right - pos
  numl <- pos - left 
  return(c("left" = left, "right" = right))
}

lorentzian <- function(x, x0, gamma, A) {
  A * (1/(pi*gamma)) * ((gamma^2) / ((x-x0)^2 + (gamma)^2))
}


fitting <- function(estimated_peak, d, l, r, pos, ymaxi){
  ms <- sum(
    (d$y[l:r] - (estimated_peak[l:r] + ymaxi - estimated_peak[pos]))**2
  )/(r - l + 1)
  rms <- sqrt(ms)
  norm_rmse <- rms/estimated_peak[pos]
  return(norm_rmse)
}

if (file.exists("/home/laura/IBEC/AlpsNMR-Tutorial/results/07-alignment/peak_data_bo.rds")) {
  nmr_dataset <- readRDS("/home/laura/IBEC/AlpsNMR-Tutorial/results/07-alignment/nmr_dataset.rds")
  peak_data_bo <- readRDS("/home/laura/IBEC/AlpsNMR-Tutorial/results/07-alignment/peak_data_bo.rds")
} else {
  nmr_dataset <- readRDS("/home/laura/IBEC/AlpsNMR-Tutorial/results/06-filter-samples/nmr_dataset.rds")
  peak_data <- nmr_detect_peaks(nmr_dataset,
                                nDivRange_ppm = 0.1,
                                scales = seq(1, 16, 2),
                                baselineThresh = NULL,
                                SNR.Th = 3)
  message("Choosing alignment reference...")
  NMRExp_ref <- nmr_align_find_ref(nmr_dataset, peak_data)
  message("Starting alignment...")
  nmr_dataset <- nmr_align(nmr_dataset, peak_data,
                           NMRExp_ref = NMRExp_ref,
                           maxShift_ppm = 0.0015,
                           acceptLostPeak = FALSE)
  peak_data <- nmr_detect_peaks(nmr_dataset,
                                nDivRange_ppm = 0.1,
                                scales = seq(1, 16, 2),
                                baselineThresh = NULL,
                                SNR.Th = 3)
  peak_data_bo <- peak_data
  saveRDS(nmr_dataset, "/home/laura/IBEC/AlpsNMR-Tutorial/results/07-alignment/nmr_dataset.rds")
  saveRDS(peak_data_bo, "/home/laura/IBEC/AlpsNMR-Tutorial/results/07-alignment/peak_data_bo.rds")
  rm(peak_data, NMRExp_ref)
}


peak_data <- peak_data_bo
area_estimation <- function(peak_data, error){
  peak_data$gamma <- NA_real_
  peak_data$area <- NA_real_
  peak_data$norm_rmse <- NA_real_
  peak_data$accepted <- NA
  
  sindex_anterior <- -1
  for (i in 1:nrow(peak_data)){
    xmaxi <- peak_data$ppm[i]
    ymaxi <- peak_data$intensity[i]
    posi <- peak_data$pos[i]
    sindex <- peak_data$sample_idx[i]
    if (sindex != sindex_anterior) {
      d2 <- data.frame(x = nmr_dataset$axis, y = nmr_dataset$data_1r[sindex,])
      sgf <- signal::sgolayfilt(d2$y,p=2, n=21, m = 2, ts=d2$x[2]-d2$x[1]) #segona derivada
      up <- which((diff(sign(sgf)))== 2)
      down <- which((diff(sign(sgf)))== -2)
    }
    out2 <- get_peak_bounds(up, down, length(sgf), posi)
    gamma <- sqrt(3)/2 * (d2$x[out2["right"]]-d2$x[out2["left"]])
    estimated_A <- - (pi * gamma^3 * sgf[posi])/2
    estimated_peak <- lorentzian(
      x = d2$x,
      x0 = xmaxi,
      gamma = gamma,
      A = estimated_A
    )
    norm_rmse <- fitting(estimated_peak, d2, out2["left"], out2["right"], posi, ymaxi)
    peak_data$gamma[i] <- gamma
    peak_data$area[i] <- estimated_A
    peak_data$norm_rmse[i] <- norm_rmse
    if (norm_rmse < error & norm_rmse > 0){ peak_data$accepted[i] <- TRUE }else{peak_data$accepted[i] <- FALSE}
    sindex_anterior <- sindex
  }
  return(peak_data)
}


peak_data <- area_estimation(peak_data, 0.5)

