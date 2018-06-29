#' Batman helpers
#' 
#' These are helper functions to make Batman tests easier
#' 
#' @param ppmRange ppm ranges for analys
#' @param specNo Ranges of spectra number to be included
#' @param bopts Batman options
#' @param batman_dir Batman input directorye
#' @param filename Filename to use, inside `batman_dir`
#' @param metabolite_names A character vector of the metabolite names to consider
#' 
#' @name nmr_batman
#' 
NULL

#' @rdname nmr_batman
#' @export
nmr_batman_options <- function(ppmRange = matrix(c(3.0, 3.1,
                                                   3.6, 3.7,
                                                   3.9, 4.0,
                                                   4.0, 4.1, 
                                                   6.95, 7.05,
                                                   7.6, 7.7,
                                                   7.8, 7.9), ncol = 2, byrow = TRUE),
                               specNo = "1",
                               paraProc = 4L,
                               negThresh = -0.5,
                               scaleFac = 1000000,
                               downSamp = 1,
                               hiresFlag = 1,
                               randSeed = 100025L,
                               nItBurnin = 200L,
                               nItPostBurnin = 5000L,
                               multFile = 2L,
                               thinning = 50L,
                               cfeFlag = 0,
                               nItRerun = 5000L,
                               startTemp = 1000,
                               specFreq = 600,
                               a = 0.00001,
                               b = 0.000000001,
                               muMean = 1.1,
                               muVar = 0.2,
                               muVar_prop = 0.002,
                               nuMVar = 0.0025,
                               nuMVarProp = 0.1,
                               tauMean = -0.05,
                               tauPrec = 2,
                               rdelta = 0.02,
                               csFlag = 0
                               ) {

  # ppmRange:
  ppmRange <- as.matrix(ppmRange)
  stopifnot(ncol(ppmRange) == 2)
  stopifnot(!anyNA(ppmRange))
  
  
  out <- list(ppmRange = ppmRange,
              specNo = specNo, paraProc = paraProc, negThresh = negThresh,
              scaleFac = scaleFac,
              downSamp = downSamp, hiresFlag = hiresFlag,
              randSeed = randSeed, nItBurnin = nItBurnin, nItPostBurnin = nItPostBurnin,
              multFile = multFile,
              thinning = thinning, cfeFlag = cfeFlag, nItRerun = nItRerun,
              startTemp = startTemp, specFreq = specFreq, a = a, b = b,
              muMean = muMean, muVar = muVar, muVar_prop = muVar_prop, nuMVar = nuMVar,
              nuMVarProp = nuMVarProp, tauMean = tauMean, tauPrec = tauPrec, rdelta = rdelta, csFlag = csFlag)
  class(out) <- "batman_options"
  out
}

#' @rdname nmr_batman
#' @export
nmr_batman_write_options <- function(bopts, batman_dir = "BatmanInput", filename = "batmanOptions.txt") {
  # ppmRange:
  full_filename <- batman_get_full_filename(batman_dir, filename)
  bopts2 <- bopts
  ppmRange <- as.matrix(bopts2$ppmRange)
  stopifnot(ncol(ppmRange) == 2)
  stopifnot(!anyNA(ppmRange))
  ppmRange <- paste0("(", apply(ppmRange, 1, function(x) paste(x, collapse = ",")), ")", collapse = " ")
  bopts2$ppmRange <- ppmRange
  # specNo
  specNo <- as.numeric(bopts2$specNo)
  specNo <- paste0(specNo, collapse = ",")
  bopts2$specNo <- specNo
  
  
  bopts_txt <- glue::glue_data(bopts2,
    "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Options file for BATMAN %%%%%%%%%%%%%%%%%%%%%%%%%",
    "%% There should be no empty line in the file, if line                               %%   ",
    "%% is not used, please comment the start of the line with %.                        %%",
    "%% Please do not comment or remove any of the parameter input lines.                %% ",
    "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%",
    "%%",
    "%% General parameters:",
    "%% ppm ranges for analysis: Please put each ppm range in a pair of parentheses in the same line, ",
    "%% separate start and end ppm values with a comma, and ",
    "%% separate each set of ppm range with space, e.g. (-0.05,0.05) (1.44,1.52)   ",
    "ppmRange - ppm ranges for analysis: {ppmRange}",
    "specNo - Ranges of spectra number to be included (e.g. 1,3-4 etc.):  {specNo}",
    "paraProc - No of parallel processes (multicores) (only 1 core will be used for single spectrum): {paraProc}",
    "negThresh - Truncation threshold for negative intensities: {negThresh}",
    "scaleFac - Intensity scale factor: {scaleFac}",
    "downSamp - Down sampling factor: {downSamp}",
    "hiresFlag - Save metabolite fit at resolution of original spectrum? (Yes - 1 / No - 0): {hiresFlag}",
    "randSeed  - Random number seed: {randSeed}  ",
    "nItBurnin - Number of burn-in iterations: {nItBurnin}",
    "nItPostBurnin - Number of post-burn-in iterations: {nItPostBurnin}",
    "multFile - Choose template of multiplets file from options below: {multFile}",
    "%% 1, The default template of multiplets in multi_data.csv file,",
    "%% 2, The user input template of multiplets in multi_data_user.csv file,",
    "%% 3, Both the default and user input template of multiplets files.",
    "%%",
    "thinning - Save MCMC state in every ? iterations: {thinning}",
    "cfeFlag - Same concentration for all spectra (fixed effect)? (Yes - 1 / No - 0): {cfeFlag}",
    "nItRerun - Number of iterations for batmanrerun: {nItRerun}",
    "startTemp - Start temperature: {startTemp}",
    "specFreq - Spectrometer frequency (MHz): {specFreq}",
    "%% ",
    "%% Uncatalogued (wavelet) component: ",
    "%% Hyper parameters for the global precision priors, lambda is inversely proportional to error variance. ",
    "%% lambda ~ Gamma(a,b/2) on wavelet coefficients:",
    "a - Gamma-distributed with shape a: {a}",
    "b - Gamma-distributed with scale b: {b}",
    "%% ",
    "%% Catalogued metabolite component: gamma is the peak width (full width at half max) of metabolite m. ",
    "%% The model for gamma is ln(gamma)= mu + nu_m, where mu is the spectrum-wide average log-peak width  ",
    "%% and nu_m is a random effect for each metabolite deviating from mu.  ",
    "muMean - Mean of prior on global peak width (mu) in ln(Hz): {muMean} ",
    "muVar - Variance of prior on global peak width (mu) in ln(Hz): {muVar}",
    "muVar_prop - Variance of proposal distribution for mu in ln(Hz): {muVar_prop}",
    "%% The mean of each prior on nu_m is 0. ",
    "%% Set variance of prior on peak width offset (nu_m) to 0 to turn off the random effect on peak width:  ",
    "nuMVar - Variance of prior on peak width offset (nu_m) in ln(Hz): {nuMVar}",
    "nuMVarProp - Variance of proposal distribution for nu_m in ln(Hz): {nuMVarProp}",
    "%%",
    "%% Wavelet truncation: tau is a vector of truncation limits, which sets a lower bound on the wavelets",
    "tauMean - mean of the prior on tau: {tauMean}",
    "tauPrec - inverse of variance of prior on tau: {tauPrec}",
    "%%",
    "%% Global prior on peak shift (truncated Gaussian)",
    "%% Note: individual priors for each multiplet can be changed in the multi_data(_user).csv file",
    "rdelta - Truncation of the prior on peak shift (ppm): {rdelta}",
    "csFlag - Specify chemical shift for each multiplet in each spectrum? (chemShiftperSpectra.csv file) (Yes - 1 / No - 0): {csFlag}",
    .sep = "\n")
  writeLines(text = bopts_txt, con = full_filename)
  bopts
}


batman_get_full_filename <- function(batman_dir, filename) {
  if (!dir.exists(batman_dir)) {
    dir.create(batman_dir, recursive = TRUE)
  }
  full_filename <- file.path(batman_dir, filename)
  if (file.exists(full_filename)) {
    stop("File ", full_filename, " already exists. Remove it and try again.")
  }
  full_filename
}

#' @rdname nmr_batman
#' @export
nmr_batman_export_dataset <- function(nmr_dataset, batman_dir = "BatmanInput", filename = "NMRdata.txt") {
  full_filename <- batman_get_full_filename(batman_dir, filename)
  nmrdata <- t(as.matrix(nmr_dataset[["data_1r"]]))
  colnames(nmrdata) <- paste0("injection_id_", nmr_dataset$metadata$injection_id)
  nmrdata <- cbind(ppm = nmr_dataset[["axis"]][[1]], nmrdata)
  utils::write.table(nmrdata, full_filename, row.names = FALSE, sep = "\t")
}

#' @rdname nmr_batman
#' @export
nmr_batman_multi_data_user_hmdb <- function(batman_dir = "BatmanInput", filename = "multi_data_user.csv") {
  full_filename <- batman_get_full_filename(batman_dir, filename)
  utils::data("hmdb", envir = environment())
  hmdb[["overwrite_pos"]] <- -50
  hmdb[["overwrite_truncation"]] <- -50
  hmdb[["Include_multiplet"]] <- 1
  cols <- c("Metabolite", "pos_in_ppm", "couple_code", "J_constant",
            "relative_intensity", "overwrite_pos", "overwrite_truncation",
            "Include_multiplet")
  
  hmdb[["J_constant"]][is.na(hmdb[["J_constant"]])] <- 0 # batman does not accept NA...
  hmdb[["Metabolite"]] <- gsub(pattern = ",", replacement = "_", hmdb[["Metabolite"]])
  utils::write.csv(hmdb[,cols], full_filename, row.names = FALSE)
  hmdb
}

#' @rdname nmr_batman
#' @export
nmr_batman_metabolites_list <- function(metabolite_names,
                                        batman_dir = "BatmanInput",
                                        filename = "metabolitesList.csv") {
  full_filename <- batman_get_full_filename(batman_dir, filename)
  write(unique(metabolite_names), file = full_filename)
  metabolite_names
}






