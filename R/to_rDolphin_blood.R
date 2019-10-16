#' to rDolphin
#' 
#' Import an spectra object from the rDolphin package using the information provided by the `parameters` file.
#' The obtained rDolphin object can be used to extract the metabolite profiling using `Automatic_targeted_profiling` function.
#' @importFrom rDolphin import_data
#' @param parameters the path in which the parameters CSV file is stored to create an rDolphin object
#' @return an `rDolphin_object`
#' @family import/export functions
#' @family nmr_dataset_1D functions
#' @examples 
#' \dontrun{
#' ## After running AlpsNMR and get aligned spectra along with the corresponding metadata we run:
#' library(AlpsNMR)
#' setwd("your_directory_with_parameters.csv")
#' rDolphin_object = to_rDolphin(parameters.csv)
#' targeted_profiling = Automatic_targeted_profiling(rDolphin_object)
#' save_profiling_output(targeted_profiling, output_directory)
#' }
#' @export
to_rDolphin <- function(...) {
UseMethod("to_rDolphin")
}
to_rDolphin <- function (parameters){
  rDolphin::import_data(parameters)
}

NULL
#' Automatic targeted profiling
#'
#' Automatic quantification of metabolites for all experiments using the information located in the ROI patterns file
#' @importFrom rDolphin automatic_profiling
#' @param imported_data an `rDolphin_object` created with `to_rDolphin` function
#' @param ROI a ROI file containing a targeted list of metabolites to fit
#' @param optimization By default TRUE. If TRUE, signals parameters are optimized for profiling quality
#' @param spectra_to_profile By default NA. If NA, all spectra are considered. Otherwise, a vector of selected spectra.
#' @return A Spectra object from the rDolphin package
#' @family import/export functions
#' @family nmr_dataset_1D functions
#' @return A list with `final_output` containing the ROI-template-metabolites, their intensities, their relative quantification and quality indicators,
#'  and `reproducibility_data` (to reproduce the automatic profiling). Use `save_profiling_output` to save the these files into your computer.
#' @examples 
#' \dontrun{
#' ## Run after apply [to_rDolphin] and get an rDolphin spectra object:
#' rDolphin_object = to_rDolphin(parameters)
#' targeted_profiling = Automatic_targeted_profiling(rDolphin_object)
#' save_profiling_output(targeted_profiling, output_directory)
#' }
#' @export
Automatic_targeted_profiling <- function(...) {
UseMethod("Automatic_targeted_profiling")
}
Automatic_targeted_profiling= function (imported_data, ROI=imported_data$ROI, optimization = TRUE, spectra_to_profile = NULL){
  rDolphin::automatic_profiling(imported_data, ROI=imported_data$ROI, optimization = TRUE, spectra_to_profile = NULL)
}

NULL
#' Parameters for blood (plasma/serum) samples profiling
#' 
#' The template `Parameters_blood` contains the chosen normalization approach (by default, PQN), the Spectometer Frequency (by default, 600.04MHz),
#' alignment (by default, TSP 0.00 ppm), bucket resolution (by default, 0.00023)
#' @name Parameters_blood
#' @docType data
#' @references \url{github.com/danielcanueto/rDolphin}
#' @keywords data
NULL

#' ROIs for blood (plasma/serum) samples
#' 
#' The template ROI_blood contains the targeted list of metabolites to be quantified (blood samples)
#'
#' @name ROI_blood
#' @docType data
#' @references \url{github.com/danielcanueto/rDolphin}
#' @keywords data

NULL
#' Files to rDoplhin (blood)
#'
#' The function generates the files required by [to_rDolphin] function.
#' @param nmr_dataset An [nmr_dataset] object
#' @return a list containing: 
#'  - `meta_rDolphin`: metadata in rDolphin format, 
#'  - `NMR_spectra`: spectra matrix
#'  - `ROI_blood`: ROI template
#'  - `Parameters_blood`: parameters file
#' @family import/export functions
#' @family nmr_dataset_1D functions
#' @family to_rDolphin_blood functions
#' @examples 
#' \dontrun{
#' ## library(AlpsNMR)
#' files_to_rDolphin_blood(nmr_dataset)
#' }
#' @export
files_to_rDolphin_blood = function (nmr_dataset){
  message("you can edit obtained files for better performance in rDolphin")
  meta_D_3col= c("NMRExperiment", "SubjectID", "Group")
  meta_rDolphin= AlpsNMR::nmr_meta_get(nmr_dataset)[meta_D_3col]
  newcolnames=c("sample", "individual", "type")
  colnames(meta_rDolphin)=newcolnames
  meta_rDolphin$type=as.numeric(as.factor(meta_rDolphin$type))
  
  NMR_spectra = nmr_data(nmr_dataset)
  Parameters_blood = NULL
  utils::data("Parameters_blood", package = "AlpsNMR", envir = environment())
  ROI_blood = NULL
  utils::data("ROI_blood", package = "AlpsNMR", envir = environment())
  
  files_rDolphin=list(Parameters_blood, meta_rDolphin, NMR_spectra, ROI_blood)
  names(files_rDolphin)=c("Parameters", "meta_rDolphin", "NMR_spectra", "ROI")
  return(files_rDolphin)
  
}
NULL

#' Files to rDoplhin (cell)
#'
#' The function generates the files required by [to_rDolphin] function.
#' @param nmr_dataset An [nmr_dataset] object
#' @return a list containing: 
#'  - `meta_rDolphin`: metadata in rDolphin format, 
#'  - `NMR_spectra`: spectra matrix
#'  - `ROI_cell`: ROI template
#'  - `Parameters_cell`: parameters file
#' @family import/export functions
#' @family nmr_dataset_1D functions
#' @family to_rDolphin_cell functions
#' @examples 
#' \dontrun{
#' ## library(AlpsNMR)
#' files_to_rDolphin_cell(nmr_dataset)
#' }
#' @export
files_to_rDolphin_cell = function (nmr_dataset){
  message("you can edit obtained files for better performance in rDolphin")
  meta_D_3col= c("NMRExperiment", "SubjectID", "Group")
  meta_rDolphin= AlpsNMR::nmr_meta_get(nmr_dataset)[meta_D_3col]
  newcolnames=c("sample", "individual", "type")
  colnames(meta_rDolphin)=newcolnames
  meta_rDolphin$type=as.numeric(as.factor(meta_rDolphin$type))
  
  NMR_spectra = nmr_data(nmr_dataset)
  Parameters_cell = NULL
  utils::data("Parameters_cell", package = "AlpsNMR", envir = environment())
  ROI_cell = NULL
  utils::data("ROI_cell", package = "AlpsNMR", envir = environment())
  
  files_rDolphin=list(Parameters_cell, meta_rDolphin, NMR_spectra, ROI_cell)
  names(files_rDolphin)=c("Parameters", "meta_rDolphin", "NMR_spectra", "ROI")
  return(files_rDolphin)
  
}
NULL

#' ROIs for cell samples
#' 
#' The template ROI_cell contains the targeted list of metabolites to be quantified (cell samples)
#'
#' @name ROI_cell
#' @docType data
#' @references \url{github.com/danielcanueto/rDolphin}
#' @keywords data

NULL

#' Parameters for cell samples profiling
#' 
#' The template `Parameters_cell` contains the chosen normalization approach (by default, PQN), the Spectometer Frequency (by default, 600.04MHz),
#' alignment (by default, TSP 0.00 ppm), bucket resolution (by default, 0.00023)
#' @name Parameters_cell
#' @docType data
#' @references \url{github.com/danielcanueto/rDolphin}
#' @keywords data
NULL

NULL

#' Files to rDoplhin (urine)
#'
#' The function generates the files required by [to_rDolphin] function.
#' @param nmr_dataset An [nmr_dataset] object
#' @return a list containing: 
#'  - `meta_rDolphin`: metadata in rDolphin format, 
#'  - `NMR_spectra`: spectra matrix
#'  - `URINE_cell`: ROI template
#'  - `Parameters_cell`: parameters file
#' @family import/export functions
#' @family nmr_dataset_1D functions
#' @family to_rDolphin_cell functions
#' @examples 
#' \dontrun{
#' ## library(AlpsNMR)
#' files_to_rDolphin_urine(nmr_dataset)
#' }
#' @export
files_to_rDolphin_urine = function (nmr_dataset){
  message("you can edit obtained files for better performance in rDolphin")
  meta_D_3col= c("NMRExperiment", "SubjectID", "Group")
  meta_rDolphin= AlpsNMR::nmr_meta_get(nmr_dataset)[meta_D_3col]
  newcolnames=c("sample", "individual", "type")
  colnames(meta_rDolphin)=newcolnames
  meta_rDolphin$type=as.numeric(as.factor(meta_rDolphin$type))
  
  NMR_spectra = nmr_data(nmr_dataset)
  Parameters_cell = NULL
  utils::data("Parameters_cell", package = "AlpsNMR", envir = environment())
  ROI_urine = NULL
  utils::data("ROI_urine", package = "AlpsNMR", envir = environment())
  
  files_rDolphin=list(Parameters_cell, meta_rDolphin, NMR_spectra, ROI_urine)
  names(files_rDolphin)=c("Parameters", "meta_rDolphin", "NMR_spectra", "ROI")
  return(files_rDolphin)
  
}
NULL

#' ROIs for urine samples
#' 
#' The template ROI_urine contains the targeted list of metabolites to be quantified (urine samples)
#'
#' @name ROI_urine
#' @docType data
#' @references \url{github.com/danielcanueto/rDolphin}
#' @keywords data
NULL

#' Parameters for urine samples profiling
#' 
#' The template `Parameters_urine` contains the chosen normalization approach (by default, PQN), the Spectometer Frequency (by default, 600.04MHz),
#' alignment (by default, TSP 0.00 ppm), bucket resolution (by default, 0.00023)
#' @name Parameters_urine
#' @docType data
#' @references \url{github.com/danielcanueto/rDolphin}
#' @keywords data
NULL

#' Save files to rDoplhin
#'
#' The function saves the CSV files required by [to_rDolphin] and [Automatic_targeted_profiling] functions for metabolite profiling.
#' @param files_rDolphin a list containing 4 elements from `files_to_rDolphin`
#'  - `meta_rDolphin`: metadata in rDolphin format, 
#'  - `NMR_spectra`: spectra matrix
#'  - `ROI`: ROI template
#'  - `Parameters_blood`: parameters file
#' @param output_directory a directory in which the CSV files are saved
#' @family import/export functions
#' @family nmr_dataset_1D functions
#' @family to_rDolphin_blood functions
#' @return CSV files containing: 
#' @examples 
#' \dontrun{
#' ## library(AlpsNMR)
#' files_rDolphin = files_to_rDolphin_blood(nmr_dataset)
#' output_directory = "C:/directory..."
#' save_files_to_rDolphin(files_rDolphin, output_directory)
#' }
#' @export
save_files_to_rDolphin = function (files_rDolphin, output_directory){
  
  output_dir_10_rDolphin <- file.path(output_directory)
  fs::dir_create(output_dir_10_rDolphin)
  
  Parameters <- file.path(output_dir_10_rDolphin, "Parameters.csv")
  meta_rDolphin <- file.path(output_dir_10_rDolphin, "meta_rDolphin.csv")
  NMR_spectra <- file.path(output_dir_10_rDolphin, "NMR_spectra.csv")
  ROI <- file.path(output_dir_10_rDolphin, "ROI.csv")
  
  utils::write.csv(files_rDolphin$Parameters, Parameters, row.names = FALSE)
  utils::write.csv(files_rDolphin$meta_rDolphin, meta_rDolphin, row.names = FALSE)
  utils::write.csv(files_rDolphin$NMR_spectra, NMR_spectra, row.names = FALSE)
  utils::write.csv(files_rDolphin$ROI, ROI, row.names = FALSE)
  
}

NULL

#' Save rDoplhin output
#'
#' The function saves the output from [Automatic_targeted_profiling] function in CSV format.
#' @param targeted_profiling A list from [Automatic_targeted_profiling] function
#' @param output_directory a directory in which the CSV files are saved
#' @family import/export functions
#' @family to_rDolphin_blood functions
#' @return rDolphin output from [Automatic_targeted_profiling] function:
#' - metabolites_intensity
#' - metabolites_quantification
#' - ROI_profiles_used
#' - chemical_shift
#' - fitting_error
#' - half_bandwidth
#' - signal_area_ratio
#' 
#' @examples 
#' \dontrun{
#' ## library(AlpsNMR)
#' rDolphin_object = to_rDolphin(parameters)
#' targeted_profiling = Automatic_targeted_profiling(rDolphin)
#' save_profiling_output(targeted_profiling, output_directory)
#' }
#' @export
save_profiling_output = function (targeted_profiling, output_directory){
  
  output_dir_10_rDolphin <- file.path(output_directory, "rDolphin_output")
  fs::dir_create(output_dir_10_rDolphin)
  
  intensity_fn <- file.path(output_dir_10_rDolphin, "metabolites_intensity.csv")
  quantification_fn <- file.path(output_dir_10_rDolphin, "metabolites_quantification.csv")
  ROI_profiles_used_fn <- file.path(output_dir_10_rDolphin, "ROI_profiles_used.csv")
  chemical_shift_fn <- file.path(output_dir_10_rDolphin, "chemical_shift.csv")
  fitting_error_fn <- file.path(output_dir_10_rDolphin, "fitting_error.csv")
  half_bandwidth_fn <- file.path(output_dir_10_rDolphin, "half_bandwidth.csv")
  signal_area_ratio_fn <- file.path(output_dir_10_rDolphin, "signal_area_ratio.csv")
  
  intensity <- tibble::as_tibble(targeted_profiling[["final_output"]][["intensity"]], rownames = "NMRExperiment")
  quantification <- tibble::as_tibble(targeted_profiling[["final_output"]][["quantification"]], rownames = "NMRExperiment")
  ROI_profiles_used <- tibble::as_tibble(targeted_profiling[["final_output"]][["ROI_profiles_used"]], rownames = "NMRExperiment")
  chemical_shift <- tibble::as_tibble(targeted_profiling[["final_output"]][["chemical_shift"]], rownames = "NMRExperiment")
  fitting_error <- tibble::as_tibble(targeted_profiling[["final_output"]][["fitting_error"]], rownames = "NMRExperiment")
  half_bandwidth <- tibble::as_tibble(targeted_profiling[["final_output"]][["half_bandwidth"]], rownames = "NMRExperiment")
  signal_area_ratio <- tibble::as_tibble(targeted_profiling[["final_output"]][["signal_area_ratio"]], rownames = "NMRExperiment")
  
  utils::write.csv(intensity, intensity_fn, row.names = FALSE)
  utils::write.csv(quantification, quantification_fn, row.names = FALSE)
  utils::write.csv(ROI_profiles_used, ROI_profiles_used_fn, row.names = FALSE)
  utils::write.csv(chemical_shift, chemical_shift_fn, row.names = FALSE)
  utils::write.csv(fitting_error, fitting_error_fn, row.names = FALSE)
  utils::write.csv(half_bandwidth, half_bandwidth_fn, row.names = FALSE)
  utils::write.csv(signal_area_ratio, signal_area_ratio_fn, row.names = FALSE)
}

NULL
#' rDolphin plot
#'
#' The function provides an interactive plot with the spectra of random samples.
#' @param rDolphin_object an `rDolphin_object` created with `to_rDolphin` function
#' @return rDolphin plot
#' @examples 
#' \dontrun{
#' ## library(AlpsNMR)
#' rDolphin_object = to_rDolphin(parameters)
#' rDolphin_plot(rDolphin_object)
#' }
#' @export
rDolphin_plot = function (rDolphin_object){
rDolphin::exemplars_plot(rDolphin_object)
}

NULL

#' Save rDolphin output and profiling plots
#' 
#' The function writes the plots of the adjusted compounds.
#' 
#' @inheritParams rDolphin::write_plots
#'
#' @return Plots of the fitted compounds in pdf files
#' @export
#'
#' @examples
#' \dontrun{
#' ## library(AlpsNMR)
#' rDolphin_object = to_rDolphin(parameters)
#' targeted_profiling = Automatic_targeted_profiling(rDolphin_object)
#' write_plots(dir_10_rDolphin, targeted_profiling$final_output,
#' targeted_profiling$reproducibility_data)
#' }
#' 
save_profiling_plots = function(...){
rDolphin::write_plots(...)
}












