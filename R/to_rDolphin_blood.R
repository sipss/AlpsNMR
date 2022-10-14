#' to rDolphin
#'
#' Parameters for blood (plasma/serum) samples profiling
#'
#' The template `Parameters_blood` contains the chosen normalization approach (by default, PQN), the Spectometer Frequency (by default, 600.04MHz),
#' alignment (by default, TSP 0.00 ppm), bucket resolution (by default, 0.00023)
#' @name Parameters_blood
#' @docType data
#' @references \url{github.com/danielcanueto/rDolphin}
#' @keywords data
#' @examples
#' data("Parameters_blood")
#' Parameters_blood
NULL

#' ROIs for blood (plasma/serum) samples
#'
#' The template ROI_blood contains the targeted list of metabolites to be quantified (blood samples)
#'
#' @name ROI_blood
#' @docType data
#' @references \url{github.com/danielcanueto/rDolphin}
#' @keywords data
#' @examples
#' data("ROI_blood")
#' ROI_blood[ROI_blood$Metabolite == "Valine", ]
NULL

#' Files to rDoplhin
#'
#' The rDolphin family functions are introduced to perform automatic targeted
#' metabolite profiling. Therefore, ensure that you interpolated from -0.1 ppm
#' in order to consider the TSP/DSS signal at 0.0 ppm. The function generates a
#' list with the files required by to_rDolphin function. Then, it is required
#' to save them with the `save_files_to_rDolphin`. to_rDolphin function will
#' read the generated "parameters.csv" file.
#' function.
#' @param nmr_dataset An [nmr_dataset] object
#' @param biological_origin String specify the type of sample (blood, urine, cell)
#' @return a list containing:
#'    - `meta_rDolphin`: metadata in rDolphin format,
#'    - `NMR_spectra`: spectra matrix
#'    - `ROI`: ROI template
#'    - `Parameters`: parameters file
#' @family import/export functions
#' @examples
#' \dontrun{
#' # Set the directory in which rDolphin files will be saved
#' output_dir_10_rDolphin <- file.path(your_path, "10-rDolphin")
#' fs::dir_create(output_dir_10_rDolphin)
#'
#' # Generate the files (for plasma/serum)
#' files_rDolphin <- files_to_rDolphin(nmr_dataset_0_10_ppm, blood)
#'
#' # Save the files
#' save_files_to_rDolphin(files_rDolphin, output_dir_10_rDolphin)
#'
#' # Build the rDolphin object. Do not forget to set the directory
#' setwd(output_dir_10_rDolphin)
#' rDolphin_object <- to_rDolphin("Parameters.csv")
#'
#' # Visualize your spectra
#' rDolphin_plot(rDolphin_object)
#'
#' # Run the main profiling function (it takes a while)
#' targeted_profiling <- Automatic_targeted_profiling(rDolphin_object)
#'
#' # Save results
#' save_profiling_output(targeted_profiling, output_dir_10_rDolphin)
#'
#' save_profiling_plots(
#'     output_dir_10_rDolphin, targeted_profiling$final_output,
#'     targeted_profiling$reproducibility_data
#' )
#'
#' # Additionally, you can run some stats
#' intensities <- targeted_profiling$final_output$intensity
#' group <- as.factor(rDolphin_object$Metadata$type)
#' model_PLS <- rdCV_PLS_RF(X = intensities, Y = group)
#' }
#' @export
files_to_rDolphin <- function(nmr_dataset, biological_origin) {
    message("you can edit obtained files for better performance in rDolphin")
    meta_D_3col <- c("NMRExperiment", "SubjectID", "Group")
    meta_rDolphin <- AlpsNMR::nmr_meta_get(nmr_dataset)[meta_D_3col]
    newcolnames <- c("sample", "individual", "type")
    colnames(meta_rDolphin) <- newcolnames
    meta_rDolphin$type <- as.numeric(as.factor(meta_rDolphin$type))

    NMR_spectra <- nmr_data(nmr_dataset)
    ROI <- NULL
    Parameters <- NULL
    Parameters_blood <- NULL
    ROI_blood <- NULL
    Parameters_cell <- NULL
    ROI_cell <- NULL
    Parameters_urine <- NULL
    ROI_urine <- NULL
    if (biological_origin == "blood") {
        utils::data("Parameters_blood",
            package = "AlpsNMR",
            envir = environment()
        )
        utils::data("ROI_blood", package = "AlpsNMR", envir = environment())
        Parameters <- Parameters_blood
        ROI <- ROI_blood
    } else if (biological_origin == "cell") {
        utils::data("Parameters_cell",
            package = "AlpsNMR",
            envir = environment()
        )
        utils::data("ROI_cell", package = "AlpsNMR", envir = environment())
        Parameters <- Parameters_cell
        ROI <- ROI_cell
    } else if (biological_origin == "urine") {
        NMR_spectra <- nmr_data(nmr_dataset)
        utils::data("Parameters_urine",
            package = "AlpsNMR",
            envir = environment()
        )
        utils::data("ROI_urine", package = "AlpsNMR", envir = environment())
        Parameters <- Parameters_urine
        ROI <- ROI_urine
    } else {
        message("Unknown biological origin")
        return(NULL)
    }

    files_rDolphin <- list(Parameters, meta_rDolphin, NMR_spectra, ROI)
    names(files_rDolphin) <- c("Parameters", "meta_rDolphin", "NMR_spectra", "ROI")
    return(files_rDolphin)
}

#' ROIs for cell samples
#'
#' The template ROI_cell contains the targeted list of metabolites to be quantified (cell samples)
#'
#' @name ROI_cell
#' @docType data
#' @references \url{github.com/danielcanueto/rDolphin}
#' @keywords data
#' @examples
#' data("ROI_cell")
#' ROI_cell[ROI_cell$Metabolite == "Valine", ]
#'
NULL

#' Parameters for cell samples profiling
#'
#' The template `Parameters_cell` contains the chosen normalization approach (by default, PQN), the Spectometer Frequency (by default, 600.04MHz),
#' alignment (by default, TSP 0.00 ppm), bucket resolution (by default, 0.00023)
#' @name Parameters_cell
#' @docType data
#' @references \url{github.com/danielcanueto/rDolphin}
#' @keywords data
#' @examples
#' data("Parameters_cell")
#' Parameters_cell
NULL

#' ROIs for urine samples
#'
#' The template ROI_urine contains the targeted list of metabolites to be quantified (urine samples)
#'
#' @name ROI_urine
#' @docType data
#' @references \url{github.com/danielcanueto/rDolphin}
#' @keywords data
#' @examples
#' data("ROI_urine")
#' ROI_urine[ROI_urine$Metabolite == "Valine", ]
NULL

#' Parameters for urine samples profiling
#'
#' The template `Parameters_urine` contains the chosen normalization approach (by default, PQN), the Spectometer Frequency (by default, 600.04MHz),
#' alignment (by default, TSP 0.00 ppm), bucket resolution (by default, 0.00023)
#' @name Parameters_urine
#' @docType data
#' @references \url{github.com/danielcanueto/rDolphin}
#' @keywords data
#' @examples
#' data("Parameters_urine")
#' Parameters_urine
NULL

#' Save files to rDoplhin
#'
#' The function saves the CSV files required by to_rDolphin and Automatic_targeted_profiling functions for metabolite profiling.
#' @param files_rDolphin a list containing 4 elements from `files_to_rDolphin`
#'    - `meta_rDolphin`: metadata in rDolphin format,
#'    - `NMR_spectra`: spectra matrix
#'    - `ROI`: ROI template
#'    - `Parameters_blood`: parameters file
#' @param output_directory a directory in which the CSV files are saved
#' @family import/export functions
#' @return CSV files containing:
#' @examples
#' \dontrun{
#' dataset <- system.file("dataset-demo", package = "AlpsNMR")
#' excel_file <- system.file("dataset-demo", "dummy_metadata.xlsx", package = "AlpsNMR")
#' nmr_dataset <- nmr_read_samples_dir(dataset)
#' files_rDolphin <- files_to_rDolphin_blood(nmr_dataset)
#' save_files_to_rDolphin(files_rDolphin, output_directory = ".")
#' }
#' @export
save_files_to_rDolphin <- function(files_rDolphin, output_directory) {
    output_dir_10_rDolphin <- file.path(output_directory)
    fs::dir_create(output_dir_10_rDolphin)

    Parameters <-
        file.path(output_dir_10_rDolphin, "Parameters.csv")
    meta_rDolphin <-
        file.path(output_dir_10_rDolphin, "meta_rDolphin.csv")
    NMR_spectra <-
        file.path(output_dir_10_rDolphin, "NMR_spectra.csv")
    ROI <- file.path(output_dir_10_rDolphin, "ROI.csv")

    utils::write.csv(files_rDolphin$Parameters, Parameters, row.names = FALSE)
    utils::write.csv(files_rDolphin$meta_rDolphin, meta_rDolphin, row.names = FALSE)
    utils::write.csv(files_rDolphin$NMR_spectra, NMR_spectra, row.names = FALSE)
    utils::write.csv(files_rDolphin$ROI, ROI, row.names = FALSE)
}

#' Save rDoplhin output
#'
#' The function saves the output from Automatic_targeted_profiling function in CSV format.
#' @param targeted_profiling A list from Automatic_targeted_profiling function
#' @param output_directory a directory in which the CSV files are saved
#' @family import/export functions
#' @return rDolphin output from Automatic_targeted_profiling function:
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
#' rDolphin_object <- to_rDolphin(parameters)
#' targeted_profiling <- Automatic_targeted_profiling(rDolphin)
#' save_profiling_output(targeted_profiling, output_directory)
#' }
#' @export
save_profiling_output <- function(targeted_profiling, output_directory) {
    output_dir_10_rDolphin <-
        file.path(output_directory, "rDolphin_output")
    fs::dir_create(output_dir_10_rDolphin)

    intensity_fn <-
        file.path(output_dir_10_rDolphin, "metabolites_intensity.csv")
    quantification_fn <-
        file.path(
            output_dir_10_rDolphin,
            "metabolites_quantification.csv"
        )
    ROI_profiles_used_fn <-
        file.path(output_dir_10_rDolphin, "ROI_profiles_used.csv")
    chemical_shift_fn <-
        file.path(output_dir_10_rDolphin, "chemical_shift.csv")
    fitting_error_fn <-
        file.path(output_dir_10_rDolphin, "fitting_error.csv")
    half_bandwidth_fn <-
        file.path(output_dir_10_rDolphin, "half_bandwidth.csv")
    signal_area_ratio_fn <-
        file.path(output_dir_10_rDolphin, "signal_area_ratio.csv")

    intensity <-
        tibble::as_tibble(targeted_profiling[["final_output"]][["intensity"]], rownames = "NMRExperiment")
    quantification <-
        tibble::as_tibble(targeted_profiling[["final_output"]][["quantification"]], rownames = "NMRExperiment")
    ROI_profiles_used <-
        tibble::as_tibble(targeted_profiling[["final_output"]][["ROI_profiles_used"]], rownames = "NMRExperiment")
    chemical_shift <-
        tibble::as_tibble(targeted_profiling[["final_output"]][["chemical_shift"]], rownames = "NMRExperiment")
    fitting_error <-
        tibble::as_tibble(targeted_profiling[["final_output"]][["fitting_error"]], rownames = "NMRExperiment")
    half_bandwidth <-
        tibble::as_tibble(targeted_profiling[["final_output"]][["half_bandwidth"]], rownames = "NMRExperiment")
    signal_area_ratio <-
        tibble::as_tibble(targeted_profiling[["final_output"]][["signal_area_ratio"]], rownames = "NMRExperiment")

    utils::write.csv(intensity, intensity_fn, row.names = FALSE)
    utils::write.csv(quantification, quantification_fn, row.names = FALSE)
    utils::write.csv(ROI_profiles_used, ROI_profiles_used_fn, row.names = FALSE)
    utils::write.csv(chemical_shift, chemical_shift_fn, row.names = FALSE)
    utils::write.csv(fitting_error, fitting_error_fn, row.names = FALSE)
    utils::write.csv(half_bandwidth, half_bandwidth_fn, row.names = FALSE)
    utils::write.csv(signal_area_ratio, signal_area_ratio_fn, row.names = FALSE)
}
