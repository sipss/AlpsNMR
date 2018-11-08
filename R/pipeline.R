#' Pipeline: Load NMR samples
#'
#' @family pipeline functions
#' @param samples_dir The directory where the samples are
#' @param output_dir Directory where the nmr_dataset and the excel files will be saved
#' @inheritParams fs::dir_ls
#'
#' @return This function saves the result to the output directory
#' @export
#'
#' @examples
#' \dontrun{
#' pipe_load_samples("/dir/with/nmr/samples/", "/dir/to/save/output", glob = "*0")
#' }
pipe_load_samples <- function(samples_dir, output_dir, glob = "*0") {
  fs::dir_create(output_dir)
  NMRExperiments <- as.character(fs::dir_ls(samples_dir, glob = glob))
  nmr_dataset <- nmr_read_samples(NMRExperiments)
  nmr_dataset_rds <- fs::path(output_dir, "nmr_dataset.rds")
  nmr_dataset_save(nmr_dataset, nmr_dataset_rds)
  nmr_export_metadata(nmr_dataset, fs::path(output_dir, "nmr_dataset_metadata.xlsx"))
  message(nmr_dataset$num_samples, " samples loaded.")
}


#' Pipeline: Add Metadata
#'
#' @family pipeline functions
#' @param nmr_dataset_rds The nmr_dataset.rds file name coming from previous nodes
#' @param excel_file An excel file name. See details for the requirements
#' 
#' 
#' The excel file can have one or more sheets. The excel sheets need to be as
#' simple as possible: One header column on the first row and values below.
#' 
#' Each of the sheets contain metadata that has to be integrated. The merge
#' (technically a left join) is done using the first column of each sheet as key.
#' 
#' In practical terms this means that the first sheet of the excel file MUST 
#' start with an "NMRExperiment" column, and as many additional columns to add
#' (e.g. FluidXBarcode, SampleCollectionDate,  TimePoint and SubjectID).
#' 
#' The second sheet can have as the first column any of the already added columns,
#' for instance the "SubjectID", and any additional columns (e.g. Gender, Age).
#' 
#' The first column on each sheet, named the key column, MUST have unique values.
#' For instance, a sheet starting with "SubjectID" MUST specify each subject ID
#' only once (without repetitions).
#' @param output_dir The output directory for this pipe element
#' 
#' @return This function saves the result to the output directory
#' @export
#'
#' @examples
#' dataset <- system.file("dataset-demo", package = "NIHSnmr")
#' excel_file <- system.file("dataset-demo", "dummy_metadata.xlsx", package = "NIHSnmr")
#' nmr_dataset <- nmr_read_samples_dir(dataset)
#' nmr_dataset_rds <- tempfile(fileext = ".rds")
#' nmr_dataset_save(nmr_dataset, nmr_dataset_rds)
#' output_dir <- tempdir()
#' pipe_add_metadata(nmr_dataset_rds = nmr_dataset_rds, output_dir = output_dir,
#'                   excel_file = excel_file)
#' # Check out: output_dir
#' 
pipe_add_metadata <- function(nmr_dataset_rds, excel_file, output_dir) {
  fs::dir_create(output_dir)
  env <- new.env()
  env$nmr_dataset <- nmr_dataset_load(nmr_dataset_rds)
  env$excel_file <- excel_file
  env$xlsx_file <- as.character(fs::path(output_dir, "nmr_metadata_added.xlsx"))
  env$nmr_dataset_outfile <- as.character(fs::path(output_dir, "nmr_dataset.rds"))
  rmd_file <- system.file("pipeline-rmd", "add-metadata.Rmd", package = "NIHSnmr")
  rmarkdown::render(input = rmd_file, output_dir = output_dir, envir = env)
  message("Add metadata completed")
}


#' Pipeline: Interpolate 1D samples
#'
#' @family pipeline functions
#' @inheritParams pipe_add_metadata
#' @inheritParams nmr_interpolate_1D 
#'
#' @return This function saves the result to the output directory
#' @export
pipe_interpolate_1D <- function(nmr_dataset_rds, axis1, output_dir) {
  fs::dir_create(output_dir)
  nmr_dataset <- nmr_dataset_load(nmr_dataset_rds)
  
  metadata_fn <- file.path(output_dir, "metadata.xlsx")
  raw_data_matrix_fn <- file.path(output_dir, "raw_data.csv")
  nmr_dataset_outfile <- file.path(output_dir, "nmr_dataset.rds")
  plot_html <- file.path(output_dir, "plot-samples.html")
  
  nmr_dataset <- nmr_interpolate_1D(nmr_dataset, axis1 = axis1)
  
  nmr_export_data_1r(nmr_dataset, raw_data_matrix_fn)
  nmr_export_metadata(nmr_dataset, metadata_fn, groups = "external")
  nmr_dataset_save(nmr_dataset, nmr_dataset_outfile)
  plot_webgl(nmr_dataset, html_filename = plot_html)
  
  message("Interpolation finished")
}




#' Pipeline: Exclude regions
#'
#' @family pipeline functions
#' @inheritParams pipe_add_metadata
#' @inheritParams nmr_exclude_region
#'
#' @return This function saves the result to the output directory
#' @export
#'
pipe_exclude_regions <- function(nmr_dataset_rds,
                                 exclude,
                                 output_dir) {
  fs::dir_create(output_dir)
  
  metadata_fn <- file.path(output_dir, "metadata.xlsx")
  raw_data_matrix_fn <- file.path(output_dir, "raw_data.csv")
  nmr_dataset_outfile <- file.path(output_dir, "nmr_dataset.rds")
  plot_html <- file.path(output_dir, "plot-samples.html")
  
  
  nmr_dataset <- nmr_dataset_load(nmr_dataset_rds)
  nmr_dataset <- nmr_exclude_region(nmr_dataset, exclude = exclude)
  
  nmr_export_data_1r(nmr_dataset, raw_data_matrix_fn)
  nmr_export_metadata(nmr_dataset, metadata_fn, groups = "external")
  nmr_dataset_save(nmr_dataset, nmr_dataset_outfile)
  plot_webgl(nmr_dataset, html_filename = plot_html)
  
  message("Regions excluded")
}

#' Pipeline: Peak detection and Alignment
#' @inheritParams pipe_add_metadata
#' @inheritParams nmr_detect_peaks
#' @inheritParams nmr_align
#'
#' @export
#'
pipe_peakdet_align <- function(nmr_dataset_rds,
                         nDivRange = 128, scales = seq(1, 16, 2),
                         baselineThresh = 0.01, SNR.Th = -1,
                         maxShift = 3, acceptLostPeak = FALSE,
                         output_dir = NULL) {
  if (is.null(output_dir)) {
    stop("An output directory must be specified")
  }
  
  fs::dir_create(output_dir)
  
  metadata_fn <- file.path(output_dir, "metadata.xlsx")
  raw_data_matrix_fn <- file.path(output_dir, "raw_data.csv")
  nmr_dataset_outfile <- file.path(output_dir, "nmr_dataset.rds")
  plot_html <- file.path(output_dir, "plot-samples.html")
  peak_data_fn <- file.path(output_dir, "peak_data.csv")
  
  
  nmr_dataset <- nmr_dataset_load(nmr_dataset_rds)
  peak_data <- nmr_detect_peaks(nmr_dataset,
                                nDivRange = nDivRange,
                                scales = scales,
                                baselineThresh = baselineThresh,
                                SNR.Th = SNR.Th)
  nmr_dataset <- nmr_align(nmr_dataset, peak_data, maxShift = maxShift,
                           acceptLostPeak = acceptLostPeak)
  
  nmr_export_data_1r(nmr_dataset, raw_data_matrix_fn)
  nmr_export_metadata(nmr_dataset, metadata_fn, groups = "external")
  nmr_dataset_save(nmr_dataset, nmr_dataset_outfile)
  plot_webgl(nmr_dataset, html_filename = plot_html)
  utils::write.csv(peak_data, peak_data_fn)
  
  message("Peaks detected and spectra aligned")
}
