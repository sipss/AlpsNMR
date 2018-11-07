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
  nmr_dataset_rds
}


#' Pipeline: Add Metadata
#'
#' @family pipeline functions
#' @param nmr_dataset_rds The nmr_dataset.rds coming from
#' @param output_dir The output directory for this pipe element
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
#' # Check out output_dir
#' 
pipe_add_metadata <- function(nmr_dataset_rds,
                              output_dir,
                              excel_file) {
  env <- new.env()
  env$nmr_dataset <- nmr_dataset_load(nmr_dataset_rds)
  env$excel_file <- excel_file
  env$xlsx_file <- as.character(fs::path(output_dir, "nmr_metadata_added.xlsx"))
  env$nmr_dataset_outfile <- as.character(fs::path(output_dir, "nmr_dataset.rds"))
  rmd_file <- system.file("pipeline-rmd", "load-metadata.Rmd", package = "NIHSnmr")
  rmarkdown::render(input = rmd_file, output_dir = output_dir, envir = env)
  message("Done")
}


