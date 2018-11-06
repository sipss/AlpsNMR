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
  NMRExperiments <- as.character(fs::dir_ls(samples_dir, type = "directory", glob = glob))
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
#' @param excel_files A character vector with excel file names
#' @param key_columns A character vector with the key column to merge. See Details
#' 
#' In the simplest case we will want to add an excel file with metadata to our
#' dataset. The excel file needs to be tidy and this means: A single excel sheet
#' with a single header column at the first row. This excel file will have
#' an "NMRExperiment" column (please respect the exact name) with values like
#' "10", "20", "30"... referring to the NMR samples and as many other columns
#' as desired. The `key_columns` will be `"NMRExperiment"` as it is the common
#' column present both in the dataset and in the excel file.
#' 
#' A more useful scenario is when you provide several excel files: After the
#' first file has been added as described above, the dataset will have
#' additional metadata columns, for instance a `SubjectID` column. The second
#' excel file may use the `SubjectID` as `key_column` so the criteria to add
#' the metadata will be based on matching `SubjectID`
#'
#' @return This function saves the result to the output directory
#' @export
#'
#' @examples
#' 
#' # Load a demo dataset with four samples:
#' dataset <- system.file("dataset-demo", package = "NIHSnmr")
#' nmr_dataset <- nmr_read_samples_dir(dataset)
#' 
#' # At first we just have the NMRExperiment column
#' print(nmr_dataset$metadata$external)
#' # The first excel file contains an "NMRExperiment" column:
#' first_excel_fn <- system.file("dataset-demo", "first_excel.xlsx", package = "NIHSnmr")
#' first_excel <- readxl::read_excel(first_excel_fn)
#' print(first_excel)
#' # We can link the SubjectID column of the first excel into the dataset
#' nmr_dataset <- nmr_add_metadata(nmr_dataset, first_excel, by = "NMRExperiment")
#' print(nmr_dataset$metadata$external)
#' # The second excel can use the SubjectID:
#' second_excel_fn <- system.file("dataset-demo", "second_excel.xlsx", package = "NIHSnmr")
#' second_excel <- readxl::read_excel(second_excel_fn)
#' print(second_excel)
#' # Add the metadata by its SubjectID:
#' nmr_dataset <- nmr_add_metadata(nmr_dataset, second_excel, by = "SubjectID")
#' # The final loaded metadata:
#' print(nmr_dataset$metadata$external)
#' 
#' # This example is equivalent to:
#' \dontrun{
#' nmr_dataset_rds <- tempfile(fileext = ".rds")
#' nmr_dataset_save(nmr_dataset, nmr_dataset_rds)
#' output_dir <- tempdir()
#' excel_files <-  c(first_excel_fn, second_excel_fn)
#' key_columns <- c("NMRExperiment", "SubjectID")
#' pipe_add_metadata(nmr_dataset_rds = a_tmp_file, output_dir = output_dir,
#'    excel_files = excel_files, key_columns = key_columns)
#' }
pipe_add_metadata <- function(nmr_dataset_rds,
                              output_dir,
                              excel_files = character(0),
                              key_columns = character(0)) {
  num_excel_files <- length(excel_files)
  assertthat::assert_that(num_excel_files == length(key_columns),
                          msg = "Each excel file needs to have a merge column name")
  message("Loading all the excel files")
  env <- new.env()
  env$nmr_dataset <- nmr_dataset_load(nmr_dataset_rds)
  env$excel_files <- excel_files
  env$key_columns <- key_columns
  env$xlsx_file <- as.character(fs::path(output_dir, "nmr_metadata_added.xlsx"))
  rmd_file <- system.file("pipeline-rmd", "load-metadata.Rmd", package = "NIHSnmr")
  rmarkdown::render(input = rmd_file, output_dir = output_dir, envir = env)
  message("Done")
}


