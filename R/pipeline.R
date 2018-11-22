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
#' pipe_load_samples("/dir/with/nmr/samples/",
#'                   glob = "*0",
#'                   output_dir = "/dir/to/save/output")
#' }
pipe_load_samples <- function(samples_dir, glob = "*0", output_dir = NULL) {
  message("Starting pipe_load_samples at ", Sys.time())
  
  if (is.null(output_dir)) {
    stop("An output directory must be specified")
  }
  
  fs::dir_create(output_dir)
  NMRExperiments <- as.character(fs::dir_ls(samples_dir, glob = glob))
  nmr_dataset <- nmr_read_samples(NMRExperiments)
  nmr_dataset_rds <- fs::path(output_dir, "nmr_dataset.rds")
  nmr_dataset_save(nmr_dataset, nmr_dataset_rds)
  nmr_meta_export(nmr_dataset, fs::path(output_dir, "nmr_dataset_metadata.xlsx"))
  message(nmr_dataset$num_samples, " samples loaded.")
  message("Ending pipe_load_samples at ", Sys.time())
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
  message("Starting pipe_add_metadata at ", Sys.time())

  if (is.null(output_dir)) {
    stop("An output directory must be specified")
  }
  
  fs::dir_create(output_dir)
  env <- new.env()
  env$nmr_dataset <- nmr_dataset_load(nmr_dataset_rds)
  env$excel_file <- excel_file
  env$xlsx_file <- as.character(fs::path(output_dir, "nmr_metadata_added.xlsx"))
  env$nmr_dataset_outfile <- as.character(fs::path(output_dir, "nmr_dataset.rds"))
  rmd_file <- system.file("pipeline-rmd", "add-metadata.Rmd", package = "NIHSnmr")
  rmarkdown::render(input = rmd_file, output_dir = output_dir, envir = env)
  message("Add metadata completed")
  message("Ending pipe_add_metadata at ", Sys.time())
}


#' Pipeline: Interpolate 1D samples
#'
#' @family pipeline functions
#' @inheritParams pipe_add_metadata
#' @inheritParams nmr_interpolate_1D 
#'
#' @return This function saves the result to the output directory
#' @export
pipe_interpolate_1D <- function(nmr_dataset_rds, axis, output_dir) {
  message("Starting pipe_interpolate_1D at ", Sys.time())
  
  if (is.null(output_dir)) {
    stop("An output directory must be specified")
  }
  
  fs::dir_create(output_dir)
  nmr_dataset <- nmr_dataset_load(nmr_dataset_rds)
  
  metadata_fn <- file.path(output_dir, "metadata.xlsx")
  raw_data_matrix_fn <- file.path(output_dir, "raw_data.csv")
  nmr_dataset_outfile <- file.path(output_dir, "nmr_dataset.rds")
  plot_html <- file.path(output_dir, "plot-samples.html")
  
  nmr_dataset <- nmr_interpolate_1D(nmr_dataset, axis = axis)
  
  nmr_export_data_1r(nmr_dataset, raw_data_matrix_fn)
  nmr_meta_export(nmr_dataset, metadata_fn, groups = "external")
  nmr_dataset_save(nmr_dataset, nmr_dataset_outfile)
  plot_webgl(nmr_dataset, html_filename = plot_html)
  
  message("Ending pipe_interpolate_1D at ", Sys.time())
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
  message("Starting pipe_exclude_regions at ", Sys.time())

  if (is.null(output_dir)) {
    stop("An output directory must be specified")
  }
  
  fs::dir_create(output_dir)
  
  metadata_fn <- file.path(output_dir, "metadata.xlsx")
  raw_data_matrix_fn <- file.path(output_dir, "raw_data.csv")
  nmr_dataset_outfile <- file.path(output_dir, "nmr_dataset.rds")
  plot_html <- file.path(output_dir, "plot-samples.html")
  
  
  nmr_dataset <- nmr_dataset_load(nmr_dataset_rds)
  nmr_dataset <- nmr_exclude_region(nmr_dataset, exclude = exclude)
  
  nmr_export_data_1r(nmr_dataset, raw_data_matrix_fn)
  nmr_meta_export(nmr_dataset, metadata_fn, groups = "external")
  nmr_dataset_save(nmr_dataset, nmr_dataset_outfile)
  plot_webgl(nmr_dataset, html_filename = plot_html)
  
  message("Ending pipe_exclude_regions at ", Sys.time())
}


#' Pipeline: Remove blatant outliers
#'
#' Uses [nmr_pca_outliers_robust] to perform the detection of outliers
#' 
#' @inheritParams pipe_add_metadata
#'
#' @return This function saves the result to the output directory
#' @export
#'
pipe_outlier_detection <- function(nmr_dataset_rds, output_dir)  {
  message("Starting pipe_outlier_detection at ", Sys.time())
  if (is.null(output_dir)) {
    stop("An output directory must be specified")
  }
  
  fs::dir_create(output_dir)
  
  metadata_fn <- file.path(output_dir, "metadata.xlsx")
  full_spectra_matrix_fn <- file.path(output_dir, "full_spectra_matrix.csv")
  nmr_dataset_outfile <- file.path(output_dir, "nmr_dataset.rds")
  plot_outlier_html <- file.path(output_dir, "plot-outlier-samples.html")
  plot_outlier_QT2 <- file.path(output_dir, "plot-Qresiduals_vs_Tscores.png")
  
  
  nmr_dataset <- nmr_dataset_load(nmr_dataset_rds)

  pca_outliers <- nmr_pca_outliers_robust(nmr_dataset)
  gplt <- nmr_pca_outliers_plot(nmr_dataset, pca_outliers)
  
  nmr_dataset_no_out <- nmr_pca_outliers_filter(nmr_dataset, pca_outliers)
  nmr_exp_all <- nmr_meta_get_column(nmr_dataset, "NMRExperiment")
  nmr_exp_noout <- nmr_meta_get_column(nmr_dataset_no_out, "NMRExperiment")
  nmr_exp_out <- setdiff(nmr_exp_all, nmr_exp_noout)
  
  nmr_export_data_1r(nmr_dataset_no_out, full_spectra_matrix_fn)
  nmr_meta_export(nmr_dataset_no_out, metadata_fn, groups = "external")
  nmr_dataset_save(nmr_dataset_no_out, nmr_dataset_outfile)
  ggplot2::ggsave(filename = plot_outlier_QT2, plot = gplt)
  
  if (length(nmr_exp_out) > 0) {
    message("The following NMRExperiments have been flagged and excluded as outliers:\n",
            glue::glue_collapse(nmr_exp_out, sep = ", ", last = " and "))
    plot_webgl(nmr_dataset, NMRExperiment = nmr_exp_out, 
               quantile_plot = TRUE, html_filename = plot_outlier_html)
  } else {
    message("No outlier detected on a first unscaled PCA (further outliers may be detected later)")
  }

  message("Ending pipe_outlier_detection at ", Sys.time())
}

#' Pipeline: Filter samples according to metadata conditions
#'
#' @inheritParams pipe_add_metadata
#' @param conditions A character vector with conditions to filter metadata.
#' 
#' The `conditions` parameter should be a character vector of valid R logical conditions.
#' Some examples:
#' 
#' - conditions <- 'Gender == "Female"'
#' - conditions <- 'Cohort == "Chuv"'
#' - conditions <- 'TimePoint %in% c("T0", "T31")'
#' - conditions <- c(Cohort == "Chuv", 'TimePoint %in% c("T0", "T31")')
#'
#' Only samples fullfilling all the given conditions are kept in further analysis.
#' 
#' @export
#'
pipe_filter_samples <- function(nmr_dataset_rds, conditions, output_dir) {
  message("Starting pipe_filter_samples at ", Sys.time())
  
  if (is.null(output_dir)) {
    stop("An output directory must be specified")
  }
  
  fs::dir_create(output_dir)
  
  metadata_fn <- file.path(output_dir, "metadata.xlsx")
  raw_data_matrix_fn <- file.path(output_dir, "raw_data.csv")
  nmr_dataset_outfile <- file.path(output_dir, "nmr_dataset.rds")
  
  nmr_dataset <- nmr_dataset_load(nmr_dataset_rds)

  conditions_expr <- rlang::parse_exprs(conditions)
  
  nmr_dataset <- NIHSnmr::filter(nmr_dataset, !!!conditions_expr)
  
  message("Saving results...")
  nmr_export_data_1r(nmr_dataset, raw_data_matrix_fn)
  nmr_meta_export(nmr_dataset, metadata_fn, groups = "external")
  nmr_dataset_save(nmr_dataset, nmr_dataset_outfile)
  pasted_conditions <- glue::glue_collapse(conditions, sep = ", ", width = 80, last = ", ")
  message("Dataset filtered by: ", pasted_conditions)
  
  message("Ending pipe_filter_samples at ", Sys.time())
}

#' Pipeline: Peak detection and Alignment
#' @inheritParams pipe_add_metadata
#' @inheritParams nmr_detect_peaks
#' @inheritParams nmr_align
#'
#' @export
#'
pipe_peakdet_align <- function(nmr_dataset_rds,
                         nDivRange_ppm = 0.1,
                         scales = seq(1, 16, 2),
                         baselineThresh = 0.01,
                         SNR.Th = -1,
                         maxShift_ppm = 0.0015, acceptLostPeak = FALSE,
                         output_dir = NULL) {
  message("Starting pipe_peakdet_align at ", Sys.time())
  
  if (is.null(output_dir)) {
    stop("An output directory must be specified")
  }
  
  fs::dir_create(output_dir)
  
  metadata_fn <- file.path(output_dir, "metadata.xlsx")
  raw_data_matrix_fn <- file.path(output_dir, "raw_data.csv")
  nmr_dataset_outfile <- file.path(output_dir, "nmr_dataset.rds")
  plot_peak_detection_html <- file.path(output_dir, "peak-detection-diagnostic.html")
  
  plot_html <- file.path(output_dir, "plot-samples.html")
  peak_data_fn <- file.path(output_dir, "peak_data.csv")
  NMRExp_ref_fn <- file.path(output_dir, "NMRExperiment_align_ref.txt")
  
  nmr_dataset <- nmr_dataset_load(nmr_dataset_rds)
  message("Detecting peaks...")
  peak_data <- nmr_detect_peaks(nmr_dataset,
                                nDivRange_ppm = nDivRange_ppm,
                                scales = scales,
                                baselineThresh = baselineThresh,
                                SNR.Th = SNR.Th)
  message("Choosing alignment reference...")
  NMRExp_ref <- nmr_align_find_ref(nmr_dataset, peak_data)
  message("Starting alignment...")
  nmr_dataset <- nmr_align(nmr_dataset, peak_data,
                           NMRExp_ref = NMRExp_ref,
                           maxShift_ppm = maxShift_ppm,
                           acceptLostPeak = acceptLostPeak)
  
  message("Saving alignment results...")

  
  gplt <- nmr_detect_peaks_plot(nmr_dataset, peak_data, NMRExperiment = NMRExp_ref)
  plotly::toWebGL(gplt) %>%
    htmltools::as.tags(standalone = TRUE) %>%
    htmltools::save_html(html = plot_peak_detection_html)
  
  nmr_export_data_1r(nmr_dataset, raw_data_matrix_fn)
  nmr_meta_export(nmr_dataset, metadata_fn, groups = "external")
  nmr_dataset_save(nmr_dataset, nmr_dataset_outfile)
  plot_webgl(nmr_dataset, html_filename = plot_html)
  utils::write.csv(peak_data, peak_data_fn, row.names = FALSE)
  write(NMRExp_ref, NMRExp_ref_fn)
  
  message("Ending pipe_peakdet_align at ", Sys.time())
}


#' Pipeline: Peak integration
#'
#' @inheritParams pipe_add_metadata
#' @param peak_det_align_dir Output directory from [pipe_peakdet_align]
#' @param peak_width_ppm A peak width in ppm
#'
#' @importFrom rlang .data
#' @export
#'
pipe_peak_integration <- function(nmr_dataset_rds, peak_det_align_dir, peak_width_ppm, output_dir) {
  message("Starting pipe_peak_integration at ", Sys.time())
  
  if (is.null(output_dir)) {
    stop("An output directory must be specified")
  }
  
  fs::dir_create(output_dir)
  
  # Input files from previous node:
  peak_data_fn <- file.path(peak_det_align_dir, "peak_data.csv")
  NMRExp_ref_fn <- file.path(peak_det_align_dir, "NMRExperiment_align_ref.txt")
  
  # Output files:
  metadata_fn <- file.path(output_dir, "metadata.xlsx")
  peak_table_fn <- file.path(output_dir, "peak_table.csv")
  
  nmr_dataset <- nmr_dataset_load(nmr_dataset_rds)
  peak_data <- utils::read.csv(file = peak_data_fn)
  NMRExperimentRef <- readLines(NMRExp_ref_fn)
  peak_data_integ <- dplyr::filter(peak_data, .data$NMRExperiment == !!NMRExperimentRef)
  peak_table <- nmr_integrate_peak_positions(samples = nmr_dataset,
                                             peak_pos_ppm = peak_data_integ$ppm,
                                             peak_width_ppm = peak_width_ppm)
  
  nmr_meta_export(nmr_dataset, metadata_fn, groups = "external")
  utils::write.csv(peak_table, peak_table_fn, row.names = FALSE)
  message("Ending pipe_peak_integration at ", Sys.time())
}


#' Pipe: Full spectra normalization
#' 
#' Normalize the full spectra to the internal calibrant region, then exclude
#' that region and finally perform PQN normalization.
#' 
#' If there is no internal calibrant, only the PQN normalization is done.
#'
#' @inheritParams pipe_add_metadata
#' @param internal_calibrant A ppm range where the internal calibrant is, or `NULL`.
#'
#' @export
pipe_normalization <- function(nmr_dataset_rds, internal_calibrant = NULL, output_dir = NULL) {
  message("Starting pipe_normalization at ", Sys.time())
  
  if (is.null(output_dir)) {
    stop("An output directory must be specified")
  }
  
  fs::dir_create(output_dir)
  
  metadata_fn <- file.path(output_dir, "metadata.xlsx")
  full_spectra_matrix_fn <- file.path(output_dir, "full_spectra_matrix.csv")
  nmr_dataset_outfile <- file.path(output_dir, "nmr_dataset.rds")
  plot_html <- file.path(output_dir, "plot-samples.html")
  
  nmr_dataset <- nmr_dataset_load(nmr_dataset_rds)
  if (!is.null(internal_calibrant)) {
    nmr_dataset <- nmr_dataset %>%
      nmr_normalize(method = "region", values = internal_calibrant) %>%
      nmr_exclude_region(exclude = list(ic = internal_calibrant))
  }
  nmr_dataset <- nmr_normalize(nmr_dataset, method = "pqn")
  
  nmr_export_data_1r(nmr_dataset, full_spectra_matrix_fn)
  nmr_meta_export(nmr_dataset, metadata_fn, groups = "external")
  nmr_dataset_save(nmr_dataset, nmr_dataset_outfile)
  plot_webgl(nmr_dataset, html_filename = plot_html)
  message("Ending pipe_normalization at ", Sys.time())
}
