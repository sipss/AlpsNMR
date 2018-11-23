#' Get metadata
#' 
#' @family metadata functions
#' @param samples a [nmr_dataset] object
#' @param columns Columns to get. By default gets all the columns.
#' @param groups Groups to get. Groups are predefined of columns. Typically 
#' `"external"` for metadata added with [nmr_meta_add].
#' 
#' Both `groups` and `columns` can't be given simultaneously.
#' 
#' @return a data frame with the injection metadata
#' @family metadata management functions
#' @family nmr_dataset functions
#' @family nmr_dataset_1D functions
#' @family nmr_dataset_peak_table functions
#' @export
nmr_meta_get <- function(samples, columns = NULL, groups = NULL) {
  metadata_list <- samples[["metadata"]]
  metadata <- metadata_list[[1]]
  for (i in utils::tail(seq_along(metadata_list), -1)) {
    metadata <- dplyr::left_join(metadata,
                                 metadata_list[[i]],
                                 by = "NMRExperiment")
  }

  # Default columns means all columns
  if (!is.null(columns) && !is.null(groups)) {
    stop("groups and columns can't be given simultaneously")
  }
  if (is.null(columns) && !is.null(groups)) {
    columns <- metadata_list[groups] %>%
      purrr::map(colnames) %>%
      purrr::flatten_chr() %>%
      unique()
  } else if (is.null(columns)) {
    columns <- colnames(metadata)
  }
  
  # NMRExperiment is always present in the output
  if (!"NMRExperiment" %in% columns) {
    columns <- c("NMRExperiment", columns)
  }
  
  # Report if user wants columns that are not present:
  if (!all(columns %in% colnames(metadata))) {
    cols_miss <- columns[!columns %in% colnames(metadata)]
    # If there are less than 10 missing columns warn all the missing names,
    # otherwise warn the first 5
    if (length(cols_miss) < 10) {
      show_cols <- length(cols_miss)
    } else {
      show_cols <- 5
    }
    warning("Missing columns: ",
            paste(utils::head(cols_miss, n = show_cols), collapse = ", "),
            " and ", length(cols_miss) - show_cols, " columns more.")
    rm(cols_miss, show_cols)
  }
  
  columns <- columns[columns %in% colnames(metadata)]
  # drop = FALSE ensures we never return a vector (always a data frame/tibble)
  metadata <- metadata[, columns, drop = FALSE]
  return(metadata)
}

#' Get a single metadata column
#' 
#' @family metadata functions
#' @param samples a [nmr_dataset] object
#' @param column A column to get
#' @return A vector with the column
#' @family metadata management functions
#' @family nmr_dataset functions
#' @family nmr_dataset_1D functions
#' @family nmr_dataset_peak_table functions
#' @export
nmr_meta_get_column <- function(samples, column = "NMRExperiment") {
  nmr_meta_get(samples, columns = column)[[column]]
}

#' Add metadata to an nmr_dataset object
#' 
#' This is useful to add metadata to datasets that can be later used for
#' plotting spectra or further analysis (PCA...).
#' 
#' @family metadata functions
#' @param nmr_data an [nmr_dataset] object
#' @param metadata A data frame with metadata to add
#' @param by A column name of both the `nmr_dataset$metadata$external` and the metadata
#' data.frame. If you want to merge two columns with different headers you can
#' use a named character vector `c("NMRExperiment" = "ExperimentNMR")` where
#' the left side is the column name of the `nmr_dataset$metadata$external` and the right side is
#' the column name of the metadata data frame.
#' 
#' @return
#' The nmr_dataset object with the added metadata
#' @export
#' @examples 
#' # Load a demo dataset with four samples:
#' dataset <- system.file("dataset-demo", package = "NIHSnmr")
#' nmr_dataset <- nmr_read_samples_dir(dataset)
#' 
#' # At first we just have the NMRExperiment column
#' print(nmr_meta_get(nmr_dataset, groups = "external"))
#' # Get a table with NMRExperiment -> SubjectID
#' dummy_metadata <- system.file("dataset-demo", "dummy_metadata.xlsx", package = "NIHSnmr")
#' NMRExp_SubjID <- readxl::read_excel(dummy_metadata, sheet = 1)
#' 
#' print(NMRExp_SubjID)
#' # We can link the SubjectID column of the first excel into the dataset
#' nmr_dataset <- nmr_meta_add(nmr_dataset, NMRExp_SubjID, by = "NMRExperiment")
#' print(nmr_meta_get(nmr_dataset, groups = "external"))
#' # The second excel can use the SubjectID:
#' SubjID_Age <- readxl::read_excel(dummy_metadata, sheet = 2)
#' print(SubjID_Age)
#' # Add the metadata by its SubjectID:
#' nmr_dataset <- nmr_meta_add(nmr_dataset, SubjID_Age, by = "SubjectID")
#' # The final loaded metadata:
#' print(nmr_meta_get(nmr_dataset, groups = "external"))
#' 
#' @family metadata management functions
#' @family nmr_dataset functions
#' @family nmr_dataset_1D functions
#' @family nmr_dataset_peak_table functions
#' @family import/export functions
nmr_meta_add <- function(nmr_data, metadata, by = "NMRExperiment") {
  nmr_meta <- nmr_meta_get(nmr_data, groups = "external")
  by_left <- ifelse(is.null(names(by)), by, names(by))
  existing_vars <- base::setdiff(colnames(nmr_meta), by_left)
  conflict <- base::intersect(existing_vars, colnames(metadata))
  # We must ensure metadata[[by]] is unique:
  metadata <- dplyr::distinct(metadata, !!!rlang::syms(by), .keep_all = TRUE)
  nmr_meta_new <- dplyr::left_join(nmr_meta, metadata, by = by, suffix = c("", "__REMOVE__"))
  are_identical <- purrr::map_lgl(conflict, function(col) {
    col1 <- col
    col2 <- paste0(col, "__REMOVE__")
    identical(nmr_meta_new[[col1]], nmr_meta_new[[col2]])
  })
  if (!all(are_identical)) {
    stop("Can't add metadata because of column conflict at: ", paste(conflict[!are_identical], sep = ", ", collapse = ", "))
  }
  nmr_meta_new <- dplyr::select(nmr_meta_new, -dplyr::ends_with("__REMOVE__"))
  nmr_data$metadata$external <- nmr_meta_new
  nmr_data
}

#' Export Metadata to an Excel file
#'
#' @family metadata functions
#' @param nmr_dataset An nmr_dataset object
#' @param xlsx_file "The .xlsx excel file"
#' @param groups A character vector. Use `"external"` for the external metadata or
#'  the default for a more generic solution
#' @return The Excel file name
#' @export
#' @family metadata management functions
#' @family nmr_dataset functions
#' @family nmr_dataset_1D functions
#' @family nmr_dataset_peak_table functions
#' @family import/export functions
nmr_meta_export <- function(nmr_dataset, 
                            xlsx_file,
                            groups = c("info", "orig", "title", "external")) {
  groups_present <- groups %in% names(nmr_dataset$metadata)
  if (!all(groups_present)) {
    warning("These metadata groups are missing and will be ignored: \n", paste(groups[!groups_present]), collapse = ", ")
    groups <- groups[groups_present]
  }
  writexl::write_xlsx(x = nmr_dataset$metadata[groups], path = xlsx_file)
}
