#' Get metadata
#'
#' @param samples a [nmr_dataset_family] object
#' @param columns Columns to get. By default gets all the columns.
#' @param groups Groups to get. Groups are predefined of columns. Typically
#' `"external"` for metadata added with [nmr_meta_add].
#'
#' Both `groups` and `columns` can't be given simultaneously.
#'
#' @return a data frame with the injection metadata
#' @family metadata functions
#' @family nmr_dataset functions
#' @family nmr_dataset_1D functions
#' @family nmr_dataset_peak_table functions
#' @export
#' @examples
#' dir_to_demo_dataset <- system.file("dataset-demo", package = "AlpsNMR")
#' dataset <- nmr_read_samples_dir(dir_to_demo_dataset)
#' metadata <- nmr_meta_get(dataset)
#' 
nmr_meta_get <- function(samples,
                         columns = NULL,
                         groups = NULL) {
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
        warning(
            "Missing columns: ",
            paste(utils::head(cols_miss, n = show_cols), collapse = ", "),
            " and ",
            length(cols_miss) - show_cols,
            " columns more."
        )
        rm(cols_miss, show_cols)
    }
    
    columns <- columns[columns %in% colnames(metadata)]
    # drop = FALSE ensures we never return a vector (always a data frame/tibble)
    metadata <- metadata[, columns, drop = FALSE]
    return(metadata)
}

#' Get a single metadata column
#'
#' @param samples a [nmr_dataset_family] object
#' @param column A column to get
#' @return A vector with the column
#' @family metadata functions
#' @family nmr_dataset functions
#' @family nmr_dataset_1D functions
#' @family nmr_dataset_peak_table functions
#' @export
#' @examples 
#' dir_to_demo_dataset <- system.file("dataset-demo", package = "AlpsNMR")
#' dataset <- nmr_read_samples_dir(dir_to_demo_dataset)
#' metadata_column <- nmr_meta_get_column(dataset)
#' 
nmr_meta_get_column <- function(samples, column = "NMRExperiment") {
    nmr_meta_get(samples, columns = column)[[column]]
}

#' Add metadata to an nmr_dataset object
#'
#' This is useful to add metadata to datasets that can be later used for
#' plotting spectra or further analysis (PCA...).
#'
#' @param nmr_data an [nmr_dataset_family] object
#' @param metadata A data frame with metadata to add
#' @param by A column name of both the `nmr_data$metadata$external` and the metadata
#' data.frame. If you want to merge two columns with different headers you can
#' use a named character vector `c("NMRExperiment" = "ExperimentNMR")` where
#' the left side is the column name of the `nmr_data$metadata$external` and the right side is
#' the column name of the metadata data frame.
#'
#' @return
#' The nmr_dataset_family object with the added metadata
#' @export
#' @examples
#' # Load a demo dataset with four samples:
#' dataset <- system.file("dataset-demo", package = "AlpsNMR")
#' nmr_dataset <- nmr_read_samples_dir(dataset)
#'
#' # At first we just have the NMRExperiment column
#' nmr_meta_get(nmr_dataset, groups = "external")
#' # Get a table with NMRExperiment -> SubjectID
#' dummy_metadata <- system.file("dataset-demo", "dummy_metadata.xlsx", package = "AlpsNMR")
#' NMRExp_SubjID <- readxl::read_excel(dummy_metadata, sheet = 1)
#'
#' NMRExp_SubjID
#' # We can link the SubjectID column of the first excel into the dataset
#' nmr_dataset <- nmr_meta_add(nmr_dataset, NMRExp_SubjID, by = "NMRExperiment")
#' nmr_meta_get(nmr_dataset, groups = "external")
#' # The second excel can use the SubjectID:
#' SubjID_Age <- readxl::read_excel(dummy_metadata, sheet = 2)
#' SubjID_Age
#' # Add the metadata by its SubjectID:
#' nmr_dataset <- nmr_meta_add(nmr_dataset, SubjID_Age, by = "SubjectID")
#' # The final loaded metadata:
#' nmr_meta_get(nmr_dataset, groups = "external")
#'
#' @family metadata functions
#' @family nmr_dataset functions
#' @family nmr_dataset_1D functions
#' @family nmr_dataset_peak_table functions
nmr_meta_add <- function(nmr_data, metadata, by = "NMRExperiment") {
    nmr_meta <- nmr_meta_get(nmr_data, groups = "external")
    by_left <- ifelse(is.null(names(by)), by, names(by))
    existing_vars <- base::setdiff(colnames(nmr_meta), by_left)
    conflict <- base::intersect(existing_vars, colnames(metadata))
    # We must ensure metadata[[by]] is unique:
    metadata <- do.call(dplyr::distinct, c(
        list(.data = metadata),
        rlang::syms(by),
        list(.keep_all = TRUE)
    ))
    nmr_meta_new <-
        dplyr::left_join(nmr_meta,
                         metadata,
                         by = by,
                         suffix = c("", "__REMOVE__"))
    are_identical <- purrr::map_lgl(conflict, function(col) {
        col1 <- col
        col2 <- paste0(col, "__REMOVE__")
        identical(nmr_meta_new[[col1]], nmr_meta_new[[col2]])
    })
    if (!all(are_identical)) {
        stop(
            "Can't add metadata because of column conflict at: ",
            paste(conflict[!are_identical], sep = ", ", collapse = ", ")
        )
    }
    nmr_meta_new <-
        dplyr::select(nmr_meta_new, -dplyr::ends_with("__REMOVE__"))
    nmr_data$metadata$external <- nmr_meta_new
    nmr_data
}

#' @rdname nmr_meta_add
#'
#' @param excel_file Path to a tidy Excel file name. The Excel can consist
#' of multiple sheets, that are added sequentially. The first column of the first
#' sheet MUST be named as one of the metadata already present in the dataset,
#' typically will be "NMRExperiment". The rest of the columns of the first sheet
#' can be named at will. Similary, the first column of the second sheet must be
#' named as one of the metadata already present in the dataset, typically
#' "NMRExperiment" or any of the columns of the first sheet. The rest of the columns
#' of the second sheet can be named at will. See the package vignette for an
#' example.
#'
#' @export
#' @examples
#' # Read a tidy excel file:
#'
#' dataset <- system.file("dataset-demo", package = "AlpsNMR")
#' nmr_dataset <- nmr_read_samples_dir(dataset)
#'
#' # At first we just have the NMRExperiment column
#' nmr_meta_get(nmr_dataset, groups = "external")
#' # Get a table with NMRExperiment -> SubjectID
#' dummy_metadata <- system.file("dataset-demo", "dummy_metadata.xlsx", package = "AlpsNMR")
#'
#' nmr_dataset <-nmr_meta_add_tidy_excel(nmr_dataset, dummy_metadata)
#' # Updated Metadata:
#' nmr_meta_get(nmr_dataset, groups = "external")
nmr_meta_add_tidy_excel <- function(nmr_data, excel_file) {
    excel_sheets <- readxl::excel_sheets(excel_file)
    excel_dfs <-
        purrr::map(excel_sheets,
                   ~ readxl::read_excel(path = excel_file, sheet = .))
    
    for (excel_df in excel_dfs) {
        # If the sheet is empty, skip
        if (nrow(excel_df) == 0 || ncol(excel_df) == 0) {
            # Empty sheet
            next
        }
        # Check that the key column is already in the metadata:
        key_column <- colnames(excel_df)[1]
        avail_columns <-
            colnames(nmr_meta_get(nmr_data, groups = "external"))
        if (!key_column %in% avail_columns) {
            stop(
                glue::glue(
                    "Key Column {key_column} is not present in the metadata.",
                    "Can't add metadata.\n Available columns:",
                    glue::glue_collapse(avail_columns, sep = ", ", last = " and ")
                )
            )
        }
        # The NMRExperiment column should of type character, although Excel may use numeric
        # because usually all the identifiers are numeric. We convert it manually:
        if ("NMRExperiment" %in% colnames(excel_df)) {
            excel_df$NMRExperiment <- as.character(excel_df$NMRExperiment)
        }
        # Add the metadata
        nmr_data <-
            nmr_meta_add(nmr_data, metadata = excel_df, by = key_column)
    }
    nmr_data
}

#' Export Metadata to an Excel file
#'
#' @param nmr_dataset An [nmr_dataset_family] object
#' @param xlsx_file "The .xlsx excel file"
#' @param groups A character vector. Use `"external"` for the external metadata or
#'    the default for a more generic solution
#' @return The Excel file name
#' @export
#' @family metadata functions
#' @family nmr_dataset functions
#' @family nmr_dataset_1D functions
#' @family nmr_dataset_peak_table functions
#' @family import/export functions
#' @examples 
#' dir_to_demo_dataset <- system.file("dataset-demo", package = "AlpsNMR")
#' dataset <- nmr_read_samples_dir(dir_to_demo_dataset)
#' #nmr_meta_export(dataset, "metadata.xlsx")
#' 
nmr_meta_export <- function(nmr_dataset,
                            xlsx_file,
                            groups = c("info", "orig", "title", "external")) {
    require_pkgs("writexl")
    groups_present <- groups %in% names(nmr_dataset$metadata)
    if (!all(groups_present)) {
        warning(
            "These metadata groups are missing and will be ignored: \n",
            paste(groups[!groups_present]),
            collapse = ", "
        )
        groups <- groups[groups_present]
    }
    writexl::write_xlsx(x = nmr_dataset$metadata[groups], path = xlsx_file)
}
