# This functions follow NIHS internal conventions and therefore are not
# suitable for general use.

#' Reads an exported xlsx slims spreadsheet
#'
#' This function reads an xlsx slims spreadsheet and converts it into a data frame
#'
#' @param file_name The exported slims spreadsheet in excel format
#' @param sheet passed on to read_excel
#' @param all_character Read all columns as text
#' @return a data frame with the slims information
#'
#' Key columns on the output data frame are "cntn_cf_fluidx_barcode" and "cntn_cf_subject_id" "cntn_cf_fk_masterSample"
#'
#' @export
read_exported_slims <- function(file_name, sheet = NULL, all_character = FALSE) {
  # Row 1: Slims logo
  # Row 2: Internal column names (this column is sometimes hidden)
  # Row 3: Human readable column names
  # Then the rest of the excel sheet is the actual data
  # We read the excel file twice, on the first time we skip the logo (skip = 1)
  slims <- readxl::read_excel(file_name, skip = 1, sheet = sheet)
  # This allows us to have the internal column names that are more
  # reliable (they should never change, slims uses them internally)
  standard_slims_names <- colnames(slims)
  if (all_character) {
    # Remove the first row (human readable title)
    slims <- slims[-1, , drop = FALSE]
    return(slims)
  }
  # However, when we do this read_excel assumes that all the spreadsheet is of
  # type character, because row 3 is text.
  # Therefore we read again the spreadsheet skipping the first two rows (skip = 2)
  slims <- readxl::read_excel(file_name, skip = 2, sheet = sheet)
  # human_readable_colnames <- colnames(slims)
  # And finally we set the standard slims names as column names
  colnames(slims) <- standard_slims_names
  # Check the column names are not completely wrong
  return(slims)
}


#' Get metadata from an nmr_dataset for irods
#' @param nmr_data The \code{\link{nmr_dataset}} object
#' @return A data frame (tibble) with the columns used by NIHSrods
#' @export
nmr_get_irods_meta <- function(nmr_data) {
  meta_irods <- tibble::tibble(injection_id = nmr_data$metadata$injection_id)

  check_and_add <- function(meta_irods, irods_column, injection_column) {
    if (injection_column %in% colnames(nmr_data$metadata)) {
      meta_irods[[irods_column]] <- nmr_data$metadata[[injection_column]]
    } else {
      meta_irods[[irods_column]] <- NA
    }
    return(meta_irods)
  }

  meta_irods <- meta_irods %>%
    check_and_add(irods_column = "sample_path", injection_column = "info_sample_path") %>%
    dplyr::mutate(src_path = NA) %>%
    check_and_add(irods_column = "project_NPDI_ID", injection_column = "orig_NPDI_ID") %>%
    check_and_add(irods_column = "study_nickname", injection_column = "orig_StudyNickname") %>%
    check_and_add(irods_column = "assay_ID", injection_column = "orig_AssayID") %>%
    dplyr::mutate(file_type = "Raw data file") %>%
    check_and_add(irods_column = "file_format", injection_column = "info_file_format") %>%
    check_and_add(irods_column = "master_sample_accession", injection_column = "orig_MasterRBarcode") %>%
    check_and_add(irods_column = "operational_sample_accession", injection_column = "orig_RBarcode") %>%
    check_and_add(irods_column = "run_ID", injection_column = "orig_RunID") %>%
    dplyr::mutate(assay_platform = "Metabolomics") %>%
    check_and_add(irods_column = "software_platform", injection_column = "acqus_TITLE") %>%
    check_and_add(irods_column = "taxonomy", injection_column = "orig_Taxonomy") %>%
    check_and_add(irods_column = "hardware_platform", injection_column = "orig_Spectrometer") %>%
    check_and_add(irods_column = "taxonomy", injection_column = "orig_Taxonomy")


  # stringr::str_c is like paste0 BUT:
  # paste0("Bruker ", c("RoomTemp", NA)) is c("Bruker RoomTemp", "Bruker NA")
  # and str_c does "the right thing"
  # stringr::str_c("Bruker ", c("RoomTemp", NA)) is c("Bruker RoomTemp", NA)
  meta_irods$assay_technique <- stringr::str_c(
    "Nuclear Magnetic Resonance-",
    nmr_data$metadata$info_dimension,
    "D", "-", nmr_data$metadata$info_nuclei,
    "-", nmr_data$metadata$info_pulse_sequence)

  # Remove "Parameter file, "
  meta_irods$software_platform <- gsub(pattern = "Parameter file, ",
                                       replacement = "",
                                       meta_irods$software_platform)

  # Change Topspin from all caps to something standard:
  meta_irods$software_platform <- gsub("TOPSPIN", "TopSpin",
                                       meta_irods$software_platform,
                                       ignore.case = TRUE)

  # Remove Version (we know TopSpin Version3.2 means TopSpin 3.2)
  meta_irods$software_platform <- gsub("Version", " ",
                                       meta_irods$software_platform,
                                       ignore.case = TRUE)

  # Remove spaces (as was defined on the ontology)
  meta_irods$software_platform <-  gsub(pattern = "[\t ]",
                                        replacement = "",
                                        meta_irods$software_platform)

  meta_irods$hardware_platform <- stringr::str_c("Bruker 600MHz NMR Spectrometer ",
                                                 meta_irods$hardware_platform)
  meta_irods$author <- NA # Not mandatory for raw data
  meta_irods$relates_to <- NA # Unknown
  return(meta_irods)
}

#' Create one zip file for each sample
#'
#' @param meta_irods Data frame given by \code{\link{nmr_get_irods_meta}}
#' @param workdir Directory to store zip files
#' @param overwrite Should existing zip files be overwritten?
#' @param ... Passed to \code{\link[utils]{zip}}
#' @return A meta_irods data frame with an additional \code{src_path} column
#'         containing the zip file path location.
#' @export
nmr_prepare_zip_files <- function(meta_irods, workdir, overwrite = FALSE, ...) {
  current_wd <- getwd()
  dir.create(workdir, recursive = TRUE)
  workdir <- normalizePath(workdir)
  new_meta_irods <- meta_irods
  new_meta_irods$src_path <- NA
  # tryCatch to make sure the working directory is restored
  # The zip utility in R is very basic. On Windows it relies on the zip.exe
  # program being installed and in the path. zip.exe is found on RTools.
  tryCatch({
    for (row_idx in seq_len(nrow(meta_irods))) {
      # The directory where the NMR sample is: (e.g. "/nihs/Instrument/.../DUND-plasma/10")
      dir_to_compress <- normalizePath(meta_irods$sample_path[row_idx])
      # dir_name: (e.g. "10")
      dir_name <- basename(dir_to_compress)
      # parent_dir: "/nihs/Instrument/.../DUND-plasma"
      parent_dir <- normalizePath(dirname(dir_to_compress))
      # destination_file: /tmp/my_temp_dir/10.zip
      destination_file <- file.path(workdir, paste0(dir_name, ".zip"))
      new_meta_irods$src_path[row_idx] <- destination_file
      if (file.exists(destination_file)) {
        if (overwrite) {
          # deletes destination zip file if overwrite is true
          unlink(destination_file)
        } else {
          # Skips otherwise
          warning("File ", destination_file, " already exists. Skipping")
          next
        }
      }
      # Change to the parent_dir so the zip file does not include the whole
      # directory tree folders
      setwd(parent_dir)
      # Zip it!
      utils::zip(zipfile = destination_file, files = dir_name, ...)
      # Go back to the initial directory, so the next file can be
      # found if relative paths are used
      setwd(current_wd)
    }
  }, finally = {
    # Always restore the working directory, even if there are errors somewhere
    setwd(current_wd)
  })
  return(new_meta_irods)
}

nmr_push_one_to_irods <- function(element, dest_path, extra_columns = NULL) {
  if (!requireNamespace("NIHSrods", quietly = TRUE)) {
    stop("NIHSrods needed for this function to work. Please install it.",
         call. = FALSE)
  }
  full_irods_path <-
    NIHSrods::iput(src_path = element$src_path,
                 dest_path = dest_path,
                 project_NPDI_ID = element$project_NPDI_ID,
                 study_nickname = element$study_nickname,
                 assay_ID = element$assay_ID,
                 file_type = element$file_type,
                 file_format = element$file_format,
                 master_sample_accession = element$master_sample_accession,
                 operational_sample_accession = element$operational_sample_accession,
                 run_ID = element$run_ID,
                 assay_platform = element$assay_platform,
                 software_platform = element$software_platform,
                 taxonomy = element$taxonomy,
                 hardware_platform = element$hardware_platform,
                 assay_technique = element$assay_technique,
                 author = element$author,
                 relates_to = element$relates_to)

  if (is.null(full_irods_path)) {
    stop("Could not upload ", element$src_path, " to ", dest_path)
  }
  if (!is.null(extra_columns) && is.null(names(extra_columns))) {
    names(extra_columns) <- extra_columns
  }
  for (extra_column_name in names(extra_columns)) {
    extra_column <- unname(extra_columns[extra_column_name])
    if (is.null(element[[extra_column]]) | is.na(element[[extra_column]])) {
      next
    }

    NIHSrods::imeta_set("d", full_irods_path,
                        attribute = extra_column_name,
                        value = element[[extra_column]])
  }
  return(full_irods_path)
}

#' Push NMR data to irods
#' @param meta_irods a dataframe as obtained from \code{\link{nmr_get_irods_meta}} with
#'                   these columns: \code{injection_id}, \code{sample_path}, \code{"src_path",
#'                   "project_NPDI_ID", "study_nickname", "assay_ID", "file_type",
#'                   "file_format", "master_sample_accession", "operational_sample_accession",
#'                   "run_ID", "assay_platform", "software_platform", "taxonomy",
#'                   "hardware_platform", "assay_technique", "author", "relates_to"}
#' @param dest_path The irods destination directory
#' @param check_metadata_issues if \code{FALSE} missing values in mandatory
#'                              \code{meta_irods} columns are allowed.
#' @param extra_columns Column names from \code{meta_irods} that will
#'                      be stored into irods as additional metadata values.
#'                      If \code{extra_columns} has names, then those names will be
#'                      used as attribute names.
#' @return The meta_irods data frame with an additional column named \code{full_irods_path}.
#' @export
#' @examples
#' \dontrun{
#'  nmrdata <- nmr_read_samples_dir("/dir/DUND-100713-Constitutive Thinness-plasma/",
#'                                  metadata_only = TRUE)
#'  irods_metadata <- nmr_get_irods_meta(nmrdata)
#'  # explore and confirm everything is in there
#'  View(irods_metadata)
#'  # optionally add extra columns (but try to adhere to a standard if possible):
#'  irods_metadata$looks_weird <- NA
#'  irods_metadata$looks_weird[1] <- TRUE
#'  # Now, let's compress the samples to zip files for irods
#'  irods_metadata <- nmr_prepare_zip_files(meta_irods = irods_metadata,
#'                                          workdir = 'ct_plasma_zip_files')
#'  # Finally push the data to the irods directory
#'  nmr_push_to_irods(irods_metadata, "/NIHSData/DUND-100713/ct/Data/Data_raw/Metabolomics",
#'                    extra_columns = "looks_weird")
#' }
nmr_push_to_irods <- function(meta_irods, dest_path, check_metadata_issues = TRUE,
                              extra_columns = NULL) {
  if (!requireNamespace("NIHSrods", quietly = TRUE)) {
    stop("NIHSrods needed for this function to work. Please install it.",
         call. = FALSE)
  }
  if (is.nmr_dataset(meta_irods)) {
    stop("To push to irods you need to use first nmr_get_irods_meta and nmr_prepare_zip_files")
  }

  mandatory_cols <- c("injection_id", "src_path", "project_NPDI_ID",
                      "study_nickname", "assay_ID", "file_type", "file_format", "run_ID",
                      "assay_platform", "software_platform", "taxonomy", "hardware_platform",
                      "assay_technique", "operational_sample_accession",
                      "master_sample_accession")

  have_mandatory_cols <- mandatory_cols %in% colnames(meta_irods)
  if (!all(have_mandatory_cols)) {
    if (check_metadata_issues) {
      stop("meta_irods has missing columns: ",
           paste(colnames(meta_irods)[have_mandatory_cols]), collapse = ", ")
    } else {
      warning("meta_irods has missing columns: ",
           paste(colnames(meta_irods)[have_mandatory_cols]), collapse = ", ")
    }
  }

  irods_metadata_plasma_problem <- meta_irods[rowSums(is.na(meta_irods[,mandatory_cols])) > 0,]

  if (nrow(irods_metadata_plasma_problem) > 0) {
    if (check_metadata_issues) {
      stop("The following rows in meta_irods have missing values:\n",
           paste(utils::capture.output(print(as.data.frame(irods_metadata_plasma_problem))), collapse = "\n"))
    } else {
      warning("The following rows in meta_irods have missing values:\n",
              paste(utils::capture.output(print(as.data.frame(irods_metadata_plasma_problem))), collapse = "\n"))
    }
  }

  NIHSrods::imkdir(dest_path, parents = TRUE)
  # First choice: it is a data frame
  if (is.data.frame(meta_irods)) {
    res <- vector('list', nrow(meta_irods))
    tryCatch({
      pb = NULL
      if (show_progress_bar(1, 0)) {
        pb <- utils::txtProgressBar(min = 1, max = nrow(meta_irods), style = 3)
      }
      for (i in seq_len(nrow(meta_irods))) {
        res[[i]] <- nmr_push_one_to_irods(meta_irods[i,], dest_path, extra_columns = extra_columns)
        utils::setTxtProgressBar(pb, value = i)
      }
      return(unlist(res, recursive = FALSE))
    }, finally = {
      close(pb)
    })
  } else if (is.list(meta_irods)) {
    res <- nmr_push_one_to_irods(meta_irods, dest_path)
  } else {
    stop("meta_irods should be a list or a data frame")
  }
  meta_irods$full_irods_path <- res
  return(meta_irods)
}


#' Search NMR samples in irods
#'
#' Search for NMR samples inside irods.
#'
#' This function is slow because NIHSrods::imeta_ls does not accept multiple files
#' at once.
#'
#' @param project_NPDI_ID The NPDI project ID (e.g. 'DUND-700713')
#' @param study_nickname The study nickname, typically defined by a Project Manager
#'                       (e.g. 'st_etienne'). A study nickname belongs to a
#'                       project ID.
#' @param assay_ID An identifier which must be unique within a study and which corresponds
#'                 to a piece of work in the project management sense. It may be defined by
#'                 the project manager or the platform, but should be known to the project manager.
#'                 For an analytical platform it corresponds to a request to perform analysis on a
#'                 set of samples from a Project. e.g. "Plasma NMR metabolites"
#'                 or "Urine NMR metabolites".
#' @param file_type By default: \code{"Raw data file"}.
#'                  The full list of accepted file types is available here:
#'                  \code{NIHSrods::getValidTerms('File type')}
#' @param file_format By default \code{"Zipped Bruker NMR directory"}. The full list of
#'                    accepted file formats is available here:
#'                    \code{NIHSrods::getValidTerms('File format')}
#' @param master_sample_accession The master sample accession. It is a code like \code{'R0000123456'}.
#' @param operational_sample_accession The operational sample accession. It is a code like \code{'R0000123456'}.
#' @param run_ID Optional field for operational use to separate distinct batches
#'               of samples within an Assay ID.
#' @param assay_platform By default: \code{'Metabolomics'}.
#'                       The full list of defined assay platforms is available here:
#'                       \code{NIHSrods::getValidTerms('Assay platform')}.
#' @param software_platform The software platform used to acquire or process the data.
#'                          One of \code{NIHSrods::getValidTerms('Software platform')}
#' @param taxonomy The biological taxonomy of the sample. For instance \code{'Homo sapiens (Human)'}.
#'                 The list of defined biological taxonomies at NIHS is available here:
#'                 \code{NIHSrods::getValidTerms('Taxonomy')}.
#' @param hardware_platform Machine which created the raw data. One of
#'                          \code{NIHSrods::getValidTerms('Hardware platform')}
#' @param assay_technique The analytical technique used for creating the sample.
#'                        One of \code{NIHSrods::getValidTerms('Assay technique')}
#' @param author This will be the person who performed the processing or analysis and made the decisions.
#'               Required for processed data, the author who did the processing.
#' @param exclude_trash Exclude samples in \code{/NIHSData/trash} (Default: \code{TRUE})
#' @return a data frame with the irods metadata that results from the search
#' @export
nmr_irods_search <- function(project_NPDI_ID = NA, study_nickname = NA,
                             assay_ID = NA, master_sample_accession = NA,
                             operational_sample_accession = NA,
                             run_ID = NA, software_platform = NA,
                             taxonomy = NA, hardware_platform = NA,
                             assay_technique = NA, file_type = "Raw data file",
                             file_format = "Zipped Bruker NMR directory",
                             assay_platform = "Metabolomics", author = NA,
                             exclude_trash = TRUE) {
  if (!requireNamespace("NIHSrods", quietly = TRUE)) {
    stop("NIHSrods needed for this function to work. Please install it.",
         call. = FALSE)
  }
  samples <- NIHSrods::isearch(
    project_NPDI_ID = project_NPDI_ID, study_nickname = study_nickname,
    assay_ID = assay_ID, file_type = file_type, file_format = file_format,
    master_sample_accession = master_sample_accession,
    operational_sample_accession = operational_sample_accession,
    run_ID = run_ID, assay_platform = assay_platform,
    software_platform = software_platform, taxonomy = taxonomy,
    hardware_platform = hardware_platform, assay_technique = assay_technique,
    author = author)

  samples$full_irods_path <- paste(samples$Collection_name, samples$Data_name, sep = "/")
  # Exclude the trash:
  if (exclude_trash) {
    samples <- samples[!grepl("^/NIHSData/trash", samples$full_irods_path),]
  }
  if (nrow(samples) == 0) {
    return(data.frame())
  }

  # Progress bar needs to be closed on error, so we use tryCatch
  metadata <- NULL
  tryCatch({
    pb <- NULL
    if (show_progress_bar(nrow(samples), 3)) {
      pb <- utils::txtProgressBar(min = 1, max = nrow(samples), style = 3)
    }

    for (i in seq_len(nrow(samples))) {
      metadata <- rbind(metadata,
                        cbind(Data_path = samples$Data_path[i],
                              full_irods_path = samples$full_irods_path[i],
                              NIHSrods::imeta_ls(type = samples$Data_type[i],
                                                 name = samples$full_irods_path[i]),
                              stringsAsFactors = FALSE))
      if (!is.null(pb)) {
        utils::setTxtProgressBar(pb, value = i)
      }
    }
  }, finally = {if (!is.null(pb)) {close(pb)}})

  # Check if for the same attribute all units are equal
  units <- by(metadata[["Unit"]],
              metadata[["Attribute"]],
              function(x) length(unique(x)))

  if (all(units == 1)) {
    # Rename Attribute to Attribute_Unit
    names(metadata)[names(metadata) == 'Attribute'] <- 'Attribute_Unit'
    # Remove Unit column:
    metadata <- metadata[, !colnames(metadata) %in% "Unit", drop = FALSE]
  } else {
    metadata <- tidyr::unite_(metadata, "Attribute_Unit", from = c("Attribute", "Unit"))
  }
  remove_rows <- metadata$Attribute_Unit == "" & metadata$Value == ""
  metadata <- metadata[!remove_rows,]
  metadata_table <- tidyr::spread(metadata, "Attribute_Unit", "Value")

  # Rename the columns to R-friendly names
  # (matching the arguments of NIHSrods::iput):
  new_column_names <- c(
    "Data_path" = "data_path", "full_irods_path" = "full_irods_path",
    "Assay ID" = "assay_ID", "Assay platform" = "assay_platform",
    "Assay technique" = "assay_technique", "File format" = "file_format",
    "File type" = "file_type", "Hardware platform" = "hardware_platform",
    "Master sample accession" = "master_sample_accession",
    "Operational sample accession" = "operational_sample_accession",
    "Project NPDI ID" = "project_NPDI_ID", "Run ID" = "run_ID",
    "Software platform" = "software_platform",
    "Study nickname" = "study_nickname", "Taxonomy" = "taxonomy",
    "Author" = "author", "Relates to" = "relates_to")

  new_col_names_present <- new_column_names[names(new_column_names) %in% colnames(metadata_table)]
  idx <- match(names(new_col_names_present), colnames(metadata_table))
  colnames(metadata_table)[idx] <- new_col_names_present

  return(metadata_table)
}

#' Read samples from irods
#'
#' @param irods_search_results The data frame resulting of nmr_irods_search that at least
#'                             must have the \code{full_irods_path} column with the
#'                             irods paths that look like
#'                             \code{"/NIHSData/home/rdollermse/DUND-100713/plasma/10.zip"}
#' @param ... passed to \code{\link{nmr_read_samples}}
#' @return a \code{\link{nmr_dataset}} object
#' @examples
#' \dontrun{
#' irods_results <-
#'   nmr_irods_search(project_NPDI_ID = "DUND-100713",
#'                    study_nickname = "st_etienne",
#'                    assay_ID = "Plasma NMR metabolites",
#'                    assay_technique = "Nuclear Magnetic Resonance-1D-1H-NOESY")
#' plasmadata <- nmr_read_samples_irods(irods_results)
#' }
#' @export
nmr_read_samples_irods <- function(irods_search_results, ...) {
  if (!requireNamespace("NIHSrods", quietly = TRUE)) {
    stop("NIHSrods needed for this function to work. Please install it.",
         call. = FALSE)
  }
  full_irods_path <- irods_search_results$full_irods_path
  sample_names_local <- character(length = length(full_irods_path))
  message("Fetching samples from irods...")
  tryCatch({
    pb <- NULL
    if (show_progress_bar(length(full_irods_path))) {
      pb <- utils::txtProgressBar(min = 0, max = length(full_irods_path), style = 3)
    }
    for (i in seq_along(full_irods_path)) {
      sample_names_local[i] <-
        NIHSrods::iget(src_path = full_irods_path[i],
                       dest_path = tempfile(fileext = paste0("_", basename(full_irods_path[i]))))
      if (!is.null(pb)) {
        utils::setTxtProgressBar(pb, value = i)
      }
    }
  }, finally = {
    if (!is.null(pb)) {
      close(pb)
    }
  })
  message("Loading samples...")
  dataset <- nmr_read_samples(sample_names = sample_names_local,
                              overwrite_sample_names = full_irods_path, ...)
  message("Appending irods metadata to nmr_dataset")
  colnames(irods_search_results) <- paste0("irods_", colnames(irods_search_results))
  dataset$metadata <- dplyr::left_join(dataset$metadata,
                                       irods_search_results,
                                       by = c("info_sample_path" = "irods_full_irods_path"))
  return(dataset)
}


# Validate metadata:

#' Validate irods_meta
#' @param irods_meta A data frame with irods metadata
#' @return a logical. \code{TRUE} if checked columns are in the ontology, \code{FALSE} otherwise.
#' @export
nmr_irods_validate_meta <- function(irods_meta) {
  file_types <- unique(irods_meta$file_type)
  output <- TRUE
  if (!(all(NIHSrods::isValidTermName(file_types, "File type")))) {
    warning("File type invalid")
    output <- FALSE
  }

  file_formats <- unique(irods_meta$file_format)
  if (!(all(NIHSrods::isValidTermName(file_formats, "File format")))) {
    warning("File format invalid")
    output <- FALSE
  }

  assay_platforms <- unique(irods_meta$assay_platform)
  if (!(all(NIHSrods::isValidTermName(assay_platforms, "Assay platform")))) {
    warning("Assay platform invalid")
    output <- FALSE
  }

  soft_platforms <- unique(irods_meta$software_platform)
  if (!(all(NIHSrods::isValidTermName(soft_platforms, "Software platform")))) {
    warning("Software platform invalid")
    output <- FALSE
  }

  taxonomies <- unique(irods_meta$taxonomy)
  if (!(all(NIHSrods::isValidTermName(taxonomies, "Taxonomy")))) {
    warning("Taxonomy invalid")
    output <- FALSE
  }

  hw_plat <- unique(irods_meta$hardware_platform)
  if (!(all(NIHSrods::isValidTermName(hw_plat, "Hardware platform")))) {
    warning("Hardware platform invalid")
    output <- FALSE
  }

  assay_techniques <- unique(irods_meta$assay_technique)
  if (!(all(NIHSrods::isValidTermName(assay_techniques, "Assay technique")))) {
    warning("Assay technique invalid")
    output <- FALSE
  }
  return(output)
}
