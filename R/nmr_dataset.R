#' nmr_dataset (S3 class)
#'
#' The `nmr_dataset` represents a set of NMR samples.
#' It is defined as an S3 class, and it can be treated as a regular list.
#'
#' It currently has the following elements:
#'
#' - `metadata`: The metadata fields. A data frame with one row per sample
#' 
#' - `axis`: A list with length equal to the dimensionality of the data.
#' For 1D spectra it is a list with a numeric vector
#' 
#' - `data_*`: Data arrays with the actual spectra. The first index represents
#' the sample, the rest of the indices match the length of each `axis`.
#' Typically `data_1r` is a matrix with one sample on each row and the chemical 
#' shifts in the columns.
#' 
#' - `processing`: The processing steps performed on the object, such as
#' interpolation, exclusion or normalization.
#' 
#' @name nmr_dataset
NULL


#' Read NMR samples
#'
#' These functions load samples from files and return a [nmr_dataset].
#'
#' @name nmr_read_samples
#' @param sample_names A character vector with file or directory names.
#' @param samples_dir A directory that contains multiple samples
#' @param metadata_only A logical, to load only metadata (default: `FALSE`)
#' @param pulse_sequence If it is set to a pulse sequence
#'                       ("NOESY", "JRES", "CPMG"...) it will only load
#'                       the samples that match that pulse sequence.
#' @param overwrite_sample_names (Internal Use)
#' @param ... Arguments passed to [read_bruker_sample()] for data loading
#' @return a [nmr_dataset] object
NULL

#' @rdname nmr_read_samples
#' @export
nmr_read_samples_dir <- function(samples_dir, pulse_sequence = NULL,
                                 metadata_only = FALSE,
                                 overwrite_sample_names = NULL,
                                 ...) {
  if (!dir.exists(samples_dir)) {
    stop("Invalid directory: ", samples_dir)
  }
  all_samples <- c(list.dirs(path = samples_dir, full.names = TRUE, recursive = FALSE),
                   list.files(path = samples_dir, full.names = TRUE, pattern = ".*zip$"),
                   list.files(path = samples_dir, full.names = TRUE, pattern = ".*jdx$"))

  dataset <- nmr_read_samples(sample_names = all_samples,
                              pulse_sequence = pulse_sequence,
                              metadata_only = metadata_only,
                              overwrite_sample_names = overwrite_sample_names,
                              ...)
  return(dataset)
}


#' @rdname nmr_read_samples
#' @export
nmr_read_samples <- function(sample_names, pulse_sequence = NULL,
                             metadata_only = FALSE,
                             overwrite_sample_names = NULL, ...) {
  # If all samples are directories or zips, use the bruker directory format:
  if (all(dir.exists(sample_names) | grepl('\\.zip$', sample_names))) {
    samples <- nmr_read_samples_bruker(sample_names = sample_names,
                                       metadata_only = metadata_only,
                                       pulse_sequence = pulse_sequence,
                                       overwrite_sample_names = overwrite_sample_names,
                                       ...)
  } else if (all(file.exists(sample_names))) {
    # otherwise the jdx format
    samples <- nmr_read_samples_jdx(sample_names = sample_names,
                                    metadata_only = metadata_only)
  } else {
    stop("Are all samples in the same format?")
  }
  return(samples)
}

# @rdname nmr_read_samples
# @keywords internal
nmr_read_samples_bruker <- function(sample_names, pulse_sequence = NULL,
                                    metadata_only = FALSE,
                                    overwrite_sample_names = NULL,
                                     ...) {
  if (length(sample_names) == 0) {
    stop("No samples to load")
  }
  # overwrite_sample_names is used when downloading a temporary file from
  # irods. In that case we want to preserve the irods path if possible
  if (!is.null(overwrite_sample_names)) {
    stopifnot(length(sample_names) == length(overwrite_sample_names))
  } else {
    # We overwrite with the same name:
    overwrite_sample_names <- sample_names
  }
  if (show_progress_bar()) {
    prgrs <- "text"
  } else {
    prgrs <- "none"
  }
  list_of_samples <-
    plyr::llply(seq_len(length(sample_names)),
                function(sampl_idx, ...) {
                  sampl <- sample_names[sampl_idx]
                  overwr <- overwrite_sample_names[sampl_idx]
                  is_zip <- NULL
                  loaded_sample <- tryCatch({
                    sampl <- normalizePath(sampl)
                    if (grepl("\\.zip$", sampl)) {
                      is_zip <- TRUE
                      NMRExperiment <- gsub(pattern = "\\.zip$", replacement = "", basename(overwr))
                      sampl_temp_dir <- tempfile(pattern = paste0("nmr_sample_", NMRExperiment, "_"))
                      utils::unzip(sampl, exdir = sampl_temp_dir)
                      sampl_dir <- normalizePath(file.path(sampl_temp_dir, NMRExperiment))
                    } else {
                      is_zip <- FALSE
                      sampl_dir <- sampl
                    }
                    # Ignore internal TopSpin directory used for sample processing
                    if (basename(sampl_dir) == "98888") {
                      return(NULL)
                    }
                    meta <- read_bruker_metadata(sampl_dir)
                    meta$info$sample_path <- overwr
                    if (!is.null(pulse_sequence) &&
                        toupper(meta$info$pulse_sequence) != toupper(pulse_sequence)) {
                      return(NULL)
                    }
                    if (metadata_only) {
                      pdata <- NULL
                    } else {
                      pdata <- read_bruker_pdata(sample_path = sampl_dir, ...)
                    }
                    output <- bruker_merge_meta_pdata(meta, pdata)
                    return(output)
                  }, error = function(err) {
                    warning("Error loading sample: ", sampl)
                    msg <- conditionMessage(err)
                    message(msg)
                    return(NULL)
                  }, finally = {
                    if (is_zip) {
                      unlink(sampl_temp_dir, recursive = TRUE)
                    }
                  })
                  return(loaded_sample)
                },
                ...,
                .progress = prgrs)

  # Remove samples that could not be loaded:
  any_error <- vapply(X = list_of_samples, FUN = is.null, FUN.VALUE = logical(1))
  list_of_samples <- list_of_samples[!any_error]
  # The number of samples
  num_samples <- length(list_of_samples)

  if (length(list_of_samples) == 0) {
    stop("No samples loaded")
  }

  # merge the sample information:
  all_fields <- unique(do.call(c, lapply(list_of_samples, function(x) names(x))))

  axis_fields <- "axis"
  data_fields <- all_fields[grepl(pattern = "^data_.*", x = all_fields)]
  metadata_fields <- setdiff(all_fields, c(axis_fields, data_fields))

  sample_meta <- list()
  for (meta_field in metadata_fields) {
    sample_meta[[meta_field]] <- list_of_lists_to_tibble(lapply(list_of_samples,
                                                                function(x) x[[meta_field]]))
    if (ncol(sample_meta[[meta_field]]) > 0) {
      colnames(sample_meta[[meta_field]]) <- paste(meta_field, colnames(sample_meta[[meta_field]]), sep = "_")
    }
  }

  samples <- list()
  metadata <- dplyr::bind_cols(sample_meta)
  names(metadata)[names(metadata) == 'info_NMRExperiment'] <- 'NMRExperiment'
  samples[["metadata"]] <- metadata

  if (!metadata_only) {
    for (data_field in data_fields) {
      samples[[data_field]] <- lapply(list_of_samples, function(x) x[[data_field]])
    }
    samples[["axis"]] <- lapply(list_of_samples, function(x) x[["axis"]])
  }
  samples[["num_samples"]] <- num_samples
  samples[["processing"]] <- list(data_loaded = !metadata_only,
                                  interpolation = FALSE,
                                  exclusion = FALSE,
                                  normalization = FALSE)
  class(samples) <- "nmr_dataset"
  return(samples)
}

# @rdname nmr_read_samples
nmr_read_samples_jdx <- function(sample_names, metadata_only = FALSE) {
  sample_names <- normalizePath(sample_names, mustWork = FALSE)
  raw_samples <- read_jdx(sample_names, metadata_only = metadata_only)
  # Assume 1-D
  samples <- list()
  if (!metadata_only) {
    block_with_data_per_sample <-
      vapply(raw_samples,
             FUN = function(sample) {
               block_with_xydata <-
                 vapply(sample$block,
                        function(block) "XYDATA" %in% names(block), logical(1))
               if (sum(block_with_xydata) == 1) {
                 return(which(block_with_xydata))
               } else {
                 return(-1)
               }
             }, numeric(1))
    if (any(block_with_data_per_sample == -1)) {
      stop("Error loading: ", sample_names[block_with_data_per_sample == -1])
    }
  }
  num_samples <- length(raw_samples)

  # Metadata:
  metadata <- dplyr::bind_rows(lapply(raw_samples, create_df_from_jdx_sample))
  metadata$file_name <- sample_names
  # Make a reasonable NMRExperiment:
  if (!("NMRExperiment" %in% colnames(metadata))) {
    NMRExperiments <- basename(sample_names)
    if (any(duplicated(NMRExperiments))) {
      NMRExperiments <- sample_names
    }
    metadata$NMRExperiment <- NMRExperiments
  }
  samples[["metadata"]] <- metadata
  if (!metadata_only) {
    samples[["data_1r"]] <- vector(mode = "list", length = num_samples)
    samples[["axis"]] <- vector(mode = "list", length(num_samples))
    for (sample_idx in seq_along(raw_samples)) {
      xydata <- raw_samples[[sample_idx]]$blocks[[block_with_data_per_sample[sample_idx]]][["XYDATA"]]
      samples[["data_1r"]][[sample_idx]] <- xydata$y
      samples[["axis"]][[sample_idx]] <- list(x = xydata$x)
    }
  }
  samples[["num_samples"]] <- num_samples
  samples[["processing"]] <- list(data_loaded = !metadata_only,
                                  interpolation = FALSE,
                                  exclusion = FALSE,
                                  normalization = FALSE)
  class(samples) <- "nmr_dataset"
  return(samples)
}


#' Object is of nmr_dataset class
#' @param x An object
#' @return \code{TRUE} if the object is an \code{\link{nmr_dataset}}, \code{FALSE} otherwise
#' @export
is.nmr_dataset <- function(x) inherits(x, "nmr_dataset")

# Needed for filter.nmr_dataset:
#' @importFrom dplyr filter
#' @export
dplyr::filter

#' Extract samples from the [nmr_dataset] based on metadata column criteria
#'
#' @param .data An [nmr_dataset] object
#' @param ... conditions, as in [dplyr]
#' @return an [nmr_dataset] object, with the matching rows
#' @importFrom dplyr filter
#' @export
filter.nmr_dataset <- function(.data, ...) {
  dots <- rlang::quos(...)
  meta <- nmr_get_metadata(.data)
  meta$tmp_row_idx <- seq_len(nrow(meta))
  indices_to_keep <- dplyr::filter(meta, !!! dots)$tmp_row_idx
  return(.data[indices_to_keep])
}

# From rlang::have_name
has_names <- function(x) {
  nms <- names(x)
  if (is.null(nms)) {
    rep(FALSE, length(x))
  }
  else {
    !(is.na(nms) | nms == "")
  }
}


#' Extract parts of an nmr_dataset
#' @param x an [nmr_dataset] object
#' @param i indices of the samples to keep
#' @return an nmr_dataset with the extracted samples
#' @examples
#' \dontrun{
#' data <- nm_read_samples_dir("your_dir")
#' data2 <- data[1:3] # get the first 3 samples
#' }
#' @export
`[.nmr_dataset` <- function(x, i) {
  output <- x
  output$metadata <- output$metadata[i, , drop = FALSE]
  data_fields <- names(output)[grepl(pattern = "^data_.*", x = names(output))]

  if (!output[["processing"]][["interpolation"]]) {
    # Without interpolation axis and data field contain a list, with the axis
    # or the data field for each sample
    output[["axis"]] <- output[["axis"]][i]
    for (data_field in data_fields) {
      output[[data_field]] <- output[[data_field]][i]
    }
  } else {
    # With interpolation the axis is shared among all samples
    # and the data_field is a matrix, with the first dimension of the matrix
    # being the sample index.
    for (data_field in data_fields) {
      dimensionality <- length(dim(output[[data_field]]))
      if (dimensionality == 2) {
        output[[data_field]] <- output[[data_field]][i, , drop = FALSE]
      } else if (dimensionality == 3) {
        output[[data_field]] <- output[[data_field]][i, , drop = FALSE]
      } else {
        stop("[.nmr_dataset not implemented for dimensionality", dimensionality, ".")
      }
    }
  }
  output$num_samples <- nrow(output$metadata)
  return(output)
}

#' Add metadata to an nmr_dataset object
#' 
#' This is useful to add metadata to datasets that can be later used for
#' plotting spectra or further analysis (PCA...).
#' 
#' @param nmr_data an [nmr_dataset] object
#' @param metadata A data frame with metadata to add
#' @param by A column name of both the `nmr_dataset$metadata` and the metadata
#' data.frame. If you want to merge two columns with different headers you can
#' use a named character vector `c("NMRExperiment" = "ExperimentNMR")` where
#' the left side is the column name of the `nmr_dataset$metadata` and the right side is
#' the column name of the metadata data frame.
#' 
#' @return
#' The nmr_dataset object with the added metadata
#' @export
nmr_add_metadata <- function(nmr_data, metadata, by = "NMRExperiment") {
  nmr_meta <- nmr_get_metadata(nmr_data)
  by_left <- ifelse(is.null(names(by)), by, names(by))
  existing_vars <- base::setdiff(colnames(nmr_meta), by_left)
  conflict <- base::intersect(existing_vars, colnames(metadata))
  # We must ensure metadata[[by]] is unique:
  metadata <- dplyr::distinct(metadata, !!!by, .keep_all = TRUE)
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
  nmr_data$metadata <- nmr_meta_new
  nmr_data
}
