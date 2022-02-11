#' nmr_dataset (S3 class)
#'
#' An `nmr_dataset` represents a set of NMR samples.
#' It is defined as an S3 class, and it can be treated as a regular list.
#'
#' It currently has the following elements:
#'
#'
#' - `metadata`: A list of data frames. Each data frame contains metadata of
#' a given area (acquisition parameters, preprocessing parameters, general sample information...)
#'
#' - `axis`: A list with length equal to the dimensionality of the data.
#' For 1D spectra it is a list with a numeric vector
#'
#' - `data_*`: Data arrays with the actual spectra. The first index represents
#' the sample, the rest of the indices match the length of each `axis`.
#' Typically `data_1r` is a matrix with one sample on each row and the chemical
#' shifts in the columns.
#'
#' - `num_samples`: The number of samples in the dataset
#'
#' @name nmr_dataset
#' @family AlpsNMR dataset objects
#' @seealso [Functions to save and load these objects][load_and_save_functions]
NULL


#' Read NMR samples
#'
#' These functions load samples from files and return a [nmr_dataset].
#'
#' @name nmr_read_samples
#' @param sample_names A character vector with file or directory names.
#' @param samples_dir A directory that contains multiple samples
#' @param format Either "bruker" or "jdx"
#' @param metadata_only A logical, to load only metadata (default: `FALSE`)
#' @param pulse_sequence If it is set to a pulse sequence
#'                                             ("NOESY", "JRES", "CPMG"...) it will only load
#'                                             the samples that match that pulse sequence.
#' @param ... Arguments passed to [read_bruker_sample()] for data loading
#' @return a [nmr_dataset] object
NULL

#' @rdname nmr_read_samples
#' @family nmr_dataset functions
#' @family import/export functions
#' @export
#' @examples 
#' dir_to_demo_dataset <- system.file("dataset-demo", package = "AlpsNMR")
#' dataset <- nmr_read_samples_dir(dir_to_demo_dataset)
#' 
nmr_read_samples_dir <- function(samples_dir,
                                 format = "bruker",
                                 pulse_sequence = NULL,
                                 metadata_only = FALSE,
                                 ...) {
    nmr_read_samples_dir_internal(
        samples_dir = samples_dir,
        format = format,
        pulse_sequence = pulse_sequence,
        metadata_only = metadata_only,
        ...
    )
}

#' @noRd
#' @inheritParams nmr_read_samples_dir
#' @inheritParams nmr_read_samples_internal
nmr_read_samples_dir_internal <- function(samples_dir,
                                          format = "bruker",
                                          pulse_sequence = NULL,
                                          metadata_only = FALSE,
                                          overwrite_sample_names = NULL,
                                          ...) {
    samples_dir <- as.character(samples_dir)
    if (!dir.exists(samples_dir)) {
        stop("Invalid directory: ", samples_dir)
    }
    if (format == "bruker") {
        all_samples <-
            c(
                list.dirs(
                    path = samples_dir,
                    full.names = TRUE,
                    recursive = FALSE
                ),
                list.files(
                    path = samples_dir,
                    full.names = TRUE,
                    pattern = ".*zip$"
                )
            )
    } else if (format == "jdx") {
        all_samples <-
            list.files(path = samples_dir,
                       full.names = TRUE,
                       pattern = ".*jdx$")
    } else {
        stop("Unsupported sample format: ", format)
    }
    
    dataset <- nmr_read_samples_internal(
        sample_names = all_samples,
        format = format,
        pulse_sequence = pulse_sequence,
        metadata_only = metadata_only,
        overwrite_sample_names = overwrite_sample_names,
        ...
    )
    return(dataset)
}


#' @rdname nmr_read_samples
#' @export
#' @examples 
#' dir_to_demo_dataset <- system.file("dataset-demo", package = "AlpsNMR")
#' zip_files <- fs::dir_ls(dir_to_demo_dataset, glob = "*.zip")
#' dataset <- nmr_read_samples(sample_names = zip_files)
#' 
nmr_read_samples <- function(sample_names,
                             format = "bruker",
                             pulse_sequence = NULL,
                             metadata_only = FALSE,
                             ...) {
    nmr_read_samples_internal(
        sample_names = sample_names,
        format = format,
        pulse_sequence = pulse_sequence,
        metadata_only = metadata_only,
        ...
    )
}

#' @noRd
#' @inheritParams nmr_read_samples
#' @inheritParams nmr_read_samples_bruker
nmr_read_samples_internal <- function(sample_names,
                                      format = "bruker",
                                      pulse_sequence = NULL,
                                      metadata_only = FALSE,
                                      overwrite_sample_names = NULL,
                                      ...) {
    sample_names <- as.character(sample_names)
    if (format == "bruker") {
        samples <- nmr_read_samples_bruker(
            sample_names = sample_names,
            metadata_only = metadata_only,
            pulse_sequence = pulse_sequence,
            overwrite_sample_names = overwrite_sample_names,
            ...
        )
    } else if (format == "jdx") {
        # otherwise the jdx format
        samples <- nmr_read_samples_jdx(sample_names = sample_names,
                                        metadata_only = metadata_only)
    } else {
        stop("Unsupported format")
    }
    return(samples)
}

#' @param overwrite_sample_names This is only is used internally when downloading a
#'    temporary file from irods. In that case we will want to replace the temporary
#'    directory with the irods path if possible
#' @noRd
nmr_read_samples_bruker <-
    function(sample_names,
             pulse_sequence = NULL,
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

        warn_future_to_biocparallel()
        list_of_samples <- BiocParallel::bplapply(
            X = seq_along(sample_names),
            FUN = function(sampl_idx, ...) {
                sampl <- sample_names[sampl_idx]
                overwr <- overwrite_sample_names[sampl_idx]
                is_zip <- NULL
                loaded_sample <-
                    tryCatch({
                        sampl <- normalizePath(sampl)
                        if (grepl("\\.zip$", sampl)) {
                            is_zip <- TRUE
                            NMRExperiment <- gsub(
                                pattern = "\\.zip$",
                                replacement = "",
                                basename(overwr)
                            )
                            sampl_temp_dir <-
                                tempfile(pattern = paste0("nmr_sample_", NMRExperiment, "_"))
                            utils::unzip(sampl, exdir = sampl_temp_dir)
                            sampl_dir <-
                                normalizePath(file.path(sampl_temp_dir, NMRExperiment))
                        } else {
                            is_zip <- FALSE
                            sampl_dir <- sampl
                        }
                        # Ignore internal TopSpin directory used for sample processing
                        if (basename(sampl_dir) == "98888") {
                            return(NULL)
                        }
                        meta <- read_bruker_metadata(sampl_dir)
                        if (is_zip) {
                            meta$info$file_format <- "Zipped Bruker NMR directory"
                        }
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
            }
        )

        # Remove samples that could not be loaded:
        any_error <- purrr::map_lgl(list_of_samples, is.null)
        list_of_samples <- list_of_samples[!any_error]
        
        if (length(list_of_samples) == 0) {
            stop("No samples loaded")
        }
        
        # merge the sample information:
        all_fields <-
            unique(do.call(c, lapply(list_of_samples, function(x)
                names(x))))
        
        axis_fields <- "axis"
        data_fields <-
            all_fields[grepl(pattern = "^data_.*", x = all_fields)]
        metadata_fields <-
            setdiff(all_fields, c(axis_fields, data_fields))
        
        sample_meta <- list()
        for (meta_field in metadata_fields) {
            sample_meta[[meta_field]] <-
                list_of_lists_to_tibble(purrr::map(list_of_samples, meta_field))
            if (ncol(sample_meta[[meta_field]]) > 0) {
                colnames(sample_meta[[meta_field]]) <-
                    paste(meta_field, colnames(sample_meta[[meta_field]]), sep = "_")
            }
        }

        nmr_experiment_col <- sample_meta[["info"]][["info_NMRExperiment"]]
        nmr_experiment_col <- vctrs::vec_as_names(nmr_experiment_col, repair = "unique")
        sample_meta <- purrr::map(sample_meta,
                                  function(x) {
                                      x %>%
                                          dplyr::mutate(NMRExperiment = nmr_experiment_col) %>%
                                          dplyr::select(.data$NMRExperiment, dplyr::everything())
                                  })
        sample_meta[["external"]] = tibble::tibble(NMRExperiment = nmr_experiment_col)
        data_fields_full <- list()
        axis <- NULL
        if (!metadata_only) {
            for (data_field in data_fields) {
                data_fields_full[[data_field]] <- purrr::map(list_of_samples, data_field)
            }
            axis <- purrr::map(list_of_samples, "axis")
        }
        samples <- new_nmr_dataset(metadata = sample_meta,
                                   data_fields = data_fields_full,
                                   axis = axis)
        return(samples)
    }

# @rdname nmr_read_samples
nmr_read_samples_jdx <-
    function(sample_names, metadata_only = FALSE) {
        sample_names <- normalizePath(sample_names, mustWork = FALSE)
        raw_samples <-
            read_jdx(sample_names, metadata_only = metadata_only)
        # Assume 1-D
        if (!metadata_only) {
            block_with_data_per_sample <-
                vapply(
                    raw_samples,
                    FUN = function(sample) {
                        block_with_xydata <-
                            vapply(sample$block,
                                   function(block)
                                       "XYDATA" %in% names(block), logical(1))
                        if (sum(block_with_xydata) == 1) {
                            return(which(block_with_xydata))
                        } else {
                            return(-1)
                        }
                    },
                    numeric(1)
                )
            if (any(block_with_data_per_sample == -1)) {
                stop("Error loading: ", sample_names[block_with_data_per_sample == -1])
            }
        }
        num_samples <- length(raw_samples)
        
        # Metadata:
        metadata <-
            dplyr::bind_rows(lapply(raw_samples, create_df_from_jdx_sample))
        metadata$file_name <- sample_names
        # Make a reasonable NMRExperiment:
        if (!("NMRExperiment" %in% colnames(metadata))) {
            NMRExperiments <- basename(sample_names)
            if (any(duplicated(NMRExperiments))) {
                NMRExperiments <- sample_names
            }
            NMRExperiments <- vctrs::vec_as_names(NMRExperiments, repair = "unique")
            metadata$NMRExperiment <- NMRExperiments
        }
        metadata <-
            dplyr::select(metadata, .data$NMRExperiment, dplyr::everything())
        metadata_external = tibble::tibble(NMRExperiment = metadata$NMRExperiment)
        
        
        axis <- NULL
        data_fields <- list()
        if (!metadata_only) {
            data_fields[["data_1r"]] <-
                vector(mode = "list", length = num_samples)
            axis <- vector(mode = "list", length(num_samples))
            for (sample_idx in seq_along(raw_samples)) {
                xydata <-
                    raw_samples[[sample_idx]]$blocks[[block_with_data_per_sample[sample_idx]]][["XYDATA"]]
                data_fields[["data_1r"]][[sample_idx]] <- xydata$y
                axis[[sample_idx]] <- list(x = xydata$x)
            }
        }
        samples <-
            new_nmr_dataset(
                metadata = list(external = metadata_external,
                                metadata = metadata),
                data_fields = data_fields,
                axis = axis
            )
        return(samples)
    }


#' Object is of [nmr_dataset] class
#' @param x An object
#' @return `TRUE` if the object is an [nmr_dataset], `FALSE` otherwise
#' @family nmr_dataset manipulation functions
#' @export
#' @examples
#' dir_to_demo_dataset <- system.file("dataset-demo", package = "AlpsNMR")
#' dataset <- nmr_read_samples_dir(dir_to_demo_dataset)
#' is(dataset)
#' 
is.nmr_dataset <- function(x)
    inherits(x, "nmr_dataset")


#' Extract parts of an nmr_dataset
#' @param x an [nmr_dataset] object
#' @param i indices of the samples to keep
#' @return an nmr_dataset with the extracted samples
#' @family subsetting functions
#' @family nmr_dataset functions
#' @export
#' @examples
#' dir_to_demo_dataset <- system.file("dataset-demo", package = "AlpsNMR")
#' dataset <- nmr_read_samples_dir(dir_to_demo_dataset)
#' dataset2 <- dataset[1:3] # get the first 3 samples
#' 
`[.nmr_dataset` <- function(x, i) {
    output <- x
    output$metadata <- purrr::map(output$metadata, function(metad) {
        metad[i, , drop = FALSE]
    })
    data_fields <-
        names(output)[grepl(pattern = "^data_.*", x = names(output))]
    
    output[["axis"]] <- output[["axis"]][i]
    for (data_field in data_fields) {
        output[[data_field]] <- output[[data_field]][i]
    }
    output$num_samples <- nrow(output$metadata[[1]])
    validate_nmr_dataset(output)
    return(output)
}


#' Print for nmr_dataset
#' @param x an [nmr_dataset] object
#' @param ... for future use
#' @family class helper functions
#' @family nmr_dataset functions
#' @return Print for nmr_dataset
#' @export
#' @examples
#' dir_to_demo_dataset <- system.file("dataset-demo", package = "AlpsNMR")
#' dataset <- nmr_read_samples_dir(dir_to_demo_dataset)
#' print(dataset)
#' 
print.nmr_dataset <- function(x, ...) {
    cat(format(x, ...), "\n")
    invisible(x)
}

#' Format for nmr_dataset
#' @param x an [nmr_dataset] object
#' @param ... for future use
#' @family class helper functions
#' @family nmr_dataset functions
#' @return Format for nmr_dataset
#' @export
#' @examples
#' dir_to_demo_dataset <- system.file("dataset-demo", package = "AlpsNMR")
#' dataset <- nmr_read_samples_dir(dir_to_demo_dataset)
#' format(dataset)
#' 
format.nmr_dataset <- function(x, ...) {
    paste0("An nmr_dataset (", x$num_samples, " samples)")
}

#' Validate nmr_dataset objects
#' 
#' @param samples An nmr_dataset object
#' @family class helper functions
#' @family nmr_dataset functions
#' @export
#' @return Validate nmr_dataset objects
#' @name validate_nmr_dataset
#' @examples
#' dir_to_demo_dataset <- system.file("dataset-demo", package = "AlpsNMR")
#' dataset <- nmr_read_samples_dir(dir_to_demo_dataset)
#' validate_nmr_dataset(dataset)
#' 
validate_nmr_dataset <- function(samples) {
    validate_nmr_dataset_family(samples)
    abort_if_not(
        inherits(samples, "nmr_dataset"),
        message = "Not an nmr_dataset object"
    )
    samples
}

#' Create an nmr_dataset object
#' 
#' @param metadata A named list of data frames
#' @param data_fields A named list. Check the examples
#' @param axis A list. Check the examples
#' @family class helper functions
#' @family nmr_dataset functions
#' @name new_nmr_dataset 
#' @return Create an nmr_dataset object
#' @export
#' @return Create an nmr_dataset object
#' @examples
#' #
#' metadata_1D <- list(external = data.frame(NMRExperiment = c("10", "20")))
#' # Sample 10 and Sample 20 can have different lengths (due to different setups)
#' data_fields_1D <- list(data_1r = list(runif(16), runif(32)))
#' # Each sample has its own axis list, with one element (because this example is 1D)
#' axis_1D <- list(list(1:16), list(1:32))
#' my_1D_data <- new_nmr_dataset(metadata_1D, data_fields_1D, axis_1D)
#'
#' # Example for 2D samples
#' metadata_2D <- list(external = data.frame(NMRExperiment = c("11", "21")))
#' data_fields_2D <- list(data_2rr = list(matrix(runif(16*3), nrow=16, ncol=3),
#'                         runif(32*3), nrow=32, ncol=3))
#' # Each sample has its own axis list, with one element (because this example is 1D)
#' axis_2D <- list(list(1:16, 1:3), list(1:32, 1:3))
#' my_2D_data <- new_nmr_dataset(metadata_2D, data_fields_2D, axis_2D)
#'
new_nmr_dataset <- function(metadata, data_fields, axis) {
    samples <- list()
    samples[["metadata"]] <- metadata
    samples <- append(x = samples, values = data_fields)
    samples[["axis"]] <- axis
    samples[["num_samples"]] <- nrow(metadata[[1]])
    class(samples) <- c("nmr_dataset", "nmr_dataset_family")
    validate_nmr_dataset(samples)
    samples
}
