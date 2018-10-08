#' NIHSnmr: Routines for importing and processing NMR data
#'
#' NIHSnmr allows you to import NMR samples into R and applying basic data
#' processing to them. See an interactive example with
#' [NIHSnmr_interactive()]
#'
#' The following functions can be combined with the pipe. They
#' create or modify the [nmr_dataset] object.
#'
#' - [nmr_read_samples_dir()] or [nmr_read_samples()]
#' - [nmr_interpolate()]
#' - [nmr_exclude_region()]
#' - [nmr_normalize()]
#' - [plot()][plot.nmr_dataset()]
#' 
#' There are also functions to extract the metadata and submit the samples to
#' irods, see the example below.
#'
#' The [nmr_dataset] object is essentially a list, so it is easy to access
#' its components for further analysis.
#'
#' @examples
#' \dontrun{
#' library(NIHSnmr)
#' sample_names <- c("dataset/10/", "dataset/20/", "dataset/30/")
#' my_nmr_dataset <- nmr_read_samples(sample_names) %>%
#'   nmr_interpolate(axis1 = c(0.4, 10)) %>%
#'   nmr_exclude_region(exclude = list(water = c(4.7, 5))) %>%
#'   nmr_normalize(method = "area") %>%
#'   plot
#' }
#' \dontrun{
#' library(NIHSnmr)
#' sample_names <- c("dataset/10/", "dataset/20/", "dataset/30/")
#' nmrdata <- nmr_read_samples(sample_names, metadata_only = TRUE)
#' irods_metadata <- nmr_get_irods_meta(nmrdata)
#' # explore and confirm everything is in irods_metadata
#' # Now, let's compress the samples to zip files for irods
#' irods_metadata <- nmr_prepare_zip_files(meta_irods = irods_metadata,
#'                                         workdir = 'my_dataset_zip_files')
#' # Finally push the data to the irods directory
#' nmr_push_to_irods(irods_metadata, "/NIHSData/DUND-XXXXXX/study/Metabolomics")
#' }
"_PACKAGE"



#' PQN normalization
#' @noRd
#' @param spectra A matrix with one spectrum on each row
#' @return A matrix with one spectrum on each row (normalized)
norm_pqn <- function(spectra) {
  num_samples <- nrow(spectra)
  if (num_samples < 10) {
    warning("The Probabalistic Quotient Normalization requires several samples ",
            "to compute the median spectra. Your number of samples is low")
  }
  # Normalize to the area
  spectra <- spectra / rowSums(spectra)
  if (num_samples == 1) {
    # We have warned, and here there is nothing to do anymore
    warning("PQN is absurd with a single sample. We have normalized it to the area.")
    return(spectra)
  }
  # Move spectra above zero:
  if (any(spectra < 0)) {
    spectra <- spectra - min(spectra)
  }
  # Median of each ppm: (We need multiple spectra in order to get a reliable median!)
  m <- matrixStats::colMedians(spectra)
  # Divide at each ppm by its median:
  f <- spectra/m[col(spectra)]
  f <- matrixStats::rowMedians(f)
  # Divide each spectra by its f value
  spectra <- spectra / f
  spectra
}


#' Normalize NMR samples
#'
#' @param samples A [nmr_dataset] object
#' @param method The criteria to be used for normalization
#' @param values If `method == "value"` then values is a list
#'               with the normalization values. The list must be named as the
#'               data fields to normalize. Typically would be something like:
#'               `values = list(data_1r = c(val1, val2, val3))`.
#'               If `method == "region"` then values is the chemical shift
#'               region to integrate.
#' @return The [nmr_dataset] object, with the samples normalized
#' @export
nmr_normalize <- function(samples,
                          method = c("area", "max", "value", "region", "pqn", "none"),
                          values = NULL) {
  # This function does not consider >1D samples. Some things may work by chance,
  # but it needs testing and revision.

  method <- tolower(method[1])
  if (!(method %in% c("area", "max", "value", "region", "pqn", "none"))) {
    stop("Unknown method: ", method)
  }
  if (method == "none") {
    return(samples)
  }

  if (samples[["processing"]][["normalization"]]) {
    warning("Samples were already normalized. Applying further normalizations")
  }

  data_fields <- names(samples)[grepl(pattern = "^data_.*", x = names(samples))]
  if (method == "area") {
    norm_factor <- list()
    for (data_field in data_fields) {
      norm_factor[[data_field]] <- rowSums(samples[[data_field]])
    }
    samples[["processing"]][["normalization_params"]] <- list(method = method,
                                                              norm_factor = norm_factor)
  } else if (method == "max") {
    norm_factor <- list()
    for (data_field in data_fields) {
      norm_factor[[data_field]] <- apply(samples[[data_field]], 1, max)
    }
    samples[["processing"]][["normalization_params"]] <- list(method = method,
                                                              norm_factor = norm_factor)
  } else if (method == "value") {
    norm_factor <- values
    samples[["processing"]][["normalization_params"]] <- list(method = method,
                                                              norm_factor = norm_factor)
  } else if (method == "region") {
    norm_factor <- list()
    for (data_field in data_fields) {
      if (length(samples[["axis"]]) > 1) {
        stop("Region normalization not implemented for dimensionality > 1")
      }
      axis1 <- samples[["axis"]][[1]]
      region_range <- axis1 >= min(values) & axis1 <= max(values)
      norm_factor[[data_field]] <- rowSums(samples[[data_field]][,region_range])
    }
    samples[["processing"]][["normalization_params"]] <- list(method = method,
                                                              norm_factor = norm_factor)
  } else if (method == "pqn") {
    # only for 1D
    if (length(samples[["axis"]]) > 1) {
      stop("PQN normalization not implemented for dimensionality > 1")
    }
    for (data_field in data_fields) {
      samples[[data_field]] <- norm_pqn(samples[[data_field]])
    }
    samples[["processing"]][["normalization"]] <- TRUE
    return(samples)
  } else {
    stop("Unimplemented method: ", method)
  }
  for (data_field in data_fields) {
    samples[[data_field]] <- sweep(x = samples[[data_field]],
                                   MARGIN = 1,
                                   STATS = norm_factor[[data_field]],
                                   FUN = "/")
  }
  samples[["processing"]][["normalization"]] <- TRUE
  return(samples)
}


#' Exclude region from samples
#'
#' Sets to zero a given region (for instance to remove the water peak)
#' @param samples A [nmr_dataset] object
#' @param exclude A list with regions to be zeroed. Typically:
#'                `exclude = list(water = c(4.7, 5.0))`
#' @return The [nmr_dataset] object, with the regions excluded
#' @export
nmr_exclude_region <- function(samples, exclude = list(water = c(4.7, 5.0))) {
  if (is.null(exclude) || length(exclude) == 0) {
    return(samples)
  }
  if (!samples[["processing"]][["interpolation"]]) {
    stop("exclusion not implemented for non-interpolated samples")
  }

  samples[["processing"]][["exclusion_params"]] <- exclude
  data_fields <- names(samples)[grepl(pattern = "^data_.*", x = names(samples))]
  for (data_field in data_fields) {
    dimension <- length(dim(samples[[data_field]])) - 1 # first dim is samples
    for (i_region in 1:length(exclude)) {
      region <- exclude[[i_region]]
      if (dimension == 1) {
        excl_dim1 <- samples[["axis"]][[1]] >= region[1] & samples[["axis"]][[1]] <= region[2]
        samples[[data_field]][,excl_dim1] <- 0
      } else if (dimension == 2) {
        excl_dim1 <- samples[["axis"]][[1]] >= region[1] & samples[["axis"]][[1]] <= region[2]
        excl_dim2 <- samples[["axis"]][[2]] >= region[3] & samples[["axis"]][[2]] <= region[4]
        samples[[data_field]][,excl_dim1, excl_dim2] <- 0
      } else {
        stop("Not implemented error")
      }
    }
  }
  samples[["processing"]][["exclusion"]] <- TRUE
  return(samples)
}

#' Interpolate samples
#'
#' Interpolate the samples so they share the same axis and can be placed
#' in the same matrix
#'
#' For 2-D samples, we use two 1-D linear interpolations.
#'
#'
#' @param samples A [nmr_dataset] object
#' @param axis1,axis2 Axis given as `c(minimum, maximum, step)`
#' @return The [nmr_dataset] object, with interpolated samples
#' @export
nmr_interpolate <- function(samples,
                            axis1=c(min = 0.4, max = 10, by = 0.0008),
                            axis2=NULL) {

  # If samples have been interpolated don't do it again.
  if (samples[["processing"]][["interpolation"]]) {
    warning("Samples already interpolated. Skipping interpolation.")
    return(samples)
  }


  # Check if we can interpolate:

  # 1. Check that each of the loaded samples has the same dimensionality
  dimensions_per_sample <- vapply(samples$axis, length, numeric(1))
  if (any(dimensions_per_sample != dimensions_per_sample[1])) {
    # Have different dimensionality, can't merge
    warning("Samples have different dimensionality. Cant't interpolate.")
    return(samples)
  }

  dimensions_per_sample <- dimensions_per_sample[1]

  # 2. Check that we have the interpolation axis
  verify_axisn <- function(axisn, one_sample_axis) {
    if (is.null(axisn)) {
      # axisn will be built on sample_axis:
      axisn <- c(min = min(one_sample_axis),
                 max = max(one_sample_axis),
                 by = abs(one_sample_axis[2] - one_sample_axis[1]))
    } else {
      if (length(axisn) == 2) {
        axisn <- c(axisn, abs(one_sample_axis[2] - one_sample_axis[1]))
        names(axisn) <- c("min", "max", "by")
      }
      if (is.null(names(axisn))) {
        names(axisn) <- c("min", "max", "by")
      }
    }
    return(axisn)
  }
  if (dimensions_per_sample >= 1) {
    # axis1 will be based on sample 1, dimension 1
    axis1 <- verify_axisn(axis1, samples[["axis"]][[1]][[1]])
  }
  if (dimensions_per_sample >= 2) {
    # axis2 will be based on sample 1, dimension 2
    axis2 <- verify_axisn(axis2, samples[["axis"]][[1]][[2]])
  }
  if (dimensions_per_sample >= 3) {
    warning("Interpolation not implemented for dimensions > 2. Skipping interpolation.")
    return(samples)
  }

  # Do interpolation:
  orig_axis <- samples[["axis"]]
  data_fields <- names(samples)[grepl(pattern = "^data_.*", x = names(samples))]
  num_samples <- samples[["num_samples"]]
  if (dimensions_per_sample == 1) {
    axis1_full <- seq(from = axis1["min"], to = axis1["max"], by = axis1["by"])
    samples[["axis"]] <- list(axis1_full)
    for (data_field in data_fields) {
      if (show_progress_bar(samples$num_samples > 5)) {
        message("Interpolating ", data_field, "...")
      }
      data_matr <- matrix(NA, nrow = num_samples, ncol = length(axis1_full))
      tryCatch({
        pb <- NULL
        if (show_progress_bar(samples$num_samples > 5)) {
          pb <- utils::txtProgressBar(min = 0, max = samples$num_samples, style = 3)
        }
        for (i in seq_len(samples$num_samples)) {
          data_matr[i, ] <- signal::interp1(x = orig_axis[[i]][[1]],
                                            y = samples[[data_field]][[i]],
                                            xi = axis1_full, method = "spline")
          if (!is.null(pb)) {
            utils::setTxtProgressBar(pb, i)
          }
        }
      }, finally = {
        if (!is.null(pb)) {
          close(pb)
        }
      })
      samples[[data_field]] <- data_matr
    }
  } else if (dimensions_per_sample == 2) {
    # For 2-D samples, we use two 1-D linear interpolations instead of a bilinear one. (FIXME)
    # Akima needs similar axis for interpolation. We will scale the axes.
    # rescale <- function(x0, y0, x1, y1) {
    #   slope <- (y1 - y0)/(x1 - x0)
    #   intercept <- y0 - slope * x0
    #   rescale_fun <- function(x) {
    #     return(slope*x + intercept)
    #   }
    #   return(rescale_fun)
    # }
    axis1_full <- seq(from = axis1["min"], to = axis1["max"], by = axis1["by"])
    axis2_full <- seq(from = axis2["min"], to = axis2["max"], by = axis2["by"])
    # axis1_scalefun <- rescale(axis1["min"], 0, axis1["max"], 1)
    # axis1_scaled <- axis1_scalefun(axis1_full)
    # axis2_scalefun <- rescale(axis2["min"], 0, axis2["max"], axis1["by"]/axis2["by"])
    # axis2_scaled <- axis2_scalefun(axis2_full)

    for (data_field in data_fields) {
      if (show_progress_bar()) {
        message("Interpolating ", data_field, "...")
      }
      data_matr <- array(data = NA, dim = c(num_samples, length(axis1_full), length(axis2_full)))
      tryCatch({
        pb <- NULL
        if (show_progress_bar(samples$num_samples > 0)) {
          pb <- utils::txtProgressBar(min = 0, max = 2*samples$num_samples, style = 3)
        }
        for (i in seq_len(samples$num_samples)) {
          data_interp_1ax <- apply(samples[[data_field]][[i]], 2,
                                   function(col) {
                                     signal::interp1(x = orig_axis[[i]][[1]],
                                                     y = col,
                                                     xi = axis1_full,
                                                     method = "spline")
                                   })
          if (!is.null(pb)) {
            utils::setTxtProgressBar(pb, 2*i - 1)
          }
          data_interp_1ax <- t(apply(data_interp_1ax, 1,
                                     function(row) {
                                       signal::interp1(x = orig_axis[[i]][[2]],
                                                       y = row,
                                                       xi = axis2_full,
                                                       method = "spline")
                                     }))
          if (!is.null(pb)) {
            utils::setTxtProgressBar(pb, 2*i)
          }
          data_matr[i, , ] <- data_interp_1ax
          # I could not make this work with the time I had
          # orig_axis_i_1_scaled <- axis1_scalefun(orig_axis[[i]][[1]])
          # orig_axis_i_2_scaled <- axis2_scalefun(orig_axis[[i]][[2]])
          #
          # x <- rep(orig_axis_i_1_scaled, times = length(orig_axis_i_2_scaled))
          # y <- rep(orig_axis_i_2_scaled, each = length(orig_axis_i_1_scaled))
          # z <- samples[[data_field]][[i]]
          # z <- as.vector(z)
          # z0 <- akima::interp(x = x, y = y, z = z,
          #                     xo = axis1_scaled, yo = axis2_scaled,
          #                     linear = TRUE)
          # data_matr[i, , ] <- z0$z
        }
      }, finally = {
        if (!is.null(pb)) {
          close(pb)
        }
      })
      samples[["axis"]] <- list(axis1_full, axis2_full)
      samples[[data_field]] <- data_matr
    }
  }
  samples[["processing"]][["interpolation"]] <- TRUE
  return(samples)
}


#' Integrate regions
#' 
#' Integrate given regions and return a data frame with them
#' @param samples A [nmr_dataset] object
#' @param regions A named list. Each element of the list is a region,
#'                given as a named numeric vector of length two with the range
#'                to integrate. The name of the region will be the name of the
#'                column
#'
#' The integration is very naïve and consists of the sum of the intensities
#' in the given ppm range.
#' 
#' @return A data frame with the NMRExperiment column and one additional column
#'         for each given region.
#' @export
nmr_integrate_regions <- function(samples, regions) {
  areas <- purrr::map_dfc(regions, function(region) {
    to_sum <- samples$axis[[1]] >= min(region) & samples$axis[[1]] < max(region)
    rowSums(samples$data_1r[, to_sum])
  })
  dplyr::bind_cols(nmr_get_metadata(samples, "NMRExperiment"),
                   areas)
}


#' Get metadata
#' @param samples a [nmr_dataset] object
#' @param columns Columns to get. By default gets all the columns.
#' @param simplify Removes columns that are constant along all samples
#' @return a data frame with the injection metadata
#' @export
nmr_get_metadata <- function(samples, columns = NULL, simplify = FALSE) {
  metadata <- dplyr::left_join(samples[["metadata_ext"]],
                               samples[["metadata"]],
                               by = "NMRExperiment")
  
  # Default columns means all columns
  if (is.null(columns)) {
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
  metadata <- metadata[,columns, drop = FALSE]
  if (simplify) {
    metadata <- simplify_df(metadata)$diff
  }
  return(metadata)
}


#' Load a [nmr_dataset] from a file
#'
#' Loads a [nmr_dataset] saved with [nmr_dataset_save()]
#'
#' @param file_name The file name with the [nmr_dataset]
#' @return the [nmr_dataset] object stored in `file_name`
#' @export
nmr_dataset_load <- function(file_name) {
  return(readRDS(file_name))
}

#' Save an [nmr_dataset] object
#'
#' Wraps `save` to save the [nmr_dataset] object in `.RDS` format
#'
#' @param nmr_dataset The [nmr_dataset] object to save
#' @param file_name The output file name.
#' @param ... Optional arguments passed to [saveRDS].
#' @return the passed [nmr_dataset] object
#' @export
nmr_dataset_save <- function(nmr_dataset, file_name, ...) {
  saveRDS(nmr_dataset, file_name)
  return(nmr_dataset)
}

# This function can be removed once we know it is not needed by our users
nmr_dataset_load_old_and_save <- function(old_file_name, new_file_name) {
  nmr_dataset <- NULL
  load(old_file_name)
  if (is.null(nmr_dataset)) {
    stop("This was not saved with the NIHSnmr version 1.0 or lower")
  }
  nmr_dataset_save(nmr_dataset, new_file_name)
}


#' Read Free Induction Decay file
#' 
#' Reads an FID file. This is a very simple function.
#' 
#' @param sample_name A single sample name
#' @param endian Endianness of the fid file ("little" by default, use "big" if acqus$BYTORDA == 1)
#' @return A numeric vector with the free induction decay values
#' @export
nmr_read_bruker_fid <- function(sample_name, endian = "little") {
  fid_file <- file.path(sample_name, "fid")
  num_numbers <- file.size(fid_file)/8
  fid <- readBin(fid_file, what = "integer", n = num_numbers, size = 4, signed = TRUE, endian = endian)
  fid
}
