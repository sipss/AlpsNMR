#' nmr_dataset_1D (S3 class)
#'
#' An `nmr_dataset_1D` represents a set of 1D interpolated NMR samples.
#' It is defined as an S3 class, and it can be treated as a regular list.
#'
#' It currently has the following elements:
#'
#'
#' - `metadata`: A list of data frames. Each data frame contains metadata of
#' a given area (acquisition parameters, preprocessing parameters, general sample information...)
#'
#' - `axis`: A numeric vector with the chemical shift axis in ppm.
#'
#' - `data_1r`: A matrix with one sample on each row and the chemical
#' shifts in the columns.
#'
#'
#' @name nmr_dataset_1D
#' @family AlpsNMR dataset objects
NULL

#' Validate 1D nmr datasets
#' @name validate_nmr_dataset
#' @param nmr_dataset_1D An [nmr_dataset_1D] object
#' @return The [nmr_dataset_1D] unchanged
#' 
#' This function is useful for its side-effects. Stopping in case of error
#'
#' @family class helper functions
#' @family nmr_dataset_1D functions
#' @export
#' @examples 
#' dir_to_demo_dataset <- system.file("dataset-demo", package = "AlpsNMR")
#' dataset <- nmr_read_samples_dir(dir_to_demo_dataset)
#' dataset_1D <- nmr_interpolate_1D(dataset, axis = c(min = -0.5, max = 10, by = 2.3E-4))
#' dataset_1D_validated <- validate_nmr_dataset_1D(dataset_1D)
#' 
validate_nmr_dataset_1D <- function(nmr_dataset_1D) {
    validate_nmr_dataset_family(nmr_dataset_1D)
    abort_if_not(
        inherits(nmr_dataset_1D, "nmr_dataset_1D"),
        message = "Not an nmr_dataset_1D"
    )
    
    abort_if_not(
        "axis" %in% names(nmr_dataset_1D),
        message = "nmr_dataset_1D must have a ppm axis"
    )
    abort_if_not(
        "data_1r" %in% names(nmr_dataset_1D),
        message = "nmr_dataset_1D must have a data_1r matrix"
    )
    
    ppm_axis <- nmr_dataset_1D[["axis"]]
    data_1r <- nmr_dataset_1D[["data_1r"]]
    abort_if_not(
        is.vector(ppm_axis) && is.numeric(ppm_axis),
        message = "axis must be a numeric vector"
    )
    abort_if_not(
        is.matrix(data_1r) && is.numeric(data_1r),
        message = "data_1r must be a numeric matrix"
    )
    
    abort_if_not(
        length(ppm_axis) == ncol(data_1r),
        message = "ppm axis does not have a length equal to ncol(data_1r)"
    )
    num_samples <- nrow(data_1r)
    
    abort_if_not(
        num_samples == nmr_dataset_1D[["num_samples"]],
        message = "The num_samples value does not match nrow(data_1r)"
    )
    
    if (!"excluded_regions" %in% names(nmr_dataset_1D)) {
        rlang::warn(
            message = c(
                'The dataset should have a "excluded_regions" element with the excluded regions.',
                "i" = paste0(
                    'If you saved and restored a dataset from a previous AlpsNMR version, and you ',
                    'know the regions you excluded, you can set it manually: if `x` is your ',
                    'dataset, just use:\n  `x[["excluded_regions"]] <- list(water = c(4.7, 5.0))`\n',
                    '  where the list can be just empty or whatever you passed to nmr_exclude_regions().'
                ),
                "i" = "Otherwise, AlpsNMR will assume there are no excluded regions on some peak detection stages."
            )
        )
        nmr_dataset_1D[["excluded_regions"]] <- list()
    }
    nmr_dataset_1D
}

#' Creates a new 1D nmr_dataset object from scratch
#'
#' @name new_nmr_dataset_1D
#' @param ppm_axis A numeric vector with the ppm values for the columns of data_1r
#' @param data_1r A numeric matrix with one NMR spectrum on each row
#' @param metadata A list of data frames with at least the `NMRExperiment` column
#'
#' @importFrom glue glue
#' @family class helper functions
#' @family nmr_dataset_1D functions
#' @export
#' @return Creates a new 1D nmr_dataset object from scratch
#' 
#' @examples
#' # Create a random spectra matrix
#' nsamp <- 12
#' npoints <- 20
#' dummy_ppm_axis <- seq(from = 0.2, to = 10, length.out = npoints)
#' dummy_spectra_matrix <- matrix(runif(nsamp*npoints), nrow = nsamp, ncol = npoints)
#' metadata <- list(external = data.frame(NMRExperiment = paste0("Sample", 1:12),
#'                                        DummyClass = c("a", "b"),
#'                                        stringsAsFactors = FALSE))
#' dummy_nmr_dataset_1D <- new_nmr_dataset_1D(ppm_axis = dummy_ppm_axis,
#'                                            data_1r = dummy_spectra_matrix,
#'                                            metadata = metadata)
#'                                                  
#'                                                  
new_nmr_dataset_1D <- function(ppm_axis, data_1r, metadata) {
    samples <- list()
    samples[["metadata"]] <- metadata
    samples[["data_1r"]] <- data_1r
    samples[["axis"]] <- ppm_axis
    samples[["num_samples"]] <- nrow(data_1r)
    samples[["excluded_regions"]] <- list()
    class(samples) <- c("nmr_dataset_1D", "nmr_dataset_family")
    validate_nmr_dataset_1D(samples)
}

#' Object is of [nmr_dataset_1D] class
#' @param x an [nmr_dataset_1D] object
#' @return `TRUE` if the object is an [nmr_dataset_1D], `FALSE` otherwise
#' @export
#' @family class helper functions
#' @family nmr_dataset_1D functions
#' @examples 
#' dir_to_demo_dataset <- system.file("dataset-demo", package = "AlpsNMR")
#' dataset <- nmr_read_samples_dir(dir_to_demo_dataset)
#' dataset_1D <- nmr_interpolate_1D(dataset, axis = c(min = -0.5, max = 10, by = 2.3E-4))
#' result <- is(dataset_1D)
is.nmr_dataset_1D <- function(x)
    inherits(x, "nmr_dataset_1D")

#' print for nmr_dataset_1D
#' @param x an [nmr_dataset_1D] object
#' @param ... for future use
#' @family class helper functions
#' @family nmr_dataset_1D functions
#' @return print for nmr_dataset_1D
#' @export
#' @examples 
#' dir_to_demo_dataset <- system.file("dataset-demo", package = "AlpsNMR")
#' dataset <- nmr_read_samples_dir(dir_to_demo_dataset)
#' dataset_1D <- nmr_interpolate_1D(dataset, axis = c(min = -0.5, max = 10, by = 2.3E-4))
#' print(dataset_1D)
print.nmr_dataset_1D <- function(x, ...) {
    cat(format(x, ...), "\n")
    invisible(x)
}

#' format for nmr_dataset_1D
#' @param x an [nmr_dataset_1D] object
#' @param ... for future use
#' @family class helper functions
#' @family nmr_dataset_1D functions
#' @return format for nmr_dataset_1D
#' @export
#' @examples 
#' dir_to_demo_dataset <- system.file("dataset-demo", package = "AlpsNMR")
#' dataset <- nmr_read_samples_dir(dir_to_demo_dataset)
#' dataset_1D <- nmr_interpolate_1D(dataset, axis = c(min = -0.5, max = 10, by = 2.3E-4))
#' format(dataset_1D)
format.nmr_dataset_1D <- function(x, ...) {
    paste0("An nmr_dataset_1D (", x$num_samples, " samples)")
}

#' Extract parts of an nmr_dataset_1D
#' @param x an [nmr_dataset_1D] object
#' @param i indices of the samples to keep
#' @return an nmr_dataset_1D with the extracted samples
#' @family subsetting functions
#' @family nmr_dataset_1D functions
#' @export
#' @examples 
#' dir_to_demo_dataset <- system.file("dataset-demo", package = "AlpsNMR")
#' dataset <- nmr_read_samples_dir(dir_to_demo_dataset)
#' dataset_1D <- nmr_interpolate_1D(dataset, axis = c(min = -0.5, max = 10, by = 2.3E-4))
#' dataset_1D[0]
`[.nmr_dataset_1D` <- function(x, i) {
    output <- x
    output$metadata <- purrr::map(output$metadata, function(metad) {
        metad[i, , drop = FALSE]
    })
    output[["data_1r"]] <- output[["data_1r"]][i, , drop = FALSE]
    output$num_samples <- nrow(output$metadata[[1]])
    validate_nmr_dataset_1D(output)
}

#' Export 1D NMR data to a CSV file
#'
#' @param nmr_dataset An [nmr_dataset_1D] object
#' @param filename The csv filename
#' 
#' @return The nmr_dataset object (unmodified)
#' @export 
#' @examples 
#' dir_to_demo_dataset <- system.file("dataset-demo", package = "AlpsNMR")
#' dataset <- nmr_read_samples_dir(dir_to_demo_dataset)
#' dataset_1D <- nmr_interpolate_1D(dataset, axis = c(min = -0.5, max = 10, by = 2.3E-4))
#' #nmr_export_data_1r(dataset_1D, "exported_nmr_dataset")
nmr_export_data_1r <- function(nmr_dataset, filename) {
    # FIXME: remove me (nmr_data() covers for this)
    abort_if_not(
        is.nmr_dataset_1D(nmr_dataset),
        message = "An nmr_dataset_1D should be given"
    )
    data_1r <- nmr_data(nmr_dataset)
    utils::write.csv(data_1r, file = filename, row.names = FALSE)
    nmr_dataset
}


#' Export 1D NMR data to SummarizedExperiment
#'
#' @param nmr_dataset An [nmr_dataset_1D] object
#'
#' @return SummarizedExperiment An SummarizedExperiment object (unmodified)
#' @export 
#' @examples 
#' dir_to_demo_dataset <- system.file("dataset-demo", package = "AlpsNMR")
#' dataset <- nmr_read_samples_dir(dir_to_demo_dataset)
#' dataset_1D <- nmr_interpolate_1D(dataset, axis = c(min = -0.5, max = 10, by = 2.3E-4))
#' se <- nmr_data_1r_to_SummarizedExperiment(dataset_1D)
nmr_data_1r_to_SummarizedExperiment <- function(nmr_dataset) {
    if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
        rlang::abort("Please install SummarizedExperiment.")
    }
    abort_if_not(
        is.nmr_dataset_1D(nmr_dataset), 
        message = "An nmr_dataset_1D should be given"
    )
    data_1r <- nmr_data(nmr_dataset)
    # SummarizedExperiment work trasposed
    SummarizedExperiment::SummarizedExperiment(assays=list(data_1r=data_1r),
                         metadata = nmr_meta_get(nmr_dataset),
                         colData=t(data_1r))
}

#' Import SummarizedExperiment as 1D NMR data
#'
#' @param se An SummarizedExperiment object
#'
#' @return nmr_dataset An [nmr_dataset_1D] object (unmodified)
#' @export 
#' @examples 
#' dir_to_demo_dataset <- system.file("dataset-demo", package = "AlpsNMR")
#' dataset <- nmr_read_samples_dir(dir_to_demo_dataset)
#' dataset_1D <- nmr_interpolate_1D(dataset, axis = c(min = -0.5, max = 10, by = 2.3E-4))
#' se <- nmr_data_1r_to_SummarizedExperiment(dataset_1D)
#' dataset_1D <- SummarizedExperiment_to_nmr_data_1r(se)
SummarizedExperiment_to_nmr_data_1r <- function(se) {
    if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
        rlang::abort("Please install the SummarizedExperiment package")
    }
    if (!requireNamespace("S4Vectors", quietly = TRUE)) {
        rlang::abort("Please install the S4Vectors package")
    }
    
    meta <- S4Vectors::metadata(se)
    
    col_names <- names(meta)
    nmr_meta <- list()
    
    #Orig, NMRExperiment
    NMRE <- which(col_names[]=="NMRExperiment")
    num_NMRE <- length(unlist(meta[NMRE]))
    meta_orig <- as.data.frame(matrix(unlist(meta[NMRE]), nrow=num_NMRE))
    colnames(meta_orig) <- names(meta[NMRE])
    nmr_meta[["orig"]] <- meta_orig
    meta <- meta[-NMRE]
    col_names <- col_names[-NMRE]
    #External
    col_names_split <- stringr::str_split(col_names, "_", simplify = TRUE)
    if(dim(col_names_split)[2]<2){
        #All is considered external
        meta_external <- as.data.frame(matrix(unlist(meta), nrow=num_NMRE))
        meta_external <- cbind(meta_orig, meta_external)
        colnames(meta_external) <- c(colnames(meta_orig), names(meta))
        nmr_meta[["external"]] <- meta_external
    } else {
        external_col <- meta[col_names_split[,2]==""]
        if(length(external_col)>0){
            meta_external <- as.data.frame(
                matrix(unlist(meta[col_names_split[,2]==""]), nrow=num_NMRE))
            meta_external <- cbind(meta_orig, meta_external)
            colnames(meta_external) <- c(colnames(meta_orig),
                                         names(meta[col_names_split[,2]==""]))
            nmr_meta[["external"]] <- meta_external
            meta <- meta[!col_names_split[,2]==""]
            col_names <- col_names[!col_names_split[,2]==""]
        } else {
            meta_external <- meta_orig
            nmr_meta[["external"]] <- meta_external
        }
        
        col_names_split <- stringr::str_split(col_names, "_", simplify = TRUE)
        df_names <- unique(col_names_split[,1])
        for(i in seq_len(length(df_names))){
            num <- which(col_names_split[,1]==df_names[i])
            meta_df <- as.data.frame(matrix(unlist(meta[num[1]:max(num)],
                                                   recursive = FALSE), nrow=num_NMRE))
            colnames(meta_df) <- names(meta[num[1]:max(num)])
            meta_df <- cbind(meta_orig, meta_df)
            nmr_meta[[df_names[i]]] <- meta_df
        }
    }
    new_nmr_dataset_1D(as.numeric(colnames(se)), 
                       as.matrix(t(SummarizedExperiment::colData(se))), nmr_meta)
}