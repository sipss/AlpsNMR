#' nmr_dataset_peak_table (S3 class)
#'
#' An `nmr_dataset_peak_table` represents a peak table with metadata.
#' It is defined as an S3 class, and it can be treated as a regular list.
#'
#' - `metadata`: A list of data frames. Each data frame contains metadata. Usually
#' the list only has one data frame named "external".
#'
#' - `peak_table`: A matrix with one sample on each row and the peaks in the
#' columns
#'
#' @name nmr_dataset_peak_table
NULL

#' Validate nmr_dataset_peak_table objects
#' @param nmr_dataset_peak_table An [nmr_dataset_peak_table] object
#' @return The [nmr_dataset_peak_table] unchanged
#'
#' @family nmr_dataset_peak_table functions
#' @family class helper functions
#' @export
validate_nmr_dataset_peak_table <- function(nmr_dataset_peak_table) {
    validate_nmr_dataset_family(nmr_dataset_peak_table)
    abort_if_not(
        inherits(nmr_dataset_peak_table, "nmr_dataset_peak_table"),
        message = "Not an nmr_dataset_peak_table"
    )
    
    abort_if_not(
        "peak_table" %in% names(nmr_dataset_peak_table),
        message = "nmr_dataset_peak_table must have a peak_table matrix"
    )
    
    peak_table <- nmr_dataset_peak_table[["peak_table"]]
    abort_if_not(
        is.matrix(peak_table) && is.numeric(peak_table),
        message = "peak_table must be a numeric matrix"
    )
    
    num_samples <- nrow(peak_table)
    
    abort_if_not(
        num_samples == nmr_dataset_peak_table[["num_samples"]],
        message = "The num_samples value does not match nrow(peak_table)"
    )
    
    nmr_dataset_peak_table
}

#' Creates a new nmr_dataset_peak_table object from scratch
#' 
#' @param peak_table A numeric matrix with one NMR spectrum on each row
#' @param metadata A list of data frames with at least the `NMRExperiment` column
#' @return Creates a new nmr_dataset_peak_table object from scratch
#' @name new_nmr_dataset_peak_table 
#' @importFrom glue glue
#' @family nmr_dataset_peak_table functions
#' @family class helper functions
#' @export
#' @examples
#' dir_to_demo_dataset <- system.file("dataset-demo", package = "AlpsNMR")
#' dataset <- nmr_read_samples_dir(dir_to_demo_dataset)
#' dataset_1D <- nmr_interpolate_1D(dataset, axis = c(min = -0.5, max = 10, by = 2.3E-4))
#' meta <- file.path(dir_to_demo_dataset, "dummy_metadata.xlsx")
#' metadata <- readxl::read_excel(meta, sheet = 1)
#' dataset_1D <- nmr_meta_add(dataset_1D, metadata = metadata, by = "NMRExperiment")
#' metadata <- list(external = dataset_1D[["metadata"]][["external"]])
#' peak_table <- nmr_data(dataset_1D)
#' new <- new_nmr_dataset_peak_table(peak_table, metadata)
#' 
new_nmr_dataset_peak_table <- function(peak_table, metadata) {
    samples <- list()
    samples[["metadata"]] <- metadata
    samples[["peak_table"]] <- as.matrix(peak_table)
    samples[["num_samples"]] <- nrow(peak_table)
    class(samples) <-
        c("nmr_dataset_peak_table", "nmr_dataset_family")
    validate_nmr_dataset_peak_table(samples)
    samples
}

#' Object is of [nmr_dataset_peak_table] class
#' @param x an [nmr_dataset_peak_table] object
#' @return `TRUE` if the object is an `nmr_dataset_peak_table`, `FALSE` otherwise
#' @export
#' @family nmr_dataset_peak_table functions
#' @family class helper functions
#' @examples
#' dir_to_demo_dataset <- system.file("dataset-demo", package = "AlpsNMR")
#' dataset <- nmr_read_samples_dir(dir_to_demo_dataset)
#' dataset_1D <- nmr_interpolate_1D(dataset, axis = c(min = -0.5, max = 10, by = 2.3E-4))
#' meta <- file.path(dir_to_demo_dataset, "dummy_metadata.xlsx")
#' metadata <- readxl::read_excel(meta, sheet = 1)
#' dataset_1D <- nmr_meta_add(dataset_1D, metadata = metadata, by = "NMRExperiment")
#' metadata <- list(external = dataset_1D[["metadata"]][["external"]])
#' peak_table <- nmr_data(dataset_1D)
#' new <- new_nmr_dataset_peak_table(peak_table, metadata)
#' is(new)
#' 
is.nmr_dataset_peak_table <-
    function(x)
        inherits(x, "nmr_dataset_peak_table")

#' print for nmr_dataset_peak_table
#' @param x an [nmr_dataset_peak_table] object
#' @param ... for future use
#' @export
#' @return print for nmr_dataset_peak_table
#' @family nmr_dataset_peak_table functions
#' @family class helper functions
#' @examples
#' dir_to_demo_dataset <- system.file("dataset-demo", package = "AlpsNMR")
#' dataset <- nmr_read_samples_dir(dir_to_demo_dataset)
#' dataset_1D <- nmr_interpolate_1D(dataset, axis = c(min = -0.5, max = 10, by = 2.3E-4))
#' meta <- file.path(dir_to_demo_dataset, "dummy_metadata.xlsx")
#' metadata <- readxl::read_excel(meta, sheet = 1)
#' dataset_1D <- nmr_meta_add(dataset_1D, metadata = metadata, by = "NMRExperiment")
#' metadata <- list(external = dataset_1D[["metadata"]][["external"]])
#' peak_table <- nmr_data(dataset_1D)
#' new <- new_nmr_dataset_peak_table(peak_table, metadata)
#' new
print.nmr_dataset_peak_table <- function(x, ...) {
    cat(format(x, ...), "\n")
    invisible(x)
}

#' Format for nmr_dataset_peak_table
#' @param x an [nmr_dataset_peak_table] object
#' @param ... for future use
#' @export
#' @return Format for nmr_dataset_peak_table
#' @family nmr_dataset_peak_table functions
#' @family class helper functions
#' @examples
#' dir_to_demo_dataset <- system.file("dataset-demo", package = "AlpsNMR")
#' dataset <- nmr_read_samples_dir(dir_to_demo_dataset)
#' dataset_1D <- nmr_interpolate_1D(dataset, axis = c(min = -0.5, max = 10, by = 2.3E-4))
#' meta <- file.path(dir_to_demo_dataset, "dummy_metadata.xlsx")
#' metadata <- readxl::read_excel(meta, sheet = 1)
#' dataset_1D <- nmr_meta_add(dataset_1D, metadata = metadata, by = "NMRExperiment")
#' metadata <- list(external = dataset_1D[["metadata"]][["external"]])
#' peak_table <- nmr_data(dataset_1D)
#' new <- new_nmr_dataset_peak_table(peak_table, metadata)
#' format(new)
format.nmr_dataset_peak_table <- function(x, ...) {
    paste0(
        "An nmr_dataset_peak_table (",
        x$num_samples,
        " samples, and ",
        ncol(x$peak_table),
        " peaks)"
    )
}

#' Extract parts of an nmr_dataset_peak_table
#' @param x an [nmr_dataset_peak_table] object
#' @param i indices of the samples to keep
#' @return an nmr_dataset_peak_table with the extracted samples
#' @family subsetting functions
#' @family nmr_dataset_peak_table functions
#' @export
#' @examples
#' dir_to_demo_dataset <- system.file("dataset-demo", package = "AlpsNMR")
#' dataset <- nmr_read_samples_dir(dir_to_demo_dataset)
#' dataset_1D <- nmr_interpolate_1D(dataset, axis = c(min = -0.5, max = 10, by = 2.3E-4))
#' meta <- file.path(dir_to_demo_dataset, "dummy_metadata.xlsx")
#' metadata <- readxl::read_excel(meta, sheet = 1)
#' dataset_1D <- nmr_meta_add(dataset_1D, metadata = metadata, by = "NMRExperiment")
#' metadata <- list(external = dataset_1D[["metadata"]][["external"]])
#' peak_table <- nmr_data(dataset_1D)
#' new <- new_nmr_dataset_peak_table(peak_table, metadata)
#' new[0]
`[.nmr_dataset_peak_table` <- function(x, i) {
    output <- x
    output$metadata <- purrr::map(output$metadata, function(metad) {
        metad[i, , drop = FALSE]
    })
    output[["peak_table"]] <-
        output[["peak_table"]][i, , drop = FALSE]
    output$num_samples <- nrow(output$metadata[[1]])
    validate_nmr_dataset_peak_table(output)
    return(output)
}

#' Export nmr_dataset_peak_table to SummarizedExperiment
#'
#' @param nmr_peak_table An [nmr_dataset_peak_table] object
#'
#' @return SummarizedExperiment object (unmodified)
#' @export 
#' @examples 
#' dir_to_demo_dataset <- system.file("dataset-demo", package = "AlpsNMR")
#' dataset <- nmr_read_samples_dir(dir_to_demo_dataset)
#' dataset_1D <- nmr_interpolate_1D(dataset, axis = c(min = -0.5, max = 10, by = 2.3E-4))
#' meta <- file.path(dir_to_demo_dataset, "dummy_metadata.xlsx")
#' metadata <- readxl::read_excel(meta, sheet = 1)
#' dataset_1D <- nmr_meta_add(dataset_1D, metadata = metadata, by = "NMRExperiment")
#' metadata <- list(external = dataset_1D[["metadata"]][["external"]])
#' peak_table <- nmr_data(dataset_1D)
#' nmr_peak_table <- new_nmr_dataset_peak_table(peak_table, metadata)
#' se <- nmr_dataset_peak_table_to_SummarizedExperiment(nmr_peak_table)
nmr_dataset_peak_table_to_SummarizedExperiment <- function(nmr_peak_table) {
    abort_if_not(
        inherits(nmr_peak_table, "nmr_dataset_peak_table"),
        message = "Not an nmr_dataset_peak_table"
    )
    peak_table <- nmr_peak_table[["peak_table"]]
    # SummarizedExperiment work trasposed
    SummarizedExperiment::SummarizedExperiment(assays=list(peak_table=peak_table),
                         metadata = nmr_meta_get(nmr_peak_table),
                         colData=t(peak_table))
}

#' Import SummarizedExperiment as mr_dataset_peak_table
#'
#' @param se An SummarizedExperiment object
#'
#' @return nmr_dataset_peak_table An [nmr_dataset_peak_table] object (unmodified)
#' @export 
#' @examples
#' dir_to_demo_dataset <- system.file("dataset-demo", package = "AlpsNMR")
#' dataset <- nmr_read_samples_dir(dir_to_demo_dataset)
#' dataset_1D <- nmr_interpolate_1D(dataset, axis = c(min = -0.5, max = 10, by = 2.3E-4))
#' meta <- file.path(dir_to_demo_dataset, "dummy_metadata.xlsx")
#' metadata <- readxl::read_excel(meta, sheet = 1)
#' dataset_1D <- nmr_meta_add(dataset_1D, metadata = metadata, by = "NMRExperiment")
#' metadata <- list(external = dataset_1D[["metadata"]][["external"]])
#' peak_table <- nmr_data(dataset_1D)
#' nmr_peak_table <- new_nmr_dataset_peak_table(peak_table, metadata)
#' se <- nmr_dataset_peak_table_to_SummarizedExperiment(nmr_peak_table)
#' nmr_peak_table <- SummarizedExperiment_to_nmr_dataset_peak_table(se)
SummarizedExperiment_to_nmr_dataset_peak_table <- function(se) {
    require_pkgs(pkg = c("SummarizedExperiment", "S4Vectors"))
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
    new_nmr_dataset_peak_table(as.matrix(t(SummarizedExperiment::colData(se))), nmr_meta)
}
