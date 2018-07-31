# These functions deal with reading Bruker NMR samples.


# Authors: Sergio Oller, Ivan Montoliu
# This code is a port of MATLAB routines written by Ivan Montoliu and Ziad Ramadan.
# The NMR import core functions are based on MATLAB code written by Nils Nyberg
# and the nmrglue python package.
#
# The nmrglue python package does seem to have a more comprenhensive set of
# import routines for NMR data (and in particular for Bruker data).
# http://www.nmrglue.com/


#' Reads a BRUKER parameter file to a list
#'
#' @param file_name File name of a BRUKER parameter file
#' @return list of parameter-value pairs
#' @keywords internal
#' @noRd
read_bruker_param <- function(file_name) {
  lines <- readLines(file_name)
  # There are vectors, values, stamps... many types of values
  # We will use regular expressions to detect each type of value and to
  # extract the relevant field name / field value.
  type_of_row <- NA*numeric(length(lines))
  all_matches <- vector("list", length(lines))

  # Definition of the regular expressions and the type name.
  # The order is important. If for a given line there is a match on the first
  # row, the rest of the regular expressions for that row are ignored.
  R <- as.data.frame(matrix(c(
    '^##\\$*(.+)=\\s?\\(\\d\\.\\.\\d+\\)(.+)$'       , 'ParVecVal',
    '^##\\$*(.+)=\\s?\\(\\d\\.\\.\\d+\\)$'           , 'ParVec',
    '^##\\$*(.+)=\\s?(.+)$'                          , 'ParVal',
    '^([^\\$#].*)$'                                  , 'Val',
    '^\\$\\$(.+)$'                                   , 'Stamp',
    '^##\\$*(.+)=$'                                  , 'ParEmpty',
    '^\\s*$'                                         , 'Empty'),
    ncol = 2, byrow = TRUE), stringsAsFactors = FALSE)
  colnames(R) <- c("pattern", "name")

  # Match each pattern to all lines. If there is a match and there was no
  # previous match then save the matching result.
  for (pat_idx in seq_len(nrow(R))) {
    lines_matched <- stringr::str_match_all(lines, R[["pattern"]][pat_idx])
    lines_with_match <- sapply(lines_matched, nrow) > 0
    new_lines_with_match <- is.na(type_of_row) & lines_with_match
    type_of_row[new_lines_with_match] <- pat_idx
    # Update the new lines that have match:
    all_matches[new_lines_with_match] <- lines_matched[new_lines_with_match]
  }
  # Use "ParVec"-like names instead of numbers, it is easier to read:
  if (anyNA(type_of_row)) {
    lines_no_match <- which(is.na(type_of_row))
    stop("Line without a match:", paste(lines_no_match, collapse = " "))
  }
  type_of_row <- R[["name"]][type_of_row]
  # Generate the output list. We will use the matched information of each line
  output <- list()
  for (i in seq_len(length(type_of_row))) {
    if (type_of_row[i] == "ParVecVal") {
      field_name <- all_matches[[i]][1,2]
      output[[field_name]] <- all_matches[[i]][1,4]
    } else if (type_of_row[i] == "ParVec") {
      field_name <- all_matches[[i]][1,2]
      output[[field_name]] <- NULL
    } else if (type_of_row[i] == "ParVal") {
      field_name <- all_matches[[i]][1,2]
      output[[field_name]] <- all_matches[[i]][1,3]
    } else if (type_of_row[i] == "Val") {
      if (is.null(output[[field_name]])) {
        output[[field_name]] <- all_matches[[i]][1,2]
      } else {
        output[[field_name]] <- paste(output[[field_name]], all_matches[[i]][1,2])
      }
    } else if (type_of_row[i] == "Stamp") {
      output[["Stamp"]] <- c(output[["Stamp"]], all_matches[[i]][1,2])
    } else if (type_of_row[i] == "ParEmpty") {
      field_name <- all_matches[[i]][1,2]
      output[[field_name]] <- NULL
    } else if (type_of_row[i] == "Empty") {
      # do nothing
    } else {
      # If this happens it means that the R matrix on top has a problem or that
      # there is a not-considered case here.
      stop("This should not have happened!")
    }
  }

  for (field in names(output)) {
    output[[field]] <- convert_field(output[[field]])
  }
  return(output)
}

convert_field <- function(element) {
  # Convert values to numeric if possible
  converted_val <- element
  tryCatch({
    # First remove spaces
    cleaned_spaces <- stringr::str_trim(element)
    if (nchar(cleaned_spaces) == 0) {
      return(cleaned_spaces)
    }
    # Then make sure there is a space between numbers (it's apparently optional):
    cleaned_spaces <- gsub("([0-9.])([+-])", "\\1 \\2", cleaned_spaces)
    # Then split by one or more spaces
    converted_val2 <- strsplit(cleaned_spaces, split = "[[:blank:]]+")[[1]]
    converted_val <- as.numeric(converted_val2)
  }, warning = function(warn) {
    # Do nothing
  }, error = function(err) {
    # Do nothing
  })
  # If it is still a string, remove leading < and trailing >
  if (is.character(converted_val) && length(converted_val) == 1) {
    if (grepl(pattern = "^\\s*<.*>\\s*$", x = converted_val)) {
      converted_val <- gsub(pattern = "^\\s*<(.*)>\\s*$", replacement = "\\1",
                            x = converted_val)
    }
  }
  return(converted_val)
}


#' Reads a Bruker procs file
#'
#' @param full_pdata_path The path to the procs files
#' @param procs_files a character vector with the procs file names. If
#'                    procs_files is NULL default sensible names are tested.
#' @return a list with a list of parameters for each procs file given
#' @keywords internal
#' @noRd
read_procs_file <- function(full_pdata_path, procs_files = NULL) {
  if (is.null(procs_files)) {
    procs_files <- c("procs", "proc2s", "proc3s", "proc4s")
    procs_files <- procs_files[file.exists(file.path(full_pdata_path, procs_files))]
  }
  output <- list()
  for (file in procs_files) {
    output[[file]] <- read_bruker_param(file_name = file.path(full_pdata_path, file))
  }
  return(output)
}

#' Reads all common acqus files
#'
#' @param sample_path A directory that corresponds to a sample
#' @param acqus_files The acqus files that you want to load. if \code{NULL} the
#'                    default set of acqus files is used.
#' @return a list with a list of parameters for each acqus file given
#' @keywords internal
#' @noRd
read_acqus_file <- function(sample_path, acqus_files = NULL) {
  if (is.null(acqus_files)) {
    acqus_files <- c("acqus", "acqu2s", "acqu3s", "acqu4s")
    acqus_files <- acqus_files[file.exists(file.path(sample_path, acqus_files))]
  }
  output <- list()
  for (file in acqus_files) {
    output[[file]] <- read_bruker_param(file_name = file.path(sample_path, file))
  }
  return(output)
}


read_levels <- function(full_pdata_path, endian = NULL, NC_proc = NULL) {
  levels_vec <- NULL
  # Read the level file if it exists
  # The old version (level) is a binary
  file_level <- file.path(full_pdata_path, "level")
  if (file.exists(file_level)) {
    if (is.null(endian) || is.null(NC_proc)) {
      warning("Can't read old level file without endian information and NC_proc")
    }
    lev <- read_bin_data(file_name = file_level, endian = endian)
    # The first two figures is the number of pos. and neg. levels
    levels_vec <- lev[3:length(lev)]
    # Adjust for NC-parameter
    levels_vec <- levels_vec/(2 ^ -NC_proc)
  } else {
    levels_vec <- NULL
  }

  # Read the clevel file if it exists
  # The new version clevel is a text file
  file_clevel <- file.path(full_pdata_path, "clevels")
  if (file.exists(file_clevel)) {
    clev <- read_bruker_param(file_name = file_clevel)
    if (clev[["LEVSIGN"]] == 1) {
      levels_vec <- clev[["LEVELS"]][clev[["LEVELS"]] > 0]
    } else if (clev[["LEVSIGN"]] == 2) {
      levels_vec <- clev[["LEVELS"]][clev[["LEVELS"]] < 0]
    } else if (clev[["LEVSIGN"]] == 3) {
      levels_vec <- clev$LEVELS[seq_len(clev[["MAXLEV"]])*2]
    } else {
      # this case was ommited in the MATLAB function, so maybe the
      # stop can be ignored. If you reach this point please check carefully.
      stop("Unexpected clevels case")
    }
  }
  # Check that levels is not one (large) scalar.
  if (length(levels_vec) == 1) {
    levels_vec <- c(levels_vec, levels_vec)
  }
  return(list(levels = levels_vec))
}

#' Read binary Bruker NMR data
#'
#' @param file_name The file name to read
#' @param endian Passed to \code{readBin}.
#' @keywords internal
#' @noRd
read_bin_data <- function(file_name, endian) {
  tryCatch({
    con <- file(file_name, "rb")
    # if file size is 512kb, or 131072 integers
    # By using a read size of 131073 integers in one read I know if I have
    # reached the end.
    chunksize <- 131073
    num_reads = 1
    data <- readBin(con, what = "integer", n = chunksize, size = 4,
                    signed = TRUE, endian = endian)
    while (length(data) == num_reads*chunksize) {
      data <- c(data,
                readBin(con, what = "integer", n = chunksize, size = 4,
                        signed = TRUE, endian = endian))
      num_reads = num_reads + 1
    }
  }, finally = {
    close(con)
  })
  return(data)
}



guess_shape_and_submatrix_shape <- function(sample) {
  output <- list(dimension = NULL, shape = NULL, submatrix_shape = NULL)
  if (!"procs" %in% names(sample)) {
    return(output)
  }
  procs <- sample[["procs"]]
  if (!all(c("SI", "XDIM") %in% names(procs))) {
    return(output)
  }
  si_0 <- procs[["SI"]]
  xdim_0 <- procs[['XDIM']]

  if (!"proc2s" %in% names(sample)) { # 1D data
    output[["dimension"]] <- 1
    output[["shape"]] <- si_0
    output[["submatrix_shape"]] <- xdim_0
    return(output)
  }

  proc2s <- sample[["proc2s"]]
  if (!all(c("SI", "XDIM") %in% names(proc2s))) {
    return(output)
  }

  si_1 <- proc2s[["SI"]]
  xdim_1 <- proc2s[['XDIM']]


  if (!"proc3s" %in% names(sample)) { # 2D data
    output[["dimension"]] <- 2
    output[["shape"]] <- c(si_0, si_1)
    output[["submatrix_shape"]] <- c(xdim_0, xdim_1)
    return(output)
  }


  proc3s <- sample[["proc3s"]]
  if (!all(c("SI", "XDIM") %in% names(proc3s))) {
    return(output)
  }

  si_2 <- proc3s[["SI"]]
  xdim_2 <- proc3s[['XDIM']]

  if (!"proc4s" %in% names(sample)) { # 3D data
    output[["dimension"]] <- 3
    output[["shape"]] <- c(si_0, si_1, si_2)
    output[["submatrix_shape"]] <- c(xdim_0, xdim_1, xdim_2)
    return(output)
  }


  proc4s <- sample[["proc4s"]]
  if (!all(c("SI", "XDIM") %in% names(proc4s))) {
    return(output)
  }

  si_3 <- proc4s[["SI"]]
  xdim_3 <- proc4s[['XDIM']]
  # Assume 4D
  output[["dimension"]] <- 4
  output[["shape"]] <- c(si_0, si_1, si_2, si_3)
  output[["submatrix_shape"]] <- c(xdim_0, xdim_1, xdim_2, xdim_3)
  return(output)
}

#' Read orig Bruker NMR file
#' @param sample_path A character path of the sample directory
#' @return a list with name-value pairs
#' @keywords internal
#' @noRd
read_orig_file <- function(sample_path) {
  orig_file <- file.path(sample_path, "orig")
  if (!file.exists(orig_file)) {
    return(NULL)
  }
  orig_lines <- readLines(orig_file)
  lines_split <- strsplit(orig_lines, split = " ", fixed = TRUE)
  all_names <- sapply(lines_split, function(line) line[[1]])
  all_vals <- sapply(lines_split, function(line) paste0(line[2:length(line)], collapse = " "))
  output <- list()
  output[all_names] <- all_vals
  return(output)
}

parse_title_file <- function(title_lines) {
  # Two options:
  #   (a) On each line you find: "Name Value" (and optionally a " ;" in the end)
  #   (b) On each line you find  "Value" (no names)
  # What option do we have? Lets guess:
  lines_split <- strsplit(title_lines, split = " ", fixed = TRUE)
  number_of_fields_per_line <- vapply(lines_split, length, numeric(1))
  
  # Case (a): (and empty lines are accepted)
  if (all(number_of_fields_per_line >= 2)) {
    all_names <- sapply(lines_split, function(line) line[[1]])
    all_vals <- sapply(lines_split, function(line) {
      value <- paste0(line[2:length(line)], collapse = " ")
      # Remove spaces and ";" at the end of the value, if they are present:
      gsub(pattern = "[\\s;]+$", replacement = "", x = value)
    })
    output <- list()
    output[all_names] <- all_vals
    return(output)
  }
  # Case (b):
  all_names <- paste0("V", seq_len(length(title_lines)))
  output <- as.list(title_lines)
  names(output) <- all_names
  return(output)
}

#' Read pdata title Bruker NMR file
#' @param sample_path A character path of the sample directory
#' @param pdata_path Path from `sample_path` to the preprocessed data
#' @return a list with name-value pairs. If the title file has no field names,
#'         then fields are named V1, V2...
#' @keywords internal
#' @noRd
read_pdata_title_file <- function(sample_path, pdata_path = "pdata/1") {
  title_file <- file.path(sample_path, pdata_path, "title")
  if (!file.exists(title_file)) {
    return(NULL)
  }
  title_lines <- readLines(title_file, warn = FALSE)
  return(parse_title_file(title_lines))
}

#' Read processed Bruker NMR data
#'
#' @param pdata_file File name of the binary NMR data to load. Usually "1r".
#'                   If it is null it is autodetected and all files are loaded.
#' @param sample_path A character path of the sample directory
#' @param pdata_path Path from `sample_path` to the preprocessed data
#' @param all_components If `FALSE` load only the real component. Otherwise load all of them
#' @param read_pdata_title If `TRUE` also reads metadata from pdata title file.
#' @return an NMR sample
#' @keywords internal
#' @noRd
read_bruker_pdata <- function(sample_path,
                              pdata_file = NULL, pdata_path  = "pdata/1",
                              all_components = FALSE, read_pdata_title = TRUE) {
  full_pdata_path <- file.path(sample_path, pdata_path)
  if (is.null(pdata_file)) {
    # determine the dimension
    if (file.exists(file.path(full_pdata_path, "1r"))) {
      if (all_components) {
        pdata_file <- c("1r", "1i")
      } else {
        pdata_file <- "1r"
      }
    } else if (file.exists(file.path(full_pdata_path, "2rr"))) {
      if (all_components) {
        pdata_file <- c("2rr", "2ri", "2ir", "2ii")
      } else {
        pdata_file <- "2rr"
      }
    } else if (file.exists(file.path(full_pdata_path, "3rrr"))) {
      if (all_components) {
        pdata_file <- c('3rrr', '3rri', '3rir', '3rii',
                        '3irr', '3iri', '3iir', '3iii')
      } else {
        pdata_file <- "3rrr"
      }
    }
    pdata_file <- pdata_file[file.exists(file.path(full_pdata_path, pdata_file))]
  }

  full_pdata_file <- file.path(full_pdata_path, pdata_file)
  if (length(full_pdata_file) == 0) {
    stop("No pdata files found to be loaded for sample ", sample_path)
  }
  if (!all(file.exists(full_pdata_file))) {
    stop("File does not exist: ",
         paste(full_pdata_file[!file.exists(full_pdata_file)], collapse = ", "))
  }

  # Read parameters
  procs_list <- read_procs_file(full_pdata_path)

  output <- list()
  output <- c(output, procs_list)

  if (!("procs" %in% names(output))) {
    stop("Missing procs file")
  }

  # Open and read file
  if (output$procs$BYTORDP == 0) {
    endian = 'little'
  } else {
    endian = 'big'
  }

  if (isTRUE(read_pdata_title)) {
    title <- read_pdata_title_file(sample_path = sample_path, pdata_path = pdata_path)
    names(title) <- paste0("title_", names(title))
    output <- c(output, title)
  }

  for (filename in pdata_file) {
    field_name <- paste0("data_", filename)
    full_filename <- file.path(full_pdata_path, filename)
    output[[field_name]] <- read_bin_data(full_filename, endian = endian)
    output[[field_name]] <- output[[field_name]]/(2 ^ -output$procs$NC_proc)
  }

  output$levels <- read_levels(full_pdata_path, endian = endian,
                               NC_proc = output$procs$NC_proc)

  data_shapes <- guess_shape_and_submatrix_shape(output)

  dimension <- data_shapes$dimension
  if (is.null(dimension)) {
    warning("Can't determine dimension. Are procs files loaded? Assuming dimension = 1")
    dimension <- 1
  }

  output$axis <- vector("list", length = dimension)

  if (dimension >= 1) {
    # Calculate x-axis
    output$axis[[1]] <- seq(from = output$procs$OFFSET,
                            to = output$procs$OFFSET - output$procs$SW_p/output$procs$SF,
                            length.out = output$procs$SI)
  }
  if (dimension >= 2) {
    # Additional axis and reordering if 2D-file
    output$axis[[2]] <- seq(from = output$proc2s$OFFSET,
                            to = output$proc2s$OFFSET - output$proc2s$SW_p/output$proc2s$SF,
                            length.out = output$proc2s$SI)
  }
  if (dimension >= 3) {
    output$axis[[3]] <- seq(from = output$proc3s$OFFSET,
                            to = output$proc3s$OFFSET - output$proc3s$SW_p/output$proc3s$SF,
                            length.out = output$proc3s$SI)
  }
  if (dimension >= 4) {
    output$axis[[4]] <- seq(from = output$proc4s$OFFSET,
                            to = output$proc4s$OFFSET - output$proc4s$SW_p/output$proc4s$SF,
                            length.out = output$proc4s$SI)
  }

  # Reorder submatrices (se XWinNMR-manual, chapter 17.5 (95.3))
  if (data_shapes$dimension == 2) {
    # TODO: Generalize this for dimensions > 2

    SI1 = data_shapes$shape[1]
    #SI2 = data_shapes$shape[2]
    #XDIM1 = data_shapes$submatrix_shape[1]
    XDIM2 = data_shapes$submatrix_shape[2]
    NoSM_along_dimensions <- data_shapes$shape/data_shapes$submatrix_shape
    NoSM2 = NoSM_along_dimensions[length(NoSM_along_dimensions)]   # No of SM along F1
    NoSM <- cumprod(NoSM_along_dimensions) # cummulative total number of Submatrices along dimensions
    num_submatrices <- NoSM[length(NoSM)]
    for (filename in pdata_file) {
      field_name <- paste0("data_", filename)
      dim(output[[field_name]]) <- c(data_shapes$submatrix_shape, num_submatrices)
      output[[field_name]] <- aperm(output[[field_name]], c(2, 1, 3))
      dim(output[[field_name]]) <- c(XDIM2, SI1, NoSM2)
      output[[field_name]] <- aperm(output[[field_name]], c(2, 1, 3))
      dim(output[[field_name]]) <- data_shapes$shape
    }
  }
  if (data_shapes$dimension > 2) {
    stop("Can't deal with dimensions larger than two")
  }

  return(output)
}

#' Infer dimension, pulse sequence and nuclei given the acqus file information
#'
#' This function determines the dimensionality of the sample, infers the pulse
#' sequence (based on the experiment name, that can be read on the acqus file
#' -in the EXP- field) and sets the sample nuclei field.
#'
#' The dimension is based on the number of acqus, acqu2s, acqu3s existant files
#'
#' The pulse sequence assignment is based on parsing the experiment name:
#' experiment name contains -> pulse sequence assigned
#' \itemize{
#' \item{NOESY -> NOESY}
#' \item{NOEZY -> NOESY (alias)}
#' \item{CPMG -> CPMG}
#' \item{DIFF -> DIFF}
#' \item{JRES -> JRES}
#' \item{COSY -> COSY}
#' \item{TOCSY -> TOCSY}
#' \item{DIPSI -> TOCSY (alias)}
#' \item{MLEV -> TOCSY (alias)}
#' \item{HSQC -> HSQC}
#' \item{HMBC -> HMBC}
#' }
#'
#' The nuclei in the sample are set according to the pulse sequence and the
#' acqus parameters, following Sofia's advice.
#'
#' @param acqus_list The list of information retrieved by \code{read_acqus_file}
#' @return A list with dimension, pulse_sequence and nuclei
#' @keywords internal
#' @noRd
#'
infer_dim_pulse_nuclei <- function(acqus_list) {
  output <- list(dimension = NULL,
                 pulse_sequence = NULL,
                 nuclei = NULL)
  # The dimension is easy
  output$dimension <- length(acqus_list)

  # The pulse sequence is not that obvious
  experiment_name = acqus_list$acqus$EXP
  NUCLEI <- paste0("NUC", 1:8) # NUC1... NUC8 help to tell us the nuclei present

  if (grepl(pattern = "NOESY", x = experiment_name, ignore.case = TRUE) ||
      grepl(pattern = "NOEZY", x = experiment_name, ignore.case = TRUE)) {
    output$pulse_sequence <- "NOESY"
    output$nuclei <- acqus_list$acqus[["NUC1"]]

  } else if (grepl(pattern = "CPMG", x = experiment_name, ignore.case = TRUE)) {
    output$pulse_sequence <- "CPMG"
    output$nuclei <- acqus_list$acqus[["NUC1"]]

  } else if (grepl(pattern = "DIFF", x = experiment_name, ignore.case = TRUE)) {
    output$pulse_sequence <- "DIFFUSION"
    output$nuclei <- acqus_list$acqus[["NUC1"]]

  } else if (grepl(pattern = "JRES", x = experiment_name, ignore.case = TRUE)) {
    output$pulse_sequence <- "JRES"
    output$nuclei <- acqus_list$acqus[["NUC1"]]

  } else if (grepl(pattern = "COSY", x = experiment_name, ignore.case = TRUE)) {
    output$pulse_sequence <- "COSY"
    output$nuclei <- paste0(acqus_list$acqus[["NUC1"]], "-", acqus_list$acqu2s[["NUC1"]])

  } else if (grepl(pattern = "TOCSY", x = experiment_name, ignore.case = TRUE) ||
             grepl(pattern = "MLEV", x = experiment_name, ignore.case = TRUE) ||
             grepl(pattern = "DIPSI", x = experiment_name, ignore.case = TRUE)) {
    output$pulse_sequence <- "TOCSY"
    output$nuclei <- paste0(acqus_list$acqus[["NUC1"]], "-", acqus_list$acqu2s[["NUC1"]])

  } else if (grepl(pattern = "HSQC", x = experiment_name, ignore.case = TRUE)) {
    output$pulse_sequence <- "HSQC"
    output$nuclei <- paste(acqus_list$acqus[NUCLEI][acqus_list$acqus[NUCLEI] != "off"], collapse = "-")

  } else if (grepl(pattern = "HMBC", x = experiment_name, ignore.case = TRUE)) {
    output$pulse_sequence <- "HMBC"
    output$nuclei <- paste(acqus_list$acqus[NUCLEI][acqus_list$acqus[NUCLEI] != "off"], collapse = "-")

  } else {
    warning("infer_dim_pulse_nuclei: Unknown Pulse Sequence in acqus$EXP field.\n",
            "Please add the pulse sequence '",experiment_name, "' to infer_dim_pulse_nuclei in R/bruker.R\n")
  }

  return(output)
}

#' Read Bruker NMR metadata
#'
#' @param sample_path A character path of the sample directory
#' @param pdata_path Optional character path of the processed data (to read pdata title file)
#' @param read_pdata_title logical. If \code{TRUE} reads the \code{pdata/1/title}
#'                         file if it exists.
#' @return an NMR sample
#' @export
read_bruker_metadata <- function(sample_path, pdata_path = "pdata/1",
                                 read_pdata_title = TRUE) {
  # Read parameters
  acqus_list <- read_acqus_file(sample_path)

  if (length(acqus_list) == 0) {
    stop("No acqus file found in ", sample_path)
  }
  info <- list()
  info$NMRExperiment <- basename(sample_path)
  info$file_format <- "Zipped Bruker NMR directory"
  info$sample_path <- sample_path
  info$import_time <- Sys.time()
  info <- c(info, infer_dim_pulse_nuclei(acqus_list))

  orig <- read_orig_file(sample_path = sample_path)
  # Store acquisition time information
  info$acq_datetime <- as.POSIXct(acqus_list$acqus$Stamp[1])



  output <- list(info = info,
                 orig = orig)

  if (isTRUE(read_pdata_title)) {
    pdata_1_title <- read_pdata_title_file(sample_path = sample_path,
                                           pdata_path = pdata_path)
    output$title <- pdata_1_title
  }

  output <- c(output, acqus_list)

  if (!("acqus" %in% names(output))) {
    stop("Missing acqus file")
  }


  if (!any(c("fid", "ser") %in% list.files(sample_path))) {
    stop("No raw data available in ", sample_path)
  }
  return(output)
}

#' Read a Bruker sample directory
#'
#' @param pdata_file File name of the binary NMR data to load. Usually "1r".
#'                   If it is null it is autodetected and all files are loaded.
#' @param sample_path A character path of the sample directory
#' @param pdata_path Path from `sample_path` to the preprocessed data
#' @param all_components If `FALSE` load only the real component. Otherwise load all of them
#' @param read_pdata_title If `TRUE` also reads metadata from pdata title file.
#' @return a list with all the bruker sample information
#' @export
read_bruker_sample <- function(sample_path,
                               pdata_file = NULL, pdata_path  = "pdata/1",
                               all_components = FALSE) {
  # This is the equivalent to the rbnmr function. It has been splitted into
  # two parts, the one reading the sample metadata and the one reading the
  # processed data (pdata) directory. The idea is that we may want to load
  # the metadata only, which should be faster so it is convenient to have a
  # simple read_bruker_metadata function

  meta <- read_bruker_metadata(sample_path, pdata_path)

  pdata <- read_bruker_pdata(sample_path = sample_path,
                             pdata_file = pdata_file,
                             pdata_path  = pdata_path,
                             all_components = all_components,
                             read_pdata_title = FALSE) # pdata title already read

  output <- bruker_merge_meta_pdata(meta, pdata)
  return(output)
}

bruker_merge_meta_pdata <- function(meta, pdata) {
  # This function concatenates the metadata and the processed data
  # This is a trivial function, but I prefer to define it because if something
  # more complex needs to be done in the merging in the future we can make sure
  # that any function trying to merge metadata and processed data in the future
  # will continue working if it calls this.
  return(c(meta, pdata))
}
