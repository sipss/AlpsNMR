# These functions deal with importing JDX files created with:
# ACD/Spectrus 2014 v 14.04
#
# Author: Sergio Oller


#' Remove comments from JDX files
#'
#' JDX comments start with \code{$$}.
#'
#' @param lines A character vector with the lines of the jdx file
#' @return A character vector where the lines have been removed
#' @keywords internal
#' @noRd
#'
strip_comments <- function(lines) {
    gsub("(.*)?\\$\\$.*", "\\1", lines)
}


DATA_FIELDS <- c("XYDATA", "DATA TABLE", "PEAK TABLE")
DATA_RELATED_FIELDS <-
    c(
        "FIRSTX",
        "LASTX",
        "MAXX",
        "MINX",
        "MAXY",
        "MINY",
        "XFACTOR",
        "YFACTOR",
        "FIRSTY",
        "LASTY",
        "DELTAX",
        "DELTAY",
        "NPOINTS"
    )

#' Process a JDX block from a JDX file
#'
#' This is a very limited and partial implementation of the JDX specification,
#' in order to support ACD/Spectrus JCAMP-DX files.
#'
#' JDX blocks should start with \code{##TITLE} and end with \code{##END}. Blocks
#' can be nested if they are of the LINK type.
#'
#' One data field is allowed per block.
#' @param lines A character vector with the JDX file.
#' @param metadata_only ignore data (do not parse it)
#' @return A list with two elements: the new block and the number of lines read from lines
#' @keywords internal
#' @noRd
#'
process_block <- function(lines, metadata_only = FALSE) {
    # Only one DATA_FIELD per block.
    
    data_type <- "none" # we haven't found any data field yet here
    data_field <-
        NULL # we don't know the name of the data_field yet
    block <- list()    # The block that we will return:
    
    # All blocks must start with a TITLE:
    if (!grepl("^\\s*##TITLE\\s*=.*", lines[1])) {
        stop("Block should start with ##TITLE")
    }
    # Get the title:
    title <-
        sub(pattern = "^\\s*##TITLE\\s*=\\s*(.*?)\\s*$", replacement = "\\1", lines[1])
    block$TITLE <- title
    # Scan each line and parse it:
    i <- 1 # Line 1 was the title
    while (i < length(lines)) {
        i <- i + 1
        line <- lines[i]
        if (line == "") {
            # Blank lines are ignored
            next
        }
        # Try to match a ##Name = Value field
        name_val <-
            stringr::str_match(string = line, pattern = "^\\s*##(.*?)\\s*=\\s*(.*?)\\s*$")
        if (is.na(name_val[1, 2])) {
            # We can't, we skip the line
            next
        }
        # The name and the value
        field_name <- name_val[1, 2]
        field_value <- name_val[1, 3]
        # If field_name is TITLE we have a new block
        # (The title of this block has been parsed already)
        if (field_name == "TITLE") {
            # Recursively process the block and append it to the list of subblocks:
            new_block <-
                process_block(lines[i:length(lines)], metadata_only = metadata_only)
            block[["blocks"]] <-
                c(block[["blocks"]], list(new_block$block))
            i <- i + new_block$num_lines - 1
            next
        }
        # If field_name is END then we finish reading new information
        if (field_name == "END") {
            break
        }
        # DATA fields expand multiple lines so we must deal with them:
        if (field_name %in% DATA_FIELDS) {
            # But only one DATA FIELD is allowed per block:
            if (data_type != "none") {
                stop("We already have a data field in this block")
            }
            # Now we know the name of the data_field
            data_field <- field_name
            
            # We need to handle multiple data formats. We just implement the ones we need
            #################### (XY..XY) format ####################
            type_of_table <-
                stringr::str_match(string = field_value,
                                   pattern = "^\\s*\\(XY\\.\\.XY\\)\\s*(.*)\\s*")
            if (!is.na(type_of_table[1, 1])) {
                # "(XY..XY)"
                data_type <- "(XY..XY)"
                if (type_of_table[1, 2] != "") {
                    # I do not know how to treat this (probably is the first value of the array)
                    stop("Unknown way to handle value from:", line)
                }
                i_data_start <- i + 1
                while (i < length(lines)) {
                    i <- i + 1
                    line <- lines[i]
                    # Read lines until we find a new field
                    if (startsWith(line, "##")) {
                        i <- i - 1
                        break
                    }
                }
                # Convert the lines to a data frame
                if (!metadata_only) {
                    data <- strsplit(lines[i_data_start:i], "[,;[:blank:]]+")
                    data <- unlist(data)
                    data <- matrix(as.numeric(data), nrow = 2)
                    block[[field_name]] <-
                        data.frame(
                            x = data[1,],
                            y = data[2,],
                            stringsAsFactors = FALSE
                        )
                }
                next
            }
            
            ######################### (X++(Y..Y)) format ##########################
            type_of_table <-
                stringr::str_match(string = field_value,
                                   pattern = "^\\s*\\(X\\+\\+\\(Y\\.\\.Y\\)\\)\\s*(.*)\\s*")
            if (!is.na(type_of_table[1, 1])) {
                # "(X++(Y..Y))"
                data_type <- "(X++(Y..Y))"
                if (type_of_table[1, 2] != "") {
                    # I do not know how to treat this (probably is the first value of the array)
                    stop("Unknown way to handle value from:", line)
                }
                i_data_start <- i + 1
                while (i < length(lines)) {
                    i <- i + 1
                    line <- lines[i]
                    if (startsWith(line, "##")) {
                        i <- i - 1
                        break
                    }
                }
                if (!metadata_only) {
                    # Remove meaningless leading spaces (lead to extra NA)
                    cleaned_spaces <-
                        stringr::str_trim(lines[i_data_start:i])
                    # Then make sure there is a space between numbers (it's apparently optional):
                    cleaned_spaces <-
                        gsub("([0-9.])([+-])",
                             "\\1 \\2",
                             cleaned_spaces)
                    # Ignore x (first column). Will be filled later
                    y <-
                        unlist(lapply(strsplit(cleaned_spaces, split = "[[:blank:]]+"),
                                      function(x)
                                          as.numeric(x[2:length(x)])))
                    block[[field_name]] <-
                        data.frame(
                            x = NA,
                            y = y,
                            stringsAsFactors = FALSE
                        )
                }
                next
            }
            # Unknown data format:
            stop("Error: Format not implemented")
        } # The field was not multiline
        # Default field-value pair
        block[[field_name]] <- field_value
    }
    # Convert fields:
    fields_to_convert <- setdiff(names(block), DATA_FIELDS)
    block[fields_to_convert] <-
        lapply(block[fields_to_convert], convert_field)
    # Now that we have converted the fields we can create the X axis
    # for the (X++(Y..Y)) data format
    if (data_type == "(X++(Y..Y))") {
        x <-
            seq(from = block[["FIRSTX"]],
                to = block[["LASTX"]],
                length.out = block[["NPOINTS"]])
        block[[data_field]][["x"]] <- x
    }
    
    # Apply X and Y Factors to the data:
    if (data_type != "none") {
        if ("x" %in% names(block[[data_field]]) &&
            "XFACTOR" %in% names(block)) {
            block[[data_field]][["x"]] <-
                block[[data_field]][["x"]] * block[["XFACTOR"]]
        }
        if ("y" %in% names(block[[data_field]]) &&
            "YFACTOR" %in% names(block)) {
            block[[data_field]][["y"]] <-
                block[[data_field]][["y"]] * block[["YFACTOR"]]
        }
        # Convert x axis from frequency to chemshift
        if (tolower(block[["XUNITS"]]) == "hz") {
            block[[data_field]][["x"]] <-
                block[[data_field]][["x"]] / block[[".OBSERVE FREQUENCY"]]
        }
    }
    return(list(block = block, num_lines = i))
}

#' Read a ACD/Spectrus JDX file
#' @param file_names A character vector with JDX file names
#' @param metadata_only A logical to load only metadata and skip the actual
#'                                            spectrum. (default: \code{FALSE})
#' @return A list with the JDX information
#' @keywords internal
#' @noRd
read_jdx <- function(file_names, metadata_only = FALSE) {
    warn_future_to_biocparallel()
    output <- BiocParallel::bplapply(
        X = file_names,
        FUN = function(file_name) {
            lines <- readLines(file_name)
            lines <- strip_comments(lines)
            sampl <- process_block(lines, metadata_only = metadata_only)$block
            info <- list(file_format = "JDX-JCAMP NMR sample")
            sampl$info <- info
            sampl
        }
    )
    return(output)
}

create_df_from_jdx_sample <- function(sampl,
                                      exclude = c("blocks", DATA_FIELDS,
                                                  DATA_RELATED_FIELDS)) {
    global_names <- setdiff(names(sampl), exclude)
    global_df <- do.call(tibble::tibble, sampl[global_names])
    
    # block:
    all_df <- global_df
    for (block_idx in seq_along(sampl[["blocks"]])) {
        block <- sampl[["blocks"]][[block_idx]]
        block_names <- setdiff(names(block), exclude)
        block_df <- do.call(tibble::tibble, block[block_names])
        colnames(block_df) <-
            paste0("b", block_idx, "_", colnames(block_df))
        all_df <- cbind(all_df, block_df)
    }
    return(all_df)
}
