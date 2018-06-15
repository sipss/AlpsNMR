# This script has utilities not related to NMR.

#' @importFrom magrittr %>%
#' @export
magrittr::`%>%`

#' Split data frame into constant and non-constant columns
#'
#' Columns that do not changes across observations are moved to the \code{const}
#' data frame. The rest are returned in the \code{diff} data frame.
#'
#' The \code{const} data frame only has one row as all rows would be identical.
#'
#' @param data a data frame
#' @keywords internal
simplify_df <- function(data) {
  if (is.null(data)) {
    return(NULL)
  }
  # equal_cols is a logical array of length equal to the columns of data, with
  # TRUE if the column is constant
  # It is more complex than expected because of R "list columns" (storing a list
  # inside a dataframe column)
  equal_cols <- sapply(seq_len(ncol(data)),
                       function(col) {
                         are_na <- is.na(data[[col]])
                         if (all(are_na)) {
                           # They are all equal
                           return(TRUE)
                         }
                         if (any(are_na)) {
                           # Some are NA, but not all of them (already considered)
                           # so they are not equal
                           return(FALSE)
                         }
                         # There are no NA, so let's check if they are all equal or not
                         all(sapply(seq_len(nrow(data)),
                                    function(row) {
                                      identical(data[[row, col]], data[[1, col]])
                                    }))
                       })
  data_diff <- data[,!equal_cols, drop = FALSE]
  data_const <- data[1,equal_cols, drop = FALSE]
  return(list(const = data_const, diff = data_diff))
}

#' Convert a list of lists to a tibble
#'
#' Each element of \code{ls} corresponds to 1 row of the output tibble. This
#' element is a list, and each subelement corresponds to a column of the
#' output tibble.
#'
#' One nice thing about this function is that if there are elements with missing
#' columns, those values will be filled with \code{NA}.
#'
#' The other nice thing of this function is that if any of those elements is
#' itself a list or a vector, it will be nicely preserved inside the tibble cell.
#'
#' @param ls A list
#' @return A tibble
#' @keywords internal
#' @examples
#' \dontrun{
#' data_as_list <- list(list(Gender = "Male", Height = 170),
#'                      list(Gender = "Female", Height = 160),
#'                      list(Gender = "Male", Weight = 80))
#' mydata <- list_of_lists_to_tibble(data_as_list)
#' }
list_of_lists_to_tibble <- function(ls) {
  # Get all the column names per sample:
  all_columns_per_sample <- lapply(ls, function(x) names(x))
  # all_columns_per_sample
  #  list(c("Gender", "Height), c("Gender", "Height"), c("Gender", "Weight"))
  # Concatenate and get a unique list of all column names:
  all_names <- unique(do.call(c, all_columns_per_sample))
  # all_names
  #  c("Gender", "Height", "Weight")

  # For each sample:
  ls_ordered <- lapply(ls,
                       function(sampl) {
                         # If the sample is NULL, then return a list of
                         # missing values with the column names:
                         if (is.null(sampl)) {
                           sampl <- as.list(rep(NA, length(all_names)))
                           names(sampl) <- all_names
                           return(sampl)
                         }
                         # Otherwise reorder list elements so they match the names
                         sampl <- sampl[all_names]
                         # If the sample does not have all the columns,
                         # the missing columns are named `NA` with value NULL.
                         # Make sure the names are properly set:
                         names(sampl) <- all_names
                         # Finally replace NULL values in the list with NA:
                         sampl <- lapply(sampl, function(value) {
                           if (is.null(value)) {
                             return(NA)
                           } else {
                             return(value)
                           }
                         })
                         # Return the sample:
                         return(sampl)
                       })

  data <- tibble::as_tibble(do.call(rbind, ls_ordered))
  # dataframe with lists on each column. Let's convert them if possible
  data <- tibble_lists_columns_to_vector_columns(data)
  return(data)
}

#' Simplifies a tibble with lists columns of length 1 of the same type
#'
#'
#' @param data a tibble
#' @return a tibble with converted columns
#' @keywords internal
tibble_lists_columns_to_vector_columns <- function(data) {
  # based on http://stackoverflow.com/questions/40046603/tibble-with-list-columns-convert-to-array-if-possible/

  ### Step 1: Find which columns have to be converted:

  # 1.1 Convert only columns of type "list"
  to_simplify_cols <- which(vapply(data,
                                   FUN = function(x) "list" %in% class(x),
                                   FUN.VALUE = logical(1)))

  # 1.1b If there are none, return:
  if (length(to_simplify_cols) == 0) {
    # Restore original colnames
    return(data)
  }

  # 1.2 Convert only list columns that have all elements of length 1:

  # Max length of the lists columns:
  length_column_elements <- apply(data[,to_simplify_cols], 2,
                                  function(x) max(vapply(x, length, numeric(1))))
  # We just simplify list columns of length 1
  to_simplify_cols <- to_simplify_cols[length_column_elements == 1]

  # 1.2.b No list columns can be simplified:
  if (length(to_simplify_cols) == 0) {
    return(data)
  }

  # 1.3 Convert only list columns that have all elements of length 1 and belong
  #     to the same class (allowing for NA values)

  # For each list column of length 1:
  types <- apply(data[,to_simplify_cols], 2,
                 function(data_column) {
                   # For each value in this column get the class, missing values are given
                   # their own class because they are always allowed
                   value_classes <- lapply(data_column,
                                           function(value) {
                                             if (is.na(value)) {
                                               "__NAVALUE__"
                                             } else {
                                               class(value)
                                             }
                                           })
                   # Remove repeated types in the column:
                   value_classes <- unique(value_classes)
                   # Remove the missing value placeholders:
                   idx <- vapply(value_classes,
                                 function(value_class) !identical(value_class, "__NAVALUE__"),
                                 logical(1))
                   value_classes <- value_classes[idx]
                   if (length(value_classes) == 0) {
                     value_classes <- list("logical")
                   }
                   value_classes
                 })
  # types is a list of the same length than to_simplify_cols.
  # types[[1]] is a list corresponding to the column to_simplify_cols[1]
  # types[[1]] contains the classes of the column to_simplify_cols[1]
  # We simplify only columns with just one class, so we check which types have length 1

  number_of_types <- vapply(types, length, numeric(1))
  # filter columns of a single class
  number_of_types <- number_of_types[number_of_types == 1]
  # Keep those columns only
  to_simplify_cols <- to_simplify_cols[names(number_of_types)]

  # No list columns can be simplified:
  if (length(to_simplify_cols) == 0) {
    return(data)
  }

  # Get all column names:
  data_col_names <- colnames(data)
  # Get the column names of the columns to simplify
  to_simplify <- data_col_names[to_simplify_cols]

  # Do the conversion
  data2 <- tidyr::unnest_(data, unnest_cols = lapply(to_simplify, as.name))
  data2 <- data2[, colnames(data)] # Preserve original column order
  return(data2)
}


# Returns a logical that determines if there is a need to show a progress bar
# More than 5 iterations, interactive and not knitting
show_progress_bar <- function(num_iterations, num_iterations_threshold = 5) {
  return(num_iterations > num_iterations_threshold &&
           interactive() &&
           is.null(getOption("knitr.in.progress")))
}
