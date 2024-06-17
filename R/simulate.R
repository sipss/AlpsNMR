split_multiplet_peaks <- function(multiplet_type) {
    splits <- strsplit(multiplet_type, split = " of ", fixed = TRUE)[[1]]
    splits <- purrr::map(splits, function(x) {
        if (grepl("^[sdtq]+$", x)) {
            strsplit(x, split="")[[1]]
        } else {
            x
        }
    })
    purrr::list_c(splits)
}

multiplet_name_to_peak_number <- function(multiplet) {
    peak_types <- c(
        "s" = 1L,
        "singlet" = 1L,
        "singlets" = 1L,
        "d" = 2L,
        "doublet" = 2L,
        "doublets" = 2L,
        "t" = 3L, 
        "triplet" = 3L,
        "triplets" = 3L,
        "q" = 4L,
        "quartet" = 4L,
        "quartets" = 4L,
        "quadruplet" = 4L,
        "quadruplets" = 4L,
        "quintet" = 5L,
        "quintets" = 5L,
        "quintuplet" = 5L,
        "quintuplets" = 5L,
        "sextet" = 6L,
        "sextets" = 6L,
        "sextuplet" = 6L,
        "sextuplets" = 6L,
        "septuplet" = 7L,
        "septuplets" = 7L,
        "octuplet" = 8L,
        "octuplets" = 8L,
        "nonetuplet" = 9L,
        "nonetuplets" = 9L
    )
    peak_nums <- peak_types[multiplet]
    as_numbers <- grepl("^[0-9]+$", multiplet)
    peak_nums[as_numbers] <- as.integer(multiplet[as_numbers])
    if (anyNA(peak_nums)) {
        cli::cli_abort("Could not parse {multiplet}. Best effort: {peak_nums}")
    }
    peak_nums
}

#' Parse a multiplet structure returning the number of peaks to generate
#' @param multiplet_types A character vector describing the multiplets, see examples
#' @return A vector of the same length as `multiplet_types` with the cardinality of the peaks in the multiplet
#' @examples
#' multiplet_to_cardinality(c("ddd", "doublet of 3", "singlet"))
#' 
#' @export
multiplet_to_cardinality <- function(multiplet_types) {
    multiplet_types |>
        purrr::set_names() |>
        purrr::map(
            function(multiplet_type) {
                multiplet_type |>
                    split_multiplet_peaks() |>
                    purrr::map(multiplet_name_to_peak_number) |>
                    purrr::list_c(ptype=integer(1)) |>
                    purrr::set_names(NULL)
            }
        )
}

pascal_triangle_row_n <- function(n) {
    out <- 1L
    if (n == 1) {
        return(out)
    }
    prev <- 1
    for (i in seq_len(n)) {
        curr <- prev * (n -i + 1) / i
        out <- c(out, curr)
        prev <- curr
    }
    out
}


multiplet_to_amplitudes <- function(peak_cardinality) {
    peaks_to_convolve <- purrr::map(peak_cardinality, pascal_triangle_row_n)
    peaks_to_convolve
}
