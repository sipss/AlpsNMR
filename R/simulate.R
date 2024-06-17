split_multiplet_orders <- function(multiplet_type) {
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
    unname(peak_nums)
}

#' Parse a multiplet structure returning the number of peaks to generate
#' @param multiplet_types A character vector describing the multiplets, see examples
#' @return A vector of the same length as `multiplet_types` with the cardinality of the peaks in the multiplet
#' @examples
#' multiplet_to_cardinality("ddd")
#' multiplet_to_cardinality("doublet of 3")
#' multiplet_to_cardinality("singlet")
#' 
#' @noRd
multiplet_to_cardinality <- function(multiplet_type) {
  multiplet_orders_str <- split_multiplet_orders(multiplet_type)
  multiplet_orders_int <- purrr::map(
    multiplet_orders_str,
    multiplet_name_to_peak_number
  )
  purrr::list_c(
    multiplet_orders_int,
    ptype=integer(1)
  )
}

multiplet_to_peaks <- function(multiplet_type, multiplet_center, coupling_constants) {
    peak_cardinality <- multiplet_to_cardinality(multiplet_type)
    peaks_to_convolve <- purrr::map(
      peak_cardinality,
      function(plet) {
        choose(plet-1, seq(0, plet-1))
      }
    )
    peak_amplitudes <- purrr::reduce(
      peaks_to_convolve,
      function(x, y) {
        as.numeric(t(outer(x, y)))
      }
    )
    peak_deltas <- purrr::map2(
      peak_cardinality,
      coupling_constants,
      function(plet, J) {
        n <- (plet-1)/2
        J*seq(from = -n, to = n, length.out=plet)
      }
    )
    peak_positions <- purrr::reduce(
      peak_deltas,
      function(x, y) {
        as.numeric(t(outer(x, y, `+`)))
      },
      .init = multiplet_center
    )
    data.frame(
      position = peak_positions,
      amplitude = peak_amplitudes
    )
}
