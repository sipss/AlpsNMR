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
    # "dt" to c(2, 3)
    peak_cardinality <- multiplet_to_cardinality(multiplet_type)
    names(peak_cardinality) <- paste0("Order", seq_along(peak_cardinality))
    
    # c(2, 3) to list(c(1,1), c(1,2,1))
    peaks_to_convolve <- purrr::map(
      peak_cardinality,
      function(plet) {
        choose(plet-1, seq(0, plet-1))
      }
    )
    # c(2,3) to list(c(1,2), c(1,2,3))
    peak_src <- purrr::map(
      peaks_to_convolve,
      seq_len
    )

    # data.frame with one row per peak and one column per order
    peak_src_df <- tidyr::expand_grid(!!!peak_src)
    # Add amplitude contribution from each order:
    peak_src_df <- mutate(
      peak_src_df,
      across(
        starts_with("Order"), 
        function(col) peaks_to_convolve[[cur_column()]][col],
        .names = "Amp_{.col}"
      )
    )
    # Multiply amplitude contributions
    to_mult <- grepl("^Amp_", colnames(peak_src_df))
    peak_src_df$amplitude <- do.call(prod, peak_src_df[to_mult])

    # Compute shifts based on coupling constants
    peak_deltas <- purrr::map2(
      peak_cardinality,
      coupling_constants,
      function(plet, J) {
        n <- (plet-1)/2
        J*seq(from = -n, to = n, length.out=plet)
      }
    )
    
    peak_src_df <- mutate(
      peak_src_df,
      across(
        starts_with("Order"), 
        function(col) peak_deltas[[cur_column()]][col],
        .names = "Delta_{.col}"
      )
    )

    # Sum delta contributions
    to_sum <- grepl("^Delta_", colnames(peak_src_df))
    peak_src_df$Delta <- do.call(sum, peak_src_df[to_sum])
    peak_src_df$position <- peak_src_df$Delta + multiplet_center
    peak_amplitudes <- purrr::reduce(
      peaks_to_convolve,
      function(x, y) {
        as.numeric(t(outer(x, y)))
      }
    )
    peak_positions <- purrr::reduce(
      peak_deltas,
      function(x, y) {
        as.numeric(t(outer(x, y, `+`)))
      },
      .init = multiplet_center
    )
    list(
         data.frame(
           position = peak_positions,
           amplitude = peak_amplitudes
         ),
        dplyr::select(peak_src_df, position, amplitudes, everything())
    )
}
