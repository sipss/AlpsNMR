#' Build a PCA on for an nmr_dataset
#' 
#' This function builds a PCA model with all the NMR spectra. Regions with
#' zero values (excluded regions) or near-zero variance regions are automatically
#' excluded from the analysis.
#' 
#' @param nmr_data a \code{\link{nmr_dataset}} object
#' @inheritParams mixOmics::pca
#' @param ... Additional arguments passed on to [mixOmics::pca]
#' 
#' @return
#' A PCA model as given by [mixOmics::pca] with two additional attributes:
#'  - `nmr_data_axis` containing the full ppm axis
#'  - `nmr_included` with the data points included in the model
#' These attributes are used internally by NIHSnmr to create loading plots
#' 
#' @export
nmr_pca_build_model <- function(nmr_data, ncomp = NULL, center = TRUE, scale = FALSE, ...) {
  UseMethod("nmr_pca_build_model")
}

#' @rdname nmr_pca_build_model
#' @export
nmr_pca_build_model.nmr_dataset <- function(nmr_data, ncomp = NULL, center = TRUE, scale = FALSE, ...) {
  zero_var_cols <- mixOmics::nearZeroVar(nmr_data$data_1r) # excluded
  data_1r <- nmr_data$data_1r[,-zero_var_cols$Position]
  rownames(data_1r) <- nmr_meta_get(nmr_data, "NMRExperiment")$NMRExperiment
  pca_model <- mixOmics::pca(X = data_1r, ncomp = ncomp, center = center, scale = scale, ...)
  # These attributes are used by nmr_pca_loadingplot:
  attr(pca_model, "nmr_data_axis") <- nmr_data$axis[[1]]
  attr(pca_model, "nmr_included") <- base::setdiff(seq_along(nmr_data$axis[[1]]), zero_var_cols$Position)
  pca_model
}

#' @rdname nmr_pca_build_model
#' @export
nmr_pca_build_model.nmr_dataset_1D <- function(nmr_data, ncomp = NULL, center = TRUE, scale = FALSE, ...) {
  data_1r <- nmr_data$data_1r
  rownames(data_1r) <- nmr_meta_get(nmr_data, "NMRExperiment")$NMRExperiment
  pca_model <- mixOmics::pca(X = data_1r, ncomp = ncomp, center = center, scale = scale, ...)
  # These attributes are used by nmr_pca_loadingplot:
  attr(pca_model, "nmr_data_axis") <- nmr_data$axis
  attr(pca_model, "nmr_included") <- seq_along(nmr_data$axis)
  pca_model
}


#' Plotting functions for PCA
#' 
#' @param nmr_data a \code{\link{nmr_dataset}} object
#' @param pca_model A PCA model trained with [nmr_pca_build_model]
#' @param comp Components to represent
#' @param ... Additional aesthetics passed on to [ggplot2::aes] (use bare unquoted names)
#' 
#' @name nmr_pca_plots
#' 
NULL

#' @rdname nmr_pca_plots
#' @export
nmr_pca_plot_variance <- function(pca_model) {
  cum_var_percent <- 100*cumsum(pca_model$sdev^2/pca_model$var.tot)
  ggplot2::qplot(x = seq_along(cum_var_percent), y = cum_var_percent, geom = "line") +
    ggplot2::xlab("Number of Principal Components") +
    ggplot2::ylab("Cummulated Explained Variance (%)")
}

#' @rdname nmr_pca_plots
#' @export
nmr_pca_scoreplot <- function(nmr_data, pca_model, comp = 1:2, ...) {
  nmr_metadata <- nmr_meta_get(nmr_data)
  scores <- tibble::as_tibble(pca_model$x, rownames = "NMRExperiment") %>%
    dplyr::left_join(nmr_metadata, by = "NMRExperiment")
  var_percent <- 100*pca_model$sdev^2/pca_model$var.tot
  axis_labels <- paste0("PC", seq_along(var_percent), " (", round(var_percent, 2),"%)")
  names(axis_labels) <- paste0("PC", seq_along(var_percent))
  if (length(comp) == 2) {
    xy <- paste0("PC", comp)
    gplt <- ggplot2::ggplot(
      data = scores, 
      mapping = ggplot2::aes(x = !!rlang::sym(xy[1]), y = !! rlang::sym(xy[2]), ...)
    ) +
      ggplot2::geom_point() +
      ggplot2::xlab(axis_labels[xy[1]]) +
      ggplot2::ylab(axis_labels[xy[2]])
    } else {
      gplt <- GGally::ggpairs(data = scores, mapping = ggplot2::aes(...),
                              columns = comp + 1, # +1 because the NMRExperiment becomes the first column
                              progress = FALSE,
                              labeller = ggplot2::as_labeller(axis_labels))
  }
  gplt
}

#' @rdname nmr_pca_plots
#' @export
nmr_pca_loadingplot <- function(pca_model, comp) {
  ppm_axis <- attr(pca_model, "nmr_data_axis")
  loadings <- matrix(0, nrow = length(ppm_axis), ncol = length(comp))
  loadings[attr(pca_model, "nmr_included"),] <- pca_model$rotation[, comp, drop = FALSE]
  # loadings[,1] # first loading
  loadings <- as.data.frame(loadings)
  loadings$ppm <- ppm_axis
  loadings_long <- reshape2::melt(loadings, id.vars = "ppm", variable.name = "component", value.name = "loading")
  ggplot2::ggplot(
    loadings_long,
    ggplot2::aes_string(x = "ppm", y = "loading", group = "component")
  ) +
    ggplot2::geom_line() +
    ggplot2::scale_x_reverse()
}

