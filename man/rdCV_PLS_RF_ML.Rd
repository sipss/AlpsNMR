% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Classifier_MUVR.R
\name{rdCV_PLS_RF_ML}
\alias{rdCV_PLS_RF_ML}
\title{Feature selection and validation in MULTILEVEL analysis}
\usage{
rdCV_PLS_RF_ML(nmr_peak_table, label, scale = TRUE, nRep = 10,
  nOuter = 5, nInner, varRatio = 0.75, DA = FALSE,
  fitness = "MISS", method = "PLS", ML = TRUE, modReturn = FALSE,
  logg = FALSE, parallel = TRUE)
}
\arguments{
\item{nmr_peak_table}{an AlpsNMR integration object (2 classes)}

\item{label}{the name of the variable to test (e.g. "Timepoint")}

\item{scale}{If TRUE, the predictor variable matrix is scaled to unit variance for PLS modelling.}

\item{nRep}{Number of repetitions of double CV. (Defaults to 5)}

\item{nOuter}{Number of outer CV loop segments. (Defaults to 6)}

\item{nInner}{Number of inner CV loop segments. (Defaults to nOuter-1)}

\item{varRatio}{Ratio of variables to include in subsequent inner loop iteration. (Defaults to 0.75)}

\item{DA}{Boolean for Classification (discriminant analysis) (By default, if Y is numeric -> DA=FALSE. If Y is factor (or character) -> DA=TRUE)}

\item{fitness}{Fitness function for model tuning (choose either 'AUROC' or 'MISS' (default) for classification; or 'RMSEP' (default) for regression.)}

\item{method}{Multivariate method. Supports 'PLS' and 'RF' (default)}

\item{ML}{Boolean for multilevel analysis (defaults to FALSE)}

\item{modReturn}{Boolean for returning outer segment models (defaults to FALSE). Setting modReturn=TRUE is required for making MUVR predictions using predMV().}

\item{logg}{Boolean for whether to sink model progressions to `log.txt`}

\item{parallel}{Boolean for whether to perform `foreach` parallel processing (Requires a registered parallel backend; Defaults to `TRUE`)}
}
\value{
a MUVR model containing selection parameters, validation and fitness
}
\description{
Statistical analysis and feature selection in a repeated double
cross-validation frame based on the partial least squares
(PLS) or random forest (RF) analyses using an algorithm
for multivariate modelling with minimally biased variable
selection (MUVR) from the \code{MUVR} package. The function
\code{rdCV_PLS_RF_ML} allows the multilevel comparison,
especially useful in crossover or longitudinal studies
(2 timepoints) considering the same individual (it
requires 2 samples of the same observation).
}
\examples{
\dontrun{
model = rdCV_PLS_RF_ML(nmr_peak_table, label = "Timepoint", ML = TRUE)
MUVR_model_plot(model)
}
}
\references{
Shi,L. et al. (2018) Variable selection and validation in multivariate modelling. Bioinformatics.
}
\seealso{
Other nmr_dataset_1D functions: \code{\link{Automatic_targeted_profiling}},
  \code{\link{[.nmr_dataset_1D}},
  \code{\link{computes_peak_width_ppm}},
  \code{\link{files_to_rDolphin_blood}},
  \code{\link{files_to_rDolphin_cell}},
  \code{\link{files_to_rDolphin_urine}},
  \code{\link{get_integration_with_metadata}},
  \code{\link{is.nmr_dataset_1D}},
  \code{\link{load_and_save_functions}},
  \code{\link{new_nmr_dataset_1D}},
  \code{\link{nmr_align_find_ref}},
  \code{\link{nmr_align}},
  \code{\link{nmr_baseline_removal}},
  \code{\link{nmr_baseline_threshold}},
  \code{\link{nmr_detect_peaks_plot}},
  \code{\link{nmr_detect_peaks_tune_snr}},
  \code{\link{nmr_detect_peaks}},
  \code{\link{nmr_exclude_region}},
  \code{\link{nmr_integrate_peak_positions}},
  \code{\link{nmr_integrate_regions}},
  \code{\link{nmr_interpolate_1D}},
  \code{\link{nmr_meta_add}},
  \code{\link{nmr_meta_export}},
  \code{\link{nmr_meta_get_column}},
  \code{\link{nmr_meta_get}}, \code{\link{nmr_normalize}},
  \code{\link{nmr_pca_build_model}},
  \code{\link{nmr_pca_outliers_filter}},
  \code{\link{nmr_pca_outliers_plot}},
  \code{\link{nmr_pca_outliers_robust}},
  \code{\link{nmr_pca_outliers}},
  \code{\link{nmr_ppm_resolution}},
  \code{\link{plot.nmr_dataset_1D}},
  \code{\link{plot_webgl}}, \code{\link{rdCV_PLS_RF}},
  \code{\link{save_files_to_rDolphin}},
  \code{\link{to_ASICS}}, \code{\link{to_ChemoSpec}},
  \code{\link{to_rDolphin}},
  \code{\link{validate_nmr_dataset_1D}}
}
\concept{nmr_dataset_1D functions}