#' Feature selection and validation in multivariate analysis
#' 
#' Numeric VIPs vector
#'
#' The function extracts the VIPs vector (numeric) from the model_VIP(MVObj)
#' function. It is not necessary if you have the ppm values in a numeric vector.
#' This is needed in case that an automated pipeline is applied, connecting the
#' output from model_VIP(MVObj) to `nmr_identify_regions` family functions.
#'
#' @param VIPs a dataframe from the model_VIP(MVObj) function. It requires a
#'     "ppms" variable
#'
#' @return a numeric ppm vector ready to be identified with
#'     `nmr_identify_regions_blood`, `nmr_identify_cell` or
#'     `nmr_identify_regions_urine`
#' @export
#'
#' @examples 
#' message("Deprecated. MUVR is not compatible with Bioconductor, 
#' use bp_kfold_VIP_analysis method instead")
#'
#'
ppm_VIP_vector <- function(VIPs) {
    ppm_to_assign = tidyr::separate(VIPs,
                                    col = "name",
                                    into = c("x1", "ppms"),
                                    sep = "_")
    ppm_to_assign = as.numeric(ppm_to_assign$ppms)
    return(ppm_to_assign)
}

#' Deprecated function
#' @seealso nmr_data_analysis
#' Feature selection and validation in multivariate analysis
#'
#' Statistical analysis and feature selection in a repeated double
#' cross-validation frame based on the partial least squares (PLS) or random
#' forest (RF) analyses using an algorithm for multivariate modelling with
#' minimally biased variable selection (MUVR) from the `MUVR` package. If your
#' work with a `nmr_peak_table` object from AlpsNMR, first you need to extract
#' the X data from the main nmr_dataset object (e.g. your peak table) with the
#' `nmr_data` function, otherwise you would try to set a list on the X. You also
#' need to set the class from this object, or just set it from another Y vector.
#'
#' @family nmr_dataset_1D functions
#' @param nmr_peak_table an AlpsNMR integration object (2 classes)
#' @export
#' @examples 
#' message("Deprecated. MUVR is not compatible with Bioconductor, 
#' use bp_kfold_VIP_analysis method instead")
#' @return a MUVR model containing selection parameters, validation and fitness
#' @references Shi,L. et al. (2018) Variable selection and validation in multivariate modelling. Bioinformatics.
rdCV_PLS_RF <- function(nmr_peak_table){
    .Defunct("nmr_data_analysis")
}


#' Deprecated function
#' Model plot
#'
#' Plot the model (a `MUVR` object) obtained from `rdCV_PLS_RF` function. For
#' more information about the multivariate modelling with minimally biased
#' variable selection (MUVR) from the `MUVR` package, see Shi et al., 2018 (DOI:
#' 10.1093/bioinformatics/bty710).
#' 
#' @param MVObj a MUVR model
#' @param model What type of model to plot ('min', 'mid' or 'max'). Defaults to 'mid'.
#' @param factCols An optional vector with colors for the factor levels (in the same order as the levels).
#' @param sampLabels Sample labels (optional; implemented for classification)
#' @param ylim Optional for imposing y-limits for regression and classification analysis
#'
#' @return A plot with the model performance
#' @export
#' @examples 
#' message("MUVR is not compatible with Bioconductor, 
#' use bp_kfold_VIP_analysis method instead")
MUVR_model_plot = function (MVObj,
                            model = "mid",
                            factCols,
                            sampLabels,
                            ylim = NULL) {
    .Defunct()
}

#' Deprecated function
#' Model VIP values
#'
#' Once, the MVObj is created and validated, this function extracts autoselected
#' ranked variables from the model (MUVR object). See `rdCV_PLS_RF` function.
#'
#' @param MVObj a MUVR model
#' @return a data frame with the order, name and average rank of selected variables
#' @export
#' @examples 
#' message("Deprecated. MUVR is not compatible with Bioconductor, 
#' use bp_kfold_VIP_analysis method instead")
model_VIP = function(MVObj) {
    .Defunct()
}   

#' Deprecated function
#' p-Value from permutation test
#'
#' The fucntion calculates the cumulative (1-tailed) probability of 'actual'
#' belonging to 'h0' (`permutation_object` from the `permutation_test_model` function).
#'
#' @param model_actual The actual model performance    (e.g. misclassifications or Q2)
#' @param permutation_object Null hypothesis distribution from permutation test
#'     from `permutation_test_model` function
#'
#' @return The p-value indicating if there is significant differences between
#'     the model performance and the null hypothesis distribution from permutation
#'     test test
#' @export
#' @examples 
#' message("MUVR is not compatible with Bioconductor, 
#' use bp_kfold_VIP_analysis method instead")
p_value_perm = function (model_actual, permutation_object) {
    .Defunct()
}   

#' Deprecated function
#' Confusion matrix of the MUVR model
#'
#' After creating a model with the `rdCV_PLS_RF` function, you can run
#' `confusion_matrix` on the model to make a confusion matrix from MUVR. This
#' gives information about the model performance (e.g. classification rate).
#'
#' @param MVObj a MUVR model
#' @param model What type of model to plot ('min', 'mid' or 'max'). Defaults to 'mid'.
#'
#' @return A confusion matrix of the model comparing actual vs predicted class
#' @export
#' @examples 
#' message("Deprecated. MUVR is not compatible with Bioconductor, 
#' use bp_kfold_VIP_analysis method instead")
confusion_matrix = function(MVObj, model = "mid") {
    .Defunct()
}   

#' Deprecated function
#' Feature selection and validation in MULTILEVEL analysis
#'
#' Statistical analysis and feature selection in a repeated double
#' cross-validation frame based on the partial least squares
#' (PLS) or random forest (RF) analyses using an algorithm
#' for multivariate modelling with minimally biased variable
#' selection (MUVR) from the `MUVR` package. The function
#' `rdCV_PLS_RF_ML` allows the multilevel comparison,
#' especially useful in crossover or longitudinal studies
#' (2 timepoints) considering the same individual (it
#' requires 2 samples of the same observation).
#'
#' @family nmr_dataset_1D functions
#' @param nmr_peak_table an AlpsNMR integration object (2 classes)
#'
#' @export
#' @examples 
#' message("Deprecated. MUVR is not compatible with Bioconductor, 
#' use bp_kfold_VIP_analysis method instead")
#' @return a MUVR model containing selection parameters, validation and fitness
#' @references Shi,L. et al. (2018) Variable selection and validation in multivariate modelling. Bioinformatics.
rdCV_PLS_RF_ML = function (nmr_peak_table){
    .Defunct()
}   

#' Deprecated function
#' Extracts AUC value
#'
#' The function extracts the AUC value from the middle `MUVR` model
#'
#' @param MVObj a MUVR model
#'
#' @return the AUC value of the middle model
#' @export
#' @examples 
#' message("Deprecated. MUVR is not compatible with Bioconductor, 
#' use bp_kfold_VIP_analysis method instead")
AUC_model <- function (MVObj) {
    .Defunct()
}   
