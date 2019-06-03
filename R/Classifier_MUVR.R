#' Feature selection and validation in multivariate analysis
#'
#' Statistical analysis and feature selection in a repeated double 
#' cross-validation frame based on the partial least squares
#' (PLS) or random forest (RF) analyses using an algorithm 
#' for multivariate modelling with minimally biased variable 
#' selection (MUVR) from the `MUVR` package.
#' 
#' @family nmr_dataset_1D functions
#' @inheritParams MUVR::MUVR
#'
#' @export
#' @examples 
#' \dontrun{
#' model = rdCV_PLS_RF(X = nmr_data(nmr_peak_table),
#'                     Y = nmr_meta_get(nmr_peak_table, groups = "external")$Timepoint,
#'                     ML = F, method = "PLS", 
#'                     fitness = "AUROC", 
#'                     nRep = 12, 
#'                     nOuter = 4, 
#'                     varRatio = 0.8,
#'                     scale = T)
#' }
#' @return a MUVR model containing selection parameters, validation and fitness
#' @references Shi,L. et al. (2018) Variable selection and validation in multivariate modelling. Bioinformatics.
rdCV_PLS_RF = function (X, Y, ID, scale = TRUE, nRep = 10, nOuter = 5, 
                        nInner, varRatio = 0.75, DA = FALSE, 
                        fitness = "MISS", 
                        method = "PLS", nCompMax, 
                        methParam, ML = FALSE, modReturn = FALSE, 
                        logg = FALSE, parallel = TRUE)
{
  cl=parallel::makeCluster(parallel::detectCores()-1)
  doParallel::registerDoParallel(cl)
  model = MUVR::MUVR(X, Y, ID, scale, nRep, nOuter, 
                     nInner, varRatio, DA, 
                     fitness, 
                     method, nCompMax, 
                     methParam, ML = FALSE , modReturn, 
                     logg, parallel)
  parallel::stopCluster(cl)# Stop parallel processing
  return(model)
}

#' Permutation test on MUVR
#' 
#' Make permutations with data and default settings from an actual MUVR object
#' 
#' @inheritParams MUVR::permutations
#' @return A permutation matrix with permuted values
#' @export
#' @examples
#' \dontrun{
#' P = permutation_test_model(model)
#' }

permutation_test_model = function (MVObj, nPerm = 50, nRep, nOuter, varRatio, parallel)
{
cl=parallel::makeCluster(parallel::detectCores()-1)
doParallel::registerDoParallel(cl)
permMatrix = MUVR::permutations(MVObj, nPerm, nRep, nOuter, varRatio, parallel)
parallel::stopCluster(cl)
return (permMatrix)
}

#' Model plot
#' 
#' Plot the model obtained from `rdCV_PLS_RF` function (a `MUVR` object)
#'
#' @inheritParams MUVR::plotMV
#'
#' @return A plot with the model performance
#' @export
#' @examples 
#' \dontrun{
#' model = rdCV_PLS_RF(X = nmr_data(nmr_peak_table),
#'                     Y = nmr_meta_get(nmr_peak_table, groups = "external")$Timepoint,
#'                     ML = F, method = "PLS", 
#'                     fitness = "AUROC", 
#'                     nRep = 12, 
#'                     nOuter = 4, 
#'                     varRatio = 0.8,
#'                     scale = T)
#'
#' MUVR_model_plot(model)
#' }
#' 
MUVR_model_plot = function (MVObj, model = "mid", factCols, sampLabels, ylim = NULL)
{
  MUVR::plotMV(MVObj, model, factCols, sampLabels, ylim)
}

#' Permutation test plot
#' 
#' Plot permutation test using actual model and permutated models
#'
#' @inheritParams MUVR::permutationPlot
#'
#' @return A plot with the comparison between the actual model versus the permuted models
#' @export
#' @examples
#' \dontrun{
#' P = permutation_test_model(model)
#' permutation_test_plot (model, P)
#' }
#' 
permutation_test_plot = function (MVObj, permMatrix, model = "mid")
{
MUVR::permutationPlot(MVObj, permMatrix, model)
}

#' model VIP values
#' 
#' The function extracts autoselected ranked variables from the model (MUVR object)
#'
#' @inheritParams MUVR::getVIP
#' @param MVObj a MUVR model
#' @return a data frame with the order, name and average rank of selected variables
#' @export
#'
#' @examples
#' \dontrun{
#' VIPs = model_VIP(MVObj)
#' }
#' 
model_VIP = function(MVObj, model = "mid"){
  MUVR::getVIP(MVObj, model)
}


#' Confusion matrix of the MUVR model
#' 
#' The function makes a confusion matrix from a MUVR model
#'
#' @inheritParams MUVR::confusionMatrix
#'
#' @return A confusion matrix of the model comparing actual vs predicted class
#' @export
#'
#' @examples
#' \dontrun{
#' confusion_matrix(MVObj)
#' }
#' 
confusion_matrix = function(MVObj, model = "mid"){
  MUVR::confusionMatrix(MVObj, model)
}

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
#' @param label the name of the variable to test (e.g. "Timepoint")
#' @inheritParams MUVR::MUVR
#'
#' @export
#' @examples 
#' \dontrun{
#' model = rdCV_PLS_RF_ML(nmr_peak_table, label = "Timepoint", ML = TRUE)
#' MUVR_model_plot(model)
#' }
#' @return a MUVR model containing selection parameters, validation and fitness
#' @references Shi,L. et al. (2018) Variable selection and validation in multivariate modelling. Bioinformatics.
rdCV_PLS_RF_ML = function (nmr_peak_table, label, scale = TRUE, nRep = 10, nOuter = 5, 
                        nInner, varRatio = 0.75, DA = FALSE, 
                        fitness = "MISS", 
                        method = "PLS", ML = TRUE, modReturn = FALSE, 
                        logg = FALSE, parallel = TRUE)
{
  ordered = AlpsNMR::get_integration_with_metadata(nmr_peak_table)
  ordered = ordered[order(ordered[[label]],ordered$NMRExperiment),]
  levelAB = levels(as.factor(ordered[[label]]))
  A = ordered[ordered[[label]]==levelAB[1],]
  B = ordered[ordered[[label]]==levelAB[2],]
  
  coln = colnames(nmr_peak_table[["metadata"]][["external"]])
  a = A %>% dplyr::select(-coln)
  b = B %>% dplyr::select(-coln)
  X = a[1:nrow(a),] - b[1:nrow(a),]
  X = X[, colSums(is.na(X)) != nrow(X)]
  X = na.omit(X)
  
  cl=parallel::makeCluster(parallel::detectCores()-1)
  doParallel::registerDoParallel(cl)
  model = MUVR::MUVR(X,  ML = TRUE)
  parallel::stopCluster(cl)# Stop parallel processing
  return(model)
}



#' #' Machine learning (delete if everything works properly, kept as example)
#' #'
#' #' Machine learning on a integration object `nmr_peak_table`
#' #' @param nmr_peak_table a dataset with integration values
#' #' @param groups a vector with the class groups
#' #' @inheritParams MUVR::MUVR
#' #' @return a MUVR model
#' #'
#' machine_learning = function (nmr_peak_table, groups, scale = T){
#'   cl=parallel::makeCluster(2)
#'   doParallel::registerDoParallel(cl)
#'   model=MUVR::MUVR(nmr_data(nmr_peak_table), groups, ML=F, 
#'                    method ="PLS", fitness = "MISS", nRep=10, 
#'                    nOuter=5, varRatio=0.8, scale = T)
#'   parallel::stopCluster(cl)# Stop parallel processing
#'   return(model)
#' }
