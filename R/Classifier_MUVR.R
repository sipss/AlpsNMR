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
#' Plot the model (a `MUVR` object) obtained from `rdCV_PLS_RF` function. For
#' more information about the multivariate modelling with minimally biased
#' variable selection (MUVR) from the `MUVR` package, see Shi et al., 2018 (DOI:
#' 10.1093/bioinformatics/bty710).
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
  if (!any(class(MVObj) == "MVObject")) {
    cat("\nWrong object class: Return NULL")
    return(NULL)
  }
  modNum = ifelse(model == "min", 1, ifelse(model == "mid", 
                                            2, 3))
  Y = MVObj$inData$Y
  nSamp = length(Y)
  if (missing(sampLabels)) 
    sampLabels = Y
  if (length(sampLabels) != nSamp) {
    warning("Length of sampLabels not equal to number of samples in Y. \n   Autonumbers used instead.")
    sampLabels = 1:nSamp
  }
  if (class(MVObj)[3] == "Regression") {
    YP = MVObj$yPred[, modNum]
    YPR = MVObj$yPredPerRep[[modNum]]
    if (is.null(ylim)) 
      ylim <- range(YPR)
    graphics::matplot(Y, YPR, pch = 20, xlab = "Original Y", ylab = "Predicted Y", 
            col = "grey", bty = "l", cex = 0.5, ylim = ylim)
    graphics::points(Y, YP, pch = 20)
    reg = stats::lm(YP ~ Y)
    graphics::abline(reg)
    graphics::legend("topleft", legend = c(paste("Model R2 =", signif(MVObj$fitMetric$R2[modNum], 
                                                            3)), paste("Model Q2 =", signif(MVObj$fitMetric$Q2[modNum], 
                                                                                            3))), bty = "n")
  }
  else if (class(MVObj)[3] == "Classification") {
    YP = MVObj$yPred[[modNum]]
    YPR = MVObj$yPredPerRep[[modNum]]
    if (is.null(ylim)) 
      ylim <- range(YPR)
    classes = 1:length(levels(Y))
    if (missing(factCols)) 
      factCols = classes + 1
    if (length(factCols) != length(classes)) {
      warning("Length of factCols not equal to number of levels in Y. \n   Autocolors used instead.")
      factCols = classes + 1
    }
    classNudge = 0.2 * ((classes - mean(classes))/(mean(classes) - 
                                                     1))
    plot(1:nSamp, Y, type = "n", ylim = ylim, xlab = "", 
         ylab = "Class prediction score", xaxt = "n")
    graphics::axis(1, at = 1:length(Y), labels = sampLabels, las = 3)
    for (cl in classes) {
      graphics::matpoints((1:nSamp) + classNudge[cl], YPR[, cl, 
                                                ], pch = 20, col = factCols[cl], cex = 0.5)
      graphics::points((1:nSamp) + classNudge[cl], YP[, cl], pch = 20, 
             col = factCols[cl])
    }
    for (li in 1:(nSamp + 1)) {
      graphics::abline(v = li - 0.5, lty = 3, col = "grey")
    }
    yClass = MVObj$yClass[, modNum]
    whichWrong = which(yClass != Y)
    wrongClass = as.numeric(Y[whichWrong])
    for (w in 1:length(wrongClass)) {
      graphics::points(whichWrong[w] + classNudge[wrongClass[w]], 
             YP[whichWrong[w], wrongClass[w]], cex = 2)
    }
    graphics::legend("topleft", legend = c(levels(Y), "misclassified"), 
           pch = c(rep(16, length(classes)), 1), col = c(factCols, 
                                                         1), cex = 0.8, pt.cex = c(rep(0.5, length(classes)), 
                                                                                   2), bty = "n")
  }
  else {
    YP = MVObj$yPred[, modNum]
    YPR = MVObj$yPredPerRep[[modNum]]
    graphics::matplot(YPR, 1:nSamp, pch = 20, col = "grey", cex = 0.5, 
            ylim = c(nSamp, 1), ylab = "Sample number", xlab = "Predicted Y")
    graphics::points(YP, 1:nSamp, pch = 20, col = "black")
    graphics::abline(h = nSamp/2 + 0.5, lty = 2)
    graphics::abline(v = 0, lty = 2)
  }
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
permutation_test_plot = function (MVObj, permMatrix, model = "mid", type = type,
                                 pos, xlab = NULL, xlim,
                                  ylim = NULL, breaks = "Sturges", main = NULL)
{
MUVR::permutationPlot(MVObj, permMatrix, model, type = type,
                      pos = pos, xlab = xlab, xlim = xlim,
                      ylim = NULL, breaks = breaks, main = NULL)
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


#' p-Value from permutation test
#' 
#' The fucntion calculates the cumulative (1-tailed) probability of 'actual'
#' belonging to 'h0' (`permutation_object` from the `permutation_test_model` function).
#'
#' @param model_actual The actual model performance  (e.g. misclassifications or Q2)
#' @param permutation_object Null hypothesis distribution from permutation test
#'   from `permutation_test_model` function
#'
#' @return The p-value indicating if there is significant differences between
#'   the model performance and the null hypothesis distribution from permutation
#'   test test
#' @export
#'
#' @examples
#' \dontrun{
#' P = permutation_test_model(MVObj)
#' p.value = p_value_perm(MVObj$miss[[2]], P[,2])
#' }
#' 
#' 
p_value_perm = function (model_actual, permutation_object){
  MUVR::pPerm(actual = model_actual, h0 = permutation_object)
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
#' # Build a model with the X data from your nmr object and your class
#' MVObj <- rdCV_PLS_RF(nmr_data(nmr_peak_table),
#' Y = nmr_peak_table_completed$Timepoint)
#' 
#' # Model performance
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
#' # Creating the model
#' model = rdCV_PLS_RF_ML(nmr_peak_table, label = "Timepoint", ML = TRUE)
#' 
#' # Model performance
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
  X = stats::na.omit(X)
  
  cl=parallel::makeCluster(parallel::detectCores()-1)
  doParallel::registerDoParallel(cl)
  model = MUVR::MUVR(X,  ML = TRUE)
  parallel::stopCluster(cl)# Stop parallel processing
  return(model)
}
