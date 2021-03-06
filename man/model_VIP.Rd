% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Classifier_MUVR.R
\name{model_VIP}
\alias{model_VIP}
\title{Deprecated function
Model VIP values}
\usage{
model_VIP(MVObj)
}
\arguments{
\item{MVObj}{a MUVR model}
}
\value{
a data frame with the order, name and average rank of selected variables
}
\description{
Once, the MVObj is created and validated, this function extracts autoselected
ranked variables from the model (MUVR object). See \code{rdCV_PLS_RF} function.
}
\examples{
message("Deprecated. MUVR is not compatible with Bioconductor, 
use bp_kfold_VIP_analysis method instead")
}
