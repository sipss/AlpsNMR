#' @title Perform an automatic 1D NMR phasing
#' @description
#' It is quite often the case that (e.g. Bruker's) 1D autophase isn't quite right.
#' This uses `NMRphasing::NMRphasing` to automatically rephase data in the spectral domain. 
#' A number of algorithms are available (see NMRphasing's documentation),
#' of which NLS, MPC_DANM and SPC_DANM are the most recent. 
#' 
#' Perhaps obviously, you are likely to get better results loading the complex data 
#' than having real-only spectra (i.e. pass `all_components = T` to `nmr_read_samples`). 
#' 
#' Run this before spectral interpolation. 
#' 
#' @param samples An `nmr_dataset_family` 1D object 
#' @param method The autophasing method -- see `NMRphasing::NMRphasing` for details. 
#' @param withBC `NMRphasing::NMRphasing` performs an integrated baseline correction -- this parameter enables or disables it. 
#' @return A (hopefully better phased) `nmr_dataset_family` 1D object, with updated real and imaginary parts. 
#' @examples
#' dir_to_demo_dataset <- system.file("dataset-demo", package = "AlpsNMR")
#' dataset <- nmr_read_samples_dir(dir_to_demo_dataset)
#' dataset <- nmr_dataset_autophase(dataset,method="MPC_DANM")
#' dataset <- nmr_interpolate_1D(dataset, axis = c(min = 1, max = 2, by = 0.002))
#' plot(dataset)
#' 
#' @export 
nmr_dataset_autophase <- function(samples, method= c("NLS", "MPC_DANM", "MPC_EMP", "SPC_DANM", "SPC_EMP", "SPC_AAM", "SPC_DSM"), withBC=F) {
    
    method <- match.arg(method)
    BiocParallel::bplapply(seq_len(length(samples[["data_1r"]])), function(i) {
        
        data_to_phase <- if("data_1i" %in% objects(samples)) {
            complex(real = samples[["data_1r"]][[i]], 
                    imaginary = samples[["data_1i"]][[i]], 
                    length.out = length(samples[["data_1r"]][[i]]
                                        ))
        } else {
            samples[["data_1r"]][[i]]
        }
        
        phased_data <- NMRphasing::NMRphasing(data_to_phase, method = method, 
                                              absorptionOnly = !"data_1i" %in% objects(samples), 
                                              withBC = F)
        
        if("data_1i" %in% objects(samples)) {
            samples[["data_1r"]][[i]] <- Re(phased_data)
            samples[["data_1i"]][[i]] <- Im(phased_data)
        } else {
            samples[["data_1r"]][[i]] <- phased_data
        }
        
    })
    
    return(samples)
}

#' @title Export data for the spectral quantification library ASICS  
#' @description
#' A simple helper function for mangling data in the right format for the 
#' spectral quantification library ASICS. 
#' 
#' @param samples An `nmr_dataset_family` 1D object 
#' @param ... Additional arguments are passed 
#'            directly to `ASICS::createSpectra`, which (in theory) provide an 
#'            opportunity to use distinct normalisation methods. 
#' 
#' @return An `ASICS::Spectra` object 
#' @examples
#' # forAsics <- alps_asics(dataset)
#' # ASICS(forAsics)
#' 
#' @export 

alps_asics <- function(samples, ...) {
  forAsics <- t(nmr_data(samples))
  forAsics <- ASICS::createSpectra(forAsics, ...)
  forAsics@sample.name <-samples$metadata$external$NMRExperiment
  return(forAsics)
}


