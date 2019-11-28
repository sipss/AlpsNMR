#' NMR peak identification (plasma/serum samples)
#' 
#' Identify given regions and return a data frame with plausible assignations 
#' in human plasma/serum samples.
#' 
#' @return a data frame with plausible assignations.
#'
#' @examples 
#' \dontrun{
#' # We identify regions from from the corresponding ppm storaged in a vector.
#' library(AlpsNMR)
#' ppm_to_assign <- c(4.060960203, 3.048970634,2.405935596,
#' 3.24146865,0.990616851,1.002075066,0.955325548)
#' identification <- nmr_identify_regions_blood (ppm_to_assign)
#' }
#' @export
#' @family peak detection functions
#' @family peak integration functions
#' @param ppm_to_assign A vector with the ppm regions to assign
#' @param num_proposed_compounds set the number of proposed metabolites sorted by the number times reported in the HMDB: `HMDB_blood`.
#' @param verbose Logical value. Set it to TRUE to print additional information
nmr_identify_regions_blood <- function(ppm_to_assign, num_proposed_compounds = 3, verbose = FALSE){
  HMDB_blood <- NULL
  utils::data("HMDB_blood", package = "AlpsNMR", envir = environment())
  output_assignation_list <- HMDB_blood[NULL,]
  
  for (ppm in ppm_to_assign) {
    lower_ppm_right_edge <- ppm - 0.015
    higher_ppm_left_edge <- ppm + 0.015
    
    ind <- intersect(which(HMDB_blood$Shift_ppm < higher_ppm_left_edge),
                     which(HMDB_blood$Shift_ppm > lower_ppm_right_edge))
    assignation_list <- as.data.frame(HMDB_blood[ind,])
    if (isTRUE(verbose)) {
      message("your peak at ",ppm, " probably corresponds to ",assignation_list[1,1], ", ",assignation_list[2,1], " or ", assignation_list[3,1])
      message("")
    }
    output_assignation_list <- rbind(output_assignation_list, assignation_list[1:num_proposed_compounds,])
  }
  output_assignation_list$ppm_to_assign <- rep(ppm_to_assign, each = num_proposed_compounds)
  
  # counts=output_assignation_list %>% dplyr::count(Metabolite) %>% dplyr::arrange(dplyr::desc(n))
  # colnames(counts) <- c("Metabolite", "Counts")
  # output_assignation_list <- merge(output_assignation_list,counts, by = "Metabolite")
  output_assignation_list <- output_assignation_list[order(output_assignation_list$ppm_to_assign,-output_assignation_list$Blood_concentration, -output_assignation_list$n_reported_in_Blood),]
  return(output_assignation_list)
}


#' The Human Metabolome DataBase multiplet table: blood metabolites normally found in NMR-based metabolomics
#'
#' @name HMDB_blood
#' @docType data
#' @references \url{hmdb.ca}
#' @keywords data
NULL

#' NMR peak identification (urine samples)
#' 
#' Identify given regions and return a data frame with plausible assignations 
#' in human urine samples. The data frame contains the column "Bouatra_2013" showing if
#' the proposed metabolite was reported in this publication as regular urinary metabolite.
#' 
#' @return a data frame with plausible assignations.
#'
#' @examples 
#' \dontrun{
#' # We identify regions from from the corresponding ppm storaged in a vector.
#' library(AlpsNMR)
#' ppm_to_assign <- c(4.060960203, 3.048970634,2.405935596,
#' 3.24146865,0.990616851,1.002075066,0.955325548)
#' identification <- nmr_identify_regions_urine (ppm_to_assign, num_proposed_compounds = 5)
#' }
#' @export
#' @family peak detection functions
#' @family peak integration functions
#' @param ppm_to_assign A vector with the ppm regions to assign
#' @param num_proposed_compounds set the number of proposed metabolites sorted by the number times reported in the HMDB: `HMDB_urine`.
#' @param verbose Logical value. Set it to TRUE to print additional information
nmr_identify_regions_urine <- function(ppm_to_assign, num_proposed_compounds = 5, verbose = FALSE){
  HMDB_urine <- NULL
  utils::data("HMDB_urine", package = "AlpsNMR", envir = environment())
  output_assignation_list <- HMDB_urine[NULL,]
  
  for (ppm in ppm_to_assign) {
    lower_ppm_right_edge <- ppm - 0.015
    higher_ppm_left_edge <- ppm + 0.015
    
    ind <- intersect(which(HMDB_urine$Shift_ppm < higher_ppm_left_edge),
                     which(HMDB_urine$Shift_ppm > lower_ppm_right_edge))
    assignation_list <- as.data.frame(HMDB_urine[ind,])
    if (isTRUE(verbose)) {
      message("your peak at ",ppm, " probably corresponds to ",assignation_list[1,1], ", ",assignation_list[2,1],", ", assignation_list[3,1],", ", assignation_list[4,1]," or ", assignation_list[5,1])
      message("")
    }
    output_assignation_list <- rbind(output_assignation_list,assignation_list[1:num_proposed_compounds,])
  }
  output_assignation_list$ppm_to_assign <- rep(ppm_to_assign,each = num_proposed_compounds)
  
  #counts=output_assignation_list %>% dplyr::count(Metabolite) %>% dplyr::arrange(dplyr::desc(n))
  #colnames(counts) <- c("Metabolite", "Counts")
  #output_assignation_list <- merge(output_assignation_list,counts, by = "Metabolite")
  output_assignation_list <- output_assignation_list[order(output_assignation_list$ppm_to_assign,-output_assignation_list$Urine_concentration, -output_assignation_list$n_reported_in_Urine),]
  return(output_assignation_list)
}

#' The Human Metabolome DataBase multiplet table: urine metabolites normally found in NMR-based metabolomics
#'
#' @name HMDB_urine
#' @docType data
#' @references \url{hmdb.ca}
#' @keywords data
NULL

#' NMR peak identification (cell samples)
#' 
#' Identify given regions and return a data frame with plausible assignations 
#' in cell samples.
#' 
#' @return a data frame with plausible assignations.
#'
#' @examples 
#' # We identify regions from from the corresponding ppm storaged in a vector.
#' library(AlpsNMR)
#' ppm_to_assign <- c(4.060960203, 3.048970634,2.405935596,
#' 3.24146865,0.990616851,1.002075066,0.955325548)
#' identification <- nmr_identify_regions_cell (ppm_to_assign, num_proposed_compounds = 3)
#' @export
#' @family peak detection functions
#' @family peak integration functions
#' @param ppm_to_assign A vector with the ppm regions to assign
#' @param num_proposed_compounds set the number of proposed metabolites in `HMDB_cell`.
#' @param verbose Logical value. Set it to TRUE to print additional information
nmr_identify_regions_cell <- function(ppm_to_assign, num_proposed_compounds = 3, verbose = FALSE){
  HMDB_cell <- NULL
  utils::data("HMDB_cell", package = "AlpsNMR", envir = environment())
  output_assignation_list <- HMDB_cell[NULL,]
  
  for (ppm in ppm_to_assign) {
    lower_ppm_right_edge <- ppm - 0.015
    higher_ppm_left_edge <- ppm + 0.015
    
    ind <- intersect(which(HMDB_cell$Shift_ppm < higher_ppm_left_edge),
                     which(HMDB_cell$Shift_ppm > lower_ppm_right_edge))
    assignation_list <- as.data.frame(HMDB_cell[ind,])
    if (isTRUE(verbose)) {
      message("your peak at ",ppm, " probably corresponds to ",assignation_list[1,1], ", ",assignation_list[2,1],", ", assignation_list[3,1])
      message("")
    }
    output_assignation_list <- rbind(output_assignation_list,assignation_list[1:num_proposed_compounds,])
  }
  output_assignation_list$ppm_to_assign <- rep(ppm_to_assign,each = num_proposed_compounds)
  output_assignation_list <- output_assignation_list[order(output_assignation_list$ppm_to_assign),]
  return(output_assignation_list)
}


#' The Human Metabolome DataBase multiplet table: cell metabolites normally found in NMR-based metabolomics
#'
#' @name HMDB_cell
#' @docType data
#' @references \url{hmdb.ca}
#' @keywords data
#' 
NULL
