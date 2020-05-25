#' Build a PLSDA model, optionally with multilevel
#' @param x the X training set
#' @param y the y training class to predict
#' @param identity the multilevel variable in [mixOmics::plsda]
#' @param ncomp The number of components of the model
#' @noRd
plsda_build <- function(x, y, identity, ncomp) {
    plsda_model <- NULL
    tryCatch({
        suppressMessages(
            utils::capture.output({
                plsda_model <- mixOmics::plsda(
                    X = x, Y = y, ncomp = ncomp,
                    scale = TRUE, multilevel = identity
                )
            })
        )
    }, error = function(e) {
        message("Error building PLSDA, continuing")
    })
    plsda_model
}


#' Compute the area under the ROC curve of a PLS-DA model on a test subset
#' @param plsda_model A mixOmics plsda model
#' @param x_test the x test set
#' @param y_test the y test class to predict
#' @param identity_test the multilevel variable in [mixOmics::plsda]
#' @return A list with two elements: 
#'    - `aucs`: A data frame with two columns: `ncomp` (the number of components) and
#'        `auc` the area under roc curve for that number of components. For multiclass problems
#'        the AUC returned is the mean of all the one-vs-other AUCs.
#'    - `aucs_full`: A list of matrices, as returned by [mixOmics::auroc].
#' @noRd
plsda_auroc <- function(plsda_model, x_test, y_test, identity_test) {
    aucs <- numeric(0L)
    aucs_full <- list()
    tryCatch({
        suppressMessages(
            utils::capture.output({
                roc <- mixOmics::auroc(plsda_model, newdata = x_test,
                                                             outcome.test = y_test,
                                                             multilevel = identity_test,
                                                             plot = FALSE)
            })
        )
        aucs <- purrr::map_dbl(roc, function(x) mean(x[, "AUC"]))
        aucs_full <- roc
    }, error = function(e) {
        message("Error in auroc estimation, continuing")
    })
    
    ncomps <- as.integer(gsub(pattern = "Comp(.*)", replacement = "\\1", x = names(aucs)))
    
    list(aucs = data.frame(ncomp = ncomps,
                                                 auc = aucs,
                                                 stringsAsFactors = FALSE),
             aucs_full = aucs_full)
}

#' Compute the variable importance in the projection
#' @param plsda_model A mixOmics plsda model
#' @return A matrix with the variable importance in the projection
#' @noRd
plsda_vip <- function(plsda_model) {
    vip <- NULL
    tryCatch({
        suppressMessages(
            utils::capture.output({
                vip <- mixOmics::vip(object = plsda_model)
            })
        )
    }, error = function(e) {
        message("Error in vip, continuing")
    })
    vip
}


#' Callback for building a PLSDA model, computing the AUROC and extract the VIP
#'
#' @param x_train Training data for x
#' @param y_train Training data for y
#' @param identity_train Training data for the identities
#' @param x_test Test data for x
#' @param y_test Test data for y
#' @param identity_test Test data for the identities
#' @param ncomp Number of components to use in the model
#' @param return_model A logical.
#' @param return_auroc A logical.
#' @param return_auroc_full A logical.
#' @param return_vip A logical.
#' 
#' 
#' For multiclass problems the AUC returned is the mean of all the one-vs-other AUCs.
#'
#' @return A list with the model, the area under the roc curve and the VIP items.
#' @noRd
callback_plsda_auroc_vip <- function(x_train, y_train, identity_train, x_test, y_test, identity_test,
                                                                         ncomp, return_model = FALSE, return_auroc = TRUE,
                                                                         return_auroc_full = FALSE, return_vip = FALSE) {
    plsda_model <- plsda_build(x_train, y_train, identity_train, ncomp = max(ncomp))
    out <- list(model = NULL, auroc = NULL, auroc_full = NULL, vip = NULL)
    if (isTRUE(return_model)) {
        out$model <- plsda_model
    }
    if (isTRUE(return_auroc) || isTRUE(return_auroc_full)) {
        aurocs <- plsda_auroc(plsda_model, x_test, y_test, identity_test)
        if (isTRUE(return_auroc)) {
            out$auroc <- aurocs$aucs

        }
        if (isTRUE(return_auroc_full)) {
            out$auroc_full <- aurocs$aucs_full
        }
    }
    
    if (isTRUE(return_vip)) {
        vip <- plsda_vip(plsda_model)
        out$vip <- vip
    }
    out
}


#' Callback to choose the best number of latent variables based on the AUC threshold
#' 
#' @param auc_threshold Threshold on the increment of AUC. Increasing the number of
#' latent variables must increase the AUC at least by this threshold.
#'
#' @return The actual function to compute the best number of latent variables according to a threshold on the increment of AUC
#' @noRd
fun_choose_best_ncomp_auc_threshold <- function(auc_threshold = 0.05) {
    force(auc_threshold)
    
    # Choose best number of latent variables based on a threshold on the auc increment.
    #' @param inner_cv_results A list of elements returned by [callback_plsda_auroc_vip]
    #' @return A list with:
    #'    - `train_evaluate_model_args`: A list wit one element named `ncomp` with the number of latent variables selected
    #'             for each outer cross-validation
    #'    - `num_latent_var`: A data frame with the number of latent variables chosen for each outer cross-validation
    #'    - `diagnostic_plot`: A plot showing the evolution of the AUC vs the number of latent variables for each iteration
    #'    - `model_performances`: A data frame with the AUC model performances 
    function(inner_cv_results) {
    model_performances <- inner_cv_results %>%
        purrr::map("auroc") %>%
        purrr::map_dfr(~ ., .id = "outer_inner") %>%
        tidyr::separate("outer_inner",
                                        into = c("cv_outer_iteration", "cv_inner_iteration"),
                                        convert = TRUE)
    # There is a more elegant way to do this.
    nlv <- model_performances %>%
        dplyr::group_by(.data$cv_outer_iteration, .data$cv_inner_iteration) %>%
        dplyr::arrange(.data$cv_outer_iteration, .data$cv_inner_iteration, .data$ncomp) %>%
        dplyr::mutate(auc_diff = ifelse(is.na(dplyr::lag(.data$auc)), .data$auc, .data$auc - dplyr::lag(.data$auc))) %>%
        dplyr::mutate(auc_limit_cumany = dplyr::cumall(.data$auc_diff > !!auc_threshold)) %>%
        dplyr::mutate(auc_limit_cumanyd = .data$auc_limit_cumany == TRUE & dplyr::lead(.data$auc_limit_cumany) == FALSE) %>% 
        dplyr::filter(.data$auc_limit_cumanyd == TRUE) %>%
        dplyr::select(-.data$auc_limit_cumany, -.data$auc_limit_cumanyd) %>%
        dplyr::ungroup() %>%
        dplyr::group_by(.data$cv_outer_iteration) %>%
        dplyr::summarise(ncomp = round(stats::median(.data$ncomp)))
    
    plot_to_choose_nlv <- ggplot2::ggplot(model_performances) + 
        ggplot2::geom_jitter(ggplot2::aes(x = .data$ncomp, y = .data$auc,
                                                                            group = .data$ncomp, color = as.character(.data$cv_inner_iteration)), 
                                                 width = 0.25, height = 0) +
        ggplot2::geom_vline(data = nlv, mapping = ggplot2::aes(xintercept = .data$ncomp), color = "red") +
        ggplot2::scale_x_continuous(name = "Number of latent variables", breaks = function(limits) {
            seq(from = 1, to = max(limits))
        }) +
        ggplot2::scale_y_continuous(name = "Area Under ROC") +
        ggplot2::facet_wrap(~cv_outer_iteration) + 
        ggplot2::guides(colour = "none")
    
    list(train_evaluate_model_args = list(ncomp = nlv$ncomp),
             num_latent_var = nlv,
             diagnostic_plot = plot_to_choose_nlv,
             model_performances = model_performances)
    }
}

#################### Validation #######

#' Callback to digest the results of the outer cross validation
#' @noRd
callback_outer_cv_auroc_vip <- function(outer_cv_results) {
    auroc <- outer_cv_results %>%
        purrr::map("auroc") %>%
        purrr::map_dfr(~ ., .id = "cv_outer_iteration") %>%
        dplyr::mutate(cv_outer_iteration = as.integer(.data$cv_outer_iteration)) %>%
        dplyr::group_by(.data$cv_outer_iteration) %>%
        dplyr::filter(.data$ncomp == max(.data$ncomp)) %>%
        dplyr::ungroup() %>%
        dplyr::arrange(.data$cv_outer_iteration)
    
    vip_vectors <- outer_cv_results %>%
        purrr::map("vip") %>%
        purrr::map2(auroc$ncomp, function(vip_matrix, selected_ncomp) {
            vip_vec <- as.numeric(vip_matrix[, selected_ncomp, drop = TRUE])
            names(vip_vec) <- rownames(vip_matrix)
            vip_vec
        })
    
    vip_ranks <- do.call(cbind, purrr::map(vip_vectors, ~rank(-.)))
    
    vip_rp <- apply(vip_ranks, 1, function(x) exp(mean(log(x)))) # geom mean (RankProducts)
    
    list(auroc = auroc,
             vip_vectors = vip_vectors,
             vip_rankproducts = vip_rp)
}


#' Method for nmr_data_analysis (PLSDA model with AUROC and VIP outputs)
#' @param ncomp Max. number of latent variables to explore in the PLSDA analysis
#' @param auc_increment_threshold Choose the number of latent variables when the
#' AUC does not increment more than this threshold.
#' 
#' @return Returns an object to be used with [nmr_data_analysis] to perform a (optionally 
#' multilevel) PLS-DA model, using the area under the ROC curve as figure of
#' merit to determine the optimum number of latent variables.
#' 
#' 
#' @export
#' @examples 
#' method <- plsda_auroc_vip_method(3)
#' 
plsda_auroc_vip_method <- function(ncomp, auc_increment_threshold = 0.05) {
    new_nmr_data_analysis_method(
        train_evaluate_model = callback_plsda_auroc_vip,
        train_evaluate_model_params_inner = list(
            ncomp = ncomp,
            return_model = FALSE,
            return_auroc = TRUE,
            return_auroc_full = FALSE,
            return_vip = FALSE
        ),
        choose_best_inner = fun_choose_best_ncomp_auc_threshold(auc_threshold = auc_increment_threshold),
        train_evaluate_model_params_outer = list(
            return_model = TRUE,
            return_auroc = TRUE,
            return_auroc_full = TRUE,
            return_vip = TRUE
        ),
        train_evaluate_model_digest_outer = callback_outer_cv_auroc_vip
    )
}


#' Compare PLSDA auroc VIP results
#'
#' @param ... Results of [nmr_data_analysis] to be combined. Give each result a name.
#'
#' @return A plot of the AUC for each method
#' @export
#' @examples 
#' \dontrun{
#' method <- plsda_auroc_vip_method(3)
#' #analisis missed
#' plsda_auroc_vip_compare()
#' }
plsda_auroc_vip_compare <- function(...) {
    dots <- list(...)
    class_compare <- names(dots)
    if (is.null(class_compare) || any(nchar(class_compare) == 0)) {
        stop("All arguments should be named")
    }
    
    auroc_tables <- dots %>%
        purrr::map("outer_cv_results_digested") %>%
        purrr::map("auroc") %>%
        purrr::map2(class_compare, function(auroc, group_name) {
            auroc %>% dplyr::select(.data$auc) %>% dplyr::mutate(Group = !!group_name)
        })
    
    toplot <- do.call(rbind, c(auroc_tables, list(stringsAsFactors = FALSE)))
    ggplot2::ggplot(toplot) + 
        ggplot2::geom_boxplot(ggplot2::aes(x = .data$Group, y = .data$auc, fill = .data$Group), show.legend = FALSE) +
        ggplot2::scale_x_discrete(name = "Model") +
        ggplot2::scale_y_continuous(name = "Area under ROC")
}
