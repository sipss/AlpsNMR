
# This is the server logic for a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)
library(shinyFiles)
library(NIHSnmr)
library(DT)
library(plotly)
library(tibble)

find_peaks <- function(x, y, m = 3, y_min = NA){
  # Return indices of y that are peaks.
  # A peak here is a local maximum in a 2*m+1 points neighbourhood that is larger than y_min
  shape <- diff(sign(diff(y, na.pad = FALSE)))
  pks <- sapply(which(shape < 0), FUN = function(i){
    z <- i - m + 1
    z <- ifelse(z > 0, z, 1)
    w <- i + m + 1
    w <- ifelse(w < length(y), w, length(y))
    if (all(y[c(z:i, (i + 2):w)] <= y[i + 1])) return(i + 1) else return(numeric(0))
  })
  pks <- unlist(pks)
  #threshold limit:
  if (!is.na(y_min)) {
    pks <- pks[y[pks] >= y_min]
  }
  pks
}

shinyServer(function(input, output) {
  data <- reactiveValues(samples_raw = NULL, samples_prep = NULL)
  ################# New dataset panel #################
  volumes <- shinyFiles::getVolumes()
  shinyDirChoose(input, 'new_dataset_samples_dir', roots = volumes)

  observeEvent(input$new_dataset_load_raw, {
    samples_dir <- parseDirPath(volumes, input$new_dataset_samples_dir)
    if (is.null(samples_dir) || length(samples_dir) == 0) {
      samples_dir <- "C:/Users/RDOllerMSe/projects/CT/sample_samples"
    }
    if (!dir.exists(samples_dir)) {
      return()
    }
    data$samples_raw <-
      nmr_read_samples_dir(samples_dir, pulse_sequence = input$new_dataset_type) %>%
      nmr_interpolate(axis1 = c(min = input$new_dataset_min_chemshift,
                                max = input$new_dataset_max_chemshift,
                                by = input$new_dataset_by_chemshift))
  })

  output$new_dataset_exclusion_numbers <- renderUI({
    members <- as.integer(input$new_dataset_number_exclusions)
    if (members <= 0) {
      return(list())
    }
    do.call(c, lapply(1:members, function(i) {
      list(numericInput(inputId = paste0("new_dataset_exclude_min_", i),
                        label = paste("Exclude minimum", i), value = 0),
           numericInput(inputId = paste0("new_dataset_exclude_max_", i),
                        label = paste("Exclude maximum", i), value = 0))
    }))
  })


  observeEvent(input$new_dataset_excl_and_norm, {
    # get list of exclusions:
    excl_list <- list()
    for (i in 1:input$new_dataset_number_exclusions) {
      excl_min <- input[[paste0("new_dataset_exclude_min_", i)]]
      excl_max <- input[[paste0("new_dataset_exclude_max_", i)]]
      excl_list[[paste0("excl_", i)]] <- c(excl_min, excl_max)
    }

    data$samples_prep <- data$samples_raw %>%
      nmr_exclude_region(exclude = excl_list) %>%
      nmr_normalize(method = input$new_dataset_normalize_method)
  })

  gplot <- reactive({
    gplot <- NULL
    if (!(is.null(data$samples_prep))) {
      gplot <- plot(data$samples_prep, interactive = FALSE)
    } else if (!(is.null(data$samples_raw))) {
      gplot <- plot(data$samples_raw, interactive = FALSE)
    }
  })

  output$all_samples_plot <- output$all_samples_plot2 <- renderPlotly({
    p <- gplot()
    if (is.null(p)) {
      p <- plotly::plot_ly()
    } else {
      p <- plotly::ggplotly(p)
    }
    p
  })

  shinyFileSave(input, 'new_dataset_save', roots = volumes)
  observeEvent(input$new_dataset_save, {
    save_file <- as.character(parseSavePath(volumes, input$new_dataset_save)$datapath)
    if (!is.null(data$samples_prep)) {
      nmr_dataset_save(data$samples_prep, save_file)
      showModal(modalDialog(
        title = "Data saved successfully",
        "Data saved",
        easyClose = TRUE,
        footer = NULL
      ))
    } else if (!is.null(data$samples_raw)) {
      nmr_dataset_save(data$samples_raw, save_file)
      showModal(modalDialog(
        title = "Data saved successfully",
        "Data saved",
        easyClose = TRUE,
        footer = NULL
      ))
    } else {
      message("Nothing to save in: ", save_file)
    }
    })
  #####################################################
  ############ Load dataset panel #####################
  #####################################################
  shinyFileChoose(input, 'load_dataset_existing', roots = volumes)
  observeEvent(input$load_dataset_existing, {
    load_file <- as.character(parseFilePaths(volumes, input$load_dataset_existing)$datapath)
    data$samples_prep <- data$samples_raw <- nmr_dataset_load(load_file)
  })
  ###################################################
  ############ Peak detection plot ##################
  ###################################################
  detected_peak_list <- reactive({
    #data frame with injection id, chemshift, intensity, samplerow, chemshift_idx
    inj <- nmr_get_metadata(data$samples_prep, columns = c("injection_id"))$injection_id
    dataset <- data$samples_prep$data_1r
    chemshifts <- data$samples_prep$axis[[1]]
    peak_found <- lapply(1:data$samples_prep$num_samples, function(samplerow) {
      indices <- find_peaks(chemshifts, dataset[samplerow,],
                            m = input$peak_detection_m, y_min = input$peak_detection_y_min)
      chemshift <- chemshifts[indices]
      intensity <- dataset[samplerow, indices]
      tibble(injection_id = rep(as.character(inj[samplerow]), length(indices)),
             chemshift = chemshift,
             intensity = intensity,
             samplerow = rep(samplerow, length(indices)),
             chemshift_idx = indices)
    })
    output <- dplyr::bind_rows(peak_found)
    output
  })

  output$peak_detection_plot <- renderPlotly({
    gplt <- gplot()
    gplt <- gplt + geom_point(data = detected_peak_list(), aes(x = chemshift, y = intensity), color = "red")
    if (is.null(gplt)) {
      gplt <- plotly::plot_ly()
    } else {
      gplt <- plotly::ggplotly(gplt)
    }
    gplt
  })
})
