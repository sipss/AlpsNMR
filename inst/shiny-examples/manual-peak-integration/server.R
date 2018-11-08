#
# This is the server logic of a Shiny web application. You can run the 
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

library(shiny)
library(NIHSnmr)
library(DT)
library(dplyr)

complete_peak_info <- function(xy, peak_start, peak_end, sample_id) {
  peak_range <- range(c(peak_start, peak_end))
  peak_start <- peak_range[1]
  peak_end <- peak_range[2]
  
  xy_region <- dplyr::filter(xy, chemshift >= (!!peak_start), chemshift <= (!!peak_end)) %>%
    dplyr::arrange(chemshift)
  
  xy_apex <- dplyr::top_n(xy_region, n = 1, wt = intensity)

  peak_apex <- xy_apex[["chemshift"]][1]
  peak_intensity <- xy_apex[["intensity"]][1]
  
  peak_area_raw <- sum(xy_region[["intensity"]])
  
  peak_start_int <- head(xy_region[["intensity"]], n = 1)
  peak_end_int <- tail(xy_region[["intensity"]], n = 1)
  
  # ax+b
  a <- (peak_end_int - peak_start_int)/(peak_end - peak_start)
  b <- peak_end_int - a * peak_end
  
  xy_region[["baseline"]] <- pmin(a*xy_region[["chemshift"]] + b, xy_region[["intensity"]])
  
  peak_area_clean <- peak_area_raw - sum(xy_region[["baseline"]])
  
  list(peak = data.frame(sample_id = sample_id,
                         peak_start = peak_start,
                         peak_end = peak_end,
                         peak_apex = peak_apex,
                         peak_intensity = peak_intensity,
                         peak_area = peak_area_raw,
                         peak_area_clean = peak_area_clean,
                         stringsAsFactors = FALSE),
       xy_region = xy_region)
}

plotly_event <- function(event = c("plotly_hover", "plotly_click",
                                   "plotly_selected", "plotly_relayout"),
                         source = "A",
                         session = shiny::getDefaultReactiveDomain()) {
  src <- sprintf(".clientValue-%s-%s", event[1], source)
  session$rootScope()$input[[src]]
}


# Define server logic required to draw a histogram
shinyServer(function(input, output, session) {
  
  r_sample_idx <- reactiveVal()
  clicked_peak <- reactiveValues(peak_start = NA)
  peak_table <- reactiveVal(value = NULL)
  
  
  r_dataset <- reactive({
    validate(need(!is.null(input$dataset$datapath), "Please load a data set (r_dataset)"))
    readRDS(input$dataset$datapath)
  })

  observeEvent(input$dataset, {
    r_sample_idx(1)
    updateNumericInput(session, "spec", value = 1, min = 1, max = r_dataset()$num_samples, step = 1)
  })
  
  observeEvent(input$spec, {
    dataset <- r_dataset()
    validate(need(!is.null(dataset), "Please load a data set"))
    validate(need(input$spec >= 1 && input$spec <= dataset$num_samples,
             sprintf("Invalid sample index use integer in [%d, %d]", 1, dataset$num_samples)))
    r_sample_idx(input$spec)
  })
  
  spectrumxy <- reactive({
    dataset <- r_dataset()
    validate(need(!is.null(dataset), "Please load a data set"))
    data.frame(chemshift = dataset$axis[[1]],
               intensity = dataset$data_1r[r_sample_idx(),, drop = TRUE])
  })
  
  observeEvent(input$next_win, {
    updateNumericInput(session, "window", value = input$window - 0.5*input$window_width)
  })

  observeEvent(input$prev_win, {
    updateNumericInput(session, "window", value = input$window + 0.5*input$window_width)
  })
  
  spectrum_plot <- reactiveVal()
  spectrum_plot_limits <- reactiveVal()
  spectrum_plot_limits_peaks <- reactiveVal()
  
  # Create the plot when a new sample is chosen
  observe({
    xy <- spectrumxy()
    validate(need(!is.null(xy), "Please load a data set and a sample"))
    
    #spectrum_plot(ggplot() + geom_line(data = xy, aes(x = chemshift, y = intensity)))
    
    spectrum_plot(plotly::plot_ly(data = xy, x = ~chemshift, source = "spectrum") %>%
      plotly::add_lines(y = ~intensity))
  })
  
  # Edit the plot x limits when a new window is chosen
  observe({
    plt <- spectrum_plot()
    validate(need(!is.null(plt), "Please load a data set and a sample"))
    chemshift_max <- input$window + input$window_width/2
    chemshift_min <- input$window - input$window_width/2
    #spectrum_plot_limits(plt + xlim(chemshift_max, chemshift_min))
    spectrum_plot_limits(plt %>%
              plotly::layout(xaxis = list(range = c(chemshift_max, chemshift_min)),
                             yaxis = list(autorange = TRUE)))
  })
  
  observe({
    spectrum_plot_limits_peaks(spectrum_plot_limits())
  })
  
  
  guide_msg <- reactiveVal("Load a sample and click on the left boundary of a peak")
  
  observeEvent(input$cancel_peak, {
    clicked_peak$peak_start <- NA
    guide_msg("Peak cancelled. Click on the left boundary of the next peak")
  })
  
  observeEvent(plotly_event("plotly_click", source = "spectrum"), {
    xy <- spectrumxy()
    have_click <- plotly::event_data(event = "plotly_click", source = "spectrum")
    
    if (is.null(have_click)) {
      return()
    }
    
    plt <- spectrum_plot_limits_peaks()
    if (is.na(clicked_peak$peak_start)) {
      # New peak
      clicked_peak$peak_start <- have_click[["x"]]
      guide_msg("Click on the right boundary of the peak")
      shinyjs::enable("cancel_peak")
      plt <- plt %>%
        plotly::add_trace(x = ~x, y = ~y, data = have_click,
                          showlegend = FALSE, type = "scatter",
                          mode = "markers",
                          marker = list(color = "red"))
    } else {
      # Finish the peak:
      shinyjs::disable("cancel_peak")
      clicked_peak$peak_end <- have_click[["x"]]
      full_peak_info <- complete_peak_info(xy, clicked_peak$peak_start, clicked_peak$peak_end, isolate(r_sample_idx()))
      plt <- plt %>%
        plotly::add_trace(x = ~x, y = ~y, data = have_click,
                          showlegend = FALSE, type = "scatter",
                          mode = "markers",
                          marker = list(color = "red")) %>%
        plotly::add_trace(x = ~chemshift, y = ~intensity, data = full_peak_info$xy_region,
                          showlegend = FALSE, type = "scatter",
                          mode = "lines",
                          line = list(color = 'rgba(0,100,80,1)')) %>%
        plotly::add_trace(x = ~chemshift, y = ~baseline, data = full_peak_info$xy_region,
                          showlegend = FALSE, type = "scatter",
                          fill = "tonexty",
                          fillcolor = 'rgba(0,100,80,0.2)',
                          line = list(color = 'rgba(0,100,80,1)'),
                          mode = "lines")
      print(full_peak_info)
      clicked_peak$peak_start <- NA
      peak_table(dplyr::bind_rows(isolate(peak_table()), full_peak_info[["peak"]]))
      guide_msg("Peak added. Click on the left boundary of the next peak")
    }
    spectrum_plot_limits_peaks(plt)
  })
  
  output$tests <- renderText({
    guide_msg()
  })
  # output$tests <- renderPrint({
  #   xy <- isolate(spectrumxy())
  #   validate(need(!is.null(xy), "Please load a data set and a sample"))
  #   
  #   have_click <- plotly::event_data(event = "plotly_click", source = "spectrum")
  #   validate(need(!is.null(have_click), "Click on the left boundary of a peak (double click to cancel)"))
  # 
  #   if (is.null(have_click)) {
  #     return()
  #   }
  #   
  #   message(clicked_peak$peak_start)
  #   if (is.na(clicked_peak$peak_start)) {
  #     # New peak
  #     clicked_peak$peak_start <- have_click[["x"]]
  #     return("Click the peak end (double click to cancel)")
  #   } else {
  #     # Finish the peak:
  #     clicked_peak$peak_end <- have_click[["x"]]
  #     full_peak_info <- complete_peak_info(xy, clicked_peak$peak_start, clicked_peak$peak_end, isolate(r_sample_idx()))
  #     print(full_peak_info)
  #     clicked_peak$peak_start <- NA
  #     peak_table(dplyr::bind_rows(isolate(peak_table()), full_peak_info))
  #     return("Peak added. Click on the left boundary of the next peak (double click to cancel)")
  #   }
  # })
  
  observeEvent(input$delete_rows, {
    if (!is.null(input$peak_table)) {
      pt <- isolate(peak_table())
      peak_table(pt[-as.numeric(input$peak_table_rows_selected),,drop = FALSE])
    }
  })
  
  output$peak_table <- renderDataTable({
    peak_table()
  })

  output$spectrumPlot <- plotly::renderPlotly({
    plt <- spectrum_plot_limits_peaks()
    validate(need(!is.null(plt), "Please load a data set and a sample"))
    plt
  })
})
