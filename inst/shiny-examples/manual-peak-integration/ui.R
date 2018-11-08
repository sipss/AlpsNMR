#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

library(shiny)
library(DT)
library(shinyjs)

# Define UI for application that draws a histogram
shinyUI(fluidPage(
  useShinyjs(),
  # Application title
  titlePanel("Manual peak integration"),
  
  # Sidebar with a slider input for number of bins 
  sidebarLayout(
    sidebarPanel(
      shiny::fileInput("dataset", "NMR Dataset", multiple = FALSE),
      shiny::numericInput("spec", label = "Sample ID", value = 1),
      #splitLayout(
      #  shiny::actionButton("prev_spec", label = NULL, icon = shiny::icon("step-backward")),
      #  shiny::actionButton("next_spec", label = NULL, icon = shiny::icon(name = "step-forward")),
      #  cellWidths = c(32, 32, 32)),
      splitLayout(
          shiny::actionButton("prev_win", label = NULL, icon = shiny::icon("step-backward")),
          shiny::numericInput("window", label = NULL, value = 5),
          shiny::actionButton("next_win", label = NULL, icon = shiny::icon(name = "step-forward")),
          cellWidths = c(32, 96, 32)),
      shiny::sliderInput("window_width", "Chemshift range", min = 0.005, max = 10, value = 1)
    ),
    
    # Show a plot of the generated distribution
    mainPanel(
      plotly::plotlyOutput("spectrumPlot"),
      textOutput("tests"),
      shinyjs::disabled(actionButton("cancel_peak", "Cancel Peak")),
      dataTableOutput("peak_table"),
      splitLayout(
        shiny::actionButton("save_peak_table", label = "Save", icon = shiny::icon("save")),
        shiny::actionButton("delete_rows", label = "Delete selected rows", icon = shiny::icon("delete"))
      )
    )
  )
))
