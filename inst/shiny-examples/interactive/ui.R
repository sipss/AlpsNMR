
# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)
library(shinyFiles)
library(plotly)


shinyUI(pageWithSidebar(
  headerPanel("Import and visualize Bruker NMR data"),
  sidebarPanel(
    conditionalPanel(condition = "input.conditionedPanels==1",
                     helpText("New dataset from NIHS directory"),
                     shinyDirButton(id = "new_dataset_samples_dir", label = "Samples directory",
                                    title = "Choose a directory"),
                     selectInput("new_dataset_type", label = "Sample type",
                                 choices = c("NOESY", "JRES", "CPMG", "DIFF"), selected = "NOESY"),
                     # interpolate
                     numericInput("new_dataset_min_chemshift",
                                  label = "Minimum Chemical Shift", value = 0.4),
                     numericInput("new_dataset_max_chemshift",
                                  label = "Maximum Chemical Shift", value = 10),
                     numericInput("new_dataset_by_chemshift",
                                  label = "Chemical Shift increment", value = 0.0008),
                     actionButton("new_dataset_load_raw", label = "Load raw data"),
                     tags$hr(),
                     # exclude
                     numericInput("new_dataset_number_exclusions", "Number of exclusion ranges", value = 0),
                     uiOutput("new_dataset_exclusion_numbers"),
                     # normalize
                     selectInput("new_dataset_normalize_method", "Normalization method",
                                 choices = c("area", "max", "none"), selected = "none"),
                     actionButton("new_dataset_excl_and_norm", label = "Exclude and Normalize"),
                     shinySaveButton(id = 'new_dataset_save', label = 'Save', title = 'Save as...',
                                     filetype = c("nmrnihs" = "nmrnihs"))
    ),
    conditionalPanel(condition = "input.conditionedPanels==2",
                     helpText("Load dataset"),
                     shinyFilesButton("load_dataset_existing", "Choose a dataset",
                                      "Load nmrnihs dataset", multiple = FALSE),
                     actionButton("load_dataset_load", label = "Load")
    ),
    conditionalPanel(condition = "input.conditionedPanels==3",
                     helpText("Peak detection"),
                     numericInput("peak_detection_m", "Min peak width", min = 1, value = 3),
                     numericInput("peak_detection_y_min", "Minimum peak intensity", value = 5000)
    )
  ),
  mainPanel(
    tabsetPanel(
      tabPanel("New dataset", value = 1,
               plotlyOutput("all_samples_plot")),
      tabPanel("Load dataset", value = 2,
               plotlyOutput("all_samples_plot2")),
      tabPanel("Find Peaks", value = 3,
               plotlyOutput("peak_detection_plot"))
      , id = "conditionedPanels"
    )
  )
))

