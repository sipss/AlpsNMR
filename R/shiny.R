
#' Run interactive package example
#' @export
NIHSnmr_interactive <- function() {
  appDir <- system.file("shiny-examples", "interactive", package = "NIHSnmr")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `NIHSnmr`.", call. = FALSE)
  }
  if (!requireNamespace("shiny", quietly = TRUE)) {
    stop("shiny needed for this function to work. Please install it.",
         call. = FALSE)
  }
  shiny::runApp(appDir, display.mode = "normal")
}
