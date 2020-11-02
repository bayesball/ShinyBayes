NormalNormal <- function() {
  appDir <- system.file("shiny-examples", "NormalNormal",
                        package = "ShinyBayes")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `TeachBayes`.", call. = FALSE)
  }

  shiny::runApp(appDir, display.mode = "normal")
}
