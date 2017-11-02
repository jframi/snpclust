#' @export
runmanclust <- function() {
  appDir <- system.file("shiny-apps", "manclust", package = "snpclust")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `snpclust`.", call. = FALSE)
  }

  shiny::runApp(appDir, display.mode = "normal",launch.browser = T)
}
