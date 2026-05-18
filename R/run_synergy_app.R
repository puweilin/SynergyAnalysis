# run_synergy_app.R — Launch the bundled Shiny app

#' Launch the Synergy Analysis Shiny app
#'
#' Convenience wrapper that locates the bundled Shiny app inside the installed
#' package and starts it. Useful after installing via
#' \code{devtools::install_github("puweilin/SynergyAnalysis")}.
#'
#' @param launch.browser Logical; open the app in an external browser
#'   (\code{TRUE}, default) or the RStudio viewer (\code{FALSE}).
#' @param ... Additional arguments passed to \code{\link[shiny]{runApp}}.
#'
#' @examples
#' \dontrun{
#' run_synergy_app()
#' }
#'
#' @export
run_synergy_app <- function(launch.browser = TRUE, ...) {
  if (!requireNamespace("shiny", quietly = TRUE)) {
    stop("Package 'shiny' is required to launch the app. Install with: install.packages('shiny')")
  }
  app_dir <- system.file("shiny", package = "SynergyAnalysis")
  if (!nzchar(app_dir)) {
    stop("Could not locate the bundled Shiny app. Reinstall the package.")
  }
  shiny::runApp(app_dir, launch.browser = launch.browser, ...)
}
