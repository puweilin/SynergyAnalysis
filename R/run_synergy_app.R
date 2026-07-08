# run_synergy_app.R - Launch the bundled Shiny app

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
  # The app attaches these at startup; check them all up front so a missing one
  # gives an actionable message instead of a cryptic error mid-launch.
  needed <- c("shiny", "bslib", "ggplot2", "plotly", "DT")
  missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop("The Shiny app needs these packages: ", paste(missing, collapse = ", "),
         ".\nInstall them with:\n  install.packages(c(",
         paste(sprintf('"%s"', missing), collapse = ", "), "))", call. = FALSE)
  }
  # Enrichment (ORA) additionally needs Bioconductor packages; warn but don't
  # block, since the rest of the app works without them.
  ora_pkgs <- c("clusterProfiler", "org.Hs.eg.db", "enrichplot", "DOSE", "stringr")
  ora_missing <- ora_pkgs[!vapply(ora_pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(ora_missing) > 0) {
    message("Note: the Enrichment tab needs these packages (install to enable): ",
            paste(ora_missing, collapse = ", "), ".")
  }
  app_dir <- system.file("shiny", package = "SynergyAnalysis")
  if (!nzchar(app_dir)) {
    stop("Could not locate the bundled Shiny app. Reinstall the package.")
  }
  shiny::runApp(app_dir, launch.browser = launch.browser, ...)
}

# bslib, plotly, and DT are Imports used by the bundled Shiny app
# (inst/shiny/app.R), not by this package's R code. Reference them here so
# R CMD check records the dependency instead of flagging an unused Import.
# This helper is never called.
#' @noRd
.app_pkg_refs <- function() {
  list(bslib::bs_theme, plotly::plot_ly, DT::datatable)
}
