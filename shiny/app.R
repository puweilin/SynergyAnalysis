# shiny/app.R — thin launcher.
#
# The application itself lives in inst/shiny/app.R (the copy bundled into the
# installed package). Keeping a single source of truth avoids the two files
# drifting apart. This shim just locates that canonical file and hands it to
# Shiny, so `shiny::runApp("shiny/")` runs exactly the same app as
# SynergyAnalysis::run_synergy_app().

local({
  # Search upward from the app's own location and the working directory for the
  # canonical app file, so this works regardless of where runApp() is invoked.
  find_up <- function(rel) {
    starts <- unique(c(
      tryCatch(dirname(normalizePath(sys.frame(1)$ofile, winslash = "/")),
               error = function(e) NULL),
      normalizePath(getwd(), winslash = "/", mustWork = FALSE)
    ))
    for (s in Filter(Negate(is.null), starts)) {
      d <- s
      for (i in 1:6) {
        cand <- file.path(d, rel)
        if (file.exists(cand)) return(cand)
        parent <- dirname(d)
        if (identical(parent, d)) break
        d <- parent
      }
    }
    NULL
  }

  app_file <- find_up(file.path("inst", "shiny", "app.R"))
  if (is.null(app_file)) {
    stop("Could not locate inst/shiny/app.R. Launch from within a clone of the ",
         "SynergyAnalysis repository, or install the package and call ",
         "SynergyAnalysis::run_synergy_app().")
  }
  shiny::shinyAppFile(app_file)
})
