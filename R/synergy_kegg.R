# synergy_kegg.R - Local KEGG annotation cache with online fallback

#' Get KEGG pathway annotation, preferring a local cache
#'
#' Returns the pathway-to-gene and pathway-to-name mappings for an organism.
#' On first use the data is downloaded with
#' \code{clusterProfiler::download_KEGG()} and stored in the user's cache
#' directory (\code{tools::R_user_dir("SynergyAnalysis", "cache")}); subsequent
#' calls within \code{max_age_days} use the cached copy. If the cache is older
#' than \code{max_age_days}, a fresh download is attempted; if that download
#' fails (e.g. no internet), the stale local copy is returned with a warning
#' so analysis can still proceed.
#'
#' @param organism Three-letter KEGG organism code, e.g. \code{"hsa"} for human.
#' @param max_age_days Refresh the cache once it is older than this. Default 90.
#' @param force_refresh If \code{TRUE}, ignore any local copy and re-download.
#'
#' @return A list with elements \code{KEGGPATHID2EXTID} (term-to-gene data
#'   frame) and \code{KEGGPATHID2NAME} (term-to-name data frame), suitable for
#'   passing to
#'   \code{clusterProfiler::enricher(TERM2GENE = ., TERM2NAME = .)}. The
#'   returned object has attributes \code{cache_age_days} and \code{source}
#'   (one of \code{"local"}, \code{"online"}, or \code{"local-stale"}).
#'
#' @examples
#' \dontrun{
#'   kegg <- get_kegg_data_cached("hsa")
#'   attr(kegg, "source")
#'   # Force a re-download:
#'   kegg <- get_kegg_data_cached("hsa", force_refresh = TRUE)
#' }
#'
#' @export
get_kegg_data_cached <- function(organism = "hsa",
                                 max_age_days = 90,
                                 force_refresh = FALSE) {
  cache_dir <- tools::R_user_dir("SynergyAnalysis", which = "cache")
  if (!dir.exists(cache_dir)) {
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  }
  cache_file <- file.path(cache_dir, sprintf("kegg_%s.rds", organism))

  cache_exists <- file.exists(cache_file)
  cache_age_days <- if (cache_exists) {
    as.numeric(difftime(Sys.time(), file.mtime(cache_file), units = "days"))
  } else {
    NA_real_
  }
  is_fresh <- cache_exists && !is.na(cache_age_days) &&
              cache_age_days < max_age_days

  if (!force_refresh && is_fresh) {
    data <- readRDS(cache_file)
    attr(data, "cache_age_days") <- cache_age_days
    attr(data, "source") <- "local"
    return(data)
  }

  if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
    stop("Package 'clusterProfiler' is required to download KEGG data.")
  }

  message(sprintf("Downloading KEGG annotation for organism '%s' ...",
                  organism))
  fresh <- tryCatch(
    clusterProfiler::download_KEGG(organism),
    error = function(e) {
      message(sprintf("KEGG download failed: %s", conditionMessage(e)))
      NULL
    }
  )

  if (!is.null(fresh)) {
    saveRDS(fresh, cache_file)
    attr(fresh, "cache_age_days") <- 0
    attr(fresh, "source") <- "online"
    return(fresh)
  }

  if (cache_exists) {
    warning(sprintf(
      "Could not refresh KEGG data; using stale local cache (age: %.0f days).",
      cache_age_days
    ))
    data <- readRDS(cache_file)
    attr(data, "cache_age_days") <- cache_age_days
    attr(data, "source") <- "local-stale"
    return(data)
  }

  stop("Cannot fetch KEGG data: download failed and no local cache available. ",
       "Check your internet connection and retry.")
}

#' Path to the local KEGG cache file
#'
#' @param organism KEGG organism code.
#' @return Absolute path to the RDS file used by \code{get_kegg_data_cached}.
#'   The file may not yet exist.
#' @export
kegg_cache_path <- function(organism = "hsa") {
  cache_dir <- tools::R_user_dir("SynergyAnalysis", which = "cache")
  file.path(cache_dir, sprintf("kegg_%s.rds", organism))
}
