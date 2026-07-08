# SynergyAnalysis-package.R - package-level documentation and imports

#' SynergyAnalysis: identify synergistic genes from transcriptomics data
#'
#' Framework and Shiny app for finding genes where a combination treatment A+B
#' produces an effect exceeding the sum of A alone and B alone.
#'
#' @keywords internal
#' @importFrom grDevices colorRampPalette
#' @importFrom stats setNames
#' @importFrom utils head read.delim
"_PACKAGE"

# Column names referenced without quoting inside ggplot2::aes() in the plotting
# functions. Declaring them keeps R CMD check from flagging "no visible binding
# for global variable" NOTEs.
utils::globalVariables(c(
  "log2FC_c_vs_nt_plot", "nlog10_q", "category", "label",
  "value", "gene", "comp"
))
