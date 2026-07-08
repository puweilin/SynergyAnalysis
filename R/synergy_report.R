# synergy_report.R - Automated report generation for synergy results

#' Render an HTML report for synergy analysis results
#'
#' Generates a self-contained HTML report with embedded tables and figures.
#'
#' @param synergy_res A synergy_result object
#' @param output_path Output file path (should end in .html or .pdf)
#' @param title Report title
#' @return Invisible path to the generated report
#' @export
render_synergy_report <- function(synergy_res,
                                   output_path = "synergy_report.html",
                                   title = "Synergy Analysis Report") {

  if (!requireNamespace("rmarkdown", quietly = TRUE)) {
    stop("Package 'rmarkdown' is required for report generation")
  }

  tmp_rmd <- tempfile(fileext = ".Rmd")
  rmd_content <- generate_report_rmd(synergy_res, title)
  writeLines(rmd_content, tmp_rmd)

  rmarkdown::render(tmp_rmd, output_file = basename(output_path),
                    output_dir = dirname(normalizePath(output_path, mustWork = FALSE)),
                    quiet = TRUE)

  message("Report generated: ", output_path)
  invisible(output_path)
}

# Build the R Markdown content with all values hardcoded
generate_report_rmd <- function(synergy_res, title) {
  s <- synergy_res$summary
  lbl <- s$labels

  lc  <- lbl["c"]; la <- lbl["a"]; lb <- lbl["b"]; ln <- lbl["nt"]
  lc_nt <- paste0("log2FC_", lc, " vs ", ln)
  la_nt <- paste0("log2FC_", la, " vs ", ln)
  lb_nt <- paste0("log2FC_", lb, " vs ", ln)

  crit_text <- if (identical(s$mode, "relaxed")) {
    "the relaxed synergy criteria (C vs NT significant + magnitude additivity)"
  } else {
    "all 4 strict synergy criteria"
  }

  # Significance column actually used for the call (Qvalue or Pvalue), so the
  # report table matches the filtering rather than always showing Qvalue.
  pcol     <- s$pval_column_used            # "Qvalue" or "Pvalue"
  pcol_c   <- paste0(pcol, "_c_vs_nt")      # data column, e.g. "Pvalue_c_vs_nt"
  pcol_lab <- paste0(pcol, "_CvsNT")        # display header

  c(
    "---",
    paste0('title: "', title, '"'),
    "date: '`r Sys.Date()`'",
    "output: html_document",
    "---",
    "",
    "```{r setup, include=FALSE}",
    "knitr::opts_chunk$set(echo = FALSE, warning = FALSE, message = FALSE)",
    "```",
    "",
    "## Summary",
    "",
    "| Metric | Value |",
    "|--------|-------|",
    paste0("| Total genes tested | ", s$n_total_genes, " |"),
    paste0("| Synergistic UP genes | ", s$n_synergy_up, " |"),
    paste0("| Synergistic DOWN genes | ", s$n_synergy_down, " |"),
    paste0("| Significance column | ", pcol, " |"),
    paste0("| P-value cutoff | ", s$p_cutoff, " |"),
    paste0("| |log2FC| cutoff | ", s$fc_cutoff, " |"),
    paste0("| Group labels | NT=", ln, ", A=", la, ", B=", lb, ", C=", lc, " |"),
    "",
    "## Volcano Plot",
    "",
    "```{r volcano, fig.width=9, fig.height=7}",
    "print(plot_synergy_volcano(synergy_res))",
    "```",
    "",
    "## Synergy UP Genes",
    "",
    paste0("*", s$n_synergy_up, " genes meet ", crit_text, " for up-regulation.*"),
    "",
    "```{r up_table}",
    "if (nrow(synergy_res$synergy_up) > 0) {",
    "  up_show <- synergy_res$synergy_up[, c('gene_name', 'log2FC_c_vs_nt', 'log2FC_a_vs_nt',",
    paste0("    'log2FC_b_vs_nt', 'Increase_c_vs_nt', 'Increase_sum_ab', '", pcol_c, "')]"),
    paste0("  colnames(up_show) <- c('Gene', '", lc_nt, "', '", la_nt, "', '", lb_nt, "',"),
    paste0("    'Increase_CvsNT', 'Increase_Sum_AB', '", pcol_lab, "')"),
    "  knitr::kable(up_show, digits = 3, format = 'html')",
    "} else {",
    "  cat('No synergistic UP genes found.')",
    "}",
    "```",
    "",
    "```{r up_bar, fig.width=10, fig.height=8}",
    "if (nrow(synergy_res$synergy_up) > 0) print(plot_synergy_contrib(synergy_res, 20, 'up'))",
    "```",
    "",
    "## Synergy DOWN Genes",
    "",
    paste0("*", s$n_synergy_down, " genes meet ", crit_text, " for down-regulation.*"),
    "",
    "```{r down_table}",
    "if (nrow(synergy_res$synergy_down) > 0) {",
    "  down_show <- synergy_res$synergy_down[, c('gene_name', 'log2FC_c_vs_nt', 'log2FC_a_vs_nt',",
    paste0("    'log2FC_b_vs_nt', 'Decrease_c_vs_nt', 'Decrease_sum_ab', '", pcol_c, "')]"),
    paste0("  colnames(down_show) <- c('Gene', '", lc_nt, "', '", la_nt, "', '", lb_nt, "',"),
    paste0("    'Decrease_CvsNT', 'Decrease_Sum_AB', '", pcol_lab, "')"),
    "  knitr::kable(down_show, digits = 3, format = 'html')",
    "} else {",
    "  cat('No synergistic DOWN genes found.')",
    "}",
    "```",
    "",
    "```{r down_bar, fig.width=10, fig.height=8}",
    "if (nrow(synergy_res$synergy_down) > 0) print(plot_synergy_contrib(synergy_res, 20, 'down'))",
    "```",
    "",
    "## Heatmap",
    "",
    "```{r heatmap, fig.width=10, fig.height=10}",
    "print(plot_synergy_heatmap(synergy_res, n_top = 20))",
    "```",
    "",
    "## Session Info",
    "",
    "```{r session}",
    "sessionInfo()",
    "```"
  )
}
