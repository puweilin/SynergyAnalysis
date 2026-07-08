make_diff <- function(log2FC, gene_id = "g1", gene_name = "Gene1",
                      Pvalue = 0.001, Qvalue = 0.9, fpkm = TRUE) {
  out <- data.frame(
    gene_id = gene_id,
    gene_name = gene_name,
    log2FC = log2FC,
    Pvalue = Pvalue,
    Qvalue = Qvalue,
    updown = ifelse(log2FC >= 0, "UP", "DOWN"),
    stringsAsFactors = FALSE
  )
  if (isTRUE(fpkm)) {
    out$sample1_FPKM <- 10
  }
  out
}

make_raw_p_inputs <- function(fpkm = TRUE) {
  list(
    c_vs_nt = make_diff(3, Pvalue = 0.001, Qvalue = 0.9, fpkm = fpkm),
    a_vs_nt = make_diff(0.5, Pvalue = 0.001, Qvalue = 0.9, fpkm = fpkm),
    b_vs_nt = make_diff(0.5, Pvalue = 0.001, Qvalue = 0.9, fpkm = fpkm),
    c_vs_a  = make_diff(2, Pvalue = 0.001, Qvalue = 0.9, fpkm = fpkm),
    c_vs_b  = make_diff(2, Pvalue = 0.001, Qvalue = 0.9, fpkm = fpkm)
  )
}

test_that("raw P-value mode is propagated to volcano data", {
  skip_if_not_installed("ggrepel")

  res <- suppressMessages(calculate_synergy(
    make_raw_p_inputs(),
    use_qvalue = FALSE,
    qc = NULL
  ))

  p <- plot_synergy_volcano(res, highlight_label = 1)
  built <- ggplot2::ggplot_build(p)

  expect_equal(res$summary$pval_column_used, "Pvalue")
  expect_equal(nrow(res$synergy_up), 1)
  expect_equal(built$data[[1]]$y[1], -log10(res$merged_all$Pvalue_c_vs_nt))
})

test_that("infinite log2FC values are removed by QC", {
  inputs <- make_raw_p_inputs()
  inputs$c_vs_a$log2FC <- Inf

  res <- suppressMessages(calculate_synergy(
    inputs,
    use_qvalue = FALSE,
    qc = synergy_qc_defaults()
  ))

  expect_equal(nrow(res$merged_all), 0)
  expect_equal(nrow(res$synergy_up), 0)
  expect_equal(
    res$summary$qc_log$dropped[res$summary$qc_log$check == "Unstable log2FC"],
    1
  )
})

test_that("data.frame inputs without gene_name are standardized", {
  inputs <- lapply(make_raw_p_inputs(), function(x) {
    x$gene_name <- NULL
    x
  })

  res <- suppressMessages(calculate_synergy(
    inputs,
    use_qvalue = FALSE,
    qc = NULL
  ))

  expect_equal(nrow(res$synergy_up), 1)
  expect_equal(res$merged_all$gene_name[1], "g1")
})

test_that("heatmap handles one gene and honors FPKM mode", {
  skip_if_not_installed("pheatmap")

  plot_file <- tempfile(fileext = ".pdf")
  grDevices::pdf(plot_file)
  on.exit(grDevices::dev.off(), add = TRUE)

  res <- suppressMessages(calculate_synergy(
    make_raw_p_inputs(),
    use_qvalue = FALSE,
    qc = NULL
  ))

  expect_s3_class(plot_synergy_heatmap(res), "pheatmap")

  res_no_fc <- res
  fc_cols <- grep("^log2FC_", colnames(res_no_fc$merged_all), value = TRUE)
  res_no_fc$merged_all[fc_cols] <- NA_real_
  expect_s3_class(plot_synergy_heatmap(res_no_fc, use_fpkm = TRUE), "pheatmap")

  res_no_fpkm <- res
  res_no_fpkm$fpkm_all <- NULL
  expect_warning(
    plot_synergy_heatmap(res_no_fpkm, use_fpkm = TRUE),
    "falling back to log2FC"
  )
})

test_that("HTML report renders with a single synergy gene", {
  skip_if_not_installed("rmarkdown")
  skip_if_not_installed("knitr")
  skip_if_not_installed("ggrepel")
  skip_if_not_installed("pheatmap")
  skip_if_not(rmarkdown::pandoc_available(), "Pandoc is not available")

  res <- suppressMessages(calculate_synergy(
    make_raw_p_inputs(),
    use_qvalue = FALSE,
    qc = NULL
  ))
  out <- tempfile(fileext = ".html")

  expect_message(render_synergy_report(res, out), "Report generated")
  expect_true(file.exists(out))
  expect_gt(file.info(out)$size, 0)
})
