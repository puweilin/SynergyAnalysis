# synergy_core.R - Synergy calculation engine

#' Calculate synergistic gene effects
#'
#' Two synergy modes:
#'   * "strict"  - all 4 criteria must hold (the original definition):
#'       1-3. C vs NT, C vs A, C vs B all significant at p < p_cutoff,
#'            with the correct direction in every comparison.
#'       4.   Increase(CvsNT) > Increase(AvsNT) + Increase(BvsNT)  [UP]
#'            Decrease(CvsNT) > Decrease(AvsNT) + Decrease(BvsNT)  [DOWN]
#'   * "relaxed" - only requires C vs NT to be statistically significant,
#'                 plus the magnitude-additivity criterion (4) above.
#'                 Useful when individual A or B effects are also strong:
#'                 statistical detection of C > A or C > B is hard, but
#'                 the fold-change addition test still demonstrates
#'                 supra-additivity.
#'
#' Where: Increase = FC - 1, Decrease = 1 - FC, FC = 2^log2FC
#'
#' @param results_list A named list of 5 data.frames or file paths:
#'   c_vs_nt, a_vs_nt, b_vs_nt, c_vs_a, c_vs_b
#' @param p_cutoff P-value threshold (default 0.05). Applied to Qvalue column.
#' @param fc_cutoff Minimum |log2FC| for C vs NT (default 0, no filter).
#' @param use_qvalue Use Qvalue (TRUE, default) or Pvalue (FALSE) for significance.
#' @param mode "strict" (default; all 4 criteria) or "relaxed" (only C vs NT
#'   significance + magnitude criterion).
#' @param qc Either a list of options for the gene QC pre-filter (see
#'   \code{synergy_qc_defaults()}) or \code{NULL} to disable QC entirely.
#'   The QC step removes obvious artefacts (low / sporadic expression and
#'   genes with extreme, unstable log2FC values) before the synergy criteria
#'   are evaluated.
#' @param labels Display names for groups: c(nt="NT", a="A", b="B", c="C")
#' @return An S3 object of class "synergy_result"
#' @export
calculate_synergy <- function(results_list,
                               p_cutoff = 0.05,
                               fc_cutoff = 0,
                               use_qvalue = TRUE,
                               mode = c("strict", "relaxed"),
                               qc = synergy_qc_defaults(),
                               labels = c(nt = "NT", a = "A", b = "B", c = "C")) {

  mode <- match.arg(mode)

  # Validate and read inputs
  files <- validate_synergy_inputs(results_list)

  pval_col <- if (use_qvalue) "Qvalue" else "Pvalue"

  # Select and rename key columns for each comparison
  extract_cols <- function(df, comp_name) {
    df_sub <- df[, c("gene_id", "gene_name", "log2FC", "Pvalue", "Qvalue", "updown")]
    colnames(df_sub) <- c("gene_id", "gene_name",
                          paste0("log2FC_", comp_name),
                          paste0("Pvalue_", comp_name),
                          paste0("Qvalue_", comp_name),
                          paste0("updown_", comp_name))
    df_sub
  }

  a_nt  <- extract_cols(files$a_vs_nt,  "a_vs_nt")
  b_nt  <- extract_cols(files$b_vs_nt,  "b_vs_nt")
  c_nt  <- extract_cols(files$c_vs_nt,  "c_vs_nt")
  c_a   <- extract_cols(files$c_vs_a,   "c_vs_a")
  c_b   <- extract_cols(files$c_vs_b,   "c_vs_b")

  # Collect per-sample FPKM matrix across all 5 files, deduplicating samples.
  # edgeR `_all.xls` files include `<sample>_FPKM` columns; the same sample may
  # appear in multiple comparison files but the FPKM value is identical.
  collect_fpkm <- function(df) {
    fpkm_cols <- grep("_FPKM$", colnames(df), value = TRUE)
    if (length(fpkm_cols) == 0) return(NULL)
    df[, c("gene_id", fpkm_cols), drop = FALSE]
  }
  fpkm_chunks <- Filter(Negate(is.null), lapply(files, collect_fpkm))
  fpkm_mat <- NULL
  if (length(fpkm_chunks) > 0) {
    fpkm_mat <- Reduce(function(x, y) {
      new_cols <- setdiff(colnames(y), colnames(x))
      if (length(new_cols) == 0) return(x)
      merge(x, y[, c("gene_id", new_cols), drop = FALSE],
            by = "gene_id", all = TRUE)
    }, fpkm_chunks)
  }

  # Merge all by gene_id
  merged <- a_nt
  merged$gene_name <- NULL  # drop, will re-attach from c_nt
  merged <- merge(merged, b_nt[, !colnames(b_nt) %in% "gene_name"], by = "gene_id", all = FALSE)
  merged <- merge(merged, c_nt[, !colnames(c_nt) %in% "gene_name"], by = "gene_id", all = FALSE)
  merged <- merge(merged, c_a[, !colnames(c_a) %in% "gene_name"],   by = "gene_id", all = FALSE)
  merged <- merge(merged, c_b[, !colnames(c_b) %in% "gene_name"],   by = "gene_id", all = FALSE)

  # Re-attach gene_name from c_nt (most comprehensive)
  gene_names <- files$c_vs_nt[, c("gene_id", "gene_name")]
  merged <- merge(gene_names, merged, by = "gene_id", all.y = TRUE)

  n_total <- nrow(merged)
  message(sprintf("Merged: %d genes across 5 comparisons", n_total))

  # ---------- QC pre-filter (drops obvious artefacts) ----------
  qc_out <- apply_gene_qc(merged, fpkm_mat, qc)
  merged <- qc_out$merged
  qc_log <- qc_out$log
  n_after_qc <- nrow(merged)
  if (!is.null(qc) && n_after_qc < n_total) {
    message(sprintf("QC: kept %d of %d genes (%d dropped)",
                    n_after_qc, n_total, n_total - n_after_qc))
  }

  # Compute FC and effect metrics
  merged$FC_c_vs_nt  <- 2 ^ merged$log2FC_c_vs_nt
  merged$FC_a_vs_nt  <- 2 ^ merged$log2FC_a_vs_nt
  merged$FC_b_vs_nt  <- 2 ^ merged$log2FC_b_vs_nt

  merged$Increase_c_vs_nt <- merged$FC_c_vs_nt - 1
  merged$Increase_a_vs_nt <- merged$FC_a_vs_nt - 1
  merged$Increase_b_vs_nt <- merged$FC_b_vs_nt - 1

  merged$Decrease_c_vs_nt <- 1 - merged$FC_c_vs_nt
  merged$Decrease_a_vs_nt <- 1 - merged$FC_a_vs_nt
  merged$Decrease_b_vs_nt <- 1 - merged$FC_b_vs_nt

  merged$Increase_sum_ab <- merged$Increase_a_vs_nt + merged$Increase_b_vs_nt
  merged$Decrease_sum_ab <- merged$Decrease_a_vs_nt + merged$Decrease_b_vs_nt

  # Criterion 1-3 helper: significance in a given comparison
  sig_cols <- c(c_vs_nt = paste0(pval_col, "_c_vs_nt"),
                c_vs_a  = paste0(pval_col, "_c_vs_a"),
                c_vs_b  = paste0(pval_col, "_c_vs_b"))

  # Direction is derived from the sign of log2FC, NOT from the pre-computed
  # `updown` column. The `updown` column in edgeR/DESeq2 output is itself
  # already gated by FDR, which would silently override the user's choice of
  # Pvalue vs Qvalue. log2FC sign is the unbiased direction call.
  fc_cols <- c(c_vs_nt = "log2FC_c_vs_nt",
               c_vs_a  = "log2FC_c_vs_a",
               c_vs_b  = "log2FC_c_vs_b")

  # Pre-filter: handle Inf/NA values
  merged$log2FC_c_vs_nt[is.infinite(merged$log2FC_c_vs_nt)] <- NA
  merged$log2FC_a_vs_nt[is.infinite(merged$log2FC_a_vs_nt)] <- NA
  merged$log2FC_b_vs_nt[is.infinite(merged$log2FC_b_vs_nt)] <- NA

  # ---------- Synergy UP ----------
  if (mode == "strict") {
    up_sig <- (
      merged[[sig_cols["c_vs_nt"]]] < p_cutoff &
      merged[[fc_cols["c_vs_nt"]]]  > 0 &
      merged[[sig_cols["c_vs_a"]]]  < p_cutoff &
      merged[[fc_cols["c_vs_a"]]]   > 0 &
      merged[[sig_cols["c_vs_b"]]]  < p_cutoff &
      merged[[fc_cols["c_vs_b"]]]   > 0
    )
  } else {
    up_sig <- (
      merged[[sig_cols["c_vs_nt"]]] < p_cutoff &
      merged[[fc_cols["c_vs_nt"]]]  > 0
    )
  }
  up_sig[is.na(up_sig)] <- FALSE

  up_fc <- (
    merged$log2FC_c_vs_nt >= fc_cutoff &
    merged$Increase_c_vs_nt > merged$Increase_sum_ab
  )
  up_fc[is.na(up_fc)] <- FALSE

  synergy_up <- merged[up_sig & up_fc, ]
  synergy_up <- synergy_up[order(-synergy_up$log2FC_c_vs_nt), ]

  # ---------- Synergy DOWN ----------
  if (mode == "strict") {
    down_sig <- (
      merged[[sig_cols["c_vs_nt"]]] < p_cutoff &
      merged[[fc_cols["c_vs_nt"]]]  < 0 &
      merged[[sig_cols["c_vs_a"]]]  < p_cutoff &
      merged[[fc_cols["c_vs_a"]]]   < 0 &
      merged[[sig_cols["c_vs_b"]]]  < p_cutoff &
      merged[[fc_cols["c_vs_b"]]]   < 0
    )
  } else {
    down_sig <- (
      merged[[sig_cols["c_vs_nt"]]] < p_cutoff &
      merged[[fc_cols["c_vs_nt"]]]  < 0
    )
  }
  down_sig[is.na(down_sig)] <- FALSE

  down_fc <- (
    merged$log2FC_c_vs_nt <= -fc_cutoff &
    merged$Decrease_c_vs_nt > merged$Decrease_sum_ab
  )
  down_fc[is.na(down_fc)] <- FALSE

  synergy_down <- merged[down_sig & down_fc, ]
  synergy_down <- synergy_down[order(synergy_down$log2FC_c_vs_nt), ]

  # ---------- Summary ----------
  summary <- list(
    n_total_genes    = n_total,
    n_after_qc       = n_after_qc,
    n_qc_dropped     = n_total - n_after_qc,
    n_synergy_up     = nrow(synergy_up),
    n_synergy_down   = nrow(synergy_down),
    p_cutoff         = p_cutoff,
    fc_cutoff        = fc_cutoff,
    use_qvalue       = use_qvalue,
    pval_column_used = pval_col,
    mode             = mode,
    qc               = qc,
    qc_log           = qc_log,
    labels           = labels
  )

  # ---------- Build result object ----------
  res <- list(
    synergy_up    = synergy_up,
    synergy_down  = synergy_down,
    summary       = summary,
    params        = list(p_cutoff = p_cutoff, fc_cutoff = fc_cutoff,
                         use_qvalue = use_qvalue, mode = mode,
                         qc = qc, labels = labels),
    merged_all    = merged,
    fpkm_all      = fpkm_mat
  )
  class(res) <- "synergy_result"

  message(sprintf("Synergy analysis complete: %d UP, %d DOWN synergistic genes",
                  summary$n_synergy_up, summary$n_synergy_down))

  res
}

# ---------- S3 methods ----------

#' @export
print.synergy_result <- function(x, ...) {
  s <- x$summary
  cat("-- Synergy Analysis Result -------------------------------\n")
  cat(sprintf("  Mode                      : %s\n", s$mode))
  cat(sprintf("  Total genes (after merge) : %d\n", s$n_total_genes))
  if (!is.null(s$qc)) {
    cat(sprintf("  After QC                  : %d (dropped %d)\n",
                s$n_after_qc, s$n_qc_dropped))
  }
  cat(sprintf("  P-value threshold (%s): %.3g\n",
              if (s$use_qvalue) "Qvalue" else "Pvalue", s$p_cutoff))
  cat(sprintf("  |log2FC| threshold        : %.2f\n", s$fc_cutoff))
  cat(sprintf("  Synergistic UP   genes    : %d\n", s$n_synergy_up))
  cat(sprintf("  Synergistic DOWN genes    : %d\n", s$n_synergy_down))
  cat("----------------------------------------------------------\n")
  if (s$n_synergy_up > 0) {
    cat("Top UP genes:\n")
    print(head(x$synergy_up[, c("gene_name", "log2FC_c_vs_nt",
          "Increase_c_vs_nt", "Increase_sum_ab")], 10), row.names = FALSE)
  }
  if (s$n_synergy_down > 0) {
    cat("\nTop DOWN genes:\n")
    print(head(x$synergy_down[, c("gene_name", "log2FC_c_vs_nt",
          "Decrease_c_vs_nt", "Decrease_sum_ab")], 10), row.names = FALSE)
  }
  invisible(x)
}

#' @export
summary.synergy_result <- function(object, ...) {
  s <- object$summary
  cat(sprintf("Synergy UP: %d | Synergy DOWN: %d | Total tested: %d\n",
              s$n_synergy_up, s$n_synergy_down, s$n_total_genes))
}

# ---------- Export ----------

#' Export synergy results to a multi-sheet Excel workbook
#'
#' @param synergy_res A synergy_result object
#' @param output_path Path for the output .xlsx file
#' @export
export_synergy_excel <- function(synergy_res, output_path) {
  if (!requireNamespace("openxlsx", quietly = TRUE)) {
    stop("Package 'openxlsx' is required for Excel export. Install with: install.packages('openxlsx')")
  }

  wb <- openxlsx::createWorkbook()

  # Sheet 1: Summary
  openxlsx::addWorksheet(wb, "Summary")
  s <- synergy_res$summary
  summary_rows <- list(
    c("Total genes (after merge)", s$n_total_genes),
    c("After QC",                  if (!is.null(s$qc)) s$n_after_qc else "-"),
    c("QC dropped",                if (!is.null(s$qc)) s$n_qc_dropped else "-"),
    c("Synergistic UP genes",      s$n_synergy_up),
    c("Synergistic DOWN genes",    s$n_synergy_down),
    c("P-value cutoff",            s$p_cutoff),
    c("|log2FC| cutoff",           s$fc_cutoff),
    c("P-value column used",       s$pval_column_used),
    c("Mode",                      s$mode),
    c("Group NT",                  s$labels["nt"]),
    c("Group A",                   s$labels["a"]),
    c("Group B",                   s$labels["b"]),
    c("Group C",                   s$labels["c"])
  )
  if (!is.null(s$qc)) {
    qc <- s$qc
    summary_rows <- c(summary_rows, list(
      c("QC: min FPKM (any sample)",      qc$min_fpkm),
      c("QC: min detection rate",         qc$min_detect_frac),
      c("QC: detection threshold (FPKM)", qc$detect_fpkm),
      c("QC: max |log2FC|",               qc$max_abs_log2fc)
    ))
  }
  summary_df <- do.call(rbind, lapply(summary_rows, function(r)
    data.frame(Metric = r[1], Value = as.character(r[2]),
               stringsAsFactors = FALSE)))
  openxlsx::writeData(wb, "Summary", summary_df)

  if (!is.null(s$qc_log) && nrow(s$qc_log) > 0) {
    openxlsx::addWorksheet(wb, "QC_Log")
    openxlsx::writeData(wb, "QC_Log", s$qc_log)
  }

  # Significance columns follow the column the call used (P vs Q).
  sig_cols <- paste0(if (isTRUE(synergy_res$params$use_qvalue)) "Qvalue_" else "Pvalue_",
                     c("c_vs_nt", "c_vs_a", "c_vs_b"))

  # Sheet 2: Synergy UP genes (select key columns)
  openxlsx::addWorksheet(wb, "Synergy_UP")
  if (nrow(synergy_res$synergy_up) > 0) {
    up_cols <- c("gene_id", "gene_name", "log2FC_c_vs_nt", "log2FC_a_vs_nt",
                 "log2FC_b_vs_nt", "log2FC_c_vs_a", "log2FC_c_vs_b",
                 sig_cols,
                 "FC_c_vs_nt", "FC_a_vs_nt", "FC_b_vs_nt",
                 "Increase_c_vs_nt", "Increase_a_vs_nt", "Increase_b_vs_nt",
                 "Increase_sum_ab")
    openxlsx::writeData(wb, "Synergy_UP", synergy_res$synergy_up[, up_cols])
  }

  # Sheet 3: Synergy DOWN genes
  openxlsx::addWorksheet(wb, "Synergy_DOWN")
  if (nrow(synergy_res$synergy_down) > 0) {
    down_cols <- c("gene_id", "gene_name", "log2FC_c_vs_nt", "log2FC_a_vs_nt",
                   "log2FC_b_vs_nt", "log2FC_c_vs_a", "log2FC_c_vs_b",
                   sig_cols,
                   "FC_c_vs_nt", "FC_a_vs_nt", "FC_b_vs_nt",
                   "Decrease_c_vs_nt", "Decrease_a_vs_nt", "Decrease_b_vs_nt",
                   "Decrease_sum_ab")
    openxlsx::writeData(wb, "Synergy_DOWN", synergy_res$synergy_down[, down_cols])
  }

  # Sheets 4-5: Per-sample FPKM for synergy UP / DOWN genes
  fpkm_mat <- synergy_res$fpkm_all
  if (!is.null(fpkm_mat)) {
    write_fpkm_sheet <- function(sheet_name, gene_df) {
      openxlsx::addWorksheet(wb, sheet_name)
      if (nrow(gene_df) == 0) return(invisible(NULL))
      sub <- fpkm_mat[match(gene_df$gene_id, fpkm_mat$gene_id), , drop = FALSE]
      sub <- cbind(gene_name = gene_df$gene_name, sub)
      # Reorder samples: NT first, then A, B, C - heuristic by prefix from the
      # comparison files isn't reliable here, so keep file-discovery order.
      openxlsx::writeData(wb, sheet_name, sub)
    }
    write_fpkm_sheet("FPKM_Synergy_UP",   synergy_res$synergy_up)
    write_fpkm_sheet("FPKM_Synergy_DOWN", synergy_res$synergy_down)
  }

  openxlsx::saveWorkbook(wb, output_path, overwrite = TRUE)
  message("Synergy Excel exported to: ", output_path)
  invisible(output_path)
}
