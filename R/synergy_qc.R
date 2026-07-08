# synergy_qc.R — Upstream gene-level QC for the synergy pipeline

#' Default options for the gene QC pre-filter
#'
#' Returned as a plain list so the Shiny UI can construct it field-by-field.
#' Pass the result (or a modified version) to \code{calculate_synergy(qc = ...)}.
#' Pass \code{NULL} to \code{calculate_synergy(qc = NULL)} to disable QC entirely.
#'
#' Defaults are deliberately permissive — they drop genes that are clearly
#' artefactual without trimming real signal:
#'
#' \describe{
#'   \item{min_fpkm}{Expression threshold. A gene is kept only if its FPKM
#'     reaches this value in at least one sample (i.e. the per-gene maximum
#'     across samples is >= min_fpkm). Default 1.}
#'   \item{min_detect_frac}{Fraction of samples in which FPKM must exceed
#'     \code{detect_fpkm} for the gene to be kept. Default 0.5.}
#'   \item{detect_fpkm}{FPKM value above which a sample counts as "detected".
#'     Default 0.1.}
#'   \item{max_abs_log2fc}{Genes are dropped if any of the 5 pairwise
#'     comparisons reports \code{|log2FC| > max_abs_log2fc}. These extreme
#'     values typically come from a denominator near zero and inflate
#'     synergy calls. Default 10 (i.e. ~1000-fold change in either direction).}
#' }
#'
#' Setting any field to \code{NA} disables that individual check.
#'
#' @return Named list of QC parameters.
#' @export
synergy_qc_defaults <- function() {
  list(
    min_fpkm        = 1,
    min_detect_frac = 0.5,
    detect_fpkm     = 0.1,
    max_abs_log2fc  = 10
  )
}

#' Apply the gene QC pre-filter
#'
#' Internal helper used by \code{calculate_synergy}. Removes genes that
#' (1) never reach a meaningful expression level, (2) are detected only in a
#' minority of samples, or (3) carry an extreme log2FC in any of the 5
#' pairwise comparisons.
#'
#' @param merged Merged data.frame of all 5 comparisons (one row per gene).
#' @param fpkm_mat FPKM matrix as built inside \code{calculate_synergy}, with
#'   \code{gene_id} as a column and one column per sample, or \code{NULL}.
#' @param qc Either the QC options list (see \code{synergy_qc_defaults}) or
#'   \code{NULL} to skip QC.
#'
#' @return A list with elements \code{merged} (the filtered merged table)
#'   and \code{log} (a small data.frame summarising how many genes each
#'   check would have dropped on its own).
#' @keywords internal
apply_gene_qc <- function(merged, fpkm_mat, qc) {
  log_df <- data.frame(
    check = character(),
    threshold = character(),
    dropped = integer(),
    stringsAsFactors = FALSE
  )

  if (is.null(qc)) {
    return(list(merged = merged, log = log_df))
  }

  keep <- rep(TRUE, nrow(merged))

  # Did the user ask for expression-based filtering?
  wants_fpkm_checks <-
    (isTRUE(qc$min_fpkm > 0)        && !is.na(qc$min_fpkm)) ||
    (isTRUE(qc$min_detect_frac > 0) && !is.na(qc$min_detect_frac))
  have_fpkm <- !is.null(fpkm_mat) && nrow(fpkm_mat) > 0 &&
    length(setdiff(colnames(fpkm_mat), "gene_id")) > 0

  # --- FPKM-based checks ---
  if (have_fpkm) {
    sample_cols <- setdiff(colnames(fpkm_mat), "gene_id")
    if (length(sample_cols) > 0) {
      fpkm_sub <- fpkm_mat[match(merged$gene_id, fpkm_mat$gene_id),
                           sample_cols, drop = FALSE]
      fpkm_sub <- as.matrix(fpkm_sub)
      fpkm_sub[is.na(fpkm_sub)] <- 0

      if (isTRUE(qc$min_fpkm > 0) && !is.na(qc$min_fpkm)) {
        max_fpkm <- apply(fpkm_sub, 1, max)
        fail_expr <- max_fpkm < qc$min_fpkm
        log_df <- rbind(log_df, data.frame(
          check = "Min FPKM (any sample)",
          threshold = sprintf(">= %g", qc$min_fpkm),
          dropped = sum(fail_expr & keep),
          stringsAsFactors = FALSE
        ))
        keep <- keep & !fail_expr
      }

      if (isTRUE(qc$min_detect_frac > 0) && !is.na(qc$min_detect_frac)) {
        detect_thr <- if (is.null(qc$detect_fpkm) || is.na(qc$detect_fpkm)) 0.1 else qc$detect_fpkm
        detect_frac <- rowMeans(fpkm_sub > detect_thr)
        fail_det <- detect_frac < qc$min_detect_frac
        log_df <- rbind(log_df, data.frame(
          check = "Detection rate",
          threshold = sprintf(">= %.0f%% samples with FPKM > %g",
                              100 * qc$min_detect_frac, detect_thr),
          dropped = sum(fail_det & keep),
          stringsAsFactors = FALSE
        ))
        keep <- keep & !fail_det
      }
    }
  } else if (wants_fpkm_checks) {
    # User asked for expression QC but the inputs carry no <sample>_FPKM columns
    # (e.g. plain DESeq2 output). Surface this rather than silently skipping it.
    warning("Gene QC: no FPKM columns found in the input files; ",
            "the 'Min FPKM' and 'Detection rate' checks were skipped.")
    log_df <- rbind(log_df, data.frame(
      check = "FPKM checks",
      threshold = "skipped: no _FPKM columns in input",
      dropped = 0L,
      stringsAsFactors = FALSE
    ))
  }

  # --- log2FC sanity check (always applied if requested) ---
  if (isTRUE(qc$max_abs_log2fc > 0) && !is.na(qc$max_abs_log2fc)) {
    lfc_cols <- grep("^log2FC_", colnames(merged), value = TRUE)
    if (length(lfc_cols) > 0) {
      mat <- as.matrix(merged[, lfc_cols, drop = FALSE])
      mat[!is.finite(mat)] <- NA  # treat Inf as extreme
      max_abs <- suppressWarnings(apply(abs(mat), 1, max, na.rm = TRUE))
      max_abs[!is.finite(max_abs)] <- 0  # all-NA rows: don't drop here
      fail_lfc <- max_abs > qc$max_abs_log2fc
      log_df <- rbind(log_df, data.frame(
        check = "Unstable log2FC",
        threshold = sprintf("|log2FC| <= %g in all 5 comparisons", qc$max_abs_log2fc),
        dropped = sum(fail_lfc & keep),
        stringsAsFactors = FALSE
      ))
      keep <- keep & !fail_lfc
    }
  }

  list(merged = merged[keep, , drop = FALSE], log = log_df)
}
