# synergy_io.R - Data I/O and validation for synergy analysis

#' Read an edgeR or DESeq2 result file
#'
#' Accepts tab-delimited .xls files produced by edgeR or DESeq2 pipelines.
#' Required columns: gene_id, gene_name, log2FC, Pvalue, Qvalue, updown
#'
#' @param filepath Path to the result file
#' @return A data.frame with standardized columns
#' @export
read_diff_result <- function(filepath) {
  if (!file.exists(filepath)) stop("File not found: ", filepath)

  raw <- read.delim(filepath, header = TRUE, stringsAsFactors = FALSE,
                    check.names = FALSE, quote = "", fill = TRUE)

  standardize_diff_df(raw, basename(filepath))
}

#' Standardize a differential-expression table
#'
#' Shared by \code{read_diff_result} (file inputs) and
#' \code{validate_synergy_inputs} (data.frame inputs) so both paths get the same
#' required-column check, numeric coercion, gene_name fill, and gene_id dedup.
#'
#' @param raw A raw data.frame with at least the required columns.
#' @param source_label Label used in error/warning messages.
#' @return A standardized data.frame.
#' @keywords internal
standardize_diff_df <- function(raw, source_label = "input") {
  required_cols <- c("gene_id", "log2FC", "Pvalue", "Qvalue", "updown")
  missing_cols <- setdiff(required_cols, colnames(raw))
  if (length(missing_cols) > 0) {
    stop("Missing required columns in ", source_label, ": ",
         paste(missing_cols, collapse = ", "))
  }

  out <- raw
  out$gene_id <- as.character(out$gene_id)
  out$log2FC  <- as.numeric(out$log2FC)
  out$Pvalue  <- as.numeric(out$Pvalue)
  out$Qvalue  <- as.numeric(out$Qvalue)
  out$updown  <- as.character(out$updown)

  # Ensure gene_name exists, fill from gene_id if not
  if (!"gene_name" %in% colnames(out)) {
    out$gene_name <- out$gene_id
  } else {
    out$gene_name <- as.character(out$gene_name)
    bad <- is.na(out$gene_name) | out$gene_name == "" | out$gene_name == "-"
    out$gene_name[bad] <- out$gene_id[bad]
  }

  # Drop duplicate gene_id rows (keep first). Downstream merges join on gene_id;
  # duplicates there cause a cartesian row explosion across the 5 comparisons.
  drop_duplicate_gene_ids(out, source_label)
}

#' Drop duplicate gene_id rows, warning if any were removed
#' @param df A data.frame with a gene_id column.
#' @param source_label Label used in the warning message.
#' @return \code{df} with duplicate gene_id rows removed.
#' @keywords internal
drop_duplicate_gene_ids <- function(df, source_label = "input") {
  dup <- duplicated(df$gene_id)
  if (any(dup)) {
    warning(sprintf("%s: dropped %d duplicate gene_id row(s); keeping first occurrence.",
                    source_label, sum(dup)))
    df <- df[!dup, , drop = FALSE]
  }
  df
}

#' Validate a named list of 5 diff result data.frames for synergy analysis
#'
#' @param results_list Named list with elements: c_vs_nt, a_vs_nt, b_vs_nt,
#'   c_vs_a, c_vs_b. Each element is either a file path (character) or a data.frame.
#' @return A named list of 5 data.frames with standardized columns
#' @export
validate_synergy_inputs <- function(results_list) {
  required_names <- c("c_vs_nt", "a_vs_nt", "b_vs_nt", "c_vs_a", "c_vs_b")
  missing_names <- setdiff(required_names, names(results_list))
  if (length(missing_names) > 0) {
    stop("Missing comparison(s): ", paste(missing_names, collapse = ", "))
  }

  # Read or accept data.frames
  processed <- lapply(names(results_list[required_names]), function(nm) {
    x <- results_list[[nm]]
    if (is.character(x)) {
      read_diff_result(x)                 # standardizes + dedups internally
    } else if (is.data.frame(x)) {
      standardize_diff_df(x, nm)          # same standardization for df inputs
    } else {
      stop("Each element must be a file path (character) or data.frame")
    }
  })
  names(processed) <- required_names

  # Check gene overlap
  gene_sets <- lapply(processed, function(df) unique(df$gene_id))
  common_genes <- Reduce(intersect, gene_sets)
  n_common <- length(common_genes)

  if (n_common == 0) {
    stop("No common genes found across the 5 comparisons. Check gene_id format.")
  }

  message(sprintf("Input validation: %d common genes across 5 comparisons", n_common))

  processed
}
