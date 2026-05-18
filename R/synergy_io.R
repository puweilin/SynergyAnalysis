# synergy_io.R — Data I/O and validation for synergy analysis

#' Read an edgeR or DESeq2 result file
#'
#' Accepts tab-delimited .xls files produced by edgeR or DESeq2 pipelines.
#' Required columns: gene_id, gene_name, log2FC, Pvalue, Qvalue, updown
#'
#' @param filepath Path to the result file
#' @return A data.frame with standardized columns
read_diff_result <- function(filepath) {
  if (!file.exists(filepath)) stop("File not found: ", filepath)

  raw <- read.delim(filepath, header = TRUE, stringsAsFactors = FALSE,
                    check.names = FALSE, quote = "", fill = TRUE)

  required_cols <- c("gene_id", "log2FC", "Pvalue", "Qvalue", "updown")
  missing_cols <- setdiff(required_cols, colnames(raw))
  if (length(missing_cols) > 0) {
    stop("Missing required columns in ", basename(filepath), ": ",
         paste(missing_cols, collapse = ", "))
  }

  # Standardize key columns
  out <- raw
  out$log2FC <- as.numeric(out$log2FC)
  out$Pvalue  <- as.numeric(out$Pvalue)
  out$Qvalue  <- as.numeric(out$Qvalue)
  out$updown  <- as.character(out$updown)

  # Ensure gene_name exists, fill from gene_id if not
  if (!"gene_name" %in% colnames(out)) {
    out$gene_name <- out$gene_id
  } else {
    out$gene_name <- as.character(out$gene_name)
    out$gene_name[is.na(out$gene_name) | out$gene_name == "" | out$gene_name == "-"] <-
      out$gene_id[is.na(out$gene_name) | out$gene_name == "" | out$gene_name == "-"]
  }

  out
}

#' Validate a named list of 5 diff result data.frames for synergy analysis
#'
#' @param results_list Named list with elements: c_vs_nt, a_vs_nt, b_vs_nt,
#'   c_vs_a, c_vs_b. Each element is either a file path (character) or a data.frame.
#' @return A named list of 5 data.frames with standardized columns
validate_synergy_inputs <- function(results_list) {
  required_names <- c("c_vs_nt", "a_vs_nt", "b_vs_nt", "c_vs_a", "c_vs_b")
  missing_names <- setdiff(required_names, names(results_list))
  if (length(missing_names) > 0) {
    stop("Missing comparison(s): ", paste(missing_names, collapse = ", "))
  }

  # Read or accept data.frames
  processed <- lapply(results_list[required_names], function(x) {
    if (is.character(x)) {
      read_diff_result(x)
    } else if (is.data.frame(x)) {
      x
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
