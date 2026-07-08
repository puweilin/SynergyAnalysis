# synergy_plot.R - Publication-quality visualizations for synergy results

# Color palette
SYNERGY_COLORS <- c(
  nt      = "#4C78A8",
  a       = "#76B041",
  b       = "#E45756",
  c       = "#9467BD",
  synergy = "#FF7F0E",
  up      = "#E45756",
  down    = "#4C78A8",
  nonsig  = "#BABABA"
)

#' Volcano plot highlighting synergistic genes
#'
#' @param synergy_res A synergy_result object
#' @param highlight_label Number of top synergy genes to label (default 10)
#' @return A ggplot object
#' @export
plot_synergy_volcano <- function(synergy_res, highlight_label = 10) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package 'ggplot2' is required")
  if (!requireNamespace("ggrepel", quietly = TRUE)) stop("Package 'ggrepel' is required")

  merged <- synergy_res$merged_all

  # Significance axis follows the column the synergy call actually used, so the
  # dashed p_cutoff line and the y-axis stay consistent with the gene set.
  use_q <- isTRUE(synergy_res$params$use_qvalue)
  pcol  <- if (use_q) "Qvalue_c_vs_nt" else "Pvalue_c_vs_nt"
  plab  <- if (use_q) "Q-value" else "P-value"

  merged$category <- "Non-significant"
  merged$category[merged$gene_id %in% synergy_res$synergy_up$gene_id]   <- "Synergy UP"
  merged$category[merged$gene_id %in% synergy_res$synergy_down$gene_id] <- "Synergy DOWN"

  # Remove Inf/-Inf for plotting
  merged$log2FC_c_vs_nt_plot <- merged$log2FC_c_vs_nt
  merged$log2FC_c_vs_nt_plot[is.infinite(merged$log2FC_c_vs_nt_plot)] <- NA
  max_fc <- max(abs(merged$log2FC_c_vs_nt_plot), na.rm = TRUE)
  cap <- min(max_fc + 1, 15)

  merged$log2FC_c_vs_nt_plot <- pmax(-cap, pmin(cap, merged$log2FC_c_vs_nt_plot))
  merged$nlog10_q <- -log10(merged[[pcol]])
  # Q-values that underflow to 0 give -log10 = Inf and would silently drop the
  # most significant genes from the plot. Pin them to the finite maximum.
  finite_max <- suppressWarnings(max(merged$nlog10_q[is.finite(merged$nlog10_q)]))
  if (is.finite(finite_max)) {
    merged$nlog10_q[is.infinite(merged$nlog10_q)] <- finite_max
  }

  # Label top synergy genes
  merged$label <- ""
  up_label <- head(synergy_res$synergy_up, highlight_label)
  down_label <- head(synergy_res$synergy_down, highlight_label)
  label_genes <- c(up_label$gene_id, down_label$gene_id)
  merged$label[merged$gene_id %in% label_genes] <-
    merged$gene_name[merged$gene_id %in% label_genes]

  # unname(): SYNERGY_COLORS["up"] carries the name "up", which would otherwise
  # make the map key "Synergy UP.up" and silently break scale_color_manual().
  col_map <- c("Synergy UP"      = unname(SYNERGY_COLORS["up"]),
               "Synergy DOWN"    = unname(SYNERGY_COLORS["down"]),
               "Non-significant" = unname(SYNERGY_COLORS["nonsig"]))

  p <- ggplot2::ggplot(merged, ggplot2::aes(x = log2FC_c_vs_nt_plot,
                                             y = nlog10_q,
                                             color = category)) +
    ggplot2::geom_point(alpha = 0.5, size = 1.2) +
    ggplot2::scale_color_manual(values = col_map) +
    ggrepel::geom_text_repel(
      ggplot2::aes(label = label),
      size = 3.5, fontface = "italic", max.overlaps = 30,
      box.padding = 0.5, segment.color = "grey50"
    ) +
    ggplot2::geom_hline(yintercept = -log10(synergy_res$params$p_cutoff),
                        linetype = "dashed", color = "grey40") +
    ggplot2::labs(
      x = expression(log[2] * "(Fold Change) - " * C * " vs " * NT),
      y = bquote(-log[10] * "(" * .(plab) * ")"),
      color = "Category",
      title = paste0("Synergy Volcano: ", synergy_res$params$labels["c"],
                     " vs ", synergy_res$params$labels["nt"])
    ) +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(
      legend.position = "bottom",
      panel.grid.minor = ggplot2::element_blank()
    )

  p
}

#' Bar chart showing effect contribution for top synergistic genes
#'
#' For UP genes: shows Increase (FC-1) for CvsNT, AvsNT, BvsNT
#' For DOWN genes: shows Decrease (1-FC)
#'
#' @param synergy_res A synergy_result object
#' @param n Top N genes to show (default 20)
#' @param direction "up" or "down"
#' @return A ggplot object
#' @export
plot_synergy_contrib <- function(synergy_res, n = 20, direction = "up") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package 'ggplot2' is required")

  if (direction == "up") {
    df <- synergy_res$synergy_up
    if (nrow(df) == 0) { message("No synergistic UP genes to plot."); return(NULL) }
    df <- head(df, n)
    metric_cols <- c("Increase_c_vs_nt", "Increase_a_vs_nt", "Increase_b_vs_nt")
    metric_label <- "Increase (FC - 1)"
  } else {
    df <- synergy_res$synergy_down
    if (nrow(df) == 0) { message("No synergistic DOWN genes to plot."); return(NULL) }
    df <- head(df, n)
    metric_cols <- c("Decrease_c_vs_nt", "Decrease_a_vs_nt", "Decrease_b_vs_nt")
    metric_label <- "Decrease (1 - FC)"
  }

  labels <- synergy_res$params$labels

  long <- data.frame(
    gene  = factor(rep(df$gene_name, 3),
                   levels = rev(unique(df$gene_name))),
    comp  = rep(c(paste0(labels["c"], " vs ", labels["nt"]),
                  paste0(labels["a"], " vs ", labels["nt"]),
                  paste0(labels["b"], " vs ", labels["nt"])),
                each = nrow(df)),
    value = c(df[[metric_cols[1]]], df[[metric_cols[2]]], df[[metric_cols[3]]])
  )

  comp_colors <- c(
    setNames(SYNERGY_COLORS["c"], paste0(labels["c"], " vs ", labels["nt"])),
    setNames(SYNERGY_COLORS["a"], paste0(labels["a"], " vs ", labels["nt"])),
    setNames(SYNERGY_COLORS["b"], paste0(labels["b"], " vs ", labels["nt"]))
  )

  p <- ggplot2::ggplot(long, ggplot2::aes(x = value, y = gene, fill = comp)) +
    ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.7),
                      width = 0.65) +
    ggplot2::scale_fill_manual(values = comp_colors) +
    ggplot2::labs(
      x = metric_label,
      y = "",
      fill = "Comparison",
      title = paste0("Top ", direction, " synergistic genes")
    ) +
    ggplot2::theme_minimal(base_size = 15) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 16, face = "bold"),
      axis.text.y = ggplot2::element_text(face = "italic", size = 13),
      axis.text.x = ggplot2::element_text(size = 12),
      axis.title.x = ggplot2::element_text(size = 14),
      legend.text = ggplot2::element_text(size = 12),
      legend.title = ggplot2::element_text(size = 13),
      legend.position = "bottom",
      panel.grid.major.y = ggplot2::element_blank()
    )

  p
}

#' Venn/Euler diagram of gene overlaps across key comparisons
#'
#' @param synergy_res A synergy_result object
#' @return A ggplot or euler plot object
#' @export
plot_synergy_overlap <- function(synergy_res) {
  merged <- synergy_res$merged_all
  p_cutoff <- synergy_res$params$p_cutoff

  # Use the same significance column the synergy call used (P vs Q).
  pfx <- if (isTRUE(synergy_res$params$use_qvalue)) "Qvalue_" else "Pvalue_"
  sig_c_vs_nt <- merged$gene_id[merged[[paste0(pfx, "c_vs_nt")]] < p_cutoff]
  sig_a_vs_nt <- merged$gene_id[merged[[paste0(pfx, "a_vs_nt")]] < p_cutoff]
  sig_b_vs_nt <- merged$gene_id[merged[[paste0(pfx, "b_vs_nt")]] < p_cutoff]

  synergy_ids <- c(synergy_res$synergy_up$gene_id, synergy_res$synergy_down$gene_id)

  sets <- list(
    "C vs NT"        = sig_c_vs_nt,
    "A vs NT"        = sig_a_vs_nt,
    "B vs NT"        = sig_b_vs_nt,
    "Synergy"        = synergy_ids
  )

  if (requireNamespace("ggVennDiagram", quietly = TRUE)) {
    ggVennDiagram::ggVennDiagram(sets, label = "count", edge_size = 0.5) +
      ggplot2::scale_fill_gradient(low = "white", high = SYNERGY_COLORS["synergy"]) +
      ggplot2::theme(legend.position = "none")
  } else {
    message("Install 'ggVennDiagram' for Venn diagrams. Returning NULL.")
    NULL
  }
}

#' Heatmap of top synergistic genes
#'
#' Uses log2FC values by default. If \code{use_fpkm = TRUE} and the input files
#' contained \code{*_FPKM} columns, the heatmap uses log2(FPKM + 1) values
#' instead. If FPKM data are requested but unavailable, the function warns and
#' falls back to log2FC values.
#'
#' @param synergy_res A synergy_result object
#' @param n_top Number of top genes per direction (default 20)
#' @param use_fpkm Whether to use FPKM data (requires full input columns)
#' @return A pheatmap object
#' @export
plot_synergy_heatmap <- function(synergy_res, n_top = 20, use_fpkm = FALSE) {
  if (!requireNamespace("pheatmap", quietly = TRUE)) {
    stop("Package 'pheatmap' is required for heatmaps")
  }

  up_genes   <- head(synergy_res$synergy_up$gene_id, n_top)
  down_genes <- head(synergy_res$synergy_down$gene_id, n_top)
  all_genes  <- c(up_genes, down_genes)

  if (length(all_genes) == 0) { message("No synergistic genes for heatmap."); return(NULL) }

  labels <- synergy_res$params$labels
  sel <- synergy_res$merged_all[
    match(all_genes, synergy_res$merged_all$gene_id),
    c("gene_id", "gene_name"),
    drop = FALSE
  ]
  sel <- sel[!is.na(sel$gene_id), , drop = FALSE]
  if (nrow(sel) == 0) {
    message("No selected synergistic genes found in merged results.")
    return(NULL)
  }

  # Direction is decided per gene_id (unique). Gene symbols can be duplicated
  # across paralogues, so build unique but still-readable row labels - assigning
  # duplicated names as rownames would otherwise error.
  direction  <- ifelse(sel$gene_id %in% synergy_res$synergy_up$gene_id, "UP", "DOWN")
  row_labels <- make.unique(sel$gene_name)

  fpkm_mat <- synergy_res$fpkm_all
  fpkm_cols <- if (is.null(fpkm_mat)) character() else setdiff(colnames(fpkm_mat), "gene_id")
  use_fpkm_data <- isTRUE(use_fpkm) && length(fpkm_cols) > 0

  if (isTRUE(use_fpkm) && !use_fpkm_data) {
    warning("FPKM data were requested but no *_FPKM columns are available; ",
            "falling back to log2FC values.", call. = FALSE)
  }

  if (use_fpkm_data) {
    heat_mat <- fpkm_mat[match(sel$gene_id, fpkm_mat$gene_id),
                         fpkm_cols, drop = FALSE]
    heat_mat <- as.matrix(heat_mat)
    suppressWarnings(storage.mode(heat_mat) <- "numeric")
    rownames(heat_mat) <- row_labels
    heat_mat[!is.finite(heat_mat)] <- NA
    heat_mat <- log2(heat_mat + 1)
    finite_vals <- heat_mat[is.finite(heat_mat)]
    if (length(finite_vals) == 0) {
      message("No finite FPKM values for selected synergistic genes.")
      return(NULL)
    }
    heat_colors <- colorRampPalette(c("white", "#FEC44F", "#B2182B"))(100)
    heat_breaks <- if (length(unique(finite_vals)) > 1) {
      range_vals <- range(finite_vals)
      seq(range_vals[1], range_vals[2], length.out = 101)
    } else {
      center <- finite_vals[1]
      delta <- max(abs(center) * 0.01, 1e-6)
      seq(center - delta, center + delta, length.out = 101)
    }
    heat_title <- paste0("Synergistic Gene log2(FPKM + 1) Heatmap (top ",
                         n_top, " each)")
  } else {
    # Pull the log2FC rows for the selected genes, keyed by the unique gene_id.
    fc_cols <- c("log2FC_c_vs_nt", "log2FC_a_vs_nt", "log2FC_b_vs_nt",
                 "log2FC_c_vs_a", "log2FC_c_vs_b")
    fc_sel <- synergy_res$merged_all[
      match(sel$gene_id, synergy_res$merged_all$gene_id),
      fc_cols,
      drop = FALSE
    ]
    heat_mat <- as.matrix(fc_sel)
    rownames(heat_mat) <- row_labels

    # Cap extreme values
    heat_mat[is.infinite(heat_mat)] <- NA
    finite_vals <- heat_mat[is.finite(heat_mat)]
    if (length(finite_vals) == 0) {
      message("No finite log2FC values for selected synergistic genes.")
      return(NULL)
    }
    max_val <- max(abs(finite_vals), na.rm = TRUE)
    if (max_val == 0 || is.na(max_val)) max_val <- 1
    cap_val <- min(max_val, 10)
    heat_mat[heat_mat > cap_val]  <- cap_val
    heat_mat[heat_mat < -cap_val] <- -cap_val

    colnames(heat_mat) <- c(
      paste0(labels["c"], " vs ", labels["nt"]),
      paste0(labels["a"], " vs ", labels["nt"]),
      paste0(labels["b"], " vs ", labels["nt"]),
      paste0(labels["c"], " vs ", labels["a"]),
      paste0(labels["c"], " vs ", labels["b"])
    )
    heat_colors <- colorRampPalette(c("#2166AC", "white", "#B2182B"))(100)
    heat_breaks <- seq(-cap_val, cap_val, length.out = 101)
    heat_title <- paste0("Synergistic Gene log2FC Heatmap (top ", n_top, " each)")
  }

  # Annotation for direction
  has_up   <- any(direction == "UP")
  has_down <- any(direction == "DOWN")
  dir_levels <- c(if (has_up) "UP", if (has_down) "DOWN")
  dir_colors <- c(if (has_up) stats::setNames(SYNERGY_COLORS["up"], "UP"),
                  if (has_down) stats::setNames(SYNERGY_COLORS["down"], "DOWN"))

  anno <- data.frame(
    Direction = factor(direction, levels = dir_levels),
    row.names = row_labels
  )
  anno_colors <- list(Direction = dir_colors)

  heat_args <- list(
    mat           = heat_mat,
    color         = heat_colors,
    annotation_row = anno,
    annotation_colors = anno_colors,
    cluster_rows   = nrow(heat_mat) >= 2,
    cluster_cols   = FALSE,
    show_rownames  = TRUE,
    fontsize_row   = 7,
    fontface_row   = "italic",
    border_color   = NA,
    main           = heat_title,
    angle_col      = 45
  )
  if (!is.null(heat_breaks)) {
    heat_args$breaks <- heat_breaks
  }

  do.call(pheatmap::pheatmap, heat_args)
}
