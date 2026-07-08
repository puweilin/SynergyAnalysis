# synergy_plot.R — Publication-quality visualizations for synergy results

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
  synergy_ids <- c(synergy_res$synergy_up$gene_id, synergy_res$synergy_down$gene_id)

  merged$category <- "Non-significant"
  merged$category[merged$gene_id %in% synergy_res$synergy_up$gene_id]   <- "Synergy UP"
  merged$category[merged$gene_id %in% synergy_res$synergy_down$gene_id] <- "Synergy DOWN"

  # Remove Inf/-Inf for plotting
  merged$log2FC_c_vs_nt_plot <- merged$log2FC_c_vs_nt
  merged$log2FC_c_vs_nt_plot[is.infinite(merged$log2FC_c_vs_nt_plot)] <- NA
  max_fc <- max(abs(merged$log2FC_c_vs_nt_plot), na.rm = TRUE)
  cap <- min(max_fc + 1, 15)

  merged$log2FC_c_vs_nt_plot <- pmax(-cap, pmin(cap, merged$log2FC_c_vs_nt_plot))
  merged$nlog10_q <- -log10(merged$Qvalue_c_vs_nt)
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
      x = expression(log[2] * "(Fold Change) — " * C * " vs " * NT),
      y = expression(-log[10] * "(Q-value)"),
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

  sig_c_vs_nt <- merged$gene_id[merged$Qvalue_c_vs_nt < p_cutoff]
  sig_a_vs_nt <- merged$gene_id[merged$Qvalue_a_vs_nt < p_cutoff]
  sig_b_vs_nt <- merged$gene_id[merged$Qvalue_b_vs_nt < p_cutoff]

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
  } else if (requireNamespace("eulerr", quietly = TRUE)) {
    fit <- eulerr::euler(sets[1:3])
    plot(fit, quantities = TRUE, fills = unname(SYNERGY_COLORS[c("c", "a", "b")]))
  } else {
    message("Install 'ggVennDiagram' or 'eulerr' for Venn diagrams. Returning NULL.")
    NULL
  }
}

#' Heatmap of top synergistic genes across all 4 groups
#'
#' Uses FPKM or count data if available in the input files.
#' Falls back to log2FC values otherwise.
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

  # Pull the log2FC rows for the selected genes, keyed by the unique gene_id.
  fc_cols <- c("log2FC_c_vs_nt", "log2FC_a_vs_nt", "log2FC_b_vs_nt",
               "log2FC_c_vs_a", "log2FC_c_vs_b")
  sel <- synergy_res$merged_all[
    synergy_res$merged_all$gene_id %in% all_genes,
    c("gene_id", "gene_name", fc_cols),
    drop = FALSE
  ]

  # Direction is decided per gene_id (unique). Gene symbols can be duplicated
  # across paralogues, so build unique but still-readable row labels — assigning
  # duplicated names as rownames would otherwise error.
  direction  <- ifelse(sel$gene_id %in% synergy_res$synergy_up$gene_id, "UP", "DOWN")
  row_labels <- make.unique(sel$gene_name)

  fc_mat <- as.matrix(sel[, fc_cols, drop = FALSE])
  rownames(fc_mat) <- row_labels

  # Cap extreme values
  fc_mat[is.infinite(fc_mat)] <- NA
  max_val <- max(abs(fc_mat), na.rm = TRUE)
  if (max_val == 0 || is.na(max_val)) max_val <- 1
  cap_val <- min(max_val, 10)
  fc_mat[fc_mat > cap_val]  <- cap_val
  fc_mat[fc_mat < -cap_val] <- -cap_val

  colnames(fc_mat) <- c(
    paste0(labels["c"], " vs ", labels["nt"]),
    paste0(labels["a"], " vs ", labels["nt"]),
    paste0(labels["b"], " vs ", labels["nt"]),
    paste0(labels["c"], " vs ", labels["a"]),
    paste0(labels["c"], " vs ", labels["b"])
  )

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

  pheatmap::pheatmap(fc_mat,
    color         = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
    breaks        = seq(-cap_val, cap_val, length.out = 101),
    annotation_row = anno,
    annotation_colors = anno_colors,
    cluster_rows   = TRUE,
    cluster_cols   = FALSE,
    show_rownames  = TRUE,
    fontsize_row   = 7,
    fontface_row   = "italic",
    border_color   = NA,
    main           = paste0("Synergistic Gene log2FC Heatmap (top ", n_top, " each)"),
    angle_col      = 45
  )
}
